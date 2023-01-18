"""
Path tracking simulation with LQR steering control and PID speed control.
author Atsushi Sakai (@Atsushi_twi)
"""
import sys
sys.path.append("cell_test")
import cubic_spline
import numpy as np
import math
import matplotlib.pyplot as plt
import scipy.linalg as la
from utils_3d import * 

# 配置LQR 的参数
# === Parameters =====

# LQR parameter
lqr_Q = np.array([[8,0,0,0],
                 [0,0,0,0], 
                 [0,0,10,0],
                 [0,0,0,0]])
lqr_R = 10 * np.eye(1)
L = 0.53  # Wheel base of the vehicle [m]，车辆轴距
max_steer = np.deg2rad(30.0)  # maximum steering angle[rad]
max_accel = 1  
Kp = 1.0
dt = 1/50 # time tick[s]，采样时间


k_arr_dict = {
    1/20:"cell_test/K_arr_fps_20.npy",
    1/30:"cell_test/K_arr_fps_30.npy",
    1/40:"cell_test/K_arr_fps_40.npy", 
    1/50:"cell_test/K_arr_fps_50.npy",
}

K_arr = np.load(k_arr_dict[dt])
v_arr = np.linspace(0.005,6,1200) 

# State 对象表示自车的状态，位置x、y，以及横摆角yaw、速度v
class State:

    def __init__(self, x=0.0, y=0.0, yaw=0.0, v=0.0):
        self.x = x
        self.y = y
        self.yaw = yaw
        self.v = v


    # 更新自车的状态，采样时间足够小，则认为这段时间内速度相同，加速度相同，使用匀速模型更新位置
    def update(self, a, delta):
    
        if delta >= max_steer:
            delta = max_steer
        if delta <= - max_steer:
            delta = - max_steer
    
        self.x = self.x + self.v * math.cos(self.yaw) * dt
        self.y = self.y + self.v * math.sin(self.yaw) * dt
        self.yaw = self.yaw + self.v / L * math.tan(delta) * dt
        self.v = self.v + a * dt

    def pi_2_pi(self,angle):
        return (angle + math.pi) % (2 * math.pi) - math.pi


    # 实现离散Riccati equation 的求解方法
    def solve_dare(self,A, B, Q, R):
        """
        solve a discrete time_Algebraic Riccati equation (DARE)
        """
        x = Q
        x_next = Q
        max_iter = 30
        eps = 0.05
    
        for i in range(max_iter):
            # x_next = A.T @ x @ A - A.T @ x @ B @ \
            #          la.inv(R + B.T @ x @ B) @ B.T @ x @ A + Q
    
            try:
                x_next = A.T @ x @ A - A.T @ x @ B @ \
                         la.inv(np.nan_to_num(R + B.T @ x @ B)) @ B.T @ x @ A + Q
            except:
                x_next = A.T @ x @ A - A.T @ x @ B @ \
                         la.pinv(np.nan_to_num(R + B.T @ x @ B)) @ B.T @ x @ A + Q
            if (abs(x_next - x)).max() < eps:
                break
            x = x_next
    
        return x_next

    # 返回值K 即为LQR 问题求解方法中系数K的解
    def dlqr(self,A, B, Q, R):
        """Solve the discrete time lqr controller.
        x[k+1] = A x[k] + B u[k]
        cost = sum x[k].T*Q*x[k] + u[k].T*R*u[k]
        # ref Bertsekas, p.151
        """
    
        # first, try to solve the ricatti equation
        X = self.solve_dare(A, B, Q, R)
    
        # compute the LQR gain
        try:
            K = la.inv(np.nan_to_num(B.T @ X @ B + R)) @ (B.T @ X @ A)
        except:
            K = la.pinv(np.nan_to_num(B.T @ X @ B + R)) @ (B.T @ X @ A)
    
        eig_result = la.eig(A - B @ K)
    
        return K, X, eig_result[0]

    # 计算距离自车当前位置最近的参考点
    def calc_nearest_index(self, cx, cy, cyaw):
        dx = [self.x - icx for icx in cx]
        dy = [self.y - icy for icy in cy]
    
        d = [idx ** 2 + idy ** 2 for (idx, idy) in zip(dx, dy)]
    
        mind = min(d)
    
        ind = d.index(mind)
    
        mind = math.sqrt(mind)
    
        dxl = cx[ind] - self.x
        dyl = cy[ind] - self.y
    
        angle = self.pi_2_pi(cyaw[ind] - math.atan2(dyl, dxl))
        if angle < 0:
            mind *= -1
    
        return ind, mind

    def PIDControl(self,target, current):
        a = Kp * (target - current)
        return a


    def closed_loop_prediction(self,cx, cy, cyaw, ck, speed_profile, goal,cloud_proto,pe,pth_e):
        last_di = 0
        last_v = 0
        di_delta_list = []

        max_di = 30 * np.pi / 180
        iter_list = []
        iter = 0
        # 设置起点的参数
        T = 0.05  # max simulation time
        goal_dis = 0.2
        stop_speed = 0.05

        Time = 0.0
        x =   [self.x]
        y =   [self.y]
        yaw = [self.yaw]
        v =   [self.v]
        t = [0.0]
        di_list = [0]

        real_x = []
        real_y = []


        # pe, pth_e = 0.0, 0.0
        import time
        tic = time.time()

        ind, e = self.calc_nearest_index(cx, cy, cyaw)

        # print("calc_nearest time cost:",time.time() - tic)

        sp = speed_profile
        tv = sp[ind]

        k = ck[ind]
        v_state = self.v
        th_e = self.pi_2_pi(self.yaw - cyaw[ind])

        # 构建LQR表达式，X(k+1) = A * X(k) + B * u(k), 使用Riccati equation 求解LQR问题
        #     dt表示采样周期，v表示当前自车的速度
        #     A = [1.0, dt, 0.0, 0.0,
        #          0.0, 0.0, v,  0.0]
        #          0.0, 0.0, 1.0, dt]
        #          0.0, 0.0, 0.0, 0.0]
        # A = np.zeros((4, 4))
        # A[0, 0] = 1.0
        # A[0, 1] = dt
        # A[1, 2] = v_state
        # A[2, 2] = 1.0
        # A[2, 3] = dt

        # 构建B矩阵，L是自车的轴距
        # B = [0.0
        #     0.0
        #     0.0
        #     v/L]
        # B = np.zeros((4, 1))
        # B[3, 0] = v_state / L

        # K, _, _ = self.dlqr(A, B, lqr_Q, lqr_R)

        nearest_k_index = getnearpos_index(v_arr,self.v) 
        K = K_arr[nearest_k_index,:].reshape(-1,4)

        # state vector，构建状态矩阵
        # x = [e, dot_e, th_e, dot_th_e, delta_v]
        # e: lateral distance to the path， e是自车到轨迹的距离
        # dot_e: derivative of e， dot_e是自车到轨迹的距离的变化率
        # th_e: angle difference to the path， th_e是自车与期望轨迹的角度偏差
        # dot_th_e: derivative of th_e， dot_th_e是自车与期望轨迹的角度偏差的变化率
        # delta_v: difference between current speed and target speed，delta_v是当前车速与期望车速的偏差
        X = np.zeros((4, 1))
        X[0, 0] = e
        X[1, 0] = (e - pe) / dt
        X[2, 0] = th_e
        X[3, 0] = (th_e - pth_e) / dt

        # input vector，构建输入矩阵u
        # u = [delta, accel]
        # delta: steering angle，前轮转角
        # accel: acceleration，自车加速度
        ustar = -K @ X

        # calc steering input
        ff = math.atan2(L * k, 1)  # feedforward steering angle
        fb = self.pi_2_pi(ustar[0, 0])  # feedback steering angle
        delta = ff + fb

        # calc accel input
        # accel = ustar[1, 0]

        ################################ RC1 #################################### 
        # accel = self.PIDControl(sp[ind],self.v)
        # if k > 1: 
        #     accel = min(max(accel, -max_accel), max_accel)

        # # check goal
        # dx = self.x - goal[0]
        # dy = self.y - goal[1]

        # # if math.hypot(dx, dy) <= 1.5:
        # #     if(accel >= 0 and self.v > 0):
        # #         accel = -np.power(self.v,2) / (2 * math.hypot(dx, dy))

        # di, target_ind, pe, pth_e, ai = delta, ind, e, th_e, accel
        # tic = time.time()
        # self.update(ai, di)

        # Time = Time + dt
        # di = min(max(di, -max_di), max_di)

        # if 0.8 < k and  k < 1:
        #     self.v =  self.v * 0.8
        
        # if k >= 1 and k < 1.3:
        #     self.v =  self.v / (k)

        # elif k >= 1.3: 
        #     self.v  =  self.v / (k * 1.2)


        ######################################## RC2 ##############################
        self.v = sp[ind]
        accel = 0 

        # check goal
        dx = self.x - goal[0]
        dy = self.y - goal[1]

        # if math.hypot(dx, dy) <= 1.5:
        #     if(accel >= 0 and self.v > 0):
        #         accel = -np.power(self.v,2) / (2 * math.hypot(dx, dy))

        di, target_ind, pe, pth_e, ai = delta, ind, e, th_e, accel
        tic = time.time()
        self.update(ai, di)
        Time = Time + dt
        di = min(max(di, -max_di), max_di)

        if 0.8 < k and  k < 1:
            self.v =  self.v * 0.8
        
        if k >= 1 and k < 1.3:
            self.v =  self.v / (k)

        elif k >= 1.3: 
            self.v  =  self.v / (k * 1.2)

        return self.v , di,self.x,self.y,self.yaw, pe, pth_e

    def calc_speed_profile(self,cx, cy, ck, cyaw, target_speed):
        speed_profile = [target_speed] * len(cx)

        direction = 1.0

        # Set stop point
        for i in range(len(cx) - 1):
            dyaw = abs(cyaw[i + 1] - cyaw[i])
            switch = math.pi / 4.0 <= dyaw < math.pi / 2.0

            if switch:
                direction *= -1

            if direction != 1.0:
                speed_profile[i] = - target_speed
            else:
                speed_profile[i] = target_speed

            if switch:
                speed_profile[i] = 0.0

        speed_profile[-1] = 0.0

        return speed_profile

    def run(self): 
        pass 

def main():
    print("LQR steering control tracking start!!")


if __name__ == '__main__':
    main()


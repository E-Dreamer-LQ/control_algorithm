import numpy as np 
import scipy.linalg as la
import sys 
sys.path.append("../")
from utils_3d import * 

def solve_dare(A, B, Q, R):
    """
    solve a discrete time_Algebraic Riccati equation (DARE)
    """
    x = Q
    x_next = Q
    max_iter = 150
    eps = 0.01

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

def dlqr(A, B, Q, R):
    """Solve the discrete time lqr controller.
    x[k+1] = A x[k] + B u[k]
    cost = sum x[k].T*Q*x[k] + u[k].T*R*u[k]
    # ref Bertsekas, p.151
    """

    # first, try to solve the ricatti equation
    X = solve_dare(A, B, Q, R)

    # compute the LQR gain
    try:
        K = la.inv(np.nan_to_num(B.T @ X @ B + R)) @ (B.T @ X @ A)
    except:
        K = la.pinv(np.nan_to_num(B.T @ X @ B + R)) @ (B.T @ X @ A)

    eig_result = la.eig(A - B @ K)

    return K, X, eig_result[0]

k_list = []

lqr_Q = np.array([[8,0,0,0],
                 [0,0,0,0], 
                 [0,0,10,0],
                 [0,0,0,0]])
lqr_R = 10*np.eye(1)

v_list = np.linspace(0.005,6,1200) 

print(v_list)

A = np.zeros((4, 4))
B = np.zeros((4, 1))
dt = 1/50
L = 0.53

k_list = []

for v_state in v_list: 
            
    A[0, 0] = 1.0
    A[0, 1] = dt
    A[1, 2] = v_state
    A[2, 2] = 1.0
    A[2, 3] = dt

    # A = A / v_state
            
    B[3, 0] = v_state / L
    K, _, _ = dlqr(A, B, lqr_Q, lqr_R)
    K_arr = np.array(K).reshape(1,4)
    # print("K",K_arr.shape)
    k_list.append(K) 

res = np.array(k_list).reshape(-1,4)
# print("res shape", res.shape )
print("res:", res)

np.save("K_arr_fps_50.npy",res) 

nearest_k_index = getnearpos_index(v_list,1) 
K = res[nearest_k_index,:]
print(K.shape)




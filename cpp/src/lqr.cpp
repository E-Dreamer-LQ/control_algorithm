#include "controller/lqr.hpp"
#include <cmath>

using namespace std; 
using namespace Eigen; 


Lqr_controler::Lqr_controler(){
}

void Lqr_controler::Init(double x,double y,double yaw,double v,double gyro_z){
    this->x = x;
    this->y = y; 
    this->yaw = yaw; 
    this->v = v; 
    this->gyro_z = gyro_z;
}

void Lqr_controler::update(double a,double delta){
    if(delta > max_steer || delta == max_steer){
        delta = max_steer;
    }
    if (delta < -max_steer || delta == -max_steer){
        delta = -max_steer;
    } 
    this->x = this->x + this->v * cos(this->yaw) * dt; 
    this->y = this->y + this->v * sin(this->yaw) * dt; 
    this->yaw = this->yaw + (this->v/L) * tan(delta) * dt; 
    this->v = this->v + a*dt;  
} 

double Lqr_controler::pid_controll(double target_speed, double cur_speed){
    double acc = Kp * (target_speed - cur_speed); 
    return acc; 
}

double  Lqr_controler::pi_2_pi(double angle){ 
    double res = 0; 
    res = fmod((angle + PI),(2 * PI))  - PI;
    return res;
}

Eigen::MatrixXd Lqr_controler::solve_dare(const Matrix4d & A,const Eigen::Matrix<double, 4, 1> &B,const MatrixXd &Q,const Eigen::Matrix<double, 1, 1> &R ){
    Matrix4d x ; 
    Matrix4d x_next; 
    int max_iter = 500; 
    double eps = 0.001; 
    for (int i = 0;i < max_iter;i++){
        x_next = A.transpose() * x * A - A.transpose() * x * B * (R + B.transpose() * x *B).inverse() * B.transpose() * x * A + Q; 
        Eigen::ArrayXXd x_abs = (x-x_next).array().abs();
        double max_value = x_abs.maxCoeff();
        if (max_value < eps){
            break;
        }
        x = x_next;
    }

    return x_next;
}

Eigen::MatrixXd Lqr_controler::dlqr(const Matrix4d & A,const Eigen::Matrix<double, 4, 1> &B,const MatrixXd &Q,const Eigen::Matrix<double, 1, 1> &R ){
    Matrix4d X = solve_dare(A,B,Q,R); 
    MatrixXd K = (B.transpose() * X * B + R).inverse() * (B.transpose() * X * A);

    return K;
}


pair<int,double> Lqr_controler::calc_nearest_index(Eigen::VectorXd& cx, Eigen::VectorXd& cy,Eigen::VectorXd& cyaw){
    int rows = cx.rows(); 
    MatrixXd dx = MatrixXd::Ones(rows,1) * this->x; 
    MatrixXd dy = MatrixXd::Ones(rows,1) * this->y; 
    dx -= cx; 
    dy -= cy; 

    MatrixXd d; 

    d.resize(rows,2);

    d <<  dx,dy; 

    MatrixXd::Index   minIndex;
    double minNorm = d.rowwise().norm().minCoeff(&minIndex);
    
    double dxl = cx[minIndex] - this->x; 
    double dyl = cy[minIndex] - this->y; 

    double angle = pi_2_pi(cyaw[minIndex] - atan2(dxl,dyl)); 

    if (angle < 0){
        minNorm *= -1;
    }
    return make_pair(minIndex,minNorm);

}

tuple<double,double,double,double,double,double,double> Lqr_controler::closed_loop_prediction(Eigen::VectorXd &cx, Eigen::VectorXd & cy,Eigen::VectorXd & cyaw,Eigen::VectorXd & ck, Eigen::VectorXd & sp, double pe,double pth_e ){
    double ti = 0; 
    int ind = 0; 
    double e = 0;

    pair<int,double>  index = calc_nearest_index(cx,cy,cyaw); 
    
    ind = index.first; 
    e = index.second; 
    double k = ck[ind] ; // 曲率
    double th_e = pi_2_pi(this->yaw - cyaw[ind]); 

    Matrix4d A; 
    A(0, 0) = 1.0;
    A(0, 1) = dt;
    A(1, 2) = this->v; 
    A(2, 2) = 1.0;
    A(2, 3) = dt;

    Matrix<double, 4, 1> B; 
    B(3, 0) = this->v / L;

    Q(0,0) = 0.1; 
    Q(2,2) = 16; 
    R(0,0) = 1; 
    cout << "Q(0,0):" << Q(0,0) << " "  <<  "Q(2,2):" << Q(2,2) << " " <<  "R(0,0):" <<R(0,0) << endl; 
    MatrixXd K = dlqr(A,B,Q,R); 
    // auto l_idx = std::upper_bound(v_arr.begin(), v_arr.end(), this->v) - v_arr.begin() - 1;
    // Eigen::MatrixXd K = Eigen::Map<Eigen::MatrixXd,Eigen::RowMajor>(k_arr[l_idx].data(),4,1);
    // K.resize(1,4);
    cout << "K:" << K << endl;
    
    Matrix<double,4,1> X;
    // X(0, 0) = e;
    // X(1, 0) = (e - pe) / dt;
    // X(2, 0) = th_e;
    // X(3, 0) = (th_e - pth_e) / dt;

    X(0, 0) = e;
    X(1, 0) = this->v * sin(th_e);
    X(2, 0) = th_e;
    X(3, 0) = this->gyro_z - sp[ind] * k;

    MatrixXd ustar = -K * X; 
    double ff = atan2(L*k, 1); 
    double fb = pi_2_pi(ustar(0,0)); 
    double delta = ff + fb; 

    pe = e; 
    pth_e = th_e;

    // double accel = pid_controll(sp[ind],this->v); 
    double accel = 0; 
    update(accel,delta); 

    this->v = sp[ind]; 

    delta = min(max(delta,-max_steer), max_steer);

    ti = ti + dt;

    // if( k > 0 && k < 1){
    //     this->v = this->v  * (1 - k); 
    // }
    // if (k > 1){
    //     this->v  = this->v   / (k);
    // }

    tuple res(this->v,delta,this->x,this->y,this->yaw,pe,pth_e);

    return  res;
    
}


Eigen::VectorXd Lqr_controler::calc_speed_profile(Eigen::VectorXd &cx,Eigen::VectorXd & cyaw,Eigen::VectorXd &ck,double target_speed){
    int rows = cx.rows(); 
    VectorXd speed_profile = MatrixXd::Ones(rows,1) * target_speed; 
    double direction = 1.0; 
    // for (int i=0; i < rows-1;i++){
    //     double dyaw = abs(cyaw[i+1] - cyaw[i]); 
    //     bool sw = ((PI / 4.0) <= dyaw && dyaw < (PI / 2));
    //     // cout << boolalpha << sw << endl; 
    //     if(sw){
    //         direction *= -1;
    //     }
    //     if (direction != 1.0){
    //         speed_profile[i] = -target_speed; 
    //     }
    //     else{
    //         speed_profile[i] = target_speed;
    //     }
    //     if (sw){
    //         speed_profile[i] = 0;
    //     }
    // }
    for (int i=0; i < rows-1;i++){
        double k = abs(ck[i]); 
        speed_profile[i] = target_speed;
        // if( k > 0 && k < 1){
        //     speed_profile[i] = target_speed * (1 - k); 
        // }
        // if (k > 1){
        //     speed_profile[i] = target_speed  / (k);
        // }
    }
    // speed_profile[rows-1] = 0;  //  已有达到终点停止
    return speed_profile;
}

Lqr_controler::~Lqr_controler(){
}




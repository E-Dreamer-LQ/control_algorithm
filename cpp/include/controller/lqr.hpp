#pragma once
#include <mutex>
#include <thread>
#include <functional>
#include <sys/time.h>
#include <vector>
#include <Eigen/Dense>
#include "rclcpp/rclcpp.hpp"
#include "std_msgs/msg/string.hpp"
#include "help.hpp"



using namespace std::chrono_literals;
using namespace std::chrono;
using namespace Eigen;

#define PI 3.1415926



class Lqr_controler
{
private:
    double speed_con = 0;
    double eps_con = 0;
    double x = 0;
    double y = 0;
    double yaw = 0;
    double v = 0;
    double Kp = 1;
    double T = 5;
    Eigen::MatrixXd Q = Eigen::MatrixXd::Identity(4, 4);
    Eigen::MatrixXd R = Eigen::MatrixXd::Identity(1,1);
    double dt = 0.02; // 采样时间
    double L = 0.53;   // 轴距
    double max_steer = (30 * PI) / 180.0;


    double gyro_z = 0;
    Helper* help = new Helper();
    std::vector<double> v_arr =  help->linspace(0,6,1200);
    std::vector<Eigen::Vector4d> k_arr = help->readMatrixFile("install/app_controller/lib/app_controller/K_arr_fps_50.txt"); 


public:
    Lqr_controler();
    ~Lqr_controler();
    void Init(double x,double y,double yaw,double v,double gyro_z);

    void update(double a, double delta);
    double pid_controll(double target_speed, double cur_speed);

    double pi_2_pi(double angle);

    Eigen::MatrixXd solve_dare(const Matrix4d &A,const Eigen::Matrix<double, 4, 1> &B,const Eigen::MatrixXd &Q,const Eigen::Matrix<double, 1, 1> &R);

    Eigen::MatrixXd dlqr(const Matrix4d &A, const Eigen::Matrix<double, 4, 1> &B,const Eigen::MatrixXd &Q,const Eigen::Matrix<double, 1, 1> &R);

    std::pair<int, double> calc_nearest_index(Eigen::VectorXd &cx, Eigen::VectorXd &cy, Eigen::VectorXd &cyaw);

    std::tuple<double,double,double,double,double,double,double> closed_loop_prediction(Eigen::VectorXd &cx, Eigen::VectorXd &cy, Eigen::VectorXd &cyaw, Eigen::VectorXd &ck, Eigen::VectorXd &sp, double pe, double pth_e);

    Eigen::VectorXd calc_speed_profile(Eigen::VectorXd &cx, Eigen::VectorXd &cyaw,Eigen::VectorXd &ck, double target_speed);
};

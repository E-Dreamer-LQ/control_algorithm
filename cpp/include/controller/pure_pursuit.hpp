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



class PP_controller
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
    double dt = 0.02; // 采样时间
    double L = 0.53;   // 轴距
    double max_steer = (30 * PI) / 180.0;
    double gyro_z = 0;
    double Lfc = 1.5;  // 前视距离
    double kk = 0.1;

public:
    PP_controller();
    ~PP_controller();
    void Init(double x,double y,double yaw,double v,double gyro_z);

    void update(double a, double delta);
    double pid_controll(double target_speed, double cur_speed);
    double pi_2_pi(double angle);

    std::pair<int, double> calc_nearest_index(Eigen::VectorXd &cx, Eigen::VectorXd &cy, Eigen::VectorXd &cyaw,Eigen::VectorXd& ck);

    std::tuple<double,double,double,double,double,int> pp_steer_control(Eigen::VectorXd &cx, Eigen::VectorXd &cy, Eigen::VectorXd &cyaw, Eigen::VectorXd &ck, Eigen::VectorXd &sp,int ind);

    Eigen::VectorXd calc_speed_profile(Eigen::VectorXd &cx, Eigen::VectorXd &cyaw,Eigen::VectorXd &ck, double target_speed);
};


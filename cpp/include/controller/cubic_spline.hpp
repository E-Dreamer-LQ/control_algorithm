#pragma once
#include <iostream> 
#include <vector> 
#include <thread> 
#include <opencv2/core.hpp>
#include <opencv2/imgcodecs.hpp>
#include "opencv2/opencv.hpp"
#include "opencv2/videoio.hpp"
#include <opencv2/videoio/registry.hpp>
#include "spline.h"
#include <Eigen/Dense>

using namespace std;



class Spline_process
{
private:
    std::vector<double> X; 
    std::vector<double> Y; 

public:
    // Spline_process(std::vector<double>X,std::vector<double> Y);
    Spline_process();
    ~Spline_process();
    
    Eigen::MatrixXd calc_spline_course(std::vector<double> &X,std::vector<double> &Y); 
    
    std::vector<double> calc_s(std::vector<double>& X,std::vector<double>& Y); 

    std::vector<double> arange(double s,double ds); 

    void csv_save(Eigen::MatrixXd &res,string savename); 

    double get_yaw(tk::spline &spX, tk::spline &spY,double i_s); 

    double get_curvature(tk::spline &spX, tk::spline &spY,double i_s); 

};





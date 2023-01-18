#pragma once
#include <iostream>
#include <Eigen/Dense>
#include <vector>
#include <fstream>
#include <zconf.h>
#include <string>
#include <sstream>



class Helper
{

public:
    Helper();
    ~Helper();
    std::vector<double> linspace(double start, double end, int num); 
    std::vector<Eigen::Vector4d> readMatrixFile(const char *fileName) ; 
    void remove(std::vector<int> &v);

};

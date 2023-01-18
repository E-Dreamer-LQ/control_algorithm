#include "controller/cubic_spline.hpp"
#include <numeric>
#include "controller/help.hpp"
#include <algorithm>
#include <unordered_map> 
#include "cubic.cpp" 

using namespace std; 
using namespace Eigen; 
#define PI 3.1415926

Spline_process::Spline_process(){

}

Spline_process::~Spline_process(){

}


double Spline_process::get_curvature(tk::spline &spX, tk::spline &spY,double i_s){
    double derivX = 0; 
    double second_derivX = 0; 
    double derivY = 0; 
    double second_derivY = 0;
    double k = 0; 
    derivX = spX.deriv(1,i_s); 
    second_derivX = spX.deriv(2,i_s); 
    derivY = spY.deriv(1,i_s); 
    second_derivY = spY.deriv(2,i_s); 
    k = (second_derivY * derivX - second_derivX * derivY) / (derivX * derivX + derivY * derivY);

    return k;

}

double Spline_process::get_yaw(tk::spline &spX, tk::spline &spY,double i_s){
    double yaw = 0; 
    double derivX = 0; 
    double derivY = 0; 
    derivX = spX.deriv(1,i_s); 
    derivY = spY.deriv(1,i_s); 
    yaw = atan2(derivY,derivX);
    yaw = fmod((yaw + PI),(2 * PI))  - PI;  // pi_2_pi 
    return yaw; 
}

vector<double> Spline_process::calc_s(std::vector<double>& X,std::vector<double>& Y){
    // double s = 0;
    std::vector <double> x0(X.begin(),X.end()-1); 
    std::vector <double> x1(X.begin()+1,X.end());
    std::vector <double> dx;
    std::transform(x0.begin(), x0.end(), x1.begin(),std::back_inserter(dx),std::minus<double>());
    std::vector <double> y0(Y.begin(),Y.end()-1); 
    std::vector <double> y1(Y.begin()+1,Y.end());
    std::vector <double> dy;
    std::transform(y0.begin(),y0.end(),y1.begin(),std::back_inserter(dy),std::minus<double>());

    Eigen::MatrixXd matrix_dx = Eigen::Map<Eigen::MatrixXd>(dx.data(), dx.size(),1);
    Eigen::MatrixXd matrix_dy = Eigen::Map<Eigen::MatrixXd>(dy.data(), dy.size(),1);

    MatrixXd d; 
    d.resize(dx.size(),2);
    d <<  matrix_dx,matrix_dy; 

    MatrixXd sqrt_norm = d.rowwise().norm();
    sqrt_norm.resize(dx.size(),1);
    vector<double> vec_norm(sqrt_norm.data(),sqrt_norm.data()+sqrt_norm.rows()*sqrt_norm.cols());
    // s = accumulate(vec_norm.begin(),vec_norm.end(),0.0);

    return vec_norm; 
}


vector<double> Spline_process::arange(double s,double ds){
    vector<double> res; 
    double sum = 0; 
    while(sum < s){
        res.push_back(sum); 
        sum += ds; 
    }
    return res;
}


void Spline_process::csv_save(Eigen::MatrixXd &res,string savename){
    ofstream dataFile;
	dataFile.open(savename, ios::out | ios::trunc);
    dataFile << "x_input"
        << ","
        << "y_input"
        << ","
        << "k"
        << ","
        << "yaw"
		<< std::endl;
	for (int i = 0; i < res.rows(); i++)
	{
		for (int j = 0; j < res.cols(); j++)
		{
			dataFile << res(i,j) << ",";          // 写入数据
            // cout << res(i,j) << ","; 
		}
		dataFile <<  endl;                       // 换行
        // cout << endl; 
	}
	
	dataFile.close();                            // 关闭文档
}


MatrixXd Spline_process::calc_spline_course(std::vector<double> &X,std::vector<double> &Y){

    vector<double> vec_ds = calc_s(X, Y);
    double s_length =  accumulate(vec_ds.begin(),vec_ds.end(),0.0);
    cout << "s_length: " << s_length << endl;
    double ds = 0.5;

    const int N = vec_ds.size(); 
    vector<double> S(N+1);
    S[0] = 0; 
    for (int i = 0; i < N; i++)
        S[i + 1] = S[i] + vec_ds[i];

    vector<double> s_list = arange(s_length,ds);
    
    vector<double> X_C(X); 
    vector<double> S_C(S); 
    tk::spline spX(S,X_C); 
    vector<double> Y_C(Y); 
    tk::spline spY(S_C,Y_C);

    unordered_map<int,int> hashmap;

    vector<double> new_x; 
    vector<double> new_y; 
    vector<double> new_curvature; 
    vector<double> new_cyaw; 

    int idx = 0;
    double last_yaw = 0;
    double last_x = X[0]; 
    double last_y = Y[0]; 

    for(int j = 0; j < s_list.size(); j++){ 
        double i_k = get_curvature(spX,spY,s_list[j]); 
        double i_yaw = get_yaw(spX,spY,s_list[j]);
        double i_x; 
        double i_y;
        i_x = last_x + ds * cos(i_yaw);  
        i_y = last_y + ds * sin(i_yaw);  
        last_yaw = i_yaw; 
        last_x = i_x; 
        last_y = i_y; 
        // cout << "i_x << i_y << i_k << i_yaw:" << endl;
        // cout << i_x << " " <<  i_y << " "<< i_k  << " "<< i_yaw << endl; 
        new_x.push_back(i_x);
        new_y.push_back(i_y);   
        new_curvature.push_back(i_k); 
        new_cyaw.push_back(i_yaw);
    }

    MatrixXd res; 
    Eigen::MatrixXd Vec_xx = Eigen::Map<Eigen::MatrixXd>(new_x.data(), new_x.size(),1);
    Eigen::MatrixXd Vec_yy = Eigen::Map<Eigen::MatrixXd>(new_y.data(), new_y.size(),1);
    Eigen::MatrixXd Vec_curva = Eigen::Map<Eigen::MatrixXd>(new_curvature.data(), new_curvature.size(),1);
    Eigen::MatrixXd Vec_yaw = Eigen::Map<Eigen::MatrixXd>(new_cyaw.data(), new_cyaw.size(),1);
    res.resize(Vec_xx.size(),4);
    // cout << "length = " << length << endl;
    res << Vec_xx,Vec_yy,Vec_curva,Vec_yaw; 
    csv_save(res,"sp.csv"); 

    cout << "csv file saved " << endl;
    return res; 
}










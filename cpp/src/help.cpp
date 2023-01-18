#include "controller/help.hpp"
 
using namespace Eigen;
using namespace std;



Helper::Helper(){

}

std::vector<double> Helper::linspace(double start, double end, int num)
{

  std::vector<double> linspaced;

  if (num == 0) { return linspaced; }
  if (num == 1) 
    {
      linspaced.push_back(start);
      return linspaced;
    }
  
  double delta = (end - start) / (num - 1);

  for(int i=0; i < num-1; ++i)
    {
      linspaced.push_back(start + delta * i);
    }
  linspaced.push_back(end); // I want to ensure that start and end
                            // are exactly the same as the input
  return linspaced;
}

vector<Eigen::Vector4d> Helper::readMatrixFile(const char *fileName) {
 
    //读取文件
    ifstream file;
    stringstream ss;
    
    vector<Eigen::Vector4d>  velocity;
    Vector4d vec_tmp;
    file.open(fileName, ios::in);
    if (!file.is_open())
    {
        cout << "read file failed" << endl;
        return velocity;
    }
    string line;   
    while(getline(file, line)){
        ss.str(line);
        string single;
        int  i = 0 ;
        // 按照空格分隔
        while(getline(ss, single, ' ')){
            vec_tmp(i++) = atof(single.c_str());
        }
        ss.clear(); //必须加，不然写不到string里。
        velocity.push_back(vec_tmp);
    }
    return velocity;
}


void Helper::remove(std::vector<int> &v)
{
    auto end = v.end();
    for (auto it = v.begin(); it != end; ++it) {
        end = std::remove(it + 1, end, *it);
    }

    v.erase(end, v.end());
}

Helper::~Helper(){

}

// int main() {
//     std::cout << "Hello, World!" << std::endl;
//     // std::vector<std::vector<double >> matrixALL = readMatrixFile("K_arr_fps_20.txt");
//     vector<Eigen::Vector4d> matrixALL = readMatrixFile("K_arr_fps_20.txt");
//     // cout << Vec_xx.rows() << std::endl; 

//     for(int i=0;i<1200;i++){
//         // for (int j = 0; j < 4; ++j) {
//         //     cout << matrixALL[i][j] << " ";
//         // }
//         Eigen::MatrixXd Vec_xx = Eigen::Map<Eigen::MatrixXd,Eigen::RowMajor>(matrixALL[i].data(),4,1);
//         Vec_xx.resize(4,1);
//         cout << Vec_xx.rows() << " " << Vec_xx.cols() << endl;
//     }
//     return 0;
// }





#include<cmath>
#include<eigen3/Eigen/Core>
#include<eigen3/Eigen/Dense>
#include<iostream>

int main(){
    float x = 2.0, y = 1.0;
    Eigen::Vector3f origin(x, y, 1.0);
    float pi = acos(-1);
    float sin45 = sin(45.0/180.0*pi);
    float cos45 = cos(45.0/180.0*pi);
    Eigen::Matrix3f rotate45_translate12;
    rotate45_translate12 << cos45, -sin45, 1.0, sin45, cos45, 2.0, 0.0, 0.0, 1.0;
    std::cout << "rotate45_translate12 \n";
    std::cout << rotate45_translate12 << std::endl;
    Eigen::Vector3f result;
    result = rotate45_translate12 * origin;
    std::cout << "Result: \n";
    std::cout << result << std::endl;

    return 0;
}
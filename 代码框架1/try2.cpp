#include<cmath>
#include<eigen3/Eigen/Core>
#include<eigen3/Eigen/Dense>
#include<iostream>
#include<string>

int main()
{   float rotation_angle = 45;
    Eigen::Vector3f rotation_axis(1,3,5);
    float pi = acos(-1);
    float sinrot = sin(rotation_angle/180.0f*pi);
    float cosrot = cos(rotation_angle/180.0f*pi);
    Eigen::Matrix4f model = Eigen::Matrix4f::Identity();
    //model << cosrot, -sinrot, 0.0,0.0, sinrot, cosrot, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1;
    //Eigen::Vector3f rotation_axis(1.0f,2.0f,3.0f);
    rotation_axis = rotation_axis.normalized();
    float nx = rotation_axis[0];
    float ny = rotation_axis[1];
    float nz = rotation_axis[2];
    Eigen::Vector3f n = rotation_axis;
    Eigen::Matrix3f view = Eigen::Matrix3f::Identity();
    std::cout << n << std::endl;
    std::cout << n.transpose() << std::endl;
    //std::cout << n*n.transpose() << std::endl;
    Eigen::Matrix3f N;
    N << 0, -nz, ny, nz, 0, -nx, -ny ,nx, 0;
    //std::cout << (1-cosrot)*n*n.transpose() << std::endl;
    N = cosrot*view + (1-cosrot)*n*n.transpose() + sinrot*N;
    model.block(0,0,3,3)= N.block(0,0,3,3);
    std::cout << model << std::endl;

    return 0;
}
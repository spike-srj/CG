#include<cmath>
#include<eigen3/Eigen/Core>
#include<eigen3/Eigen/Dense>
#include<iostream>

int main()
{   
    //note: through we have use float to denfine sin45/cos45, we still have to add .0 benhind 45/180 as well as matrix and so on
    // create rotation matrix,contrarotate 45 degree
    float pi = acos(-1);
    float sin45 = sin(45.0f/180.0f*pi);
    float cos45 = cos(45.0f/180.0f*pi);
    Eigen::Matrix3f i;
    i << cos45, -sin45, 1.0, sin45, cos45, 2.0, 0.0, 0.0, 1.0;
    //defiction vector
    Eigen::Vector3f v(2.0f,1.0f,1.0f);
    //compute
    std::cout << "Example of output \n";
    std::cout << i*v << std::endl;
    return 0;
}
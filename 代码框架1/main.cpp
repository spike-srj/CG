#include "Triangle.hpp"
#include "rasterizer.hpp"
#include <eigen3/Eigen/Eigen>
#include <iostream>
#include <opencv2/opencv.hpp>

constexpr double MY_PI = 3.1415926;

Eigen::Matrix4f get_view_matrix(Eigen::Vector3f eye_pos)
{
    Eigen::Matrix4f view = Eigen::Matrix4f::Identity();

    Eigen::Matrix4f translate;
    translate << 1, 0, 0, -eye_pos[0], 0, 1, 0, -eye_pos[1], 0, 0, 1,
        -eye_pos[2], 0, 0, 0, 1;

    view = translate * view;

    return view;
}

Eigen::Matrix4f get_model_matrix(float rotation_angle, Eigen::Vector3f rotation_axis)
{
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
    
    // TODO: Implement this function
    // Create the model matrix for rotating the triangle around the Z axis.
    // Then return it.

    return model;
}

Eigen::Matrix4f get_projection_matrix(float eye_fov, float aspect_ratio,
                                      float zNear, float zFar)
{
    // Students will implement this function
    float n = zNear;
    float f = zFar;
    float t = tan(eye_fov/2)*n;
    float r = t*aspect_ratio;
    float l = -r;
    float b = -t;
    Eigen::Matrix4f projection = Eigen::Matrix4f::Identity();
    //projection << 
    projection << 2*n/(r-l),0,(l+r)/(l-r),0,0,2*n/(t-b),(b+t)/(b-t),0,0,0,(f+n)/(n-f),2*f*n/(f-n),0,0,1,0;
    // TODO: Implement this function
    // Create the projection matrix for the given parameters.
    // Then return it.

    return projection;
}

int main(int argc, const char** argv)
{
    float angle = 0;
    bool command_line = false;
    std::string filename = "output.png";

    if (argc >= 3) {
        command_line = true;
        angle = std::stof(argv[2]); // -r by default
        if (argc == 4) {
            filename = std::string(argv[3]);
        }
        else
            return 0;
    }

    rst::rasterizer r(700, 700);//实例化光栅器r，分辨率为700*700.构造frame_buf和depth_buf（见hpp）

    Eigen::Vector3f eye_pos = {0, 0, 5};//look from z axis and the distance is 5 

    std::vector<Eigen::Vector3f> pos{{2, 0, -2}, {0, 2, -2}, {-2, 0, -2}};//给定三个点

    std::vector<Eigen::Vector3i> ind{{0, 1, 2}};
    
    //然后把pos 和ind 传给光栅器，得到它们的ID,这里做的相当于顶点着色器
    auto pos_id = r.load_positions(pos);//pos是三个点，返回的id就代表三个点
    auto ind_id = r.load_indices(ind);

    int key = 0;
    int frame_count = 0;

    if (command_line) {
        r.clear(rst::Buffers::Color | rst::Buffers::Depth);

        r.set_model(get_model_matrix(angle, Eigen::Vector3f(1,0,0)));
        r.set_view(get_view_matrix(eye_pos));
        r.set_projection(get_projection_matrix(45, 1, 0.1, 50));

        r.draw(pos_id, ind_id, rst::Primitive::Triangle);
        cv::Mat image(700, 700, CV_32FC3, r.frame_buffer().data());
        image.convertTo(image, CV_8UC3, 1.0f);

        cv::imwrite(filename, image);

        return 0;
    }

    while (key != 27) //除非键盘输入“Esc”结束程序，否则一直循环
    {
        //清除画布，也就是把frame_buf 和 depth_buf 初始化为全0 和全无穷大
        r.clear(rst::Buffers::Color | rst::Buffers::Depth);

        r.set_model(get_model_matrix(angle, Eigen::Vector3f(1,0,0)));
        r.set_view(get_view_matrix(eye_pos));
        r.set_projection(get_projection_matrix(45, 1, 0.1, 50));
        
        //调用绘制函数，输入顶点ID和索引ID以及绘制图元方法
        r.draw(pos_id, ind_id, rst::Primitive::Triangle);
        
        
        //第二点，OpenCV如何把一个点向量转化为一张图片的。
        //首先用Mat构造一个image，Mat(nrows,ncols,type,fillValue)
        //r.frame_buffer().data() 返回指针，指向vector的第一个元素的地址 C++11
        //然后转换一下类型，从32位3通道Float类型转换为8位3通道unsigned类型
        //然后显示出来，每隔10ms刷新。
        cv::Mat image(700, 700, CV_32FC3, r.frame_buffer().data());
        image.convertTo(image, CV_8UC3, 1.0f);
        cv::imshow("image", image);
        key = cv::waitKey(10);

        std::cout << "frame count: " << frame_count++ << '\n';

        if (key == 'a') {
            angle += 10;
        }
        else if (key == 'd') {
            angle -= 10;
        }
    }

    return 0;
}

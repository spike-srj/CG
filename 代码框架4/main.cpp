#include <chrono>
#include <iostream>
#include <opencv2/opencv.hpp>
#define AA true
std::vector<cv::Point2f> control_points;

float aa_color(float px,float py,float cx,float cy){
    if(sqrt(pow((cx-px),2)+pow((cy-py),2))>=1)
    return 1;
    else
    return 255*(1-sqrt(pow((cx-px),2)+pow((cy-py),2)));
}

void mouse_handler(int event, int x, int y, int flags, void *userdata) 
{
    if (event == cv::EVENT_LBUTTONDOWN && control_points.size() < 7) 
    {
        std::cout << "Left button of the mouse is clicked - position (" << x << ", "
        << y << ")" << '\n';
        control_points.emplace_back(x, y);
    }     
}

void naive_bezier(const std::vector<cv::Point2f> &points, cv::Mat &window) 
{
    auto &p_0 = points[0];
    auto &p_1 = points[1];
    auto &p_2 = points[2];
    auto &p_3 = points[3];

    for (double t = 0.0; t <= 1.0; t += 0.001) 
    {
        auto point = std::pow(1 - t, 3) * p_0 + 3 * t * std::pow(1 - t, 2) * p_1 +
                 3 * std::pow(t, 2) * (1 - t) * p_2 + std::pow(t, 3) * p_3;

        window.at<cv::Vec3b>(point.y, point.x)[2] = 255;
    }
}

cv::Point2f recursive_bezier(const std::vector<cv::Point2f> &control_points, float t) 
{
    // TODO: Implement de Casteljau's algorithm
    
    int n = control_points.size();
    if (control_points.size() == 2){
        return (1-t)*control_points[0] +t*control_points[1];
    }
    std::vector<cv::Point2f> vec;
    for(int i =0; i<n-1; i++){
        vec.push_back((1-t)*control_points[i]+t*control_points[i+1]);
    }
    return recursive_bezier(vec,t);  //cv::Point2f();

}

void bezier(const std::vector<cv::Point2f> &control_points, cv::Mat &window) 
{
    // TODO: Iterate through all t = 0 to t = 1 with small steps, and call de Casteljau's 
    // recursive Bezier algorithm.
    int n = sizeof(control_points);
    for(double t = 0.0; t <= 1.0; t += 0.001){
        auto point = recursive_bezier(control_points,t);
        if(AA){
            float x_min=std::floor(point.x)+0.5;
            float x_max=std::ceil(point.x)+0.5;
            float y_min=std::floor(point.y)+0.5;
            float y_max=std::ceil(point.y)+0.5;

            float d=aa_color(x_min,y_min,point.x,point.y);
            if(window.at<cv::Vec3b>(y_min, x_min)[1]<d)
            	window.at<cv::Vec3b>(y_min, x_min)[1]= d;
            d=aa_color(x_max,y_max,point.x,point.y);
            if(window.at<cv::Vec3b>(y_max, x_max)[1]<d)
            	window.at<cv::Vec3b>(y_max, x_max)[1]= d;
            d=aa_color(x_min,y_max,point.x,point.y);
            if(window.at<cv::Vec3b>(y_max, x_min)[1]<d)
            	window.at<cv::Vec3b>(y_max, x_min)[1]= d;
            d=aa_color(x_max,y_min,point.x,point.y);
            if(window.at<cv::Vec3b>(y_min, x_max)[1]<d)
            	window.at<cv::Vec3b>(y_min, x_max)[1]= d;
        }
        else{
            window.at<cv::Vec3b>(point.y, point.x)[1] = 255; 
        }
    }

}

int main() 
{
    cv::Mat window = cv::Mat(700, 700, CV_8UC3, cv::Scalar(0));
    cv::cvtColor(window, window, cv::COLOR_BGR2RGB);
    cv::namedWindow("Bezier Curve", cv::WINDOW_AUTOSIZE);

    cv::setMouseCallback("Bezier Curve", mouse_handler, nullptr);

    int key = -1;
    while (key != 27) 
    {
        for (auto &point : control_points) 
        {
            cv::circle(window, point, 3, {255, 255, 255}, 3);
        }

        if (control_points.size() == 7) 
        {
            //naive_bezier(control_points, window);
            bezier(control_points, window);

            cv::imshow("Bezier Curve", window);
            cv::imwrite("my_bezier_curve.png", window);
            key = cv::waitKey(0);

            return 0;
        }

        cv::imshow("Bezier Curve", window);
        key = cv::waitKey(20);
    }

return 0;
}

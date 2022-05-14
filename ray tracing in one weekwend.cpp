#include <iostream>

int main() {
    const int image_width = 200;
    const int image_height = 100;

    std::cout << "P3\n" << image_width << ' ' << image_height << "\n255\n";
    //对于像素来说, 每一行是从左往右写入的。行从上开始往下写入的。
    for (int j = image_height-1; j >= 0; --j) {
        for (int i = 0; i < image_width; ++i) {
            //通常我们把RGB通道的值限定在0.0到1.0。我们之后计算颜色值的时候将使用一个动态的范围, 
            //这个范围并不是0到1。但是在使用这段代码输出图像之前, 我们将把颜色映射到0到1。所以这部分输出图像代码不会改变。
            auto r = double(i) / image_width;
            auto g = double(j) / image_height;
            auto b = 0.2;
            int ir = static_cast<int>(255.999 * r);
            int ig = static_cast<int>(255.999 * g);
            int ib = static_cast<int>(255.999 * b);
            std::cout << ir << ' ' << ig << ' ' << ib << '\n';
        }
    }
    //显示进度
    std::cerr << "\nDone.\n";
}

//对于linux或mac用户，先run一边上面的代码，会生成一个同名文件main，然后在当前目录的终端输入./main > image.ppm，就可生成一张图片
./main > image.ppm

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
创建vec3.h头文件，该文件在mac和win上运行会报错，只有虚拟机上可以不知道为啥

#include <iostream>
#include <cmath>
class vec3 {
    public:
        //构造函数，第一行是无参构造函数，第二行是有参构造函数，将e0，e1，e2传入了e{}中
        vec3() : e{0,0,0} {}
        vec3(double e0, double e1, double e2) : e{e0, e1, e2} {}
        //类的函数
        double x() const { return e[0]; }
        double y() const { return e[1]; }
        double z() const { return e[2]; }
        //操作符重载，无参操作
        vec3 operator-() const { return vec3(-e[0], -e[1], -e[2]); }
        //没看懂
        double operator[](int i) const { return e[i]; }
        double& operator[](int i) { return e[i]; }
        //操作符重载，有参操作
        vec3& operator+=(const vec3 &v) {
            e[0] += v.e[0];
            e[1] += v.e[1];
            e[2] += v.e[2];
            return *this;
        }
        //操作符重载，有参操作
        vec3& operator*=(const double t) {
            e[0] *= t;
            e[1] *= t;
            e[2] *= t;
            return *this;
        }
        //操作符重载，有参操作
        vec3& operator/=(const double t) {
            return *this *= 1/t;
        }
        //函数
        double length() const {
            return sqrt(length_squared());
        }
        //函数
        double length_squared() const {
            return e[0]*e[0] + e[1]*e[1] + e[2]*e[2];
        }
        //没看懂
        void write_color(std::ostream &out) {
            // Write the translated [0,255] value of each color component.
            out << static_cast<int>(255.999 * e[0]) << ' '
                << static_cast<int>(255.999 * e[1]) << ' '
                << static_cast<int>(255.999 * e[2]) << '\n';
        }

    public:
        double e[3];  //属性
};
//第二部分
// inline 是什么，这里也是操作符重载麻
inline std::ostream& operator<<(std::ostream &out, const vec3 &v) {
    return out << v.e[0] << ' ' << v.e[1] << ' ' << v.e[2];
}

inline vec3 operator+(const vec3 &u, const vec3 &v) {
    return vec3(u.e[0] + v.e[0], u.e[1] + v.e[1], u.e[2] + v.e[2]);
}

inline vec3 operator-(const vec3 &u, const vec3 &v) {
    return vec3(u.e[0] - v.e[0], u.e[1] - v.e[1], u.e[2] - v.e[2]);
}

inline vec3 operator*(const vec3 &u, const vec3 &v) {
    return vec3(u.e[0] * v.e[0], u.e[1] * v.e[1], u.e[2] * v.e[2]);
}

inline vec3 operator*(double t, const vec3 &v) {
    return vec3(t*v.e[0], t*v.e[1], t*v.e[2]);
}

inline vec3 operator*(const vec3 &v, double t) {
    return t * v;
}

inline vec3 operator/(vec3 v, double t) {
    return (1/t) * v;
}

inline double dot(const vec3 &u, const vec3 &v) {
    return u.e[0] * v.e[0]
         + u.e[1] * v.e[1]
         + u.e[2] * v.e[2];
}

inline vec3 cross(const vec3 &u, const vec3 &v) {
    return vec3(u.e[1] * v.e[2] - u.e[2] * v.e[1],
                u.e[2] * v.e[0] - u.e[0] * v.e[2],
                u.e[0] * v.e[1] - u.e[1] * v.e[0]);
}

inline vec3 unit_vector(vec3 v) {
    return v / v.length();
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//可以使用vec3类将我们的main函数改成这样
#include "vec3.h"

#include <iostream>

int main() {
    const int image_width = 200;
    const int image_height = 100;

    std::cout << "P3\n" << image_width << ' ' << image_height << "\n255\n";

    for (int j = image_height-1; j >= 0; --j) {
        std::cerr << "\rScanlines remaining: " << j << ' ' << std::flush;
        for (int i = 0; i < image_width; ++i) {
            vec3 color(double(i)/image_width, double(j)/image_height, 0.2);
            color.write_color(std::cout);
        }
    }

    std::cerr << "\nDone.\n";
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//创建ray.h
#ifndef RAY_H
#define RAY_H

#include "vec3.h"

class ray {
    public:
        //无参构造函数
        ray() {}
        //有参构造函数，将origin传给属性orig
        ray(const vec3& origin, const vec3& direction)
            : orig(origin), dir(direction)
        {}
        //函数，用来读取两个属性
        vec3 origin() const    { return orig; }
        vec3 direction() const { return dir; }
        //函数，需要传入参数t，来计算光线t时的位置
        vec3 at(double t) const {
            return orig + t*dir;
        }

    public:
        vec3 orig;  //属性，表示光的起点
        vec3 dir;   //属性，表示光的方向
};
#endif

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//更新main.cpp
#include "ray.h"

#include <iostream>
//ray_color(ray)函数根据y值将蓝白做了个线性插值的混合, 我们这里把射线做了个单位化, 以保证y的取值范围(-1.0<y<1.0)。因为我们使用y轴做渐变, 所以你可以看到这个蓝白渐变也是竖直的。
//为什么要将y单位化为[-1,1]
vec3 ray_color(const ray& r) {
    //unit_vector出现在vec3里,是用来归一化vec的，
    vec3 unit_direction = unit_vector(r.direction());
    //接下来使用了一个标准的小技巧将y的范围从-1.0 ≤ y ≤ 1.0映射到了0 ≤ y ≤ 1.0。这样t=1.0时就是蓝色, 而t=0.0时就是白色
    auto t = 0.5*(unit_direction.y() + 1.0);
    //现在我们采用的是线性混合(linear blend)或者说线性插值(liner interpolation)。或者简称其为lerp。一个lerp一般来说会是下面的形式:(1.0-t)*startvalue + t*endvalue
}
    return (1.0-t)*vec3(1.0, 1.0, 1.0) + t*vec3(0.5, 0.7, 1.0);
}

int main() {
    //图像的宽高
    const int image_width = 200;
    const int image_height = 100;

    std::cout << "P3\n" << image_width << " " << image_height << "\n255\n";
    //发出射线的原点从图像的左下角开始沿着xy方向做增量直至遍历全图
    vec3 lower_left_corner(-2.0, -1.0, -1.0);
    //纬度上限
    vec3 horizontal(4.0, 0.0, 0.0);
    //经度上限
    vec3 vertical(0.0, 2.0, 0.0);
    //摄像机/光线的起点，第一条光线应该是（-2，-1，-1）
    vec3 origin(0.0, 0.0, 0.0);
    //从左上角开始遍历
    //左下角坐标加上高就是左上角的位置
    for (int j = image_height-1; j >= 0; --j) {
        std::cerr << "\rScanlines remaining: " << j << ' ' << std::flush;
        for (int i = 0; i < image_width; ++i) {
            auto u = double(i) / image_width;
            auto v = double(j) / image_height;
            ray r(origin, lower_left_corner + u*horizontal + v*vertical);
            //此时光线方向尚未归一化，在ray_color处会归一化
            vec3 color = ray_color(r);
            color.write_color(std::cout);
        }
    }

    std::cerr << "\nDone.\n";
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

bool hit_sphere(const vec3& center, double radius, const ray& r) {
    vec3 oc = r.origin() - center;
    auto a = dot(r.direction(), r.direction());
    auto b = 2.0 * dot(oc, r.direction());
    auto c = dot(oc, oc) - radius*radius;
    auto discriminant = b*b - 4*a*c;
    return (discriminant > 0);
}

vec3 ray_color(const ray& r) {
    if (hit_sphere(vec3(0,0,-1), 0.5, r))
        return vec3(1, 0, 0);
    vec3 unit_direction = unit_vector(r.direction());
    auto t = 0.5*(unit_direction.y() + 1.0);
    return (1.0-t)*vec3(1.0, 1.0, 1.0) + t*vec3(0.5, 0.7, 1.0);
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//main.cc 球体表面法相
double hit_sphere(const vec3& center, double radius, const ray& r) {
    vec3 oc = r.origin() - center;
    auto a = dot(r.direction(), r.direction());
    auto b = 2.0 * dot(oc, r.direction());
    auto c = dot(oc, oc) - radius*radius;
    auto discriminant = b*b - 4*a*c;
    if (discriminant < 0) {
        return -1.0;
    } else {
        //相交时小的t
        return (-b - sqrt(discriminant) ) / (2.0*a);
    }
}

vec3 ray_color(const ray& r) {
    auto t = hit_sphere(vec3(0,0,-1), 0.5, r);
    if (t > 0.0) {
        //r.at(t)得到相交位置的坐标，用该位置减去摄像机的位置，再归一化为[-1,1]。之前是将y单位化，这里是将交点和相机的向量单位化
        vec3 N = unit_vector(r.at(t) - vec3(0,0,-1));
        return 0.5*vec3(N.x()+1, N.y()+1, N.z()+1);
    }
    vec3 unit_direction = unit_vector(r.direction());
    t = 0.5*(unit_direction.y() + 1.0);
    return (1.0-t)*vec3(1.0, 1.0, 1.0) + t*vec3(0.5, 0.7, 1.0);
}


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
#ifndef VEC3_H
#define VEC3_H
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
        void write_color(std::ostream &out, int samples_per_pixel) {
            // Divide the color total by the number of samples.
            auto scale = 1.0 / samples_per_pixel;
            auto r = scale * e[0];
            auto g = scale * e[1];
            auto b = scale * e[2];

            // Write the translated [0,255] value of each color component.
            out << static_cast<int>(256 * clamp(r, 0.0, 0.999)) << ' '
                << static_cast<int>(256 * clamp(g, 0.0, 0.999)) << ' '
                << static_cast<int>(256 * clamp(b, 0.0, 0.999)) << '\n';
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
#endif

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
            //此时r传入的参数是屏幕空间的（一个4*2的矩形），x~（-2，2）；y~（-1，1），因此要在上面将分辨率空间映射到屏幕空间！！
            //先将分辨率空间通过除法转化为标准空间（0-1），再乘以屏幕空间的宽高，转化到屏幕空间
            //raster space ——> NDC space ——> screen space; 后面还有——> world space
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

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//为了在场景中渲染不止一个球，使用一个抽象类, 任何可能与光线求交的东西实现时都继承这个类, 并且让球以及球列表也都继承这个类，命名为hittable
//hittable.h
#ifndef HITTABLE_H
#define HITTABLE_H
#include "ray.h"

//计算的结果存在一个结构体里
struct hit_record {
    vec3 p;
    vec3 normal;
    double t;
};

class hittable {
    public:
        //hittable类有个接受射线为参数的函数
        //为了便利, 加入了一个区间 [t_min，t_max] 来判断相交是否有效
        virtual bool hit(const ray& r, double t_min, double t_max, hit_record& rec) const = 0;
};

#endif
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


//这是继承自它的sphere球体类:
//sphere.h
#ifndef SPHERE_H
#define SPHERE_H

#include "hittable.h"
#include "vec3.h"

class sphere: public hittable {
    public:
        sphere() {}
        sphere(vec3 cen, double r) : center(cen), radius(r) {};
        //虚函数，具体定义在下面
        virtual bool hit(const ray& r, double tmin, double tmax, hit_record& rec) const;

    public:
        vec3 center;
        double radius;
};

bool sphere::hit(const ray& r, double t_min, double t_max, hit_record& rec) const {
    vec3 oc = r.origin() - center;
    auto a = r.direction().length_squared();
    auto half_b = dot(oc, r.direction());
    auto c = oc.length_squared() - radius*radius;
    auto discriminant = half_b*half_b - a*c;

    if (discriminant > 0) {
        auto root = sqrt(discriminant);
        auto temp = (-half_b - root)/a;
        if (temp < t_max && temp > t_min) {
            rec.t = temp;  //计算的t
            rec.p = r.at(rec.t);  //交点位置
            rec.normal = (rec.p - center) / radius;  //交点法线
            return true;
        }
        //并不是要求两个交，而是，如果上一个没有满足给定的条件，就再算另一个
        temp = (-half_b + root) / a;
        if (temp < t_max && temp > t_min) {
            rec.t = temp;
            rec.p = r.at(rec.t);  
            rec.normal = (rec.p - center) / radius;
            return true;
        }
    }
    return false;
}
#endif

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool front_face;
//入射光线与发现如果是同方向的，说明光线是从球内部射到了内壁上
if (dot(ray_direction, outward_normal) > 0.0) {
    // ray is inside the sphere
    //为了统一让法线与光线反向
    normal = -outward_normal;
    front_face = false;
}
else {
    // ray is outside the sphere
    normal = outward_normal;
    front_face = true;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//hittable.h 加入时间与面朝向
ifndef HITTABLE_H
#define HITTABLE_H

#include "ray.h"

struct hit_record {
    vec3 p;
    vec3 normal;
    double t;
    bool front_face;

    inline void set_face_normal(const ray& r, const vec3& outward_normal) {
        front_face = dot(r.direction(), outward_normal) < 0;
        normal = front_face ? outward_normal :-outward_normal;
    }
};

class hittable {
    public:
        virtual bool hit(const ray& r, double t_min, double t_max, hit_record& rec) const = 0;
};
#endif

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//sphere.h 加入射入面判别
bool sphere::hit(const ray& r, double t_min, double t_max, hit_record& rec) const {
    vec3 oc = r.origin() - center;
    auto a = r.direction().length_squared();
    auto half_b = dot(oc, r.direction());
    auto c = oc.length_squared() - radius*radius;
    auto discriminant = half_b*half_b - a*c;

    if (discriminant > 0) {
        auto root = sqrt(discriminant);
        auto temp = (-half_b - root)/a;
        if (temp < t_max && temp > t_min) {
            rec.t = temp;
            rec.p = r.at(rec.t);
            vec3 outward_normal = (rec.p - center) / radius;
            rec.set_face_normal(r, outward_normal);
            return true;
        }
        temp = (-half_b + root) / a;
        if (temp < t_max && temp > t_min) {
            rec.t = temp;
            rec.p = r.at(rec.t);
            vec3 outward_normal = (rec.p - center) / radius;
            rec.set_face_normal(r, outward_normal);
            return true;
        }
    }
    return false;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//hittable_list.h
#ifndef HITTABLE_LIST_H
#define HITTABLE_LIST_H

#include "hittable.h"
#include <memory>
#include <vector>

using std::shared_ptr;
using std::make_shared;

class hittable_list: public hittable {
    public:
        hittable_list() {}
        hittable_list(shared_ptr<hittable> object) { add(object); }

        void clear() { objects.clear(); }
        void add(shared_ptr<hittable> object) { objects.push_back(object); }

        virtual bool hit(const ray& r, double tmin, double tmax, hit_record& rec) const;

    public:
        std::vector<shared_ptr<hittable>> objects;
};

bool hittable_list::hit(const ray& r, double t_min, double t_max, hit_record& rec) const {
    hit_record temp_rec;
    bool hit_anything = false;
    auto closest_so_far = t_max;

    for (const auto& object : objects) {
        //此处的object->hit调用的hit不是hittable_list::hit，而是object的hit，如果object是sphere，则调用sphere::hit
        if (object->hit(r, t_min, closest_so_far, temp_rec)) {
            hit_anything = true;
            //此时temp_rec.t已经在调用object->hit时更新了。closest_so_far一开始设为最大的时间，经理一次更新后可能是最小的t也可能不是，但上面判断object->hit中的时间上限会更新
            //因此closest_so_far只会越来越小，它代表了与光线相交的第一个物体，后面的物体就被更新掉了
            closest_so_far = temp_rec.t;
            rec = temp_rec;
        }
    }

    return hit_anything;
}

#endif
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//shared_ptr<type>是指向一些已分配内存的类型的指针。每当你将它的值赋值给另一个智能指针时, 物体的引用计数器就会+1。当智能指针离开它所在的生存范围(例如代码块或者函数外), 物体的引用计数器就会-1。一旦引用计数器为0, 即没有任何智能指针指向该物体时, 该物体就会被销毁

//一般来说, 智能指针首先由一个刚刚新分配好内存的物体来初始化:

shared_ptr<double> double_ptr = make_shared<double>(0.37);
shared_ptr<vec3>   vec3_ptr   = make_shared<vec3>(1.414214, 2.718281, 1.618034);
shared_ptr<sphere> sphere_ptr = make_shared<sphere>(vec3(0,0,0), 1.0);

//make_shared<thing>(thing_constructor_params ...)为指定的类型分配一段内存, 使用你指定的构造函数与参数来创建这个类, 并返回一个智能指针shared_ptr<thing>

//使用C++的auto类型关键字, 可以自动推断make_shared<type>返回的智能指针类型, 于是我们可以把上面的代码简化为:

auto double_ptr = make_shared<double>(0.37);
auto vec3_ptr   = make_shared<vec3>(1.414214, 2.718281, 1.618034);
auto sphere_ptr = make_shared<sphere>(vec3(0,0,0), 1.0);

//我们在代码中使用智能指针的目的是为了能让多个几何图元共享一个实例(举个栗子, 一堆不同球体使用同一个纹理材质), 并且这样内存管理比起普通的指针更加的简单方便。

//std::shared_ptr在头文件<memory>中

#include<memory>
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


//定义一些常用的常数
//rtweekend.h
#ifndef RTWEEKEND_H
#define RTWEEKEND_H

#include <cmath>
#include <cstdlib>
#include <limits>
#include <memory>


// Usings

using std::shared_ptr;
using std::make_shared;

// Constants

const double infinity = std::numeric_limits<double>::infinity();
const double pi = 3.1415926535897932385;

// Utility Functions

inline double degrees_to_radians(double degrees) {
    return degrees * pi / 180;
}

inline double ffmin(double a, double b) { return a <= b ? a : b; }
inline double ffmax(double a, double b) { return a >= b ? a : b; }
//生成随机数
//一个简单的实现方法是, 使用<cstdlib>中的rand()函数。这个函数会返回0到RAND_MAX中的一个任意整数。
//无参函数
inline double random_double() {
    // Returns a random real in [0,1).
    return rand() / (RAND_MAX + 1.0);
}
//有参函数
inline double random_double(double min, double max) {
    // Returns a random real in [min,max).
    return min + (max-min)*random_double();
}

inline double clamp(double x, double min, double max) {
    if (x < min) return min;
    if (x > max) return max;
    return x;
}

class camera {
    public:
        camera() {
            lower_left_corner = vec3(-2.0, -1.0, -1.0);
            horizontal = vec3(4.0, 0.0, 0.0);
            vertical = vec3(0.0, 2.0, 0.0);
            origin = vec3(0.0, 0.0, 0.0);
        }

        ray get_ray(double u, double v) {
            return ray(origin, lower_left_corner + u*horizontal + v*vertical - origin);
        }

    public:
        vec3 origin;
        vec3 lower_left_corner;
        vec3 horizontal;
        vec3 vertical;
};


// Common Headers

#include "ray.h"
#include "vec3.h"

#endif
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//main.cc
#include "rtweekend.h"

#include "hittable_list.h"
#include "sphere.h"

#include <iostream>
vec3 ray_color(const ray& r, const hittable& world) {
    hit_record rec;
    if (world.hit(r, 0, infinity, rec)) {
        return 0.5 * (rec.normal + vec3(1,1,1));
    }
    vec3 unit_direction = unit_vector(r.direction());
    auto t = 0.5*(unit_direction.y() + 1.0);
    return (1.0-t)*vec3(1.0, 1.0, 1.0) + t*vec3(0.5, 0.7, 1.0);
}

int main() {
    const int image_width = 200;
    const int image_height = 100;

    std::cout << "P3\n" << image_width << ' ' << image_height << "\n255\n";

    vec3 lower_left_corner(-2.0, -1.0, -1.0);
    vec3 horizontal(4.0, 0.0, 0.0);
    vec3 vertical(0.0, 2.0, 0.0);
    vec3 origin(0.0, 0.0, 0.0);

    hittable_list world;
    world.add(make_shared<sphere>(vec3(0,0,-1), 0.5));
    world.add(make_shared<sphere>(vec3(0,-100.5,-1), 100));

    for (int j = image_height-1; j >= 0; --j) {
        std::cerr << "\rScanlines remaining: " << j << ' ' << std::flush;
        for (int i = 0; i < image_width; ++i) {
            auto u = double(i) / image_width;
            auto v = double(j) / image_height;
            ray r(origin, lower_left_corner + u*horizontal + v*vertical);

            vec3 color = ray_color(r, world);

            color.write_color(std::cout);
        }
    }

    std::cerr << "\nDone.\n";
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//我们需要一个算法来生成球体内的随机点。我们会采用最简单的做法:否定法(rejection method)。
//首先, 在一个xyz取值范围为-1到+1的单位立方体中选取一个随机点, 如果这个点在球外就重新生成直到该点在球内:
//vec3.h
class vec3 {
  public:
    ...
    inline static vec3 random() {
        return vec3(random_double(), random_double(), random_double());
    }

    inline static vec3 random(double min, double max) {
        return vec3(random_double(min,max), random_double(min,max), random_double(min,max));
    }
    vec3 random_in_unit_sphere() {
        while (true) {
            auto p = vec3::random(-1,1);
            if (p.length_squared() >= 1) continue;
            return p;
        }
    }
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
//main.cpp  
vec3 ray_color(const ray& r, const hittable& world, int depth) {
    hit_record rec;

    // If we've exceeded the ray bounce limit, no more light is gathered.
    if (depth <= 0)
        return vec3(0,0,0);

    if (world.hit(r, 0, infinity, rec)) {
        vec3 target = rec.p + rec.normal + random_in_unit_sphere();
        return 0.5 * ray_color(ray(rec.p, target - rec.p), world, depth-1);
    }

    vec3 unit_direction = unit_vector(r.direction());
    auto t = 0.5*(unit_direction.y() + 1.0);
    return (1.0-t)*vec3(1.0, 1.0, 1.0) + t*vec3(0.5, 0.7, 1.0);
} 
int main() {
    const int image_width = 200;
    const int image_height = 100;
    const int samples_per_pixel = 100;
    const int max_depth = 50;
    std::cout << "P3\n" << image_width << " " << image_height << "\n255\n";

    hittable_list world;
    world.add(make_shared<sphere>(vec3(0,0,-1), 0.5));
    world.add(make_shared<sphere>(vec3(0,-100.5,-1), 100));
    camera cam;
    for (int j = image_height-1; j >= 0; --j) {
        std::cerr << "\rScanlines remaining: " << j << ' ' << std::flush;
        for (int i = 0; i < image_width; ++i) {
            vec3 color(0, 0, 0);
            for (int s = 0; s < samples_per_pixel; ++s) {
                auto u = (i + random_double()) / image_width;
                auto v = (j + random_double()) / image_height;
                ray r = cam.get_ray(u, v);
               color += ray_color(r, world, max_depth);
            }
            color.write_color(std::cout, samples_per_pixel);
        }
    }

    std::cerr << "\nDone.\n";
}   
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//拒绝法生成的点是单位球体积内的的随机点, 这样生成的向量大概率上会和法线方向相近, 并且极小概率会沿着入射方向反射回去。这个分布律的表达式有一个 [公式] 的系数, 
//其中 [公式] 是反射光线距离法向量的夹角。这样当光线从一个离表面很小的角度射入时, 也会散射到一片很大的区域, 对最终颜色值的影响也会更低。
//然而, 事实上的lambertian的分布律并不是这样的, 它的系数是 [公式] 。真正的lambertian散射后的光线距离法相比较近的概率会更高, 
//但是分布律会更加均衡。这是因为我们选取的是单位球面上的点。我们可以通过在单位球内选取一个随机点, 然后将其单位化来获得该点。
//上面的公式部分没看懂
//vec3.h
//采用的是极坐标，圆心z坐标为0
vec3 random_unit_vector() {
    auto a = random_double(0, 2*pi);
    auto z = random_double(-1, 1);
    auto r = sqrt(1 - z*z);
    return vec3(r*cos(a), r*sin(a), z);
}
    
//vec3.h
vec3 random_in_hemisphere(const vec3& normal) {
    vec3 in_unit_sphere = random_in_unit_sphere();
    if (dot(in_unit_sphere, normal) > 0.0) // In the same hemisphere as the normal
        return in_unit_sphere;
    else
        return -in_unit_sphere;
}
//不同随机数生成函数代表不同漫反射
    
//vec3.h
vec3 ray_color(const ray& r, const hittable& world, int depth) {
    hit_record rec;

    // If we've exceeded the ray bounce limit, no more light is gathered.
    if (depth <= 0)
        return vec3(0,0,0);

    if (world.hit(r, 0.001, infinity, rec)) {
        //target是光线经过漫反射后的方向
        vec3 target = rec.p + random_in_hemisphere(rec.normal);
        //光线可能会被漫反射多次才进入相机，这里只考虑反射次数对能量的衰减因此，因此反推也是一样
        return 0.5 * ray_color(ray(rec.p, target - rec.p), world, depth-1);
    }

    vec3 unit_direction = unit_vector(r.direction());
    auto t = 0.5*(unit_direction.y() + 1.0);
    return (1.0-t)*vec3(1.0, 1.0, 1.0) + t*vec3(0.5, 0.7, 1.0);
}

    
    
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    
//设计并封装一个抽象的材质类    
//material.h
class material {
    public:
        virtual bool scatter(
            const ray& r_in, const hit_record& rec, vec3& attenuation, ray& scattered
        ) const = 0;
};

    
    
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    
//hit_record存储着光线与模型相交点的信息，包括材质，这些信息是在计算hit时更新的
//hittable.h
#ifndef HITTABLE_H
#define HITTABLE_H

#include "rtweekend.h"

class material;

struct hit_record {
    vec3 p;
    vec3 normal;
    shared_ptr<material> mat_ptr;
    double t;
    bool front_face;


    inline void set_face_normal(const ray& r, const vec3& outward_normal) {
        front_face = dot(r.direction(), outward_normal) < 0;
        normal = front_face ? outward_normal :-outward_normal;
    }
};

class hittable {
    public:
        virtual bool hit(const ray& r, double t_min, double t_max, hit_record& rec) const = 0;
};

#endif
    
    
    
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    
class sphere: public hittable {
    public:
        sphere() {}
        //初始化时加上材质指针
        sphere(vec3 cen, double r, shared_ptr<material> m)
            : center(cen), radius(r), mat_ptr(m) {};

        virtual bool hit(const ray& r, double tmin, double tmax, hit_record& rec) const;

    public:
        vec3 center;
        double radius;
        shared_ptr<material> mat_ptr;
};

bool sphere::hit(const ray& r, double t_min, double t_max, hit_record& rec) const {
    vec3 oc = r.origin() - center;
    auto a = r.direction().length_squared();
    auto half_b = dot(oc, r.direction());
    auto c = oc.length_squared() - radius*radius;
    auto discriminant = half_b*half_b - a*c;

    if (discriminant > 0) {
        auto root = sqrt(discriminant);
        auto temp = (-half_b - root)/a;
        if (temp < t_max && temp > t_min) {
            rec.t = temp;
            rec.p = r.at(rec.t);
            vec3 outward_normal = (rec.p - center) / radius;
            rec.set_face_normal(r, outward_normal);
            //更新了材质指针
            rec.mat_ptr = mat_ptr;
            return true;
        }
        temp = (-half_b + root) / a;
        if (temp < t_max && temp > t_min) {
            rec.t = temp;
            rec.p = r.at(rec.t);
            vec3 outward_normal = (rec.p - center) / radius;                
            rec.set_face_normal(r, outward_normal);
            //更新了材质指针
            rec.mat_ptr = mat_ptr;
            return true;
        }
    }
    return false;
}

    
    
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    
//定义一种材质，它的散射函数是这样的    
//material.h
class lambertian : public material {
    public:
        //初始化材质，给albedo赋值
        lambertian(const vec3& a) : albedo(a) {}

        virtual bool scatter(
            const ray& r_in, const hit_record& rec, vec3& attenuation, ray& scattered
        ) const {
            vec3 scatter_direction = rec.normal + random_unit_vector();
            scattered = ray(rec.p, scatter_direction);
            //attenuation是衰减
            attenuation = albedo;
            return true;
        }

    public:
        vec3 albedo;
};

    
    
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    
//定义反射函数    
//vec3.h
vec3 reflect(const vec3& v, const vec3& n) {
    return v - 2*dot(v,n)*n;
}

    
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//定义金属材质
//material.h
class metal : public material {
    public:
        metal(const vec3& a) : albedo(a) {}

        virtual bool scatter(
            const ray& r_in, const hit_record& rec, vec3& attenuation, ray& scattered
        ) const {
            vec3 reflected = reflect(unit_vector(r_in.direction()), rec.normal);
            scattered = ray(rec.p, reflected);
            //给反照率赋值
            attenuation = albedo;
            return (dot(scattered.direction(), rec.normal) > 0);
        }

    public:
        vec3 albedo;
};

    
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//main.cc
//
vec3 ray_color(const ray& r, const hittable& world, int depth) {
    hit_record rec;

    // If we've exceeded the ray bounce limit, no more light is gathered.
    if (depth <= 0)
        return vec3(0,0,0);

    if (world.hit(r, 0.001, infinity, rec)) {
        ray scattered;
        //attenuation是反照率，这里初始化，在材质类里传值
        vec3 attenuation;
        //如果该光线与物体相交的话，根据相交处的材质信息判断光线在相交处是否有反射或散射
        //如果该光线与物体相交且不发生反射,则将这一块颜色设为vec3(0,0,0)，相当于相机看不见。因为模型本身不发光，如果最终光线落到模型上而不是落到光源，则认为相机看到的是黑色
        if (rec.mat_ptr->scatter(r, rec, attenuation, scattered))
            //如果有散射或反射，再次判断散射后的光线是否与物体相交
            //attenuation是反照率，可以理解为辐射衰减率
            return attenuation * ray_color(scattered, world, depth-1);
        return vec3(0,0,0);
    }
    //如果该光线与背景相交或者反射后与背景相交（没有与物体相交），则将背景颜色赋予这一块像素，反射后的颜色会乘上相应的衰减率
    vec3 unit_direction = unit_vector(r.direction());
    auto t = 0.5*(unit_direction.y() + 1.0);
    return (1.0-t)*vec3(1.0, 1.0, 1.0) + t*vec3(0.5, 0.7, 1.0);
}

    
    
//给场景加入金属球
//main.cc
int main() {
    const int image_width = 200;
    const int image_height = 100;
    const int samples_per_pixel = 100;
    const int max_depth = 50;

    std::cout << "P3\n" << image_width << " " << image_height << "\n255\n";

    hittable_list world;
    //两个漫反射材质的球，make_shared<lambertian>(vec3(0.7, 0.3, 0.3))是初始化材质，给albedo赋值，(0.7, 0.3, 0.3)分别是RGB的反照率
    world.add(make_shared<sphere>(
        vec3(0,0,-1), 0.5, make_shared<lambertian>(vec3(0.7, 0.3, 0.3))));

    world.add(make_shared<sphere>(
        vec3(0,-100.5,-1), 100, make_shared<lambertian>(vec3(0.8, 0.8, 0.0))));
    //两个金属材质的球
    world.add(make_shared<sphere>(vec3(1,0,-1), 0.5, make_shared<metal>(vec3(0.8, 0.6, 0.2))));
    world.add(make_shared<sphere>(vec3(-1,0,-1), 0.5, make_shared<metal>(vec3(0.8, 0.8, 0.8))));

    camera cam;
    for (int j = image_height-1; j >= 0; --j) {
        std::cerr << "\rScanlines remaining: " << j << ' ' << std::flush;
        for (int i = 0; i < image_width; ++i) {
            vec3 color(0, 0, 0);
            for (int s = 0; s < samples_per_pixel; ++s) {
                auto u = (i + random_double()) / image_width;
                auto v = (j + random_double()) / image_height;
                ray r = cam.get_ray(u, v);
                color += ray_color(r, world, max_depth);
            }
            color.write_color(std::cout, samples_per_pixel);
        }
    }

    std::cerr << "\nDone.\n";
}

    
//给金属材质添加粗糙度，由反射方向的随机球半径表示，半径越大越有可能反射到球体内部，相当于被球吸收    
//material.h
class metal : public material {
    public:
        metal(const vec3& a, double f) : albedo(a), fuzz(f < 1 ? f : 1) {}

        virtual bool scatter(
            const ray& r_in, const hit_record& rec, vec3& attenuation, ray& scattered
        ) const {
            vec3 reflected = reflect(unit_vector(r_in.direction()), rec.normal);
            //reflected是从向交点出射的向量（单位向量），向量再加上一个单位向量，可能会达到球内部，相当于被吸收了。fuzz越大，射入球的概率越大
            scattered = ray(rec.p, reflected + fuzz*random_in_unit_sphere());
            attenuation = albedo;
            return (dot(scattered.direction(), rec.normal) > 0);//dot<0我们认为吸收
        }

    public:
        vec3 albedo;
        double fuzz;
};

    
//折射光线计算    
//vec3.h
vec3 refract(const vec3& uv, const vec3& n, double etai_over_etat) {
    auto cos_theta = dot(-uv, n);
    vec3 r_out_parallel =  etai_over_etat * (uv + cos_theta*n);
    vec3 r_out_perp = -sqrt(1.0 - r_out_parallel.length_squared()) * n;
    return r_out_parallel + r_out_perp;
}

    
    
//只会发生折射的材质    
//material.h
class dielectric : public material {
    public:
        dielectric(double ri) : ref_idx(ri) {}

        virtual bool scatter(
            const ray& r_in, const hit_record& rec, vec3& attenuation, ray& scattered
        ) const {
            attenuation = vec3(1.0, 1.0, 1.0);
            double etai_over_etat;
            //如果光线从材质外面射入，就要用外部的折射率除以材质的
            if (rec.front_face) {
                etai_over_etat = 1.0 / ref_idx;
            } 
            //如果光线从材质里面射入，材质的折射率为n，外面的折射率为n‘。
            else {
                etai_over_etat = ref_idx;
            }

            vec3 unit_direction = unit_vector(r_in.direction());
            vec3 refracted = refract(unit_direction, rec.normal, etai_over_etat);
            scattered = ray(rec.p, refracted);
            return true;
        }

        double ref_idx;
};

 
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    
//一个可以在不发生折射时反射的材质    
//material.h
class dielectric : public material {
    public:
        dielectric(double ri) : ref_idx(ri) {}

        virtual bool scatter(
            const ray& r_in, const hit_record& rec, vec3& attenuation, ray& scattered
        ) const {
            attenuation = vec3(1.0, 1.0, 1.0);
            double etai_over_etat = (rec.front_face) ? (1.0 / ref_idx) : (ref_idx);

            vec3 unit_direction = unit_vector(r_in.direction());
            double cos_theta = ffmin(dot(-unit_direction, rec.normal), 1.0);
            double sin_theta = sqrt(1.0 - cos_theta*cos_theta);
            //当光线从光密介质进入光疏介质中如果入射角大于某个临界值的时候，就会发生全反射现象，此处光线的传播路径可以参考https://www.cnblogs.com/lv-anchoret/p/10217719.html
            if (etai_over_etat * sin_theta > 1.0 ) {
                vec3 reflected = reflect(unit_direction, rec.normal);
                scattered = ray(rec.p, reflected);
                return true;
            }

            vec3 refracted = refract(unit_direction, rec.normal, etai_over_etat);
            scattered = ray(rec.p, refracted);
            return true;
        }

    public:
        double ref_idx;
};

    
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    
//main 定义的材质
world.add(make_shared<sphere>(
    vec3(0,0,-1), 0.5, make_shared<lambertian>(vec3(0.1, 0.2, 0.5))));

world.add(make_shared<sphere>(
    vec3(0,-100.5,-1), 100, make_shared<lambertian>(vec3(0.8, 0.8, 0.0))));

world.add(make_shared<sphere>(vec3(1,0,-1), 0.5, make_shared<metal>(vec3(0.8, 0.6, 0.2), 0.0)));
world.add(make_shared<sphere>(vec3(-1,0,-1), 0.5, make_shared<dielectric>(1.5)));
    
    
    
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    
//折射率会随角度变化
//material.h
double schlick(double cosine, double ref_idx) {
    auto r0 = (1-ref_idx) / (1+ref_idx);
    r0 = r0*r0;
    return r0 + (1-r0)*pow((1 - cosine),5);
}
    
class dielectric : public material {
    public:
        dielectric(double ri) : ref_idx(ri) {}

        virtual bool scatter(
            const ray& r_in, const hit_record& rec, vec3& attenuation, ray& scattered
        ) const {
            attenuation = vec3(1.0, 1.0, 1.0);
            double etai_over_etat = (rec.front_face) ? (1.0 / ref_idx) : (ref_idx);

            vec3 unit_direction = unit_vector(r_in.direction());
            double cos_theta = ffmin(dot(-unit_direction, rec.normal), 1.0);
            double sin_theta = sqrt(1.0 - cos_theta*cos_theta);
            if (etai_over_etat * sin_theta > 1.0 ) {
                vec3 reflected = reflect(unit_direction, rec.normal);
                scattered = ray(rec.p, reflected);
                return true;
            }
            double reflect_prob = schlick(cos_theta, etai_over_etat);
            if (random_double() < reflect_prob)
            {
                vec3 reflected = reflect(unit_direction, rec.normal);
                scattered = ray(rec.p, reflected);
                return true;
            }
            vec3 refracted = refract(unit_direction, rec.normal, etai_over_etat);
            scattered = ray(rec.p, refracted);
            return true;
        }

    public:
        double ref_idx;
};
    
    
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    
//此时得到的图像是实心玻璃球，会将背景位置颠倒
//这里有个简单又好用的trick, 如果你将球的半径设为负值, 形状看上去并没什么变化, 但是法相全都翻转到内部去了。
//所以就可以用这个特性来做出一个空心的玻璃球    
world.add(make_shared<sphere>(vec3(0,0,-1), 0.5, make_shared<lambertian>(vec3(0.1, 0.2, 0.5))));
world.add(make_shared<sphere>(
    vec3(0,-100.5,-1), 100, make_shared<lambertian>(vec3(0.8, 0.8, 0.0))));
world.add(make_shared<sphere>(vec3(1,0,-1), 0.5, make_shared<metal>(vec3(0.8, 0.6, 0.2), 0.3)));
world.add(make_shared<sphere>(vec3(-1,0,-1), 0.5, make_shared<dielectric>(1.5)));
world.add(make_shared<sphere>(vec3(-1,0,-1), -0.45, make_shared<dielectric>(1.5)));

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//camera.h
//为相机添加fov和纵横比，fov改变会使得画面变化，最直接的就是lower_left_corner变化
class camera {
    public:
        camera(
            double vfov, // top to bottom, in degrees
            double aspect
        ) {
            origin = vec3(0.0, 0.0, 0.0);

            auto theta = degrees_to_radians(vfov);
            auto half_height = tan(theta/2);
            auto half_width = aspect * half_height;

            lower_left_corner = vec3(-half_width, -half_height, -1.0);

            horizontal = vec3(2*half_width, 0.0, 0.0);
            vertical = vec3(0.0, 2*half_height, 0.0);
        }

        ray get_ray(double u, double v) {
            return ray(origin, lower_left_corner + u*horizontal + v*vertical - origin);
        }

    public:
        vec3 origin;
        vec3 lower_left_corner;
        vec3 horizontal;
        vec3 vertical;
};
    
    
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//进一步为相机添加起始点和朝向以及vup。此时画面（lower_left_corner）会随着起点和朝向以及fov的不同而不同
class camera {
    public:
        camera(
            vec3 lookfrom, vec3 lookat, vec3 vup,
            double vfov, // top to bottom, in degrees
            double aspect
        ) {
            origin = lookfrom;
            vec3 u, v, w;

            auto theta = degrees_to_radians(vfov);
            //这里默认相机z坐标（距离成像平面）为-1，固定不变。下面加入焦距后会变
            auto half_height = tan(theta/2);
            auto half_width = aspect * half_height;
            w = unit_vector(lookfrom - lookat);
            u = unit_vector(cross(vup, w));
            v = cross(w, u);
            //half_width、half_width只是个常数，而不是三维向量，因此计算时要乘以三维坐标的单位向量
            //u、v、w不是（1，0，0）（0，1，0）这种单位向量，而是根据上面计算lookfrom - lookat得到的向量再单位化
            lower_left_corner = origin - half_width*u - half_width*v - w;
            
            horizontal = 2*half_width*u;
            vertical = 2*half_height*v;
        }

        ray get_ray(double s, double t) {
            return ray(origin, lower_left_corner + s*horizontal + t*vertical - origin);
        }

    public:
        vec3 origin;
        vec3 lower_left_corner;
        vec3 horizontal;
        vec3 vertical;
};
    
    
//main.cc
const auto aspect_ratio = double(image_width) / image_height;
...
camera cam(vec3(-2,2,1), vec3(0,0,-1), vup, 90, aspect_ratio);

    
    
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    
//vec3.h 从一个单位小圆盘射出光线
vec3 random_in_unit_disk() {
    while (true) {
        auto p = vec3(random_double(-1,1), random_double(-1,1), 0);
        if (p.length_squared() >= 1) continue;
        return p;
    }
}

    
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//加入景深
class camera {
    public:
        camera(
            vec3 lookfrom, vec3 lookat, vec3 vup,
            double vfov, // top to bottom, in degrees
            double aspect, double aperture, double focus_dist
        ) {
            origin = lookfrom;
            lens_radius = aperture / 2;

            auto theta = degrees_to_radians(vfov);
            //  half_height = tan(theta/2)*dis，这里将dis放到了后面，也就是*focus_dist
            auto half_height = tan(theta/2);
            auto half_width = aspect * half_height;

            w = unit_vector(lookfrom - lookat);
            u = unit_vector(cross(vup, w));
            v = cross(w, u);
            lower_left_corner = origin
                              - half_width * focus_dist * u
                              - half_height * focus_dist * v
                              - focus_dist * w;

            horizontal = 2*half_width*focus_dist*u;
            vertical = 2*half_height*focus_dist*v;
        }

        ray get_ray(double s, double t) {
            // 半径为lens_radius的随机圆
            vec3 rd = lens_radius * random_in_unit_disk();
            //取得点在一个二维平面内，所以要用相机坐标系修正，令这些点位于origin所处的uv平面内。之前不考虑景深时origin只是一个点，自然不需要用相机坐标系修正
            vec3 offset = u * rd.x() + v * rd.y();

            return ray(
                //在原本只是一个点的origin附近一个圆内随机取点
                origin + offset,
                lower_left_corner + s*horizontal + t*vertical - origin - offset
           );
        }

    public:
        vec3 origin;
        vec3 lower_left_corner;
        vec3 horizontal;
        vec3 vertical;
        vec3 u, v, w;
        double lens_radius;
};    
 

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    
hittable_list random_scene() {
    hittable_list world;

    world.add(make_shared<sphere>(
        vec3(0,-1000,0), 1000, make_shared<lambertian>(vec3(0.5, 0.5, 0.5))));

    int i = 1;
    for (int a = -11; a < 11; a++) {
        for (int b = -11; b < 11; b++) {
            auto choose_mat = random_double();
            vec3 center(a + 0.9*random_double(), 0.2, b + 0.9*random_double());
            if ((center - vec3(4, 0.2, 0)).length() > 0.9) {
                if (choose_mat < 0.8) {
                    // diffuse
                    auto albedo = vec3::random() * vec3::random();
                    world.add(
                        make_shared<sphere>(center, 0.2, make_shared<lambertian>(albedo)));
                } else if (choose_mat < 0.95) {
                    // metal
                    auto albedo = vec3::random(.5, 1);
                    auto fuzz = random_double(0, .5);
                    world.add(
                        make_shared<sphere>(center, 0.2, make_shared<metal>(albedo, fuzz)));
                } else {
                    // glass
                    world.add(make_shared<sphere>(center, 0.2, make_shared<dielectric>(1.5)));
                }
            }
        }
    }

    world.add(make_shared<sphere>(vec3(0, 1, 0), 1.0, make_shared<dielectric>(1.5)));

    world.add(
        make_shared<sphere>(vec3(-4, 1, 0), 1.0, make_shared<lambertian>(vec3(0.4, 0.2, 0.1))));

    world.add(
        make_shared<sphere>(vec3(4, 1, 0), 1.0, make_shared<metal>(vec3(0.7, 0.6, 0.5), 0.0)));

    return world;
}
    
int main() {
    ...
    auto world = random_scene();

    vec3 lookfrom(13,2,3);
    vec3 lookat(0,0,0);
    vec3 vup(0,1,0);
    auto dist_to_focus = 10.0;
    auto aperture = 0.1;

    camera cam(lookfrom, lookat, vup, 20, aspect_ratio, aperture, dist_to_focus);
    ...
}

    
    
//运动模糊：当计算第n个像素的颜色时，比如左上角第一个像素，我们像之前抗锯齿做法一样发射100根光线，只是这次这些光线被赋予了
//随机的时间属性t。并且计算hit函数时，我们根据此时光线的t，更新模型的中心位置。就像马路上的一块地砖，从早到晚被不同的人经过
//我们选择的像素在不同时刻有不同的颜色，最后给他取个平均值就能产生模糊的效果
//简单说就是之前的步骤不变、计算交点的算法不变、光线的随机采样不变，只在计算交点时改变模型的位置
    
    
//首先给光线加上时间属性，这一属性将在main函数的cam.get_ray()中赋值，是随机数，随机的范围在cam的初始化时赋予（t1，t2）
class ray {
    public:
        ray() {}
        ray(const vec3& origin, const vec3& direction, double time = 0.0)
            : orig(origin), dir(direction), tm(time)
        {}

        vec3 origin() const    { return orig; }
        vec3 direction() const { return dir; }
        double time() const    { return tm; }

        vec3 at(double t) const {
            return orig + t*dir;
        }

    public:
        vec3 orig;
        vec3 dir;
        double tm;
}; 
    

//给相机添加快门时间t1,t2    
class camera {
    public:
        camera(
            vec3 lookfrom, vec3 lookat, vec3 vup,
            double vfov, // top to bottom, in degrees
            double aspect, double aperture, double focus_dist, double t0 = 0, double t1 = 0
        ) {
            origin = lookfrom;
            lens_radius = aperture / 2;
            time0 = t0;
            time1 = t1;
            auto theta = degrees_to_radians(vfov);
            auto half_height = tan(theta/2);
            auto half_width = aspect * half_height;

            w = unit_vector(lookfrom - lookat);
            u = unit_vector(cross(vup, w));
            v = cross(w, u);

            lower_left_corner = origin
                              - half_width*focus_dist*u
                              - half_height*focus_dist*v
                              - focus_dist*w;

            horizontal = 2*half_width*focus_dist*u;
            vertical = 2*half_height*focus_dist*v;
        }

        ray get_ray(double s, double t) {
            vec3 rd = lens_radius * random_in_unit_disk();
            vec3 offset = u * rd.x() + v * rd.y();
            return ray(
                origin + offset,
                lower_left_corner + s*horizontal + t*vertical - origin - offset,
                random_double(time0, time1)
            );
        }

    public:
        vec3 origin;
        vec3 lower_left_corner;
        vec3 horizontal;
        vec3 vertical;
        vec3 u, v, w;
        double lens_radius;
        double time0, time1;  // shutter open/close times
};
    
    
//    
class moving_sphere : public hittable {
    public:
        moving_sphere() {}
        moving_sphere(
            vec3 cen0, vec3 cen1, double t0, double t1, double r, shared_ptr<material> m)
            : center0(cen0), center1(cen1), time0(t0), time1(t1), radius(r), mat_ptr(m)
        {};

        virtual bool hit(const ray& r, double tmin, double tmax, hit_record& rec) const;

        vec3 center(double time) const;

    public:
        vec3 center0, center1;
        double time0, time1;
        double radius;
        shared_ptr<material> mat_ptr;
};

vec3 moving_sphere::center(double time) const{
    return center0 + ((time - time0) / (time1 - time0))*(center1 - center0);
}

bool moving_sphere::hit(
    const ray& r, double t_min, double t_max, hit_record& rec) const {
    vec3 oc = r.origin() - center(r.time());
    auto a = r.direction().length_squared();
    auto half_b = dot(oc, r.direction());
    auto c = oc.length_squared() - radius*radius;

    auto discriminant = half_b*half_b - a*c;

    if (discriminant > 0) {
        auto root = sqrt(discriminant);

        auto temp = (-half_b - root)/a;
        if (temp < t_max && temp > t_min) {
            rec.t = temp;
            rec.p = r.at(rec.t);
            vec3 outward_normal = (rec.p - center(r.time())) / radius;
            rec.set_face_normal(r, outward_normal);
            rec.mat_ptr = mat_ptr;
            return true;
        }

        temp = (-half_b + root) / a;
        if (temp < t_max && temp > t_min) {
            rec.t = temp;
            rec.p = r.at(rec.t);
            vec3 outward_normal = (rec.p - center(r.time())) / radius;
            rec.set_face_normal(r, outward_normal);
            rec.mat_ptr = mat_ptr;
            return true;
        }
    }
    return false;
}

    
    
//材质
class lambertian : public material {
    public:
        lambertian(const vec3& a) : albedo(a) {}

        virtual bool scatter(
            const ray& r_in, const hit_record& rec, vec3& attenuation, ray& scattered
        ) const {
            vec3 scatter_direction = rec.normal + random_unit_vector();
            scattered = ray(rec.p, scatter_direction, r_in.time());
            attenuation = albedo;
            return true;
        }

        vec3 albedo;
};
    
    
    
hittable_list random_scene() {
    hittable_list world;

    world.add(make_shared<sphere>(
        vec3(0,-1000,0), 1000, make_shared<lambertian>(vec3(0.5, 0.5, 0.5))));

    int i = 1;
    for (int a = -10; a < 10; a++) {
        for (int b = -10; b < 10; b++) {
            auto choose_mat = random_double();
            vec3 center(a + 0.9*random_double(), 0.2, b + 0.9*random_double());
            if ((center - vec3(4, .2, 0)).length() > 0.9) {
                if (choose_mat < 0.8) {
                    // diffuse
                    auto albedo = vec3::random() * vec3::random();
                    world.add(make_shared<moving_sphere>(
                        center, center + vec3(0, random_double(0,.5), 0), 0.0, 1.0, 0.2,
                        make_shared<lambertian>(albedo)));
                } else if (choose_mat < 0.95) {
                    // metal
                    auto albedo = vec3::random(.5, 1);
                    auto fuzz = random_double(0, .5);
                    world.add(
                        make_shared<sphere>(center, 0.2, make_shared<metal>(albedo, fuzz)));
                } else {
                    // glass
                    world.add(make_shared<sphere>(center, 0.2, make_shared<dielectric>(1.5)));
                }
            }
        }
    }

    world.add(make_shared<sphere>(vec3(0, 1, 0), 1.0, make_shared<dielectric>(1.5)));
    world.add(make_shared<sphere>(
        vec3(-4, 1, 0), 1.0, make_shared<lambertian>(vec3(0.4, 0.2, 0.1))));
    world.add(make_shared<sphere>(
        vec3(4, 1, 0), 1.0, make_shared<metal>(vec3(0.7, 0.6, 0.5), 0.0)));

    return world;
}

    
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////  

//怎么分堆。还得想想怎么去检测光线和包围盒相交
//从现在开始, 我们会把轴对齐的包围盒叫成矩形平行管道：aabb
//和击中那些会在屏幕上显示出来的物体时不同, 射线与AABB求交并不需要去获取那些法向啊交点啊这些东西, AABB不需要在屏幕上渲染出来。所以aabb是单独的类而不是hittable的子类

class aabb {
    public:
        aabb() {}
        aabb(const vec3& a, const vec3& b) { _min = a; _max = b;}

        vec3 min() const {return _min; }
        vec3 max() const {return _max; }

        bool hit(const ray& r, double tmin, double tmax) const {
            for (int a = 0; a < 3; a++) {
                auto t0 = ffmin((_min[a] - r.origin()[a]) / r.direction()[a],
                                (_max[a] - r.origin()[a]) / r.direction()[a]);
                auto t1 = ffmax((_min[a] - r.origin()[a]) / r.direction()[a],
                                (_max[a] - r.origin()[a]) / r.direction()[a]);
                tmin = ffmax(t0, tmin);
                tmax = ffmin(t1, tmax);
                if (tmax <= tmin)
                    return false;
            }
            return true;
        }

        vec3 _min;
        vec3 _max;
};
    
    
//Andrew Kensler's hit method
//可以看到在上面的基础上略去了一些重复计算, 优化了不少
inline bool aabb::hit(const ray& r, double tmin, double tmax) const {
    for (int a = 0; a < 3; a++) {
        auto invD = 1.0f / r.direction()[a];
        auto t0 = (min()[a] - r.origin()[a]) * invD;
        auto t1 = (max()[a] - r.origin()[a]) * invD;
        if (invD < 0.0f)
            std::swap(t0, t1);
        tmin = t0 > tmin ? t0 : tmin;
        tmax = t1 < tmax ? t1 : tmax;
        if (tmax <= tmin)
            return false;
    }
    return true;
}
    
    
    
class hittable {
    public:
        virtual bool hit(
            const ray& r, double t_min, double t_max, hit_record& rec) const = 0;
        virtual bool bounding_box(double t0, double t1, aabb& output_box) const = 0;
};

    
//求一模型的包围盒，与求一模型的交点属于同一类别
bool sphere::bounding_box(double t0, double t1, aabb& output_box) const {
    output_box = aabb(
        center - vec3(radius, radius, radius),
        center + vec3(radius, radius, radius));
    return true;
}
    
//为什么hittable_list里还有bounding_box函数？参考hittable_list里的hit，在对单个物体求交时要先在hittable_list里遍历所有物体    
bool hittable_list::bounding_box(double t0, double t1, aabb& output_box) const {
    if (objects.empty()) return false;

    aabb temp_box;
    bool first_box = true;

    for (const auto& object : objects) {
        //目前还没有哪个物体的bounding_box函数会返回false
        if (!object->bounding_box(t0, t1, temp_box)) return false;
        //
        output_box = first_box ? temp_box : surrounding_box(output_box, temp_box);
        first_box = false;
    }

    return true;
}
    
//求两个包围盒的包围盒
aabb surrounding_box(aabb box0, aabb box1) {
    //small是一个包含三个数的数组
    vec3 small(ffmin(box0.min().x(), box1.min().x()),
               ffmin(box0.min().y(), box1.min().y()),
               ffmin(box0.min().z(), box1.min().z()));
    //big同理
    vec3 big  (ffmax(box0.max().x(), box1.max().x()),
               ffmax(box0.max().y(), box1.max().y()),
               ffmax(box0.max().z(), box1.max().z()));
    //将这两个数组带入aabb函数返回一个新的包围盒
    return aabb(small,big);
}
    
    
//bvh_node和sphere一样是hittable的子类，因此一样由hit和bounding_box
//bvh.h
class bvh_node : public hittable {
    public:
        bvh_node();

        bvh_node(hittable_list& list, double time0, double time1)
            : bvh_node(list.objects, 0, list.objects.size(), time0, time1)
        {}

        bvh_node(
            std::vector<shared_ptr<hittable>>& objects,
            size_t start, size_t end, double time0, double time1);

        virtual bool hit(const ray& r, double tmin, double tmax, hit_record& rec) const;
        virtual bool bounding_box(double t0, double t1, aabb& output_box) const;

    public:
        shared_ptr<hittable> left;
        shared_ptr<hittable> right;
        aabb box;
};

bool bvh_node::bounding_box(double t0, double t1, aabb& output_box) const {
    output_box = box;
    return true;
}
    
    
bool bvh_node::hit(const ray& r, double t_min, double t_max, hit_record& rec) const {
    if (!box.hit(r, t_min, t_max))
        return false;

    bool hit_left = left->hit(r, t_min, t_max, rec);
    bool hit_right = right->hit(r, t_min, hit_left ? rec.t : t_max, rec);

    return hit_left || hit_right;
}
    
    
//看不懂    
#include <algorithm>
...

bvh_node::bvh_node(
    std::vector<shared_ptr<hittable>>& objects,
    size_t start, size_t end, double time0, double time1
) {
    int axis = random_int(0,2);
    //嵌套条件运算符：先判断axis是否为真，如果是则返回box_x_compare，如果不是则判断axis是否等于1，如果是则返回box_y_compare，如果不是则返回box_z_compare
    auto comparator = (axis == 0) ? box_x_compare
                    : (axis == 1) ? box_y_compare
                                  : box_z_compare;

    size_t object_span = end - start;

    if (object_span == 1) {
        left = right = objects[start];
    } else if (object_span == 2) {
        if (comparator(objects[start], objects[start+1])) {
            left = objects[start];
            right = objects[start+1];
        } else {
            left = objects[start+1];
            right = objects[start];
        }
    } else {
        std::sort(objects.begin() + start, objects.begin() + end, comparator);

        auto mid = start + object_span/2;
        left = make_shared<bvh_node>(objects, start, mid, time0, time1);
        right = make_shared<bvh_node>(objects, mid, end, time0, time1);
    }

    aabb box_left, box_right;

    if (  !left->bounding_box (time0, time1, box_left)
       || !right->bounding_box(time0, time1, box_right)
    )
        std::cerr << "No bounding box in bvh_node constructor.\n";

    box = surrounding_box(box_left, box_right);
}
    
    
    
inline bool box_compare(const shared_ptr<hittable> a, const shared_ptr<hittable> b, int axis) {
    aabb box_a;
    aabb box_b;

    if (!a->bounding_box(0,0, box_a) || !b->bounding_box(0,0, box_b))
        std::cerr << "No bounding box in bvh_node constructor.\n";

    return box_a.min().e[axis] < box_b.min().e[axis];
}


bool box_x_compare (const shared_ptr<hittable> a, const shared_ptr<hittable> b) {
    return box_compare(a, b, 0);
}

bool box_y_compare (const shared_ptr<hittable> a, const shared_ptr<hittable> b) {
    return box_compare(a, b, 1);
}

bool box_z_compare (const shared_ptr<hittable> a, const shared_ptr<hittable> b) {
    return box_compare(a, b, 2);
}

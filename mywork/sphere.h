#ifndef SPHERE_H
#define SPHERE_H

#include "hittable.h"
#include "vec3.h"






class sphere: public hittable {
    public:
        sphere() {}
        sphere(vec3 cen, double r, shared_ptr<material> m) : center(cen), radius(r), mat_ptr(m) {};

        virtual bool hit(const ray& r, double tmin, double tmax, hit_record& rec) const;
        virtual bool bounding_box(double t0, double t1, aabb& output_box) const;
        void get_sphere_uv(const vec3& p, double& u, double& v) const{
            auto phi = atan2(p.z(), p.x());
            auto theta = asin(p.y());
            u = 1-(phi + pi) / (2*pi);
            v = (theta + pi/2) / pi;
        }
        float getArea() const{
            float area = 4 * M_PI * radius;
            return area;
        }
        void Sample(hit_record &pos, float &pdf){
            float theta = 2.0 * M_PI * random_double(), phi = M_PI * random_double();
            vec3 dir(std::cos(phi), std::sin(phi)*std::cos(theta), std::sin(phi)*std::sin(theta));
            pos.p = center + radius * dir;
            pos.normal = dir;
            pos.emit = mat_ptr->emitted(pos.u, pos.v, pos.p);
            
            pdf = 1.0f / area;
        }
    public:
        vec3 center;
        double radius;
        shared_ptr<material> mat_ptr;
        float area;
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
            rec.mat_ptr = mat_ptr;
            get_sphere_uv((rec.p-center)/radius,rec.u,rec.v);
            return true;
        }
        temp = (-half_b + root) / a;
        if (temp < t_max && temp > t_min) {
            rec.t = temp;
            rec.p = r.at(rec.t);
            vec3 outward_normal = (rec.p - center) / radius;
            rec.set_face_normal(r, outward_normal);
            rec.mat_ptr = mat_ptr;
            return true;
        }
    }
    return false;
}


class moving_sphere : public hittable {
    public:
        moving_sphere() {}
        moving_sphere(
            vec3 cen0, vec3 cen1, double t0, double t1, double r, shared_ptr<material> m)
            : center0(cen0), center1(cen1), time0(t0), time1(t1), radius(r), mat_ptr(m)
        {};

        virtual bool hit(const ray& r, double tmin, double tmax, hit_record& rec) const;
        virtual bool bounding_box(double t0, double t1, aabb& output_box) const;
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




bool sphere::bounding_box(double t0, double t1, aabb& output_box) const {
    output_box = aabb(
        center - vec3(radius, radius, radius),
        center + vec3(radius, radius, radius));
    return true;
}


bool moving_sphere::bounding_box(double t0, double t1, aabb& output_box) const {
    aabb box0(
        center(t0) - vec3(radius, radius, radius),
        center(t0) + vec3(radius, radius, radius));
    aabb box1(
        center(t1) - vec3(radius, radius, radius),
        center(t1) + vec3(radius, radius, radius));
    output_box = surrounding_box(box0, box1);
    return true;
}

class xy_rect: public hittable {
        public:
            xy_rect() {}

            xy_rect(double _x0, double _x1, double _y0, double _y1, double _k, shared_ptr<material> mat)
                : x0(_x0), x1(_x1), y0(_y0), y1(_y1), k(_k), mat_ptr(mat) {};

            virtual bool hit(const ray& r, double t0, double t1, hit_record& rec) const;

            virtual bool bounding_box(double t0, double t1, aabb& output_box) const {
                output_box =  aabb(vec3(x0,y0, k-0.0001), vec3(x1, y1, k+0.0001));
                return true;
            }

        public:
            shared_ptr<material> mat_ptr;
            double x0, x1, y0, y1, k;
};

bool xy_rect::hit(const ray& r, double t0, double t1, hit_record& rec) const {
        auto t = (k-r.origin().z()) / r.direction().z();
        if (t < t0 || t > t1)
            return false;
        auto x = r.origin().x() + t*r.direction().x();
        auto y = r.origin().y() + t*r.direction().y();
        if (x < x0 || x > x1 || y < y0 || y > y1)
            return false;
        rec.u = (x-x0)/(x1-x0);
        rec.v = (y-y0)/(y1-y0);
        rec.t = t;
        vec3 outward_normal = vec3(0, 0, 1);
        rec.set_face_normal(r, outward_normal);
        rec.mat_ptr = mat_ptr;
        rec.p = r.at(t);
        return true;
    }

class xz_rect: public hittable {
        public:
            xz_rect() {}

            xz_rect(double _x0, double _x1, double _z0, double _z1, double _k, shared_ptr<material> mat)
                : x0(_x0), x1(_x1), z0(_z0), z1(_z1), k(_k), mat_ptr(mat) {};

            virtual bool hit(const ray& r, double t0, double t1, hit_record& rec) const;

            virtual bool bounding_box(double t0, double t1, aabb& output_box) const {
                output_box =  aabb(vec3(x0,k-0.0001,z0), vec3(x1, k+0.0001, z1));
                return true;
            }

        public:
            shared_ptr<material> mat_ptr;
            double x0, x1, z0, z1, k;
    };

class yz_rect: public hittable {
    public:
        yz_rect() {}

        yz_rect(double _y0, double _y1, double _z0, double _z1, double _k, shared_ptr<material> mat): y0(_y0), y1(_y1), z0(_z0), z1(_z1), k(_k), mat_ptr(mat) 
        {};

        virtual bool hit(const ray& r, double t0, double t1, hit_record& rec) const;

        virtual bool bounding_box(double t0, double t1, aabb& output_box) const {
            output_box =  aabb(vec3(k-0.0001, y0, z0), vec3(k+0.0001, y1, z1));
            return true;
        }

    public:
        shared_ptr<material> mat_ptr;
        double y0, y1, z0, z1, k;
};


bool xz_rect::hit(const ray& r, double t0, double t1, hit_record& rec) const {
        auto t = (k-r.origin().y()) / r.direction().y();
        if (t < t0 || t > t1)
            return false;
        auto x = r.origin().x() + t*r.direction().x();
        auto z = r.origin().z() + t*r.direction().z();
        if (x < x0 || x > x1 || z < z0 || z > z1)
            return false;
        rec.u = (x-x0)/(x1-x0);
        rec.v = (z-z0)/(z1-z0);
        rec.t = t;
        vec3 outward_normal = vec3(0, 1, 0);
        rec.set_face_normal(r, outward_normal);
        rec.mat_ptr = mat_ptr;
        rec.p = r.at(t);
        return true;
    }

bool yz_rect::hit(const ray& r, double t0, double t1, hit_record& rec) const {
    auto t = (k-r.origin().x()) / r.direction().x();
    if (t < t0 || t > t1)
        return false;
    auto y = r.origin().y() + t*r.direction().y();
    auto z = r.origin().z() + t*r.direction().z();
    if (y < y0 || y > y1 || z < z0 || z > z1)
        return false;
    rec.u = (y-y0)/(y1-y0);
    rec.v = (z-z0)/(z1-z0);
    rec.t = t;
    vec3 outward_normal = vec3(1, 0, 0);
    rec.set_face_normal(r, outward_normal);
    rec.mat_ptr = mat_ptr;
    rec.p = r.at(t);
    return true;
}

 class flip_face : public hittable {
        public:
            flip_face(shared_ptr<hittable> p) : ptr(p) {}

            virtual bool hit(const ray& r, double t_min, double t_max, hit_record& rec) const {
                if (!ptr->hit(r, t_min, t_max, rec))
                    return false;

                rec.front_face = !rec.front_face;
                return true;
            }

            virtual bool bounding_box(double t0, double t1, aabb& output_box) const {
                return ptr->bounding_box(t0, t1, output_box);
            }

        public:
            shared_ptr<hittable> ptr;
    };

    
class box: public hittable  {
    public:
        box() {}
        box(const vec3& p0, const vec3& p1, shared_ptr<material> ptr);

        virtual bool hit(const ray& r, double t0, double t1, hit_record& rec) const;

        virtual bool bounding_box(double t0, double t1, aabb& output_box) const {
            output_box = aabb(box_min, box_max);
            return true;
        }

    public:
        vec3 box_min;
        vec3 box_max;
        hittable_list sides;
};

box::box(const vec3& p0, const vec3& p1, shared_ptr<material> ptr) {
    box_min = p0;
    box_max = p1;

    sides.add(make_shared<xy_rect>(p0.x(), p1.x(), p0.y(), p1.y(), p1.z(), ptr));
    sides.add(make_shared<flip_face>(
        make_shared<xy_rect>(p0.x(), p1.x(), p0.y(), p1.y(), p0.z(), ptr)));

    sides.add(make_shared<xz_rect>(p0.x(), p1.x(), p0.z(), p1.z(), p1.y(), ptr));
    sides.add(make_shared<flip_face>(
        make_shared<xz_rect>(p0.x(), p1.x(), p0.z(), p1.z(), p0.y(), ptr)));

    sides.add(make_shared<yz_rect>(p0.y(), p1.y(), p0.z(), p1.z(), p1.x(), ptr));
    sides.add(make_shared<flip_face>(
        make_shared<yz_rect>(p0.y(), p1.y(), p0.z(), p1.z(), p0.x(), ptr)));
}

bool box::hit(const ray& r, double t0, double t1, hit_record& rec) const {
    return sides.hit(r, t0, t1, rec);
}

class translate : public hittable {
    public:
        translate(shared_ptr<hittable> p, const vec3& displacement)
            : ptr(p), offset(displacement) {}

        virtual bool hit(const ray& r, double t_min, double t_max, hit_record& rec) const;
        virtual bool bounding_box(double t0, double t1, aabb& output_box) const;

    public:
        shared_ptr<hittable> ptr;
        vec3 offset;
};

bool translate::hit(const ray& r, double t_min, double t_max, hit_record& rec) const {
    ray moved_r(r.origin() - offset, r.direction(), r.time());
    if (!ptr->hit(moved_r, t_min, t_max, rec))
        return false;

    rec.p += offset;
    rec.set_face_normal(moved_r, rec.normal);

    return true;
}

bool translate::bounding_box(double t0, double t1, aabb& output_box) const {
    if (!ptr->bounding_box(t0, t1, output_box))
        return false;

    output_box = aabb(
        output_box.min() + offset,
        output_box.max() + offset);

    return true;
}

class rotate_y : public hittable {
    public:
        rotate_y(shared_ptr<hittable> p, double angle);

        virtual bool hit(const ray& r, double t_min, double t_max, hit_record& rec) const;
        virtual bool bounding_box(double t0, double t1, aabb& output_box) const {
            output_box = bbox;
            return hasbox;
        }

    public:
        shared_ptr<hittable> ptr;
        double sin_theta;
        double cos_theta;
        bool hasbox;
        aabb bbox;
};


rotate_y::rotate_y(shared_ptr<hittable> p, double angle) : ptr(p) {
    auto radians = degrees_to_radians(angle);
    sin_theta = sin(radians);
    cos_theta = cos(radians);
    hasbox = ptr->bounding_box(0, 1, bbox);

    vec3 min( infinity,  infinity,  infinity);
    vec3 max(-infinity, -infinity, -infinity);

    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            for (int k = 0; k < 2; k++) {
                auto x = i*bbox.max().x() + (1-i)*bbox.min().x();
                auto y = j*bbox.max().y() + (1-j)*bbox.min().y();
                auto z = k*bbox.max().z() + (1-k)*bbox.min().z();

                auto newx =  cos_theta*x + sin_theta*z;
                auto newz = -sin_theta*x + cos_theta*z;

                vec3 tester(newx, y, newz);

                for (int c = 0; c < 3; c++) {
                    min[c] = ffmin(min[c], tester[c]);
                    max[c] = ffmax(max[c], tester[c]);
                }
            }
        }
    }

    bbox = aabb(min, max);
}

bool rotate_y::hit(const ray& r, double t_min, double t_max, hit_record& rec) const {
    vec3 origin = r.origin();
    vec3 direction = r.direction();

    origin[0] = cos_theta*r.origin()[0] - sin_theta*r.origin()[2];
    origin[2] = sin_theta*r.origin()[0] + cos_theta*r.origin()[2];

    direction[0] = cos_theta*r.direction()[0] - sin_theta*r.direction()[2];
    direction[2] = sin_theta*r.direction()[0] + cos_theta*r.direction()[2];

    ray rotated_r(origin, direction, r.time());

    if (!ptr->hit(rotated_r, t_min, t_max, rec))
        return false;

    vec3 p = rec.p;
    vec3 normal = rec.normal;

    p[0] =  cos_theta*rec.p[0] + sin_theta*rec.p[2];
    p[2] = -sin_theta*rec.p[0] + cos_theta*rec.p[2];

    normal[0] =  cos_theta*rec.normal[0] + sin_theta*rec.normal[2];
    normal[2] = -sin_theta*rec.normal[0] + cos_theta*rec.normal[2];

    rec.p = p;
    rec.set_face_normal(rotated_r, normal);

    return true;
}

#endif
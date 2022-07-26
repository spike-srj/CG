#ifndef SPHERE_H
#define SPHERE_H

#include "hittable.h"
#include "vec3.h"
#include "material.h"
#include "hittable_list.h"




class sphere: public hittable {
    public:
        sphere() {}
        sphere(vec3 cen, double r, shared_ptr<material> m) : center(cen), radius(r) {mat_ptr=m;};

        virtual bool hit(const ray& r, double tmin, double tmax, hit_record& rec) const;
        virtual bool bounding_box(double t0, double t1, aabb& output_box) const;
        virtual void get_sphere_uv(const vec3& p, double& u, double& v) const;
        virtual double getArea() const;
        virtual void Sample(hit_record &pos, float &pdf);
    public:
        vec3 center;
        double radius;
        //shared_ptr<material> mat_ptr;
        float area;
};




class moving_sphere : public hittable {
    public:
        moving_sphere() {}
        moving_sphere(
            vec3 cen0, vec3 cen1, double t0, double t1, double r, shared_ptr<material> m)
            : center0(cen0), center1(cen1), time0(t0), time1(t1), radius(r)
        { mat_ptr=m;};

        virtual bool hit(const ray& r, double tmin, double tmax, hit_record& rec) const;
        virtual bool bounding_box(double t0, double t1, aabb& output_box) const;
        vec3 center(double time) const;

    public:
        vec3 center0, center1;
        double time0, time1;
        double radius;
        //shared_ptr<material> mat_ptr;
};



class xy_rect: public hittable {
        public:
            xy_rect() {}

            xy_rect(double _x0, double _x1, double _y0, double _y1, double _k, shared_ptr<material> mat)
                : x0(_x0), x1(_x1), y0(_y0), y1(_y1), k(_k), mp(mat) {};

            virtual bool hit(const ray& r, double t0, double t1, hit_record& rec) const;

            virtual bool bounding_box(double t0, double t1, aabb& output_box) const {
                output_box =  aabb(vec3(x0,y0, k-0.0001), vec3(x1, y1, k+0.0001));
                return true;
            }
            virtual double getArea() const;
            virtual void Sample(hit_record &pos, float &pdf);
            
             

        public:
            shared_ptr<material> mp;
            double x0, x1, y0, y1, k;
            double area;
};



class xz_rect: public hittable {
        public:
            xz_rect() {}

            xz_rect(double _x0, double _x1, double _z0, double _z1, double _k, shared_ptr<material> mat)
                : x0(_x0), x1(_x1), z0(_z0), z1(_z1), k(_k), mp(mat) {};

            virtual bool hit(const ray& r, double t0, double t1, hit_record& rec) const;

            virtual bool bounding_box(double t0, double t1, aabb& output_box) const {
                output_box =  aabb(vec3(x0,k-0.0001,z0), vec3(x1, k+0.0001, z1));
                return true;
            }
            virtual double getArea() const;
            virtual void Sample(hit_record &pos, float &pdf);

        public:
            shared_ptr<material> mp;
            double x0, x1, z0, z1, k;
            double area;
    };

class yz_rect: public hittable {
    public:
        yz_rect() {}

        yz_rect(double _y0, double _y1, double _z0, double _z1, double _k, shared_ptr<material> mat): y0(_y0), y1(_y1), z0(_z0), z1(_z1), k(_k), mp(mat) 
        {};

        virtual bool hit(const ray& r, double t0, double t1, hit_record& rec) const;

        virtual bool bounding_box(double t0, double t1, aabb& output_box) const {
            output_box =  aabb(vec3(k-0.0001, y0, z0), vec3(k+0.0001, y1, z1));
            return true;
        }
        virtual double getArea() const;
        virtual void Sample(hit_record &pos, float &pdf);

    public:
        shared_ptr<material> mp;
        double y0, y1, z0, z1, k;
        double area;
};




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


#endif
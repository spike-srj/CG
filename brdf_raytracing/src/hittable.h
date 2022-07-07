#ifndef HITTABLE_H
#define HITTABLE_H

#include "ray.h"
#include "material.h"
#include "hitrecord.h"
//class material;   //??????????????????

class aabb {
    public:
        aabb() {}
        aabb(const vec3& a, const vec3& b) { _min = a; _max = b;}

        vec3 min() const {return _min; }
        vec3 max() const {return _max; }
        inline bool hit(const ray& r, double tmin, double tmax) const {
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
        


        vec3 _min;
        vec3 _max;
};

extern aabb surrounding_box(aabb box0, aabb box1);




class hittable {
    public:
        virtual bool hit(
            const ray& r, double t_min, double t_max, hit_record& rec) const = 0;
        virtual bool bounding_box(double t0, double t1, aabb& output_box) const = 0;
        virtual double getArea(){}
        virtual void Sample(hit_record &pos, float &pdf){}
        hittable(){mat_ptr=NULL;}
    public:
        shared_ptr<material> mat_ptr;
};





#endif
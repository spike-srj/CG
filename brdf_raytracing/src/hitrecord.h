#ifndef _HIT_RECORD_H_
#define _HIT_RECORD_H_


#include "vec3.h"
#include "ray.h"
class material;


struct hit_record {
    vec3 p;
    vec3 normal;
    double t;
    bool front_face;
    double u;
    double v;
    shared_ptr<material> mat_ptr;
    vec3 emit;

    hit_record()
    {
        mat_ptr=NULL;
    }

    inline void set_face_normal(const ray& r, const vec3& outward_normal) {
        front_face = dot(r.direction(), outward_normal) < 0;
        normal = front_face ? outward_normal :-outward_normal;
    }
};

#endif
#ifndef HITTABLE_LIST_H
#define HITTABLE_LIST_H

#include "hittable.h"
#include <memory>
#include <vector>
#include "rtweekend.h"
//#include "sphere.h"


using std::shared_ptr;
using std::make_shared;




class hittable_list: public hittable {
    public:
        hittable_list() {}
        hittable_list(shared_ptr<hittable> object) { add(object); }

        void clear() { objects.clear(); }
        void add(shared_ptr<hittable> object) { 
            objects.push_back(object);
            size++; }

        virtual bool hit(const ray& r, double tmin, double tmax, hit_record& rec) const;
        virtual bool bounding_box(double t0, double t1, aabb& output_box) const;
        virtual void sampleLight(hit_record &pos, float &pdf);
    public:
        std::vector<shared_ptr<hittable>> objects;
        int size;
        
};



#endif

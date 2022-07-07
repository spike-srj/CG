#include "hittable.h"



aabb surrounding_box(aabb box0, aabb box1) {
    vec3 small(ffmin(box0.min().x(), box1.min().x()),
               ffmin(box0.min().y(), box1.min().y()),
               ffmin(box0.min().z(), box1.min().z()));
    vec3 big  (ffmax(box0.max().x(), box1.max().x()),
               ffmax(box0.max().y(), box1.max().y()),
               ffmax(box0.max().z(), box1.max().z()));
    return aabb(small,big);
}

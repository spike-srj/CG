
#include "material.h"
#include "hittable.h"

bool lambertian::scatter(
            const ray& r_in, const hit_record& rec, vec3& attenuation, ray& scattered
        ) const 
        {
            vec3 scatter_direction = rec.normal + random_unit_vector();
            scattered = ray(rec.p, scatter_direction, r_in.time());
            attenuation = albedo->value(rec.u, rec.v, rec.p);
            return true;
        }
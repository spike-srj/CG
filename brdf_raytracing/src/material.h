#ifndef _MATERIAL_INCLUDE_
#define _MATERIAL_INCLUDE_

#include "rtweekend.h"
#include "hitrecord.h"
//#include "sphere.h"
#include "texture.h"
#include <iostream>
#include "vec3.h"

static double schlick2(double cosine, double ref_idx) {
    auto r0 = (1-ref_idx) / (1+ref_idx);
    r0 = r0*r0;
    return r0 + (1-r0)*pow((1 - cosine),5);
}

class material {
    public:
        virtual vec3 emitted(double u, double v, const vec3& p) const {
            return vec3(0,0,0);
        }
        virtual bool hasemission(double u, double v, const vec3& p) const {
            return false;
        }
        virtual bool scatter(
            const ray& r_in, const hit_record& rec, vec3& attenuation, ray& scattered
        ) const = 0;
        virtual vec3 eval(const vec3 &wi, const vec3 &wo, const vec3 &N) const = 0;
        virtual float pdf(const vec3 &wi, const vec3 &wo, const vec3 &N) const = 0;
        vec3 Kd, Ks;
        
};


class lambertian : public material {
    public:
        lambertian(shared_ptr<texture> a) : albedo(a) {}

        virtual bool scatter(
            const ray& r_in, const hit_record& rec, vec3& attenuation, ray& scattered
        ) const ;

        virtual vec3 eval(const vec3 &wi, const vec3 &wo, const vec3 &N) const{
            // calculate the contribution of diffuse   model
            float cosalpha = dot(N, wo);
            if (cosalpha > 0.0f) {
                vec3 diffuse = Kd / M_PI;
                return diffuse;
            }
            else
                return vec3(0.0f,0.0f,0.0f);
        }
        virtual float pdf(const vec3 &wi, const vec3 &wo, const vec3 &N) const{
            if (dot(wo, N) > 0.0f)
                return 0.5f / M_PI;
            else
                return 0.0f;
        }

    public:
        shared_ptr<texture> albedo;
};



class metal : public material {
    public:
        metal(const vec3& a, double f) : albedo(a), fuzz(f < 1 ? f : 1) {}

        virtual bool scatter(
            const ray& r_in, const hit_record& rec, vec3& attenuation, ray& scattered
        ) const {
            vec3 reflected = reflect(unit_vector(r_in.direction()), rec.normal);
            scattered = ray(rec.p, reflected + fuzz*random_in_unit_sphere());
            //scattered = ray(rec.p, reflected);
            attenuation = albedo;
            return (dot(scattered.direction(), rec.normal) > 0);
        }
        virtual vec3 eval(const vec3 &wi, const vec3 &wo, const vec3 &N) const {
            // calculate the contribution of diffuse   model
            float cosalpha = dot(N, wo);
            if (cosalpha > 0.0f) {
                vec3 diffuse = Kd / M_PI;
                return diffuse;
            }
            else
                return vec3(0.0f,0.0f,0.0f);
        }
        virtual float pdf(const vec3 &wi, const vec3 &wo, const vec3 &N) const {
            if (dot(wo, N) > 0.0f)
                return 0.5f / M_PI;
            else
                return 0.0f;
        }
    public:
        vec3 albedo;
        double fuzz;
};




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
            double reflect_prob = schlick2(cos_theta, etai_over_etat);
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
        virtual vec3 eval(const vec3 &wi, const vec3 &wo, const vec3 &N) const {
            // calculate the contribution of diffuse   model
            float cosalpha = dot(N, wo);
            if (cosalpha > 0.0f) {
                vec3 diffuse = Kd / M_PI;
                return diffuse;
            }
            else
                return vec3(0.0f,0.0f,0.0f);
        }
        virtual float pdf(const vec3 &wi, const vec3 &wo, const vec3 &N) const {
            if (dot(wo, N) > 0.0f)
                return 0.5f / M_PI;
            else
                return 0.0f;
        }
    public:
        double ref_idx;
};

class diffuse_light : public material  {
    public:
        diffuse_light(shared_ptr<texture> a) : emit(a) {}

        virtual bool scatter(
            const ray& r_in, const hit_record& rec, vec3& attenuation, ray& scattered
        ) const {
            return false;
        }

        virtual vec3 emitted(double u, double v, const vec3& p) const {
            return emit->value(u, v, p);
        }
        virtual bool hasemission(double u, double v, const vec3& p) const {
            return true;
        }

        virtual vec3 eval(const vec3 &wi, const vec3 &wo, const vec3 &N) const {
            float cosalpha = dot(N, wo);
            if (cosalpha > 0.0f) {
                vec3 diffuse = Kd / M_PI;
                return diffuse;
            }
            else
                return vec3(0.0f,0.0f,0.0f);
        }
        virtual float pdf(const vec3 &wi, const vec3 &wo, const vec3 &N) const {
            if (dot(wo, N) > 0.0f)
                return 0.5f / M_PI;
            else
                return 0.0f;
        }        
    public:
        shared_ptr<texture> emit;
};

#endif
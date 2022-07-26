#include "rtweekend.h"

#include "hittable_list.h"
#include "camera.h"
#include "hittable.h"
#include "hitrecord.h"
#include "material.h"
#include "bvh.h"
#include "ray.h"
#include "sphere.h"
#include <iostream>
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"


vec3 ray_color(const ray& r, const vec3& background, hittable_list& world, int depth) {
        hit_record rec;

        // If we've exceeded the ray bounce limit, no more light is gathered.
        if (depth <= 0)
            return vec3(0,0,0);

        // If the ray hits nothing, return the background color.
        //update rec
        if (!world.hit(r, 0.001, infinity, rec))
            return background;

        // 如果打到物体
		vec3 L_dir(0, 0, 0);
		vec3 L_indir(0, 0, 0);
        //均匀采样光源物体，取光源上一点
		hit_record lightInter;
		float pdf_light = 0.0f;
		world.sampleLight(lightInter, pdf_light);

		// 物体表面法线
		auto& N = rec.normal;
		// 灯光表面法线
		auto& NN = lightInter.normal;
		//物体点坐标
		auto& objPos = rec.p;
		//光源点坐标
		auto& lightPos = lightInter.p;

		auto diff = lightPos - objPos;
		auto lightDir = diff.normalized();
		float lightDistance = diff.x() * diff.x() + diff.y() * diff.y() + diff.z() * diff.z();

		//从物体打向光源的感光线(感光线为光线传播的逆方向)
		ray light(objPos, lightDir);
		//该光线与场景求交
		hit_record light2obj;

        // 如果反射击中光源
        if (world.hit(light, 0.001, infinity, light2obj) && (light2obj.p - lightPos).norm() < 1e-2)
        {
            //获取改材质的brdf，这里的brdf为漫反射（brdf=Kd/pi）
            vec3 f_r = rec.mat_ptr->eval(r.dir, lightDir, N);
            //直接光照光 = 光源光 * brdf * 光线和物体角度衰减 * 光线和光源法线角度衰减 / 光线距离 / 该点的概率密度（1/该光源的面积）
            L_dir = lightInter.emit * f_r * dot(lightDir, N) * dot(-lightDir, NN) / lightDistance / pdf_light;
        }
        float RussianRoulette = 0.8;
        //俄罗斯轮盘赌，确定是否继续弹射光线
        if (random_double() < RussianRoulette)
        {   
            //scattered = sample
            ray scattered;
            vec3 attenuation;
            vec3 emitted = rec.mat_ptr->emitted(rec.u, rec.v, rec.p);
            if (!rec.mat_ptr->scatter(r, rec, attenuation, scattered))
                return emitted;
            float pdf = rec.mat_ptr->pdf(r.dir, scattered.dir, N);
            vec3 f_r = rec.mat_ptr->eval(r.dir, scattered.dir, N);
            L_indir = ray_color(scattered, background, world, depth-1) * f_r * dot(scattered.dir, N) / pdf / RussianRoulette;
            L_indir = emitted + L_indir;
        }
        //最后返回直接光照和间接光照
        vec3 L = L_dir + L_indir;
        return L;
    }

hittable_list random_scene() {
    hittable_list world;
    auto checker = make_shared<checker_texture>(
        make_shared<constant_texture>(vec3(0.2, 0.3, 0.1)),
        make_shared<constant_texture>(vec3(0.9, 0.9, 0.9))
    );
    world.add(make_shared<sphere>(vec3(0,-1000,0), 1000, make_shared<lambertian>(checker)));

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
                        make_shared<lambertian>(make_shared<constant_texture>(albedo))));
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
        vec3(-4, 1, 0), 1.0, make_shared<lambertian>(make_shared<constant_texture>(vec3(0.5, 0.5, 0.5)))));
    world.add(make_shared<sphere>(
        vec3(4, 1, 0), 1.0, make_shared<metal>(vec3(0.7, 0.6, 0.5), 0.0)));

    return static_cast<hittable_list>(make_shared<bvh_node>(world,0,1));
}


hittable_list two_perlin_spheres() {
    hittable_list objects;

    auto pertext = make_shared<noise_texture>(5);
    objects.add(make_shared<sphere>(vec3(0,-1000, 0), 1000, make_shared<lambertian>(pertext)));
    objects.add(make_shared<sphere>(vec3(0, 2, 0), 2, make_shared<lambertian>(pertext)));

    return objects;
}

hittable_list earth() {
    int nx, ny, nn;
    unsigned char* texture_data = stbi_load("earthmap.jpg", &nx, &ny, &nn, 0);

    auto earth_surface =
        make_shared<lambertian>(make_shared<image_texture>(texture_data, nx, ny));
    auto globe = make_shared<sphere>(vec3(0,0,0), 2, earth_surface);

    return hittable_list(globe);
}

hittable_list simple_light() {
        hittable_list objects;

        auto pertext = make_shared<noise_texture>(4);
        auto obj1 =  make_shared<lambertian>(pertext);
        obj1->Kd = vec3(0.63f, 0.065f, 0.05f);
        objects.add(make_shared<sphere>(vec3(0,-1000, 0), 1000,obj1));
        auto obj2 =  make_shared<lambertian>(pertext);
        obj2->Kd = vec3(0.14f, 0.45f, 0.091f);
        objects.add(make_shared<sphere>(vec3(0,2,0), 2, obj2));

        auto difflight = make_shared<diffuse_light>(make_shared<constant_texture>(vec3(4,4,4)));
        objects.add(make_shared<sphere>(vec3(0,7,0), 2, difflight));
        objects.add(make_shared<xy_rect>(3, 5, 1, 3, -2, difflight));
            

        return objects;
    }

hittable_list cornell_box() {
        hittable_list objects;

        auto red = make_shared<lambertian>(make_shared<constant_texture>(vec3(0.65, 0.05, 0.05)));
        
        auto white = make_shared<lambertian>(make_shared<constant_texture>(vec3(0.73, 0.73, 0.73)));
        auto green = make_shared<lambertian>(make_shared<constant_texture>(vec3(0.12, 0.45, 0.15)));
        auto light = make_shared<diffuse_light>(make_shared<constant_texture>(vec3(15, 15, 15)));
        red->Kd = vec3(0.63f, 0.065f, 0.05f);
        green->Kd = vec3(0.14f, 0.45f, 0.091f);
        white->Kd = vec3(0.725f, 0.71f, 0.68f);
        light->Kd = vec3(0.65f,0.65f,0.65f);

        objects.add(make_shared<yz_rect>(0.0, 555.0, 0.0, 555.0, 555.0, green));
        objects.add(make_shared<yz_rect>(0.0, 555.0, 0.0, 555.0, 0.0, red));
        objects.add(make_shared<xz_rect>(213, 343, 227, 332, 554, light));
        objects.add(make_shared<xz_rect>(0, 555, 0, 555, 0, white));
        objects.add(make_shared<xy_rect>(0, 555, 0, 555, 555, white));
        objects.add(make_shared<xz_rect>(0, 555, 0, 555, 555, white));
        
        shared_ptr<hittable> box1 = make_shared<box>(vec3(0, 0, 0), vec3(165, 330, 165), white);
        box1 = make_shared<rotate_y>(box1,  15);
        box1 = make_shared<translate>(box1, vec3(265,0,295));
        objects.add(box1);

        shared_ptr<hittable> box2 = make_shared<box>(vec3(0,0,0), vec3(165,165,165), white);
        box2 = make_shared<rotate_y>(box2, -18);
        box2 = make_shared<translate>(box2, vec3(130,0,65));
        objects.add(box2);


        return static_cast<hittable_list>(make_shared<bvh_node>(objects,0,0));        
    }

//main.cc
int main() {
    clock_t start,end;
    start = clock();
    const vec3 background(0,0,0);
    const int image_width = 1920;
    const int image_height = 1080;
    const int samples_per_pixel = 300;
    const int max_depth = 300;
    std::cout << "P3\n" << image_width << " " << image_height << "\n255\n";

    auto R = cos(pi/4);
    auto world = simple_light();//cornell_box();///////

    //vec3 lookfrom(278, 278, -800);
    //vec3 lookat(278,278,0);
    //vec3 lookfrom(5, 5, -15);
    //vec3 lookat(5,5,0);
    vec3 lookfrom(13, 2, 3);
    vec3 lookat(0,0,0);
    vec3 vup(0,1,0);
    auto dist_to_focus = 10.0;
    auto aperture = 0.0;
    const auto aspect_ratio = double(image_width) / image_height;
    camera cam(lookfrom, lookat, vup, 40, aspect_ratio, aperture, dist_to_focus, 0.0, 1.0);  //20-40
    for (int j = image_height-1; j >= 0; --j) {
        std::cerr << "\rScanlines remaining: " << j << ' ' << std::flush;
        for (int i = 0; i < image_width; ++i) {
            vec3 color(0, 0, 0);
            for (int s = 0; s < samples_per_pixel; ++s) {
                auto u = (i + random_double()) / image_width;
                auto v = (j + random_double()) / image_height;
                ray r = cam.get_ray(u, v);
               color += ray_color(r, background, world, max_depth);
            }
            color.write_color(std::cout, samples_per_pixel);
        }
    }

    std::cerr << "\nDone.\n";
    end = clock();
    std::cout<<"time = "<<double(end-start)/CLOCKS_PER_SEC<<"s"<<std::endl;
}










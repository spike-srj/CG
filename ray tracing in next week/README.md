整理这篇Ray tracing in next week 的思路，介绍整套程序运行的步骤

![image](https://github.com/spike-srj/CG/blob/main/ray%20tracing%20in%20next%20week/v2-460d12c611ba27f071c9977864c6a2d0_1440w.jpeg)

Features

    光线、可自定义的摄像机、背景
    物体：球体、移动的球体、矩形、包围盒    
    多重采样抗锯齿    
    材质（漫反射材质、金属材质、绝缘体玻璃球（折射）、光源）
    散焦模糊
    动态模糊
    BVH包围盒加速
    纹理，贴图
    柏林噪声
    旋转
    体积体volumes
    
步骤

    main函数部分除了初始化摄像机多了两个时间参数，其余和上一篇一样
    首先给定图像尺寸.宽200、高100、每个像素采样100次、光线最多散射50次
    调用要展示的物体，调用random_scene()函数创建hittable_list类。这一类保存了物体信息。
    设置摄像机（起点、看向方向、vup、离焦点距离、aperture、宽高比、fov=20）
    开始遍历每一个像素，在每一个像素下设置初始颜色color（0，0，0）
        对每个像素采样100次，采样是指采颜色。
            采样是根据xxxxx来采的，这里的空格可以填“whitted-style光追”；“BRDF+路径追踪”；“其他光追模型”
            教材里最开始是根据在图中的不同位置用线性插值法给了渐变的颜色值（没有用ray_color函数）
            后来教材里用whitted-style光追，首先生成光线auto u = (i + random_double()) / image_width;???ray r = cam.get_ray(u, v);???这是配合多重采样，对给定像素发射多条射线多次采样求平均
            将光线、模型列表、最大散射次数传入ray_color函数进行递归计算。一个像素的值该由从这一像素所射到相机的光线颜色决定（包含本身发射的光、从其他位置散射、反射来的）
                ray_color首先定义一个相交点结构rec
                要进行两个判断：1是散射次数大于50次返回黑色
                2是判断光线是否与物体相交（调用hittable_list类的hit函数，该hit函数会遍历hittable_list里的所有物体，对每一个物体调用该物体类下的hit函数判断是否相交）
                    如果相交，///物体类下的hit函数会更新相交点rec的信息（包括：相交时间、位置、单位法线、法线朝向、材质指针，该指针由创建物体时传入（同时传入的参数有反照率、fuzz模糊程度），指向一个材质类（可以类比vec3类））
                    使用get_sphere_uv函数更新rec的uv坐标    
                        定义散射光线scattered以及衰减系数attenuation
                        接着判断是否发生散射（调用rec的材质指针所指材质类的散射函数，传入scattered和attenuation）
                        不同材质的散射函数不同，主要是更新scattered和attenuation。以lambertian材质为例：
                            先根据法线计算散射光线，再根据交点位置和散射光线更新scattered，根据金属的反照率更新attenuation
                            此时反照率是一个贴图类，在材质的scatter函数里调用贴图来更新attenuation
                            将更新后的scattered再次传入ray_color函数迭代计算，并乘以一个衰减率。
                            经过多次迭代直到scattered不再射到物体上，执行返回背景颜色的步骤                       
                            如果没有发生散射，则返回该材质的发光项                        
                     如果传入的光线（可能是经过反射散射后的）没有与物体相交，返回背景的颜色（可能会在前面乘以衰减系数）；
            将迭代后的颜色加入color
        采样100次，对color取平均，输出写入xxx文件
    
    
     
bvh包围盒，分两步

    第一步是在创建物体列表的最后调用，用以创建bvh包围盒的，return static_cast<hittable_list>(make_shared<bvh_node>(world,0,1))
    第二部，在ray_color中判断光线与物体是否相交时，此时物体被换成了包围盒，先调用bvh的hit函数，直到最后一层包围盒时调用物体的hit函数
    想要使用包围盒需在物体类中加入构建包围盒的函数，这会在第一步用上。最主要的还是创建包围盒，求交时与之前类似
             

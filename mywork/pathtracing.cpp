// Implementation of Path Tracing
//首先判断传入的光线是否与场景有交点
//如果有交点分两种情况：1、直接打到光源上；2、打到物体上。打到物体上后进行两次采样，一次采样击中光源的反射光、一次根据俄罗斯轮盘随机采样
Vector3f Scene::castRay(const Ray& ray, int depth) const
{	
	//获取场景中某个三角形和该光线的交点信息，相当于hit_rec
	Intersection inter = intersect(ray);

	if (inter.happened)//如果光线与场景有交点
	{
		//如果打到光源，此时交点为光源
		if (inter.m->hasEmission())
		{
			// 如果射线第一次打到光源，直接返回光源光
			if (depth == 0)
			{
				return inter.m->getEmission();
			}
			//若射线不是第一次打到光源，而是经过几次弹射和才打到光源，则返回0。该问题在交点为物体时求解。（此时交点为光源）
			//这里不太理解，就是说下面交点为物体时，进行两次采样，其中使用俄罗斯轮盘采样如果目标为光源，则返回0？
			else return Vector3f(0, 0, 0);
		}

		// 如果打到物体
		Vector3f L_dir(0, 0, 0);
		Vector3f L_indir(0, 0, 0);

		//均匀采样光源物体，取光源上一点
		Intersection lightInter;
		float pdf_light = 0.0f;
		//sample是随机采样，sampleLight是在光源上随机一个位置采样
		//此时lightInter是光源上某点，inter是光线与物体的交点
		sampleLight(lightInter, pdf_light);

		// 物体表面法线
		auto& N = inter.normal;
		// 灯光表面法线
		auto& NN = lightInter.normal;
		//物体点坐标
		auto& objPos = inter.coords;
		//光源点坐标
		auto& lightPos = lightInter.coords;
		//物体上的交点指向光源上某点
		auto diff = lightPos - objPos;
		auto lightDir = diff.normalized();
		//物体上某点到光源上某点的距离平方
		float lightDistance = diff.x * diff.x + diff.y * diff.y + diff.z * diff.z;

		//从物体打向光源的感光线(感光线为光线传播的逆方向)
		Ray light(objPos, lightDir);
		//该光线与场景求交
		//light2obj是光源上的交点，由于是根据采样来的lightInter和物体上的inter得到的，因此light2obj理应等于lightInter
		Intersection light2obj = intersect(light);

		// 如果反射击中光源，light2obj理应等于lightInter，这里用他们的差小于0.01判断
		if (light2obj.happened && (light2obj.coords - lightPos).norm() < 1e-2)
		{
			//获取改材质的brdf，这里的brdf为漫反射（brdf=Kd/pi）
			//wo=ray.direction，wi=lightDir
			Vector3f f_r = inter.m->eval(ray.direction, lightDir, N);
			//直接光照光 = 光源光 * brdf * 光线和物体角度衰减 * 光线和光源法线角度衰减 / 光线距离 / 该点的概率密度（1/该光源的面积）
			L_dir = lightInter.emit * f_r * dotProduct(lightDir, N) * dotProduct(-lightDir, NN) / lightDistance / pdf_light;
		}

		//俄罗斯轮盘赌，确定是否继续弹射光线
		if (get_random_float() < RussianRoulette)
		{
			//获取半平面上的随机弹射方向。
			//sample随机采样
			Vector3f nextDir = inter.m->sample(ray.direction, N).normalized();
			//定义弹射光线
			Ray nextRay(objPos, nextDir);
			//获取相交点
			Intersection nextInter = intersect(nextRay);
			//如果有相交，且是与物体相交
			if (nextInter.happened && !nextInter.m->hasEmission())
			{
				//该点间接光= 弹射点反射光 * brdf * 角度衰减 / pdf(认为该点四面八方都接收到了该方向的光强，为1/(2*pi)) / 俄罗斯轮盘赌值(强度矫正值)
				float pdf = inter.m->pdf(ray.direction, nextDir, N);
				Vector3f f_r = inter.m->eval(ray.direction, nextDir, N);
				L_indir = castRay(nextRay, depth + 1) * f_r * dotProduct(nextDir, N) / pdf / RussianRoulette;
			}
		}
		//最后返回直接光照和间接光照
		return L_dir + L_indir;
	}
	//如果光线与场景无交点
	return Vector3f(0, 0, 0);
}
————————————————
版权声明：本文为CSDN博主「Elsa的迷弟」的原创文章，遵循CC 4.0 BY-SA版权协议，转载请附上原文出处链接及本声明。
原文链接：https://blog.csdn.net/weixin_44518102/article/details/122495106

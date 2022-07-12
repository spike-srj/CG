#ifndef __MODEL_H__
#define __MODEL_H__

#include <vector>
#include "geometry.h"

class Model {
private:
	//verts_是一个存储一系列3个数（xyz坐标表示一个顶点）的数组
	std::vector<Vec3f> verts_;
	std::vector<std::vector<int> > faces_;
public:
	Model(const char *filename);
	~Model();
	//返回顶点总数
	int nverts();
	//返回面的总数
	int nfaces();
	//返回vert_中第i个顶点坐标，如v -0.000581696 -0.734665 -0.623267
	Vec3f vert(int i);
	//返回faces_中第idx个面，如f 24/1/24 25/2/25 26/3/26
	std::vector<int> face(int idx);
};

#endif //__MODEL_H__

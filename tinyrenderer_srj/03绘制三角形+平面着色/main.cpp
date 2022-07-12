#include <vector>
#include <cmath>
#include "tgaimage.h"   
#include "model.h"      
#include "geometry.h" 


const TGAColor white = TGAColor(255, 255, 255, 255);
const TGAColor red   = TGAColor(255, 0,   0,   255);
const TGAColor green   = TGAColor(0, 255,   0,   255);
Model *model = NULL;
const int width  = 800;
const int height = 800;

void line(int x0, int y0, int x1, int y1, TGAImage& image, TGAColor color) {
    bool steep = false; //标记当前斜率的绝对值是否大于1
	if (std::abs(x0 - x1) < std::abs(y0 - y1)) {
		//斜率绝对值>1了，此时将线段端点各自的x,y坐标对调。
		std::swap(x0, y0);
		std::swap(x1, y1);
		steep = true;
	}
	if (x0 > x1) {  //x0>x1时，对线段端点坐标进行对调
		std::swap(x0, x1);
		std::swap(y0, y1);
	}

    for (int x = x0; x <= x1; x++) {
        float t = (x - x0) / (float)(x1 - x0);
        int y = y0 * (1. - t) + y1 * t;
        if (steep) {
			//如果线段是斜率大于1的，那么线段上的点原本坐标应该是(y,x)
			image.set(y, x, color);
		}
		else {
			image.set(x, y, color);
		}
    }
}



//绘制三角形(坐标1，坐标2，坐标3，tga指针，颜色)
void triangle(Vec2i t0, Vec2i t1, Vec2i t2, TGAImage &image, TGAColor color) {
    //三角形面积为0
    if (t0.y==t1.y && t0.y==t2.y) return;
    //根据y的大小对坐标进行排序
    if (t0.y>t1.y) std::swap(t0, t1);  //如果t1小，将t1赋值给t0
    if (t0.y>t2.y) std::swap(t0, t2);  //与t1比较后，不管谁大最后t0更小，再拿t0与剩下的t2比，如果t2更小，则赋值给t0
    if (t1.y>t2.y) std::swap(t1, t2);  //现在t0最小，但t1，t2大小还未比，让t2成为最大的点
    int total_height = t2.y-t0.y;
    //以高度差作为循环控制变量，此时不需要考虑斜率，因为着色完后每行都会被填充
    for (int i=0; i<total_height; i++) {
        //根据t1将三角形分割为上下两部分
        //如果是在上半部分则second_half为true
        bool second_half = i>t1.y-t0.y || t1.y==t0.y;
        //如果在上半部分segment_height为上半三角形的高，否则为下半三角形的高
        int segment_height = second_half ? t2.y-t1.y : t1.y-t0.y;
        //y-y0/y1-y0，用来计算左端点
        float alpha = (float)i/total_height;
        //当在上半部分时，beta=i-(y1-y0)/segment_height,当在下半时beta=i/segment_height
        float beta  = (float)(i-(second_half ? t1.y-t0.y : 0))/segment_height; 
        //计算A,B两点的坐标,A应该是左端点,B应该是右端点
        Vec2i A =               t0 + (t2-t0)*alpha;
        Vec2i B = second_half ? t1 + (t2-t1)*beta : t0 + (t1-t0)*beta;
        if (A.x>B.x) std::swap(A, B);
        //根据A,B和当前高度对tga着色
        for (int j=A.x; j<=B.x; j++) {
            image.set(j, t0.y+i, color);
        }
    }
}




int main(int argc, char** argv) {

    if (2==argc) {
        model = new Model(argv[1]);  //命令行控制方式构造model
    } else {
        model = new Model("obj/african_head.obj"); //代码方式构造model，african_head
    }
    TGAImage image(width, height, TGAImage::RGB);
    Vec3f light_dir(0,0,-1);//Vec3f light_dir = Vec3f(1.0,1,1);
    for (int i=0; i<model->nfaces(); i++) {
        std::vector<int> face = model->face(i); //创建face数组用于保存一个face的三个顶点坐标，[f 24/1/24 25/2/25 26/3/26 ]
        Vec2i screen_coords[3];
        Vec3f world_coords[3];   //新加入一个数组用于存放三个顶点的世界坐标
        //遍历一个三角形面的三个顶点
        for (int j=0; j<3; j++) {
            //获取其中一个顶点的坐标
            Vec3f v = model->vert(face[j]);
            //将该坐标转换成屏幕坐标（左下角为0,0,z）
            //这里直接一步转换到了屏幕坐标，实际上，中间还有几个转换比如NDC空间、视口变换、透视变换
            screen_coords[j] = Vec2i((v.x*0.8+1.)*width/2., (v.y*0.8+1.)*height/2.);
            world_coords[j]  = v;//世界坐标    即模型坐标
            
        }
        //三角形的法线向量
        Vec3f n = (world_coords[2]-world_coords[0])^(world_coords[1]-world_coords[0]);
        //归一化
        n.normalize();
        float intensity = n*light_dir;//光照强度=法向量*光照方向   即法向量和光照方向重合时，亮度最高
        //强度小于0，说明平面朝向为内  即背面裁剪
        if (intensity>0) {
            //同一个三角形内颜色相同
            triangle(screen_coords[0], screen_coords[1], screen_coords[2], image, TGAColor(intensity*255, intensity*255, intensity*255, 255));
        }
        //triangle(screen_coords[0], screen_coords[1], screen_coords[2], image, TGAColor(255, 255, 255, 255));
    }
    image.flip_vertically();
    image.write_tga_file("output.tga");
    //不删除模型会出问题，但为什么呢？
    delete model;
    return 0;
}



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
        //Vec2i screen_coords[3];
        //Vec3f world_coords[3];   //新加入一个数组用于存放三个顶点的世界坐标
        for (int j=0; j<3; j++) {
            //这里的face已经是[f 24/1/24 25/2/25 26/3/26 ]了，所以face[0]表示24/1/24，也就是第一个顶点的序号
            //v0通过vert函数得到一个顶点的坐标
            //每两个顶点画一条线，循环三次
            Vec3f v0 = model->vert(face[j]);  
            Vec3f v1 = model->vert(face[(j+1)%3]);
            //根据顶点v0和v1画线
            //先要进行模型坐标到屏幕坐标的转换。  (-1,-1)对应(0,0)   (1,1)对应(width,height)
            int x0 = (v0.x+1.)*width/2.;
            int y0 = (v0.y+1.)*height/2.;
            int x1 = (v1.x+1.)*width/2.;
            int y1 = (v1.y+1.)*height/2.;
            //画线
            line(x0,y0, x1,y1, image, white);
        }
    }
    //这句话是干什么用的？
    image.flip_vertically(); // i want to have the origin at the left bottom corner of the image
    
    image.write_tga_file("output.tga");
}



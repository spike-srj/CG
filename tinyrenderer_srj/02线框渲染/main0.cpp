#include <vector>
#include <cmath>
#include "tgaimage.h"   //tga��ͼ��
#include "model.h"      //ģ���࣬��Ҫʵ��ģ�͵Ķ�ȡ
#include "geometry.h"   //���ο⣬��Ҫ������Vec2��Vec3����

//������ɫ
const TGAColor white = TGAColor(255, 255, 255, 255);
const TGAColor red   = TGAColor(255, 0,   0,   255);
Model *model = NULL;
//������ȸ߶�
const int width  = 800;
const int height = 800;

//�����㷨(����1������2��tgaĿ��ָ�룬ָ����ɫ)
void line(int x0, int y0, int x1, int y1, TGAImage &image, TGAColor color) {
    //�ж��߶�б�ʵľ���ֵ�Ƿ����1
    bool steep = false;
    //����1��Ϊtrue������������Ե�x��y�����任Ϊ����y=x��y=-x�ԳƵĵ�
    if (std::abs(x0-x1)<std::abs(y0-y1)) {
        std::swap(x0, y0);
        std::swap(x1, y1);
        steep = true;
    }
    //��֤����2��x,y��������1��x,y��
    if (x0>x1) {
        std::swap(x0, x1);
        std::swap(y0, y1);
    }
    //��ʱx1����x0����б����-1��1֮�䣬��x��ѭ�����Ʊ���
    for (int x=x0; x<=x1; x++) {
        //����x�����߶ζ�Ӧ��y
        float t = (x-x0)/(float)(x1-x0);
        int y = y0*(1.-t) + y1*t;
        //��б�ʴ���1����ʵ����Ϊ(y,x)������Ϊ(x,y)
        if (steep) {
            image.set(y, x, color);
        } else {
            image.set(x, y, color);
        }
    }
}


int main(int argc, char** argv) {
    
    if (2==argc) {
        model = new Model(argv[1]);  //命令行控制方式构造model
    } else {
        model = new Model("obj/cube.obj"); //代码方式构造model
    }
    TGAImage image(width, height, TGAImage::RGB);
    for (int i=0; i<model->nfaces(); i++) {
        std::vector<int> face = model->face(i); //创建face数组用于保存一个face的三个顶点坐标，[f 24/1/24 25/2/25 26/3/26 ]
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
}

#include <vector>
#include <cmath>
#include <cstdlib>
#include <limits>
#include "tgaimage.h"
#include "model.h"
#include "geometry.h"
#include <algorithm>



Vec3f light_dir(0.2,0.15,-1);
Vec3f camera(0,0,3);//Vec3f light_dir = Vec3f(1.0,1,1);
Model *model = NULL;
const int width  = 800;
const int height = 800;
const int depth  = 255;


//4d-->3d
//除以最后一个分量。（当最后一个分量为0，表示向量）
//不为0，表示坐标
Vec3f m2v(Matrix m) {
    return Vec3f(m[0][0]/m[3][0], m[1][0]/m[3][0], m[2][0]/m[3][0]);
}

//3d-->4d
//添加一个1表示坐标
Matrix v2m(Vec3f v) {
    Matrix m(4, 1);
    m[0][0] = v.x;
    m[1][0] = v.y;
    m[2][0] = v.z;
    m[3][0] = 1.f;
    return m;
}

//视角矩阵
//将物体x，y坐标(-1,1)转换到屏幕坐标(100,700)    1/8width~7/8width
//zbuffer(-1,1)转换到0~255
Matrix viewport(int x, int y, int w, int h) {
    Matrix m = Matrix::identity(4);
    //第4列表示平移信息
    m[0][3] = x+w/2.f;
    m[1][3] = y+h/2.f;
    m[2][3] = depth/2.f;
    //对角线表示缩放信息
    m[0][0] = w/2.f;
    m[1][1] = h/2.f;
    m[2][2] = depth/2.f;
    return m;
}


//后面画三角形就用不到画线了
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


/*
//计算质心坐标
Vec3f barycentric(Vec3f A, Vec3f B, Vec3f C, Vec3f P) {
    Vec3f s[2];
    //计算[AB,AC,PA]的x和y分量
    for (int i=2; i--; ) {
        s[i][0] = C[i]-A[i];
        s[i][1] = B[i]-A[i];
        s[i][2] = A[i]-P[i];
    }
    //[u,v,1]和[AB,AC,PA]对应的x和y向量都垂直，所以叉乘
    Vec3f u = cross(s[0], s[1]);
    //三点共线时，会导致u[2]为0，此时返回(-1,1,1)
    if (std::abs(u[2])>1e-2)
        //若1-u-v，u，v全为大于0的数，表示点在三角形内部
        return Vec3f(1.f-(u.x+u.y)/u.z, u.y/u.z, u.x/u.z);
    return Vec3f(-1,1,1);
}
*/


//第一种绘制三角形的办法
//绘制三角形(坐标1，坐标2，坐标3，tga指针，颜色)
void triangle(Vec3i t0, Vec3i t1, Vec3i t2, Vec2i uv0, Vec2i uv1, Vec2i uv2, float *zbuffer, TGAImage &image, float intensity) {
    //三角形面积为0
    if (t0.y==t1.y && t0.y==t2.y) return;
    //根据y的大小对坐标进行排序
    if (t0.y>t1.y) {std::swap(t0, t1); std::swap(uv0, uv1);}  //如果t1小，将t1赋值给t0
    if (t0.y>t2.y) {std::swap(t0, t2); std::swap(uv0, uv2);} //与t1比较后，不管谁大最后t0更小，再拿t0与剩下的t2比，如果t2更小，则赋值给t0
    if (t1.y>t2.y) {std::swap(t1, t2); std::swap(uv1, uv2);}  //现在t0最小，但t1，t2大小还未比，让t2成为最大的点
    int total_height = t2.y-t0.y;
    //以高度差作为循环控制变量，此时不需要考虑斜率，因为着色完后每行都会被填充
    for (int i=0; i<total_height; i++) {
        //根据t1将三角形分割为上下两部分
        bool second_half = i>t1.y-t0.y || t1.y==t0.y;
        int segment_height = second_half ? t2.y-t1.y : t1.y-t0.y;
        float alpha = (float)i/total_height;
        float beta  = (float)(i-(second_half ? t1.y-t0.y : 0))/segment_height; 
        //计算A,B两点的坐标,A应该是左端点,B应该是右端点。此时的t0是屏幕坐标，整形，不能改成浮点型，但计算时乘以的比例应该是浮点型
        Vec3i A =               t0 + Vec3f(t2-t0)*alpha;
        Vec3i B = second_half ? t1 + Vec3f(t2-t1)*beta : t0 + Vec3f(t1-t0)*beta;
        //计算UV
        Vec2i uvA =               uv0 +      (uv2-uv0)*alpha;  //第一次插值
        Vec2i uvB = second_half ? uv1 +      (uv2-uv1)*beta : uv0 +      (uv1-uv0)*beta;
        if (A.x>B.x) std::swap(A, B);
        //根据A,B和当前高度对tga着色
        for (int j=A.x; j<=B.x; j++) {
            //计算当前点在AB之间的比例
            float phi = B.x==A.x ? 1. : (float)(j-A.x)/(float)(B.x-A.x);
            //计算出当前点的坐标,A，B保存了z轴信息
            Vec3f   P = Vec3f(A) + Vec3f(B-A)*phi;
            Vec2i uvP =     uvA +   (uvB-uvA)*phi;  //第二次插值
            if (P.x < width && P.y < height)
            {
                //计算当前zbuffer下标=P.x+P.y*width
                int idx = P.x+P.y*width;
                //当前点的z大于zbuffer信息，覆盖掉，并更新zbuffer
                if (zbuffer[idx]<P.z) {
                    zbuffer[idx] = P.z;
                    TGAColor color = model->diffuse(uvP);
                    image.set(P.x, P.y, TGAColor(color.r*intensity, color.g*intensity, color.b*intensity));
                }
            }
        }
    }
}

/*
//第二种绘制三角形的办法，开头的包围盒与101里的方法一致，只是后面判断是否在三角形内的方法不同。这里的方法更好，直接给出了质心坐标
//绘制三角形(坐标数组，zbuffer指针，tga指针，颜色)
//pts是三组顶点
void triangle(Vec3f *pts, float *zbuffer, TGAImage &image, TGAColor color) {
    Vec2f bboxmin( std::numeric_limits<float>::max(),  std::numeric_limits<float>::max());
    Vec2f bboxmax(-std::numeric_limits<float>::max(), -std::numeric_limits<float>::max());
    Vec2f clamp(image.get_width()-1, image.get_height()-1);
    //确定三角形的边框
    for (int i=0; i<3; i++) {
        for (int j=0; j<2; j++) {
            bboxmin[j] = std::max(0.f,      std::min(bboxmin[j], pts[i][j]));
            bboxmax[j] = std::min(clamp[j], std::max(bboxmax[j], pts[i][j]));
        }
    }
    Vec3f P;
    //遍历边框中的每一个点
    for (P.x=bboxmin.x; P.x<=bboxmax.x; P.x++) {
        for (P.y=bboxmin.y; P.y<=bboxmax.y; P.y++) {
            //计算质心
            if (P.x > 600 && P.y > 500)
            {
                P.x += 0.01;
            }
            //bc_screen就是质心坐标，pts[0]是一个顶点坐标
            Vec3f bc_screen  = barycentric(pts[0], pts[1], pts[2], P);
            //质心坐标有一个负值，说明点在三角形外
            //bc_screen是质心坐标，xyz对应三个顶点，不是真正意义上的xyz
            if (bc_screen.x<0 || bc_screen.y<0 || bc_screen.z<0) continue;
            P.z = 0;
            //计算zbuffer，并且每个顶点的z值乘上对应的质心坐标分量
            //pts[i][2]是某个顶点的z值，这里三角形内某个点P的质心坐标的z分量是由三个顶点的z分量乘以他们对应的质心坐标
            for (int i=0; i<3; i++) P.z += pts[i][2]*bc_screen[i];
            if (zbuffer[int(P.x+P.y*width)]<P.z) {
                zbuffer[int(P.x+P.y*width)] = P.z;
                image.set(P.x, P.y, color);
            }
        }
    }
}
*/


Vec3f world2screen(Vec3f v) {
    return Vec3f(int((v.x+1.)*width/2.+.5), int((v.y+1.)*height/2.+.5), v.z);
}

int main(int argc, char** argv) {
    if (2==argc) {
        model = new Model(argv[1]);  //命令行控制方式构造model
    } else {
        model = new Model("obj/african_head.obj"); //代码方式构造model，african_head
    }
    //创建zbuffer，大小为画布大小
    float *zbuffer = new float[width*height];
    //初始化zbuffer，设定一个很小的值
    for (int i=0; i<width*height; i++) {
        //��ʼ��zbuffer
        zbuffer[i] = std::numeric_limits<int>::min();
    }
    TGAImage image(width, height, TGAImage::RGB);
    

    Matrix Projection = Matrix::identity(4);
    //初始化视角矩阵
    Matrix ViewPort   = viewport(width/8, height/8, width*3/4, height*3/4);
    //投影矩阵[3][2]=-1/c，c为相机z坐标
    Projection[3][2] = -1.f/camera.z;


    for (int i=0; i<model->nfaces(); i++) {
        std::vector<int> face = model->face(i); //创建face数组用于保存一个face的三个顶点坐标，[f 24/1/24 25/2/25 26/3/26 ]
        Vec3f screen_coords[3];
        Vec3f world_coords[3];   //新加入一个数组用于存放三个顶点的世界坐标
        for (int j=0; j<3; j++) {
            Vec3f v = model->vert(face[j]);
            screen_coords[j] = m2v(ViewPort * Projection * v2m(v));
            world_coords[j]  = v;//世界坐标    即模型坐标
            
        }
        //先不急着把屏幕坐标传入triangle函数，先把颜色算出来（这里指获得uv坐标）
        Vec3f n = (world_coords[2]-world_coords[0])^(world_coords[1]-world_coords[0]);
        n.normalize();
        //计算光照
        float intensity = n*light_dir;
        intensity = std::min(std::abs(intensity),1.f);
        if (intensity>0) {
            Vec2i uv[3];
            for (int k=0; k<3; k++) {
                //i是面的id，k是一个面里顶点的id
                uv[k] = model->uv(i, k);
            }
            //绘制三角形，这里不需要再传入颜色。三角函数里会根据uv坐标设置颜色
            triangle(screen_coords[0], screen_coords[1], screen_coords[2], uv[0], uv[1], uv[2], zbuffer, image, intensity);
        }
        //triangle(screen_coords, zbuffer, image, TGAColor(intensity * 255, intensity * 255, intensity * 255, 255));
    }
    image.flip_vertically();
    image.write_tga_file("output.tga");

    //输出zbuffer
    {
        TGAImage zbimage(width, height, TGAImage::GRAYSCALE);
        for (int i=0; i<width; i++) {
            for (int j=0; j<height; j++) {
                zbimage.set(i, j, TGAColor(zbuffer[i+j*width], 1));
            }
        }
        zbimage.flip_vertically();
        zbimage.write_tga_file("zbuffer.tga");
    }
    //不删除模型会出问题，但为什么呢？
    delete model;
    delete [] zbuffer;
    return 0;
}

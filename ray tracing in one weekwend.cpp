#include <iostream>

int main() {
    const int image_width = 200;
    const int image_height = 100;

    std::cout << "P3\n" << image_width << ' ' << image_height << "\n255\n";
    //对于像素来说, 每一行是从左往右写入的。行从上开始往下写入的。
    for (int j = image_height-1; j >= 0; --j) {
        for (int i = 0; i < image_width; ++i) {
            //通常我们把RGB通道的值限定在0.0到1.0。我们之后计算颜色值的时候将使用一个动态的范围, 
            //这个范围并不是0到1。但是在使用这段代码输出图像之前, 我们将把颜色映射到0到1。所以这部分输出图像代码不会改变。
            auto r = double(i) / image_width;
            auto g = double(j) / image_height;
            auto b = 0.2;
            int ir = static_cast<int>(255.999 * r);
            int ig = static_cast<int>(255.999 * g);
            int ib = static_cast<int>(255.999 * b);
            std::cout << ir << ' ' << ig << ' ' << ib << '\n';
        }
    }
}



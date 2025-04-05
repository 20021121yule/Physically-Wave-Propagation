#include "Wave_Propagation.h"
#include <iostream>

using namespace std;

int main() {
    // 设置参数
    string image_path = "/Users/yule/Desktop/test.jpg"; // 替换成你的图片路径
    double lambda = 632.8e-9; // HeNe激光波长（米）
    double dz = 0.01; // 传播距离（米）
    double pixel_size = 10e-6; // 像素尺寸（米）

    // 创建 Wave_Propagation 对象
    Wave_Propagation wave(image_path, lambda, dz, pixel_size);

    // 执行传播
    wave.propagate();

    cout << "光波传播计算完成！" << endl;

    return 0;
}

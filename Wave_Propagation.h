//
// Created by 余乐 on 25-3-28.
//

#ifndef WAVE_PROPAGATION_H
#define WAVE_PROPAGATION_H

#include <opencv2/opencv.hpp>
#include <fftw3.h>

using namespace cv;
using namespace std;

class Wave_Propagation {
private:
    // 物理常量
    double lambda; // 波长 (单位：米)
    double dz; // 传播距离 (单位：米)
    double pixel_size; // 像素尺寸 (单位：米)
    double k; // 总波矢

    // 图片参数
    string pic_path; // 文件路径
    Mat image; // 图片像素矩阵
    int rows, cols; // 图片像素的行数、列数

    // 计算得到的输出光场
    Mat output_OptField;

public:
    // 辅助函数
    static void fftshift(fftw_complex *data, int rows, int cols);
    static void ifftshift(fftw_complex *data, int rows, int cols);

    // 构造函数
    Wave_Propagation(const string &input_pic_path, const double &input_lambda,
                     const double &input_dz, const double &input_pixel_size);

    //传播过程
    void propagate();
};


#endif //WAVE_PROPAGATION_H

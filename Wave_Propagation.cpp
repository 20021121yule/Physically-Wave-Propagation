//
// Created by 余乐 on 25-3-28.
//

#include "Wave_Propagation.h"

void Wave_Propagation::fftshift(fftw_complex *data, int rows, int cols) {
    int halfRows = rows / 2;
    int halfCols = cols / 2;

    for (int i = 0; i < halfRows; i++) {
        for (int j = 0; j < halfCols; j++) {
            // 四象限交换
            swap(data[i * cols + j], data[(i + halfRows) * cols + (j + halfCols)]);
            swap(data[i * cols + (j + halfCols)], data[(i + halfRows) * cols + j]);
        }
    }
}

void Wave_Propagation::ifftshift(fftw_complex *data, int rows, int cols) {
    int halfRows = rows / 2;
    int halfCols = cols / 2;

    for (int i = 0; i < halfRows; i++) {
        for (int j = 0; j < halfCols; j++) {
            // 交换四个象限，使其恢复原状
            swap(data[i * cols + j], data[(i + halfRows) * cols + (j + halfCols)]);
            swap(data[i * cols + (j + halfCols)], data[(i + halfRows) * cols + j]);
        }
    }
}

Wave_Propagation::Wave_Propagation(const string &input_pic_path, const double &input_lambda, const double &input_dz,
                                   const double &input_pixel_size) : pic_path(input_pic_path),
                                                                     lambda(input_lambda), dz(input_dz),
                                                                     pixel_size(input_pixel_size) {
    // 读取图片并转换为灰度图
    this->image = imread(pic_path, IMREAD_GRAYSCALE);
    if (this->image.empty()) {
        cerr << "Could not open or find the image" << endl;
        this->rows = 0;
        this->cols = 0;
        return;
    }

    this->rows = this->image.rows;
    this->cols = this->image.cols;

    this->k = 2 * M_PI / input_lambda; // 总波矢
}


void Wave_Propagation::propagate() {
    fftw_complex *in, *out;
    fftw_plan p;
    in = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * rows * cols);
    out = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * rows * cols);

    // 先将图片像素分配给in变量
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            in[i * cols + j][0] = image.at<uchar>(i, j);
            in[i * cols + j][1] = 0;
        }
    }

    p = fftw_plan_dft_2d(rows, cols, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(p);

    // 对out做shift处理
    fftshift(out, rows, cols);

    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            // 计算kz波矢
            double k_x = (2 * M_PI / (pixel_size * rows)) * (i - rows / 2);
            double k_y = (2 * M_PI / (pixel_size * cols)) * (j - cols / 2);
            double k_z = (k_x * k_x + k_y * k_y < k * k) ? sqrt(k * k - k_x * k_x - k_y * k_y) : 0;

            double H_real = cos(k_z * dz);
            double H_imag = sin(k_z * dz);

            double real = out[i * cols + j][0];
            double imag = out[i * cols + j][1];

            out[i * cols + j][0] = real * H_real - imag * H_imag;
            out[i * cols + j][1] = imag * H_real + real * H_imag;
        }
    }

    // 对out做ifftshift处理
    ifftshift(out, rows, cols);

    // 逆变换
    p = fftw_plan_dft_2d(rows, cols, out, in, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(p);

    Mat amplitude(rows, cols, CV_64F);
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            double out_real = in[i * cols + j][0];
            double out_imag = in[i * cols + j][1];

            amplitude.at<double>(i, j) = sqrt(out_real * out_real + out_imag * out_imag);
        }
    }

    this->output_OptField.release();
    normalize(amplitude, this->output_OptField, 0, 255, NORM_MINMAX, CV_8U);


    // 显示窗口
    string windowName = "output_OptField at " + to_string(dz) + " m";

    namedWindow(windowName, WINDOW_NORMAL);
    imshow(windowName, this->output_OptField);
    waitKey(0);
    destroyAllWindows(); // 关闭所有窗口

    // 释放内存
    fftw_destroy_plan(p);
    fftw_free(in);
    fftw_free(out);
}

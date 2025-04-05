pth = '/Users/yule/Desktop/test.jpg'; % 输入路径
imag_RGB = imread(pth);

% 转换为灰度图，并转换为 double 类型
imag_RGB = double(imag_RGB);
imag_GrayScale = 0.299 * imag_RGB(:,:,1) + 0.587 * imag_RGB(:,:,2) + 0.114 * imag_RGB(:,:,3);
imag_GrayScale = double(imag_GrayScale); % 确保精度

% 显示灰度图
figure;
imagesc(imag_GrayScale);
colormap('gray');
title('灰度图');
colorbar;
axis image;

% 物理参数
lambda_ = 632.8e-9;  % 波长（米）
dz = 0.01;           % 传播距离（米）
pixel_size = 10e-6;  % 像素尺寸（米）
k = 2 * pi / lambda_; % 波矢

% 计算 FFT 并中心化
inputOptField_fft = fft2(imag_GrayScale);
inputOptField_fft_shifted = fftshift(inputOptField_fft);

% 角谱定理计算
[x_pixels, y_pixels] = size(imag_GrayScale);

% 计算 kx 和 ky 轴
kx_vals = (2 * pi / (pixel_size * x_pixels)) * ((0:x_pixels-1) - x_pixels / 2);
ky_vals = (2 * pi / (pixel_size * y_pixels)) * ((0:y_pixels-1) - y_pixels / 2);

% 生成 2D 频率网格
[k_x_mesh, k_y_mesh] = meshgrid(kx_vals, ky_vals);

% 计算 kz，防止复数出现
k_z = sqrt(k^2 - k_x_mesh.^2 - k_y_mesh.^2);

% 计算传播函数
H = exp(1i * k_z * dz);

% 传播计算
outputOptField_fft_shifted = inputOptField_fft_shifted .* H;
outputOptField_fft = ifftshift(outputOptField_fft_shifted);
outputOptField = ifft2(outputOptField_fft);  % ← 这里要用 ifft2

% 取振幅，表示光场强度
outputOptField_amp = abs(outputOptField);

% 显示传播后的光场
figure;
imagesc(outputOptField_amp);
colormap('gray');
title('传播后的光场（振幅）');
axis image;
colorbar;
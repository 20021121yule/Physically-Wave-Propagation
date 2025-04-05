cmake_minimum_required(VERSION 3.30)
project(C___OE)

set(CMAKE_CXX_STANDARD 20)

# 设置 OpenCV 的安装路径
set(OpenCV_DIR "/usr/local/lib/cmake/opencv4")

# 查找 OpenCV 包
find_package(OpenCV REQUIRED)

# 包含 OpenCV 头文件
include_directories(${OpenCV_INCLUDE_DIRS})

add_executable(C___OE main.cpp
        Wave_Propagation.cpp
        Wave_Propagation.h
)

# 链接 OpenCV 库
target_link_libraries(C___OE ${OpenCV_LIBS})

# ✅ 手动指定 FFTW 的 include 和库路径（去掉 pkg-config）
set(FFTW_INCLUDE_DIR "/usr/local/include")
set(FFTW_LIBRARY "/usr/local/lib/libfftw3.dylib")

# 添加 FFTW 头文件路径
include_directories(${FFTW_INCLUDE_DIR})

# 链接 FFTW 库
target_link_libraries(C___OE ${FFTW_LIBRARY})
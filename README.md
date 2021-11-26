# TransformCoordinate
##1. 下载及编译 proj.4
for ubuntu:
'sudo apt-get install proj-bin '

****
## 第2种方案：下载源码编译
下载发布的源码：https://proj.org/download.html#current-release \
cd proj-8.2.0 \
mkdir build \
cd build \
// From the build directory you can now configure CMake, build and install the binaries:

cmake .. \
cmake --build . \
sudo cmake --build . --target install


****
编译后的库的头文件位置：
/usr/local/include

静态库位置：
/usr/local/lib

动态库位置：
/usr/local/bin

相关编译文件位置：\
/usr/local/lib/cmake/proj4

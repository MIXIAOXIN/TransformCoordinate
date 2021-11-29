#pragma once
#include <math.h>
#include <limits>
//#include <proj/>

const double PI = acos(-1.0);
const double ee = 0.0066943799013;//wgs-84 椭球第一偏心率平方;
const double ee2 = 0.00673949674227;//wgs-84 椭球第二偏心率平方;
const double a = 6378137.0000000000;//wgs-84 椭球长半轴;
const double b = 6356752.3142451795; //wgs-84 椭球短半轴

/*
// @brief: 将wgs84 地心直角坐标转为wgs84经纬度
(x, y, z) 为地心直角坐标系下的坐标
(reference_ellipsoid_a, reference_ellipsoid_b) 为参考椭球的长轴和短轴 （）
epslon 为残差系数
(latitude, longitude, height) 为转换后的纬度、经度、大地高
*/
bool ecef2lla(const double x, const double y, const double z,
	const double reference_ellipsoid_a, const double reference_ellipsoid_b,
	const double epslon,
	double& latitude, double& longitude, double& height);

/*
// @brief: 将wgs84 经纬度大地高转为Gussian 平面投影坐标
(latitude, longitude, height) 为转换后的纬度、经度、大地高
centralLon 为中央经线的经度
(x, y, z) 为高斯投影的平面坐标
*/
void WGS84ToGauss3(const double latitude, const double longitude, const double height, 
	const double &centralLon, 
	double& utm_x, double& utm_y, double& utm_z);
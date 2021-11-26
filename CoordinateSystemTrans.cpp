#include "CoordinateSystemTrans.h"

bool ecef2lla(const double x, const double y, const double z,
	const double reference_ellipsoid_a, const double reference_ellipsoid_b,
	const double epslon,
	double& latitude, double& longitude, double& height)
{
	// 1. 计算椭球偏心率e：
	double e = (reference_ellipsoid_a * reference_ellipsoid_a - reference_ellipsoid_b * reference_ellipsoid_b) / (reference_ellipsoid_a * reference_ellipsoid_a);

	// 计算经度：longtitude:
	longitude = atan(y / x);


	// 2. 计算基准椭球体的曲率半径N
	// 3.1 计算中间变量p = sqrt(x^2 + y^2)，为当前点在地心直角坐标系xoy平面上的投影长度
	// 3.2 计算中间变量angle_B:当前点在地心直角坐标系xoy平面的投影与X轴的夹角
	// 从假设PE = PD开始，不断迭代，直至 |PE' - PE| < epslon
	double p = sqrt(x * x + y * y); 
	double angle_B = atan(y / (x * (1 - e * e)));
	double W = std::sqrt(1 - e * e * sin(angle_B) * sin(angle_B));

	// 迭代计算PE，直至满足残差条件
	double PE = DBL_MAX, PE_new = z;
	double angle_phi = 0.0, curvature_N = 0.0;
	while (std::fabs(PE - PE_new) > epslon)
	{
		PE = PE_new;
		angle_phi = atan(PE/ p);		
		double curvature_N = reference_ellipsoid_a / W;
		PE_new = z + curvature_N * e * e * sin(angle_phi);
	}

	// 4. 计算纬度、大地高
	latitude = angle_phi;
	height = sqrt(p * p + PE * PE) - curvature_N;

	return true;
}



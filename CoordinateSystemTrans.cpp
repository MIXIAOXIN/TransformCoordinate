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


void WGS84ToGauss3(const double latitude, const double longitude, const double height,
	const double &centralLon,
	double& utm_x, double& utm_y, double& utm_z) {
	double N, L0;
	double W[3];
	W[0] = latitude;
	W[1] = longitude;
	W[2] = height;

	L0 = centralLon;

	N = a / sqrt(1 - ee * sin(W[0] * PI / 180.0)*sin(W[0] * PI / 180.0));
	double m0, m2, m4, m6, m8;
	double a0, a2, a4, a6, a8;
	m0 = a * (1 - ee);
	m2 = 3 * ee*m0 / 2;
	m4 = 5 * ee*m2 / 4;
	m6 = 7 * ee*m4 / 6;
	m8 = 9 * ee*m6 / 8;
	a0 = m0 + m2 / 2 + 3 * m4 / 8 + 5 * m6 / 16 + 35 * m8 / 128;
	a2 = m2 / 2 + m4 / 2 + 15 * m6 / 32 + 7 * m8 / 16;
	a4 = m4 / 8 + 3 * m6 / 16 + 7 * m8 / 32;
	a6 = m6 / 32 + m8 / 16;
	a8 = m8 / 128;
	double X;
	X = a0 * W[0] * PI / 180.0 - a2 * sin(2 * W[0] * PI / 180.0) / 2 + a4 * sin(4 * W[0] * PI / 180.0) / 4 - a6 * sin(6 * W[0] * PI / 180.0) / 6 + a8 * sin(8 * W[0] * PI / 180.0) / 8;
	double l, t, yb;
	l = (W[1] - L0)*PI / 180.0;
	t = tan(W[0] * PI / 180.0);
	yb = sqrt(ee2)*cos(W[0] * PI / 180.0);

	utm_y = X + N * t*pow(cos(W[0] * PI / 180.0)*l, 2) / 2 + N * t*(5 - t * t + 9 * yb*yb + 4 * pow(yb, 4))*pow(cos(W[0] * PI / 180.0)*l, 4) / 24 + N * t*(61 - 58 * t*t + pow(t, 4))*pow(cos(W[0] * PI / 180.0)*l, 6) / 720;
	utm_x = N * cos(W[0] * PI / 180.0)*l + N / 6.0*(1 - t * t + yb * yb)*pow((cos(W[0] * PI / 180.0)*l), 3) + N / 120.0*(5.0 - 18.0*t*t + pow(t, 4) + 14 * yb*yb - 58 * yb*yb*t*t)*pow((cos(W[0] * PI / 180.0)*l), 5) + 500000;
	utm_x += 38 * 1000000;
	utm_z = height;
}
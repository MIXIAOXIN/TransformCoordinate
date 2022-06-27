#pragma once
#include <math.h>
#include <limits>
#include<cfloat>
#include <iostream>
//#include <proj/>

const double PI = acos(-1.0);
const double ee = 0.0066943799013;//wgs-84 �����һƫ����ƽ��;
const double ee2 = 0.00673949674227;//wgs-84 
const double a = 6378137.0000000000;//wgs-84 reference_ellipsoid_a;
const double b = 6356752.3142451795; //wgs-84 reference_ellipsoid_b;
const double f = 1.0/298.257223563;  // ellipsoid ratio
const double k0 = 0.9996; // UTM scale on the central meridian

/*
// @brief: ��wgs84 ����ֱ������תΪwgs84��γ��
(x, y, z) Ϊ����ֱ������ϵ�µ�����
(reference_ellipsoid_a, reference_ellipsoid_b) Ϊ�ο�����ĳ���Ͷ��� ����
epslon Ϊ�в�ϵ��
(latitude, longitude, height) Ϊת�����γ�ȡ����ȡ���ظ�
*/
bool ecef2lla(const double x, const double y, const double z,
	const double reference_ellipsoid_a, const double reference_ellipsoid_b,
	const double epslon,
	double& latitude, double& longitude, double& height);

/*
// @brief: ��wgs84 ��γ�ȴ�ظ�תΪGussian ƽ��ͶӰ����
(latitude, longitude, height): 
centralLon: central longitude for the point: ([l / 6] + 1) * 6 -3
(x, y, z): utm coordinate 
*/
void WGS84ToGauss3(const double latitude, const double longitude, const double height, 
	const double &centralLon, 
	double& utm_x, double& utm_y, double& utm_z);

int calLongtitudinalZone(double lon);

double calCentralMeridianLongtitude(int zone);

void LongtitudeLatitude2UTM(const double latitude, const double longitude,
							double& utm_x, double& utm_y);
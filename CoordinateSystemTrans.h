#pragma once
#include <math.h>
#include <limits>
//#include <proj/>

const double PI = acos(-1.0);
const double ee = 0.0066943799013;//wgs-84 �����һƫ����ƽ��;
const double ee2 = 0.00673949674227;//wgs-84 ����ڶ�ƫ����ƽ��;
const double a = 6378137.0000000000;//wgs-84 ���򳤰���;
const double b = 6356752.3142451795; //wgs-84 ����̰���

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
(latitude, longitude, height) Ϊת�����γ�ȡ����ȡ���ظ�
centralLon Ϊ���뾭�ߵľ���
(x, y, z) Ϊ��˹ͶӰ��ƽ������
*/
void WGS84ToGauss3(const double latitude, const double longitude, const double height, 
	const double &centralLon, 
	double& utm_x, double& utm_y, double& utm_z);
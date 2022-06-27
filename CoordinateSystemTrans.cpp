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
	//utm_x += 38 * 1000000;
	utm_z = height;
}

int calLongtitudinalZone(double lon)
{
	int zon = std::floor((lon+180.0)/6) + 1;
	return zon;
}

double calCentralMeridianLongtitude(int zone)
{
	double c_lon = (zone-1) * 6 - 180.0 + 3.0;
	return c_lon;
}

void LongtitudeLatitude2UTM(const double latitude, const double longitude,
							double& utm_x, double& utm_y)
{
	// calculate central meridian longtitude
	int zone = calLongtitudinalZone(longitude);
	double cen_lon = calCentralMeridianLongtitude(zone);
	double cen_lon_r = cen_lon * PI / 180.0;
	double lat_r = latitude * PI / 180.0, lon_r = longitude * PI / 180.0 - cen_lon_r;
	 
	// 
	double falseEasting = 500e3, falseNorthing = 10000e3;
	// ---- easting, northing: Karney 2011 Eq 7-14, 29, 35:

    double e = std::sqrt(f*(2.0-f)); // eccentricity
    double n = f / (2.0 - f);        // 3rd flattening
    double n2 = n*n, n3 = n*n2, n4 = n*n3, n5 = n*n4, n6 = n*n5;
	double cos_lon = std::cos(lon_r), sin_lon = std::sin(lon_r), tan_lon = std::tan(lon_r);
	double tan_lat = std::tan(lat_r); // τ ≡ tanφ, τʹ ≡ tanφʹ; prime (ʹ) indicates angles on the conformal sphere
    double sigma = std::sinh(e*std::atanh(e*tan_lat/std::sqrt(1+tan_lat*tan_lat)));
	
	double tan_lat_inv = tan_lat * std::sqrt(1+sigma*sigma) - sigma*std::sqrt(1+tan_lat*tan_lat);
	double epxl = std::atan2(tan_lat_inv, cos_lon); // ξʹ
    double phal = std::asinh(sin_lon / std::sqrt(tan_lat_inv * tan_lat_inv + cos_lon*cos_lon)); // ηʹ

    double A = a/(1.0+n) * (1.0 + 1.0/4.0*n2 + 1.0/64.0*n4 + 1.0/256.0*n6); // 2πA is the circumference of a meridian

	// alpha α 
    double alpha[6] = {1.0/2.0*n - 2.0/3.0*n2 + 5.0/16.0*n3 + 41.0/180.0*n4 - 127.0/288.0*n5 + 7891.0/37800.0*n6,
                  13.0/48.0*n2 -  3.0/5.0*n3 + 557.0/1440.0*n4 + 281.0/630.0*n5 - 1983433.0/1935360.0*n6,
                           61.0/240.0*n3 -  103.0/140.0*n4 + 15061.0/26880.0*n5 + 167603.0/181440.0*n6,
                                   49561.0/161280.0*n4 - 179.0/168.0*n5 + 6601661.0/7257600.0*n6,
                                                     34729.0/80640.0*n5 - 3418889.0/1995840.0*n6,
                                                                  212378941.0/319334400.0*n6 };

    double epxll = epxl;
    for (int j=0; j<6; j++){
		epxll += alpha[j] * std::sin(2*j*epxl) * std::cosh(2*j*phal);
	} 


    double phall = phal;
    for (int j=0; j<6; j++)
	{
		phall += alpha[j] * std::cos(2*j*epxl) * std::sinh(2*j*phal);
	} 

	double x = k0 * A * phall;
    double y = k0 * A * epxll;

    // ---- convergence: Karney 2011 Eq 23, 24

    double p_l = 1.0;
    for (int j=0; j<6; j++)
	{
		p_l += 2*j*alpha[j] * std::cos(2*j*epxl) * std::cosh(2*j*phal);
	} 
    double q_l = 0.0;
    for (int j=0; j<6; j++)
	{
		q_l += 2*j*alpha[j] * std::sin(2*j*epxl) * std::sinh(2*j*phal);
	}

    const double gamma_l = std::atan(tan_lat_inv / std::sqrt(1+tan_lat_inv*tan_lat_inv)*tan_lon); // γ'
    const double gamma_ll  = std::atan2(q_l, p_l);  // γʺ
	const double gamma = gamma_l + gamma_ll; // γ

    // ---- scale: Karney 2011 Eq 25
	const double sin_lat = std::sin(lat_r); // φ
    const double k_l = std::sqrt(1.0 - e*e*sin_lat*sin_lat) * std::sqrt(1 + tan_lat*tan_lat) / std::sqrt(tan_lat_inv*tan_lat_inv + cos_lon*cos_lon);
    const double k_ll = A / a * std::sqrt(p_l*p_l + q_l*q_l);

    const double k = k0 * k_l * k_ll;

    // ------------
	// shift x/y to false origins
    x = x + falseEasting;             // make x relative to false easting
    if (y < 0)
	{y = y + falseNorthing;} // make y in southern hemisphere relative to false northing

    utm_x = x;
	utm_y = y;
}
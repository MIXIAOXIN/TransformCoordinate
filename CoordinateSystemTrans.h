#pragma once
#include <math.h>
#include <limits>

#include "proj_api.h"

bool ecef2lla(const double x, const double y, const double z,
	const double reference_ellipsoid_a, const double reference_ellipsoid_b,
	const double epslon,
	double& latitude, double& longitude, double& height);




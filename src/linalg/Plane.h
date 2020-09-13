#pragma once

#include "linalg.h"
#include "point3d.h"

class Plane {

private:

	double a, b, c, d;

public:

	Plane(double a, double b, double c, double d);

	Plane(double x, double y, double z, double ax, double ay);

	Plane(Point3D& a, Point3D& b, Point3D& c);

	double signum(const Point3D& p);

};

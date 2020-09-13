#include "Plane.h"

Plane::Plane(double a, double b, double c, double d) {
	this->a = a;
	this->b = b;
	this->c = c;
	this->d = d;
}

Plane::Plane(double x, double y, double z, double ax, double ay) {
	// z0 = z + (x0 - x) * tan(ax) + (y0 - y) * tan(ay)
	a = tan(ax);
	b = tan(ay);
	c = -1;
	d = z - a * x - b * y;
}

Plane::Plane(Point3D& a, Point3D& b, Point3D& c) {
	// https://en.wikipedia.org/wiki/Plane_(geometry)#Method_2

	d = -det(
		a.x, a.y, a.z,
		b.x, b.y, b.z,
		c.x, c.y, c.z);

	this->a = det(
		1, a.y, a.z,
		1, b.y, b.z,
		1, c.y, c.z);

	this->b = det(
		a.x, 1, a.z,
		b.x, 1, b.z,
		c.x, 1, c.z);

	this->c = det(
		a.x, a.y, 1,
		b.x, b.y, 1,
		c.x, c.y, 1);
}

double Plane::signum(const Point3D& p) {
	return a * p.x + b * p.y + c * p.z + d;
}

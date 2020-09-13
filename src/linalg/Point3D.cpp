#include "Point3D.h"

Point3D::Point3D(double x, double y, double z) {
	this->x = x;
	this->y = y;
	this->z = z;
}

Point3D Point3D::multiply(double k) {
	return Point3D(k * x, k * y, k * z);
}

Point3D Point3D::add(const Point3D& p) {
	return Point3D(x + p.x, y + p.y, z + p.z);
}

Point3D Point3D::subtract(const Point3D& p) {
	return Point3D(x - p.x, y - p.y, z - p.z);
}

double Point3D::scalar(const Point3D& p) {
	return x * p.x + y * p.y + z * p.z;
}

Point3D Point3D::vector(const Point3D& p) {
	return Point3D(y * p.z - z * p.y, z * p.x - x * p.z, x * p.y - y * p.x);
}

double Point3D::dist(const Point3D& p) {
	return subtract(p).norm();
}

double Point3D::norm() {
	return sqrt(x * x + y * y + z * z);
}

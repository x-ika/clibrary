#pragma once

#include <cmath>

struct Point3D {

	double x, y, z;

	Point3D(double x, double y, double z);

	Point3D multiply(double k);

	Point3D add(const Point3D& p);

	Point3D subtract(const Point3D& p);

	double scalar(const Point3D& p);

	Point3D vector(const Point3D& p);

	double dist(const Point3D& p);

	double norm();

};

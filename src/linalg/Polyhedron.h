#pragma once

#include "mylib.h"
#include "linalg.h"
#include "point3d.h"
#include "plane.h"

#include <vector>

using namespace std;

const double EPS = 1.0e-9;

class Polyhedron {

private:

	uint n;

	vector<Point3D> ps;

	vector<vector<Point3D>> surfaces;

public:

	Polyhedron(vector<Point3D>& ps);

	Point3D& getPoint(int i);

	double getVolume();

	vector<Polyhedron> divide(Plane& plane);

private:

	void check();

	void buildSurface();

	void add(vector<Point3D>& c, Point3D& p);

};

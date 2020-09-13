#pragma once

#include "Polyhedron.h"

class TankVolumeCalculator {

private:

	Point3D underground;

	Polyhedron tank;

public:

	TankVolumeCalculator(Polyhedron const& tank);

	double getVolume(double x, double y, double z, double ax, double ay);

};

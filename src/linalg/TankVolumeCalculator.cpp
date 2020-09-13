#include "TankVolumeCalculator.h"


TankVolumeCalculator::TankVolumeCalculator(Polyhedron const& tank) : underground(0, 0, -1e9), tank(tank) {
}

double TankVolumeCalculator::getVolume(double x, double y, double z, double ax, double ay) {

	Plane fuelSurface(x, y, z, ax, ay);

	double s = fuelSurface.signum(underground);

	Polyhedron fuel = tank.divide(fuelSurface)[s < 0 ? 0 : 1];

	return fuel.getVolume();
}

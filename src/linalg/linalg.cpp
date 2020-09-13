#include "linalg.h"

using namespace std;

double det(double a0, double a1, double a2, double a3, double a4, double a5, double a6, double a7, double a8) {
	double v1 = a0 * (a4 * a8 - a5 * a7);
	double v2 = a3 * (a1 * a8 - a2 * a7);
	double v3 = a6 * (a1 * a5 - a2 * a4);

	return v1 - v2 + v3;
}

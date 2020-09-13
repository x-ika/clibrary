#include "Polyhedron.h"

Polyhedron::Polyhedron(vector<Point3D>& ps) {
	this->ps = ps;
	n = ps.size();
	check();
}

Point3D& Polyhedron::getPoint(int i) {
	return ps[i];
}

double Polyhedron::getVolume() {

	// http://wwwf.imperial.ac.uk/~rn/centroid.pdf

	buildSurface();

	double volume = 0;

	for (vector<Point3D> surface : surfaces) {
		for (uint i = 2; i < surface.size(); i++) {
			Point3D a = surface[0];
			Point3D db = surface[i - 1].subtract(a);
			Point3D dc = surface[i].subtract(a);

			volume += a.scalar(db.vector(dc));
		}
	}

	return volume / 6;
}

vector<Polyhedron> Polyhedron::divide(Plane& plane) {

	buildSurface();

	vector<Point3D> c1;
	vector<Point3D> c2;

	for (Point3D p : ps) {
		double s = plane.signum(p);
		if (s < EPS) {
			c1.push_back(p);
		}
		if (s > -EPS) {
			c2.push_back(p);
		}
	}

	for (vector<Point3D> surface : surfaces) {
		for (uint i = 0; i < surface.size(); i++) {
			Point3D a = surface[i], b = surface[(i + 1) % surface.size()];
			if (plane.signum(a) * plane.signum(b) < 0) {
				Point3D& left = a, & right = b;
				while (left.dist(right) > EPS) {
					Point3D middle = right.add(left).multiply(0.5);
					double s = plane.signum(middle);
					if (s == 0) {
						left = middle;
						break;
					}
					if (s * plane.signum(left) > 0) {
						left = middle;
					}
					else {
						right = middle;
					}
				}

				add(c1, left);
				add(c2, left);
			}
		}
	}

	vector<Polyhedron> ret;
	ret.push_back(Polyhedron(c1));
	ret.push_back(Polyhedron(c2));
	return ret;
}

void Polyhedron::check() {
	for (uint i = 0; i < n; i++) {
		for (uint j = 0; j < i; j++) {
			double d = ps[i].dist(ps[j]);
			if (d < EPS) {
				cerr << "Points " << i << " and " << j << " are too closer" << endl;
			}
			for (uint k = 0; k < n; k++) {
				if (k != i && k != j && abs(ps[i].dist(ps[k]) + ps[j].dist(ps[k]) - d) < EPS) {
					cerr << "Points " << i << " " << j << " " << k << " are too closer" << endl;
				}
			}
		}
	}
}

void Polyhedron::buildSurface() {
	if (surfaces.size() > 0) {
		return;
	}
	for (uint i = 0; i < n; i++) {
		for (uint j = 0; j < i; j++) {
			for (uint k = 0; k < j; k++) {

				bool newPlane = true;

				vector<Point3D> tmp;

				Plane plane(ps[i], ps[j], ps[k]);
				int dir = 0;
				for (uint t = 0; t < n; t++) {
					double s = plane.signum(ps[t]);
					if (abs(s) < EPS) {
						tmp.push_back(ps[t]);
						uint z = tmp.size();
						if (z == 1 && t != k || z == 2 && t != j || z == 3 && t != i) {
							newPlane = false;
							break;
						}
					}
					else {
						if (dir* s < 0) {
							newPlane = false;
							break;
						}
						dir = s > 0 ? 1 : -1;
					}
				}

				if (!newPlane) {
					continue;
				}

				for (uint a = 1; a < tmp.size(); a++) {
					for (uint b = a; b < tmp.size(); b++) {
						bool ok = false;
						Point3D u = tmp[b].subtract(tmp[a - 1]);
						for (uint c = a; c < tmp.size(); c++) {
							Point3D v = tmp[c].subtract(tmp[a - 1]);
							if (dir * plane.signum(tmp[a].add(u.vector(v))) > 0) {
								ok = true;
								break;
							}
						}
						if (ok) continue;
						swap(tmp[a], tmp[b]);
						break;
					}
				}
				surfaces.push_back(tmp);

			}
		}
	}

}

void Polyhedron::add(vector<Point3D>& c, Point3D& p) {
	for (Point3D& t : c) {
		if (p.dist(t) < EPS) {
			return;
		}
	}
	c.push_back(p);
}

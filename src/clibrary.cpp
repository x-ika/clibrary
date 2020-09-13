#include "mylib.h"

#include "mmtester.h"

#include "TankVolumeCalculator.h"

//---------------------------------------------------------------------------------------

const string BASE_DIR = "D:/ika/projects/Competitions/Marathons/TCO2017_MR1C/";

struct TCase {
	int nv;
	VI edges;
};

map<int, TCase> testCases;

void readTests() {
	ifstream test_data;
	test_data.open(BASE_DIR + "testerc/data");
	int seed, nv, ne, e;
	while (test_data >> seed >> nv >> ne) {
		TCase cur;
		cur.nv = nv;
		ne *= 3;
		while (ne-- > 0) {
			test_data >> e;
			cur.edges.push_back(e);
		}
		testCases[seed] = cur;
	}
	test_data.close();
}

double test(int seed) {
	return 0;
}

void testMarathon(int argc, char* argv[]) {

	readTests();

	vector<string> args;
	for (int i = 1; i < argc; i++) {
		args.push_back(argv[i]);
	}

	MMTester::test(*test, BASE_DIR + "testerc/scores", MMTester::ABSOLUTE_SCORING, args, 0);
}

//---------------------------------------------------------------------------------------

Polyhedron getCube() {
	int n = 8;
	vector<Point3D> ps;
	for (int i = 0; i < n; i++) {
		ps.push_back(Point3D(i & 1, i >> 1 & 1, i >> 2 & 1));
	}
	return Polyhedron(ps);
}

Polyhedron getPyramid() {
	vector<Point3D> ps;
	ps.push_back(Point3D(0, 0, 0));
	ps.push_back(Point3D(1, 0, 0));
	ps.push_back(Point3D(0, 1, 0));
	ps.push_back(Point3D(-453, 123, 36));
	return Polyhedron(ps);
}

void testVolumeCalculator() {

	FastTimer::init(9);
	FastTimer::start(0);

	//TankVolumeCalculator calculator(getCube());
	TankVolumeCalculator calculator(getPyramid());

	//cout << calculator.getVolume(0, 0, 0.5, atan(0.5), 0) << endl;
	cout << calculator.getVolume(0, 0, 18, 111233, 0) << endl;

}

//---------------------------------------------------------------------------------------

int main(int argc, char* argv[]) {

	//testMarathon(argc, argv);
	testVolumeCalculator();

	system("pause");
	return 0;
}

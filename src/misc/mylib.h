#pragma once

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>

#include <cmath>
#include <ctime>
#include <cstring>
#include <string>

#include <vector>
#include <set>
#include <bitset>
#include <map>
#include <algorithm>

using namespace std;

typedef unsigned int uint;
typedef long long int64;
typedef unsigned long long uint64;
typedef vector<int> VI;
typedef vector<uint> VUI;
typedef vector<uint64> VUI64;
typedef vector<VI> VII;
typedef vector<VUI> VUII;
typedef vector<VII> VIII;
typedef vector<double> VD;
typedef vector<string> VS;

const double PI = 3.14159265358979323846;

//---------------------------------------------------------------------------------------

namespace Common {

	string dtos(double a);

	template<class T> string toString(const T& a);

	template<class T> void out(T x);

	template<class T> void out(T* a, int n);

}

//---------------------------------------------------------------------------------------

namespace FastTimer {

	inline uint64 time();

	void init(int n);

	void start(int i);

	uint64 elapsed(int i);

	void end(int i);

	uint64 getTotal(int i);

	void print(int i);

	void printAll();

}

//---------------------------------------------------------------------------------------

namespace FastRandom {

	uint fastRandom();

	void init(int64 x);

	int next(int bits);

	inline int nextInt(int n);

	inline double nextDouble();

	template<class T> inline T randomElement(T* a, int n);

}

//---------------------------------------------------------------------------------------

namespace ArrayUtils {

	void init(int n);

	//-----------------------------------------------------------------------------------

	template<class T> T* create(int n);

	template<class T> T** create(int n, int m);

	template <class T> T* clone(T* a, int n);

	template <class T> void fill(T* a, int from, int to, T val);

	template <class T> void arraycopy(T* src, int srcInd, T* dst, int dstInd, int n);

	//-----------------------------------------------------------------------------------

	template <class T> T* createStatic(int n);

	template <class T> void swap(T* a, int i, int j);

	template <class T> void swap(T* a, int i, int j, int n);

	template <class T> int median(T* x, int* p, int a, int b, int c);

	template <class T> void quicksort(T* x, int* p, int a, int b);

	//-----------------------------------------------------------------------------------

	void transformBy(int a, int b, int* p, int* x);

	template <class T> void restore(int a, int b, T** x, int n);

	template <class T> void sortBy(int a, int b, int* p, T** x, int n);

	template <class T> void quicksort(int a, int b, T** x, int n);

}

//---------------------------------------------------------------------------------------

namespace GeometryUtils {

	void init(int n);

	double* initBox(double* box);

	double getBoxSize(double* box);

	double getBoundingBoxSize(double* b1, double* b2);

	double* updateBoundingBox(int n, int* id, double* x, double* y, double* box);

	double* updateBoundingBox(int n, double* x, double* y, double* box);

	double* ifrotate(double ox, double oy, int a, int n, int* id,
		double* x, double* y, double* sin, double* cos, double* box);

	double* ifrotate(double ox, double oy, int a, int n,
		double* x, double* y, double* sin, double* cos, double* box);

	void rotate(double ox, double oy, int a, int n, int* id, double* x, double* y, double* sin, double* cos);

	void rotate(double ox, double oy, int a, int n,
		double* x, double* y, double* sin, double* cos);

	void move(double dx, double dy, int n, int* id, double* x, double* y);

	void move(double dx, double dy, int n, double* x, double* y);

	void scale(double k, int n, int* id, double* x, double* y);

	void scale(double k, int n, double* x, double* y);

	int getConvexHull(int n, double* x, double* y);

}

//---------------------------------------------------------------------------------------

struct Graph {

	int n, m, * id, * z, ** g, ** e;

	Graph();

	Graph(int nVertex);

	void edge(int i, int j, int k);

};

namespace FastGraph {

	void init(int maxN, int maxM);

	Graph buildGraph(int n);

	void saveEdge(int i, int j, int k);

	void toEdges(Graph& graph);

	int dfs(int i, int c, int* f, Graph& g);

	/**
	* Works even if some edges are negative (bridges)
	* @param graph graph
	* @return components
	*/
	vector<Graph> getComponents(Graph& graph);

	void revertEdges(Graph& g, bool all);

	int dfs(int i, int p, int* d, int* f, Graph& g);

	vector<Graph> getStronglyConnectedComponents(Graph& g);

	Graph projectGraph(Graph& g, int* gp, int n, int type);

}

//---------------------------------------------------------------------------------------

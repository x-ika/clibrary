#pragma once

#ifdef _WIN32
#define LOCAL_TESTING
#include <intrin.h>
#endif

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

	string dtos(double a) {
		ostringstream sout;
		sout << setprecision(16) << a;
		return sout.str();
	}

	template<class T> string toString(const T &a) {
		ostringstream sout;
		sout << a;
		return sout.str();
	}

	template<class T> void out(T x) {
		ostringstream s;
		s << x;
		while (s.str().length() < 9) s << ' ';
		cout << s.str();
	}

	template<class T> void out(T *a, int n) {
		for (int i = 0; i < n; i++) {
			out(a[i]);
		}
		cout << endl;
	}

}

//---------------------------------------------------------------------------------------

namespace FastTimer {

	VUI64 startTime, totalTime;

#ifdef _WIN32
	inline uint64 time() {
		return __rdtsc() / 2200000;
	}
#else
	inline uint64 time() {
		uint a, b;
		__asm__ volatile ("rdtsc" : "=A" (a), "=d" (b));
		return ((uint64)b << 32 | (uint64)a) / 2500000;
	}
#endif

	void init(int n) {
		startTime.assign(n, 0);
		totalTime.assign(n, 0);
	}

	void start(int i) {
		startTime[i] = time();
	}

	uint64 elapsed(int i) {
		return time() - startTime[i];
	}

	void end(int i) {
		totalTime[i] += elapsed(i);
	}

	uint64 getTotal(int i) {
		return totalTime[i];
	}

	void print(int i) {
		fprintf(stdout, "Time %d: %.3f\n", i, totalTime[i] / 1e3);
	}

	void printAll() {
		for (uint i = 0; i < startTime.size(); i++) {
			if (totalTime[i]) print(i);
		}
	}

}

//---------------------------------------------------------------------------------------

namespace FastRandom {

	const int64 mp = 0x5DEECE66DLL;
	const int64 ad = 0xB;
	const int64 msk = (1LL << 48) - 1;
	const double MULT = 1.0 / (1LL << 53);

	uint64 seed;

	uint fastRandom() {
		unsigned y = 2463534242;
		return y ^= (y ^= (y ^= y << 13) >> 17) << 5;
	}

	void init(int64 x) {
		seed = (x ^ mp) & msk;
	}

	int next(int bits) {
		seed = (seed * mp + ad) & msk;
		return (int)(seed >> (48 - bits));
	}

	inline int nextInt(int n) {
		return next(31) % n;
	}

	inline double nextDouble() {
		return (((int64)(next(26)) << 27) + next(27)) * MULT;
	}

	template<class T> inline T randomElement(T* a, int n) {
		return a[nextInt(n)];
	}

}

//---------------------------------------------------------------------------------------

namespace ArrayUtils {

	int N;

	// constant array i -> i
	int* ORDER;
	// P[i] element moved to i
	int* P;
	// inverse of a P
	int* PI;

	void init(int n) {
		N = n;
		ORDER = new int[n];
		for (int i = 0; i < n; i++) {
			ORDER[i] = i;
		}
		P = new int[n];
		PI = new int[n];
	}

	//-----------------------------------------------------------------------------------

	template<class T> T* create(int n) {
		T *r = new T[n];
		memset(r, 0, sizeof(T) * n);
		return r;
	}

	template<class T> T** create(int n, int m) {
		T **r = new T*[n];
		for (int i = 0; i < n; i++) {
			r[i] = create<T>(m);
		}
		return r;
	}

	template <class T> T* clone(T *a, int n) {
		T *b = new T[n];
		memcpy(b, a, sizeof(T) * n);
		return b;
	}

	template <class T> void fill(T *a, int from, int to, T val) {
		for (int i = from; i < to; i++) {
			a[i] = val;
		}
	}

	template <class T> void arraycopy(T *src, int srcInd, T *dst, int dstInd, int n) {
		memcpy(dst + dstInd, src + srcInd, sizeof(T) * n);
	}

	//-----------------------------------------------------------------------------------

	template <class T> T* createStatic(int n) {
		static T* tmp = new T[n];
		return tmp;
	}

	template <class T> void swap(T *a, int i, int j) {
		std::swap(a[i], a[j]);
	}

	template <class T> void swap(T *a, int i, int j, int n) {
		for (; n-- > 0; i++, j++) {
			std::swap(a[i], a[j]);
		}
	}

	template <class T> int median(T *x, int *p, int a, int b, int c) {
		return x[p[a]] < x[p[b]] ?
			x[p[b]] < x[p[c]] ? b : x[p[a]] < x[p[c]] ? c : a :
			x[p[b]] > x[p[c]] ? b : x[p[a]] > x[p[c]] ? c : a;
	}

	template <class T> void quicksort(T *x, int *p, int a, int b) {
		int len = b - a;
		if (len < 7) {
			for (int i = a; i < b; i++) {
				for (int j = i; j > a && x[p[j - 1]] > x[p[j]]; j--) {
					swap(p, j, j - 1);
				}
			}
			return;
		}
		int m = a + (len >> 1);
		if (len > 7) {
			int l = a;
			int n = b - 1;
			if (len > 40) {
				int s = len >> 3;
				l = median(x, p, l, l + s, l + 2 * s);
				m = median(x, p, m - s, m, m + s);
				n = median(x, p, n - 2 * s, n - s, n);
			}
			m = median(x, p, l, m, n);
		}
		int s, i = a, j = i, k = b - 1, t = k;
		T &v = x[p[m]];
		while (true) {
			for (; j <= k && x[p[j]] <= v; j++) {
				if (x[p[j]] == v) {
					swap(p, i++, j);
				}
			}
			for (; k >= j && x[p[k]] >= v; k--) {
				if (x[p[k]] == v) {
					swap(p, k, t--);
				}
			}
			if (j > k) {
				break;
			}
			swap(p, j++, k--);
		}
		s = min(i - a, j - i);
		swap(p, a, j - s, s);
		s = min(t - k, b - t - 1);
		swap(p, j, b - s, s);
		if ((s = j - i) > 1) {
			quicksort(x, p, a, a + s);
		}
		if ((s = t - k) > 1) {
			quicksort(x, p, b - s, b);
		}
	}

	//-----------------------------------------------------------------------------------

	void transformBy(int a, int b, int* p, int* x) {
		for (int i = a; i < b; i++) {
			x[i] = p[x[i]];
		}
	}

	template <class T> void restore(int a, int b, T **x, int n) {
		int *p = P;
		T *t = createStatic<T>(N);
		for (int id = 0; id < n; id++) {
			T *y = x[id];
			for (int i = a; i < b; i++) {
				t[p[i]] = y[i];
			}
			arraycopy(t, a, y, a, b - a);
		}
	}

	template <class T> void sortBy(int a, int b, int *p, T **x, int n) {
		T *t = createStatic<T>(N);
		for (int id = 0; id < n; id++) {
			T *y = x[id];
			for (int i = a; i < b; i++) {
				t[i] = y[p[i]];
			}
			arraycopy(t, a, y, a, b - a);
		}
	}

	template <class T> void quicksort(int a, int b, T **x, int n) {
		int *p = P, *pi = PI;
		arraycopy(ORDER, a, p, a, b - a);
		quicksort(x[0], p, a, b);
		for (int i = a; i < b; i++) {
			pi[p[i]] = i;
		}
		sortBy(a, b, p, x, n);
	}

}

//---------------------------------------------------------------------------------------

namespace GeometryUtils {

	double *C;

	// indexes of the points belonged to convex hull
	int *H;

	void init(int n) {
		C = new double[n];
		H = new int[n];
	}

	double* initBox(double *box) {
		box[0] = box[2] = 1e9;
		box[1] = box[3] = -1e9;
		return box;
	}

	double getBoxSize(double *box) {
		return max(box[1] - box[0], box[3] - box[2]);
	}

	double getBoundingBoxSize(double *b1, double *b2) {
		return max( max(b1[1], b2[1]) - min(b1[0], b2[0]), max(b1[3], b2[3]) - min(b1[2], b2[2]) );
	}

	double* updateBoundingBox(int n, int *id, double *x, double *y, double *box) {
		double minx = box[0], maxx = box[1];
		double miny = box[2], maxy = box[3];
		for (int j = 0; j < n; j++) {
			double tx = x[id[j]], ty = y[id[j]];
			minx = min(minx, tx);
			maxx = max(maxx, tx);
			miny = min(miny, ty);
			maxy = max(maxy, ty);
		}
		box[0] = minx;
		box[1] = maxx;
		box[2] = miny;
		box[3] = maxy;
		return box;
	}

	double* updateBoundingBox(int n, double *x, double *y, double *box) {
		return updateBoundingBox(n, ArrayUtils::ORDER, x, y, box);
	}

	double* ifrotate(double ox, double oy, int a, int n, int *id,
		             double *x, double *y, double *sin, double *cos, double *box)
	{
		double minx = box[0], maxx = box[1];
		double miny = box[2], maxy = box[3];
		for (int i = 0; i < n; i++) {
			int j = id[i];
			double tx = x[j] - ox, ty = y[j] - oy;
			double rx = tx * cos[a] - ty * sin[a] + ox;
			double ry = tx * sin[a] + ty * cos[a] + oy;
			minx = min(minx, rx);
			maxx = max(maxx, rx);
			miny = min(miny, ry);
			maxy = max(maxy, ry);
		}
		box[0] = minx;
		box[1] = maxx;
		box[2] = miny;
		box[3] = maxy;
		return box;
	}

	double* ifrotate(double ox, double oy, int a, int n,
		             double *x, double *y, double *sin, double *cos, double *box)
	{
		return ifrotate(ox, oy, a, n, ArrayUtils::ORDER, x, y, sin, cos, box);
	}

	void rotate(double ox, double oy, int a, int n, int *id, double *x, double *y, double *sin, double *cos) {
		for (int i = 0; i < n; i++) {
			int j = id[i];
			double tx = x[j] - ox, ty = y[j] - oy;
			x[j] = tx * cos[a] - ty * sin[a] + ox;
			y[j] = tx * sin[a] + ty * cos[a] + oy;
		}
	}

	void rotate(double ox, double oy, int a, int n,
		        double *x, double *y, double *sin, double *cos)
	{
		rotate(ox, oy, a, n, ArrayUtils::ORDER, x, y, sin, cos);
	}

	void move(double dx, double dy, int n, int *id, double *x, double *y) {
		for (int i = 0; i < n; i++) {
			int j = id[i];
			x[j] += dx;
			y[j] += dy;
		}
	}

	void move(double dx, double dy, int n, double *x, double *y) {
		move(dx, dy, n, ArrayUtils::ORDER, x, y);
	}

	void scale(double k, int n, int *id, double *x, double *y) {
		for (int i = 0; i < n; i++) {
			int j = id[i];
			x[j] *= k;
			y[j] *= k;
		}
	}

	void scale(double k, int n, double *x, double *y) {
		scale(k, n, ArrayUtils::ORDER, x, y);
	}

	int getConvexHull(int n, double *x, double *y) {
		double lx = 1e9, ly = 0;
		for (int i = 0; i < n; i++) {
			if (lx > x[i] || lx == x[i] && ly > y[i]) {
				lx = x[i];
				ly = y[i];
			}
		}
		for (int i = 0; i < n; i++) {
			double dx = x[i] - lx, dy = y[i] - ly;
			double d2 = dx * dx + dy * dy;
			C[i] = d2 < 1e-12 ? -1e9 : atan2(dy, dx);
		}
		ArrayUtils::quicksort(0, n, new double*[3] {C, x, y}, 3);
		int *h = H;
		int sz = 0;
		for (int i = 0; i < n; i++) {
			h[sz++] = i;
			while (sz > 2) {
				int p = h[sz - 3], q = h[sz - 2], r = h[sz - 1];
				double v1x = x[q] - x[p], v1y = y[q] - y[p];
				double v2x = x[r] - x[q], v2y = y[r] - y[q];
				if (v1x * v2y - v1y * v2x > 0) {
					break;
				}
				h[--sz - 1] = r;
			}
		}
		ArrayUtils::restore(0, n, new double*[2] {x, y}, 2);
		ArrayUtils::transformBy(0, sz, ArrayUtils::P, h);
		return sz;
	}

}

//---------------------------------------------------------------------------------------

struct Graph {

	int n, m, *id, *z, **g, **e;

	Graph() {
		n = 0;
		m = 0;
		id = 0;
		z = 0;
		g = 0;
		e = 0;
	}

	Graph(int nVertex) {
		n = nVertex;
		id = ArrayUtils::create<int>(n);
		z = ArrayUtils::create<int>(n);
		g = ArrayUtils::create<int*>(n);
		e = ArrayUtils::create<int*>(n);
		ArrayUtils::arraycopy(ArrayUtils::ORDER, 0, id, 0, n);
	}

	void edge(int i, int j, int k) {
		g[i][z[i]] = j;
		e[i][z[i]] = k;
		z[i]++;
	}

};

namespace FastGraph {

	int *T1, *T2;

	// subgraph index by vertex
	int *GP;
	// subgraph vertext index by vertex
	int *VP;
	// edges of graph in linear format
	int TZ, *TI, *TJ, *TE;

	void init(int maxN, int maxM) {
		T1 = new int[maxN];
		T2 = new int[maxN];
		GP = new int[maxN];
		VP = new int[maxN];
		TI = new int[maxM];
		TJ = new int[maxM];
		TE = new int[maxM];
	}

	Graph buildGraph(int n) {
		Graph graph(n);
		int *z = graph.z;
		int **g = graph.g;
		int **e = graph.e;
		for (int t = 0; t < TZ; t++) {
			z[TI[t]]++;
			z[TJ[t]]++;
		}
		for (int i = 0; i < n; i++) {
			g[i] = new int[z[i]];
			e[i] = new int[z[i]];
			z[i] = 0;
		}
		for (int t = 0; t < TZ; t++) {
			graph.edge(TI[t], TJ[t], TE[t]);
			graph.edge(TJ[t], TI[t], TE[t]);
		}
		graph.m = TZ;
		TZ = 0;
		return graph;
	}

	void saveEdge(int i, int j, int k) {
		TI[TZ] = i;
		TJ[TZ] = j;
		TE[TZ] = k;
		TZ++;
	}

	void toEdges(Graph &graph) {
		int n = graph.n;
		int* gz = graph.z;
		int **g = graph.g;
		int **e = graph.e;
		for (int i = 0; i < n; i++) {
			int z = gz[i], *gg = g[i], *ee = e[i];
			for (int t = 0; t < z; t++) {
				if (i < gg[t]) {
					saveEdge(i, gg[t], ee[t]);
				}
			}
		}
	}

	int dfs(int i, int c, int *f, Graph &g) {
		if (i < 0 || f[i] == c) {
			return 0;
		}
		int z = g.z[i], *gg = g.g[i];
		f[i] = c;
		int s = 1;
		for (int t = 0; t < z; t++) {
			s += dfs(gg[t], c, f, g);
		}
		return s;
	}

	/**
	* Works even if some edges are negative (bridges)
	* @param graph graph
	* @return components
	*/
	vector<Graph> getComponents(Graph &graph) {

		int n = graph.n;
		int *gp = GP, *vp = VP, *ind = T1;

		int c = 0;
		ArrayUtils::fill(gp, 0, n, -1);
		ArrayUtils::fill(vp, 0, n, -1);
		vector<Graph> list;
		for (int p = 0; p < n; p++) {
			if (gp[p] == -1) {
				int gn = dfs(p, c++, gp, graph);
				for (int i = 0, k = 0; i < n; i++) {
					if (gp[i] == c - 1) {
						vp[ind[k] = i] = k++;
					}
				}
				for (int i = 0; i < gn; i++) {
					int fi = ind[i];
					for (int t = 0; t < graph.z[fi]; t++) {
						int fj = graph.g[fi][t];
						if (fj > fi) {
							saveEdge(i, vp[fj], graph.e[fi][t]);
						}
					}
				}
				Graph g(buildGraph(gn));
				for (int i = 0; i < gn; i++) {
					g.id[i] = graph.id[ind[i]];
				}
				list.push_back(g);
			}
		}

		int *gz = new int[c];
		for (int i = 0; i < c; i++) {
			gz[i] = -list[i].n;
		}
		ArrayUtils::quicksort(0, c, new int*[1]{ gz }, 1);
		ArrayUtils::sortBy(0, c, ArrayUtils::P, new Graph*[1]{ &list[0] }, 1);
		ArrayUtils::transformBy(0, n, ArrayUtils::PI, gp);

		return list;
	}

	void revertEdges(Graph &g, bool all) {
		int n = g.n;
		for (int i = 0; i < n; i++) {
			int z = g.z[i], *gg = g.g[i];
			for (int t = 0; t < z; t++) {
				if (all || gg[t] < 0) {
					gg[t] = -gg[t] - 1;
				}
			}
		}
	}

	int dfs(int i, int p, int *d, int *f, Graph &g) {
		int z = g.z[i], *gg = g.g[i];
		f[i] = d[i];
		for (int t = 0; t < z; t++) {
			int j = gg[t];
			if (j != p) {
				if (d[j] > 0) {
					gg[t] = -gg[t] - 1;
					f[i] = min(f[i], d[j]);
				}
				else {
					d[j] = d[i] + 1;
					if (dfs(j, i, d, f, g) <= d[i]) {
						gg[t] = -gg[t] - 1;
						int zz = g.z[j], *ggg = g.g[j];
						for (int tt = 0; tt < zz; tt++) {
							if (ggg[tt] == i) {
								ggg[tt] = -ggg[tt] - 1;
								break;
							}
						}
						f[i] = min(f[i], f[j]);
					}
				}
			}
		}
		return f[i];
	}

	vector<Graph> getStronglyConnectedComponents(Graph &g) {
		int *d = T1, *f = T2;
		ArrayUtils::fill(d, 0, g.n, 0);
		d[0] = 1;
		dfs(0, -1, d, f, g);
		revertEdges(g, true);
		vector<Graph> ret = getComponents(g);
		revertEdges(g, false);
		return ret;
	}

	Graph projectGraph(Graph &g, int *gp, int n, int type) {
		for (int i = 0; i < g.n; i++) {
			int z = g.z[i], *gg = g.g[i], *ee = g.e[i];
			for (int t = 0; t < z; t++) {
				int j = gg[t];
				if (gp[i] < gp[j]) {
					saveEdge(gp[i], gp[j], type == 0 ? ee[t] : type == 1 ? i : j);
				}
			}
		}
		return buildGraph(n);
	}

}

//---------------------------------------------------------------------------------------

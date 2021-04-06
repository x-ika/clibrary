#pragma once

#include "mylib.h"
#include <string>
#include <fstream>
#include <vector>
#include <map>

using namespace std;

template<class K> class DynamicMap {

private:

	bool keepSmaller;

	string fileName;

	map<K, double> data;

public:

	DynamicMap(bool keepSmaller, string fileName) {
		this->keepSmaller = keepSmaller;
		this->fileName = fileName;

		ifstream fileStream;
		fileStream.open(fileName);
		K key;
		double value;
		while (fileStream >> key >> value) {
			data[key] = value;
		}
		fileStream.close();

	}

	void close() {

		ofstream fileStream;
		fileStream.open(fileName);
		for (pair<const int, double> e : data) {
			fileStream << e.first << setprecision(16) << " " << e.second <<  endl;
		}
		fileStream.close();

	}

	double get(K key) {
		return data.count(key) ? data[key] : keepSmaller ? 1e308 : 0;
	}

	bool update(K key, double value) {
		if (!data.count(key) || keepSmaller == (value < data[key])) {
			data[key] = value;
			return true;
		}
		return false;
	}

};

namespace MMTester {

	const int ABSOLUTE_SCORING = 0;
	const int RELATIVE_KEEP_MAX = 1;
	const int RELATIVE_KEEP_MIN = 2;

	string SPACE = "                        ";

	void printAlligned(string s) {
		fprintf(stdout, (s + SPACE.substr(s.length())).c_str());
	}

	double calc(double score, double best, int type) {
		return type == ABSOLUTE_SCORING ? score :
			type == RELATIVE_KEEP_MAX ? best == 0 ? 1 : score / best : score == 0 ? 1 : best / score;
	}

	void test(double (*tester_function)(int), string file, int type, vector<string> &args, int s) {
		vector<int> seeds;
		if (args[s] == "a") {
			for (int i = stoi(args[s + 1]); i < stoi(args[s + 2]); i++) {
				seeds.push_back(i);
			}
		} else {
			for (int i = s; i < args.size(); i++) {
				seeds.push_back(stoi(args[i]));
			}
		}
		int n = 1;
		vector<double> sum;
		sum.resize(n + 1);
		DynamicMap<int> dmap(type == RELATIVE_KEEP_MIN, file);
		for (int seed : seeds) {
			double *sc = new double[n + 1];
			sc[n] = dmap.get(seed);
			for (int i = 0; i < n; i++) {
				sc[i] = tester_function(seed);
				dmap.update(seed, sc[i]);
			}
			fprintf(stdout, "Seed: %d\n                       ", seed);
			for (int i = 0; i <= n; i++) {
				sum[i] += sc[i] = calc(sc[i], dmap.get(seed), type);
				const string em = "";
				printAlligned(Common::dtos(sc[i]));
			}
			fprintf(stdout, "\n");
		}
		dmap.close();
		fprintf(stdout, "                       ");
		for (double score : sum) {
			printAlligned(Common::dtos(score));
		}
		fprintf(stdout, "\n");
	}

}

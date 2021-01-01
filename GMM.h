#pragma once
#include "Array.h"

class GMM {
public:
	int kmeans_iters = 10;
	int em_iters = 10;
	int k, d, n;
	double log_pid;
	Array *load=nullptr;
	Array **sigs, **isigs, *log_dets;
	Array *us;
	Array *log_ws;
	//Array &H;

	GMM(const char *filename);

	~GMM();
	//GMM(int k);

	Array& predict_proba(Array &x, Array &res);

	/*void kmeans(Array xs);

	void fit(Array xs);*/

};



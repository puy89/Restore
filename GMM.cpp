#include "GMM.h"
#include<iostream>

GMM::GMM(const char* filename){
	Array &load=*(new Array(filename));
	k = load.arr[0];
	d = load.arr[1];
	
	int end = 2 + k;
	log_ws = &load.sub(2, end);
	log_ws->ilog();
	int us_size = d * k;
	int begin = end;
	end += us_size;
	Array &us_ravel = load.sub(begin, end);
	us = &us_ravel.reshape(new int[]{k, d}, 2);
	delete &us_ravel;
	begin = end;
	int mat_size = d*d;
	end += mat_size;
	sigs = new Array*[k];
	for (int i = 0; i < k; i++, begin=end, end+=mat_size) {
		Array& sub = load.sub(begin, end);
		sigs[i] = &sub.reshape(new int[] {d, d}, 2);
		delete& sub;
	}
	isigs = new Array * [k];
	for (int i = 0; i < k; i++, begin=end, end += mat_size) {
		Array& sub = load.sub(begin, end);
		isigs[i] = &sub.reshape(new int[] {d, d}, 2);
		delete& sub;
	}
	log_dets = &load.sub(begin, begin+k);;
}

GMM::~GMM(){
	if(load)
	    delete &load;
	delete sigs;
	delete isigs;
}

/*GMM::GMM(int k) {
	this->k = k;
}*/

Array& GMM::predict_proba(Array &x, Array &res) {
	Array us_x(*us);
	us_x.iminus_row(x);
	//std::cout << "after us_x=" << us_x.toString();
	Array tar(d);
	for (int i = 0; i < k; i++){
		tar = 0;
		res.arr[i] = -(us_x[i] ^ (isigs[i]->apply(tar, us_x[i]))) / 2 + log_ws->arr[i]-log_dets->arr[i];
	}
	return res;
}

/*void GMM::kmeans(Array xs) {
	n = xs.shape[0];
	d = xs.shape[1];
	us = zeros(new int[] {k, d}, 2);
	log_pid = d * Math.log(2 * Math.PI) / 2;
	int zs[] = new int[n];
	sigs = new Array[k];
	Random rand = new Random();
	for (int i = 0; i < n; i++)
		zs[i] = rand.nextInt(k);
	Array xs2 = xs.norm2_rows();
	for (int it = 0; it < kmeans_iters; it++) {
		double hist[] = new double[k];
		for (int z = 0; z < k; z++)
			for (int i = 0; i < n; i++) {
				hist[zs[i]] += 1;
				us.row_iplus(zs[i], xs.row(i));
			}
		us.idiv_col(new Array(hist));
		Array us2 = us.norm2_rows();
		//System.out.println(new Array(hist));

		Array dist = xs.mul(us.t());
		dist.iprod(-2);
		dist.iplus_col(xs2);
		dist.iplus_row(us2);
		dist.argmin_row(zs);
		//System.out.println(us);
	}
	double hist[] = new double[k];
	for (int z = 0; z < k; z++) {
		sigs[z] = Array.eye(d);
		for (int i = 0; i < n; i++)
			hist[zs[i]] += 1;
	}
	ws = new Array(hist);
	ws.idiv(n);
	System.out.println(us);

}

void GMM::fit(Array xs) {
	if (us == null)
		kmeans(xs);
	Array H = zeros(new int[] {k, n});
	for (int i = 0; i < em_iters; i++) {

		for (int z = 0; z < k; z++) {
			Array isig = sigs[z].invMat();
			Array Axs = xs.mul(isig);
			Array Axs2 = Axs.prod(xs).sum_rows();
			Array u = us.row(z);
			System.out.println("mul");
			System.out.println(isig.mul(u));
			double Aus2 = u.dot(isig.mul(u));

			Array logs = Axs.mul(u);
			logs.iprod(-2).iplus(Aus2).iplus(Axs2);
			H.row_iplus(z, logs.idiv(-2).exp().iprod(ws.arr[z]));
		}

		H.idiv_row(H.sum_cols());
		Array norm = H.sum_rows();
		ws.set(norm);
		ws.idiv(n);
		us = H.mul(xs).idiv_col(norm);
		for (int z = 0; z < k; z++) {
			Array xsu = xs.copy();
			xsu.iminus_row(us.row(z));
			sigs[z] = xsu.t().mul(xsu.iprod_col(H.row(z))).idiv(norm.arr[z]);
		}

	}
	log_ws = ws.log();
	isigs = new Array[k];
	for (int z = 0; z < k; z++)
		isigs[z] = sigs[z].invMat();

	System.out.println(us);
	for (int z = 0; z < k; z++)
		System.out.println(sigs[z]);
}*/

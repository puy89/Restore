#pragma once
#include <math.h>
#include "Array.h"
#include "GMM.h"
#include "Lsqr.h"


class Mask :public Applyble{
    int *shape=nullptr;
    Array &mask;
public:
    Mask(Array &mask, int h, int w);
    ~Mask();
    Array& apply(Array &tar, Array &vec);

    Array& applyT(Array& tar, Array& vec);

    int* get_shape();

};

class Impaint :public Applyble {
    int* shape = nullptr;
    Array& mask;
    int a1, b1, a2, b2, h, w;
public:
    Impaint(int a1, int b1, int a2, int b2, int h, int w);
    ~Impaint();
    Array& apply(Array& tar, Array& vec);

    Array& applyT(Array& tar, Array& vec);

    int* get_shape();

};


class Restorer :public Applyble {
    int h, w, *shape_src, *shape, m, n;
    Array *b_z=nullptr;
    Array *tmp=nullptr;
    Array *beta_u=nullptr;
    Array **inv_sig_eye=nullptr;
    int n_rows, n_cols;
    GMM &gmm;
    Array *icovs, &b;
    Applyble &trans;
    Lsqr lsqr;
public:
    int outer_iters = 5;
    int patch_size = 8, patch_step = 4, patch_area=64;
    double beta_lsqr=1, beta_patch=0.001;

    Restorer(Applyble& trans, Array& b, GMM& gmm, int* shape_src);

    ~Restorer();

    Array& apply(Array& tar, Array& arr);

    Array& applyT(Array& tar, Array& arr);

    int* get_shape();

    void solve_patchs(Array &tar);

    Array& restore(Array& tar);
};

#include "Restorer.h"
#include<iostream>


Restorer::Restorer(Applyble& trans, Array& b, GMM& gmm, int* shape_src) :trans(trans), b(b), gmm(gmm) {
    h = b.shape[0];
    w = b.shape[1];
    this->shape_src = shape_src;
}

Array& Restorer::apply(Array& tar, Array& arr) {    
    Array& img = arr.reshape(shape_src, 2);
    Array& subtar = tar.sub(0, b.size);
    trans.apply(subtar, arr);
    for (int k = 0, i = 0, z_i = b.size, z_end = b.size + patch_area; i < n_rows; k += patch_step, i++)
        for (int l = 0, j = 0; j < n_cols; l += patch_step, j++, z_i = z_end, z_end += patch_area) {
            Array &p=img.crop(k, k + patch_size, l, l + patch_size);
            p *= beta_lsqr;
            tar.set_iplus(z_i, z_end, p);
            delete &p;
        }
    delete& img;
    delete& subtar;
    return tar;
}

Array& Restorer::applyT(Array& tar, Array& arr) {
    Array& subarr = arr.sub(0, b.size);
    Array& img = tar.reshape(shape_src, 2);
    trans.applyT(tar, subarr);
    for (int k = 0, i = 0, z_i = b.size, z_end = b.size + patch_area; i < n_rows; k += patch_step, i++)
        for (int l = 0, j = 0; j < n_cols; l += patch_step, j++, z_i = z_end, z_end += patch_area) {
            Array& p = arr.sub(z_i, z_end);
            p *= beta_lsqr;
            img.crop_iplus(k, k + patch_size, l, l + patch_size, p);
            delete &p;
        }
    delete& img;
    delete& subarr;
    return tar;
}

int* Restorer::get_shape() {
    return shape;
}


Restorer::~Restorer() {
    if (shape != nullptr)
        delete shape;
    if (b_z != nullptr)
        delete b_z;
    if (tmp != nullptr)
        delete tmp;
    if (beta_u != nullptr)
        delete beta_u;
    if(inv_sig_eye != nullptr){
        for (int i = 0; i < gmm.k; i++)
            delete inv_sig_eye[i];
        delete[] inv_sig_eye;
    }
}

void Restorer::solve_patchs(Array &tar) {
    Array& img = tar.reshape(shape_src, 2);
    std::cout << "patchs\n";
    //h, w = shape_src[:2]
    Array res(gmm.k);
    for (int k = 0, i = 0, z_i = b.size, z_end = b.size + patch_area; i<n_rows; k += patch_step, i++)
        for (int l = 0, j = 0; j<n_cols; l += patch_step, j++, z_i = z_end, z_end += patch_area) {
            Array& p = img.crop(k, k + patch_size, l, l + patch_size).ravel();
            int arg = gmm.predict_proba(p, res).argmax();
            Array& u = (*beta_u)[arg];
            Array& z = b_z->sub(z_i, z_end);
            *tmp = 0;
            z = 0;
            gmm.sigs[arg]->apply(*tmp, p) += u;
            inv_sig_eye[arg]->apply(z, *tmp);

            z *= beta_lsqr;
            delete& u;
            delete& p;
            delete& z;
        }
    delete& img;
}

Array& Restorer::restore(Array& tar) {
    n_rows = int(ceil((h - patch_size + 1) /(double)patch_step));
    n_cols = int(ceil((w - patch_size + 1) /(double)patch_step));
    patch_area = patch_size * patch_size;
    tmp = new Array(patch_area);
    b_z = new Array(b.size + n_rows * n_cols * patch_area);
    m = b_z->size;
    n = tar.size;
    shape = new int[2]{m, n};
    std::cout << "n_rows = " << n_rows << "n_cols = " << n_cols << "\n";
    b_z->set(0, b.size, b);
    beta_u = new Array(*gmm.us);
    *beta_u *= 1/beta_patch;
    inv_sig_eye = new Array*[gmm.k];
    Array &e=eye(gmm.d);
    e *= 1/beta_patch;
    for (int i = 0; i < gmm.k; i++){
        Array M(e);
        M += *gmm.sigs[i];
        inv_sig_eye[i] = &(M.inv());
    }
    delete& e;
    for (int i = 0; i < outer_iters; i++){
        std::cout << "iter " << i << "\n";
        solve_patchs(tar);
        lsqr.solve(*this, *b_z, tar);
    }
    std::cout << "crop=" << tar.crop(0, 8, 0, 8).toString();
    //delete& bravel;
    return tar;
}


Mask::Mask(Array& mask, int h, int w): mask(mask){
    int hw = h * w;
    shape = new int[2]{hw, hw};
}

Mask::~Mask() {
    delete shape;
}

Array& Mask::apply(Array& tar, Array& vec) {
    Array masked_vec(vec);
    masked_vec *= mask;
    tar += masked_vec;
    return tar;
}

Array& Mask::applyT(Array& tar, Array& vec) {
    Array masked_vec(vec);
    masked_vec *= mask;
    tar += masked_vec;
    return tar;
}

int* Mask::get_shape() {
    return shape;
}

Impaint::Impaint(int a1, int b1, int a2, int b2, int h, int w) : mask(mask) {
    int hw = h * w;
    this->a1 = a1;
    this->b1 = b1;
    this->a2 = a2;
    this->b2 = b2;
    this->h = h;
    this->w = w;
    shape = new int[2]{ hw, hw };
}

Impaint::~Impaint() {
    delete shape;
}

Array& Impaint::apply(Array& tar, Array& vec) {
    Array black_vec(vec);
    Array& black_img = black_vec.reshape(new int[]{h, w}, 2);
    black_img.crop_set(a1, b1, a2, b2, 0);
    tar += black_img;
    delete &black_img;
    return tar;
}

Array& Impaint::applyT(Array& tar, Array& vec) {
    Array black_vec(vec);
    Array& black_img = black_vec.reshape(new int[] {h, w}, 2);
    black_img.crop_set(a1, b1, a2, b2, 0);
    tar += black_img;
    delete& black_img;
    return tar;
}

int* Impaint::get_shape() {
    return shape;
}
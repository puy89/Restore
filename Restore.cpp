// Restore.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include "Restorer.h"

int main()
{
    Array arr_ravel("lena");
    int h = 200, w = 300;
    Array& arr = arr_ravel.reshape(new int[]{h, w}, 2);
    GMM gmm("gmm");
    int shape[] = { h, w };
    Impaint trans(100, 120, 140, 160, h, w);
    arr.crop_set(100, 120, 140, 160, 100.);
    Restorer rest(trans, arr, gmm, shape);
    Array arr_rec(arr_ravel);
    
    rest.restore(arr_rec);
    arr_rec.save("rec_impaint_lena");
}

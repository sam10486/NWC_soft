#include <iostream>
#include <vector>
#include <algorithm>
#include <iterator>
#include "NWC_math.h"
#include "BitOperate.h"
#include "math.h"
#include "NWC.h"

using namespace std;

vector<long long> NWC_forward(vector<long long> a, long long degree, long long modular){
    long long t = degree;
    vector<long long> phi_array_rev = phi_array(degree, modular);

    for(int m = 1; m < degree; m = m << 1){
        t = t / 2;
        for(int i = 0; i < m; i++){
            int j1 = 2*i*t;
            int j2 = j1 + t - 1;
            long long S = phi_array_rev.at(m+i);
            for(int j = j1; j <= j2; j++){
                long long U = a.at(j);
                long long V = MulMod(a.at(j+t), S, modular);
                a.at(j) = AddMod(U, V, modular);
                a.at(j+t) = SubMod(U, V, modular);
            }
        }
    }
    return a;
}

ostream_iterator<long long> os(cout, " | ");
vector<long long> NWC_backward(vector<long long> a, long long degree, long long modular){
    long long t = 1;
    vector<long long> phi_array_inv_rev = phi_array_inv(degree, modular);
    for(int m=degree; m>1; m = m >> 1){
        int j1 = 0;
        int h = m / 2;
        for(int i = 0; i < h; i++){
            int j2 = j1 + t - 1;
            long long S = phi_array_inv_rev.at(h+i);
            for(int j = j1; j <= j2; j++){
                long long U = a.at(j);
                long long V = a.at(j+t);
                a.at(j) = AddMod(U, V, modular);
                a.at(j+t) = MulMod(SubMod(U,V,modular), S, modular);
            }
            j1 = j1 + 2*t;
        }
        t = 2*t;
    }
    long long degree_inv = InvMod(degree, modular);
    for(int j = 0; j < degree; j++){
        a.at(j) = MulMod(a.at(j), degree_inv, modular);
    }
    return a;
}


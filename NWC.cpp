#include <iostream>
#include <vector>
#include <algorithm>
#include <iterator>
#include "NWC_math.h"
#include "BitOperate.h"
#include "math.h"
#include "NWC.h"

using namespace std;

long long NWC_forward(vector<long long> a, long long degree, long long modular){
    long long t = degree;
    vector<long long> phi_array_rev = phi_array(degree, modular);

    for(int m = 1; m < degree; m = m << 1){
        cout << "m = " << m << endl;
        t = t / 2;
        for(int i = 0; i < m; i++){
            int j1 = 2*i*t;
            int j2 = j1 + t - 1;
            long long S = phi_array_rev.at(m+i);
        }
    }
}


#include <iostream>
#include <vector>
#include <algorithm>
#include <iterator>
#include <random>
#include "NWC_math.h"
#include "BitOperate.h"
#include "math.h"
#include "NWC.h"
#include <NTL/ZZ.h>

using namespace std;
using namespace NTL;

int main(){
    ZZ modular = power2_ZZ(64) - power2_ZZ(32) + (ZZ)1;
    int n = 65536;
    BitOperate rev;
    ostream_iterator<long long> os(cout, " , ");
    
    /*std::default_random_engine generator( time(NULL) );
    std::uniform_real_distribution<float> unif(0, modular);*/

    ZZ phi = find_phi(n, modular);
    cout << "phi = " << phi << endl;

    ZZ n_th_phi = PowerMod(phi, 2*n, modular);
    cout << "n_th_phi = " << n_th_phi << endl;

    return 0;
}
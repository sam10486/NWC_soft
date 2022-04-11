#include <iostream>
#include <vector>
#include <algorithm>
#include <iterator>
#include "NWC_math.h"
#include "BitOperate.h"
#include "math.h"
#include "NWC.h"

using namespace std;

int main(){
    int in_1 = 2;
    int in_2 = 5;
    int modular = 97;
    long long ans =0;
    int n = 8;
    BitOperate rev;
    ostream_iterator<long long> os(cout, " | ");
    

    vector<long long> phi_array_rev = phi_array(n, modular);
    copy(phi_array_rev.begin(), phi_array_rev.end(), os);
    cout << endl;
    NWC_forward(phi_array_rev, n, modular);
    /*cout << "phi = " << find_phi(16,97) << endl;
    cout << "prou = " << find_prou(16,97) << endl;*/
    /*for (int i = 0; i < n; i++){
        ans = rev.BitReserve(i, log2(n));
        cout << "i = " << i << " reverse index = " << ans << endl;
    }*/
    return 0;
}
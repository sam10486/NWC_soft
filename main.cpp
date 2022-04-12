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
    int modular = 40961;
    long long ans =0;
    int n = 4096;
    BitOperate rev;
    ostream_iterator<long long> os(cout, " | ");
    

    /*vector<long long> phi_array_rev = phi_array(n, modular);
    copy(phi_array_rev.begin(), phi_array_rev.end(), os);
    cout << endl;
    vector<long long> phi_array_inv_rev = phi_array_inv(n, modular);
    copy(phi_array_inv_rev.begin(), phi_array_inv_rev.end(), os);
    cout << endl;*/
    
    cout << "-------------------------input data--------------------" << endl;
    vector<long long> a(n);
    for(int i=0; i<n; i++)
        a.at(i) = i;
    copy(a.begin(), a.end(), os);
    cout << endl;
   // cout << "-------------------DFT ans -----------------------" << endl;
    vector<long long> NTT_forward_vec ;
    NTT_forward_vec = NWC_forward(a, n, modular);
    /*copy(NTT_forward_vec.begin(), NTT_forward_vec.end(), os);
    cout << endl;*/
    cout << "---------------IDFT ans-------------------" << endl;
    vector<long long> NTT_backward_vec ;
    long long rev_index;

    NTT_backward_vec = NWC_backward(NTT_forward_vec, n, modular);
    copy(NTT_backward_vec.begin(), NTT_backward_vec.end(), os);
    cout << endl;

    /*cout << "phi = " << find_phi(16,97) << endl;
    cout << "prou = " << find_prou(16,97) << endl;*/
    /*for (int i = 0; i < n; i++){
        ans = rev.BitReserve(i, log2(n));
        cout << "i = " << i << " reverse index = " << ans << endl;
    }*/

    return 0;
}


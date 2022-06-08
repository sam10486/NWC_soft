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
    long long radix_4 = 4;
    long long m = 16;
    long long modular = 257;

    long long input_vector[m];
    long long stage0_data_in[4][radix_4];
    long long stage0_data_fft[4][radix_4];
    long long stage0_data_out_tmp[4][radix_4];

    long long stage1_data_in[4][radix_4];
    long long stage1_data_fft[4][radix_4];

    long long phi_16 = find_phi(m, modular);

    BitOperate BR;

    


}
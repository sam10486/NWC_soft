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

    long long n = 16;
    long long data_in[n];
    long long data_out[n];
    long long radix = 2;

    for(int i=0; i<n; i++){
        data_in[i] = i;
    }

    mem_AE(data_out, data_in, n, radix);
}
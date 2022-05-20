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
#include <fstream>

using namespace std;
using namespace NTL;

int main(){
    long long n = 64;
    ZZ modular = find_prime((ZZ)1,7);
    cout << "modular = " << modular << endl;
    long long bit_width = 9;
    long long alpha = bit_width + 1;

    ZZ value = precompute_value(modular, bit_width, alpha);
    cout << "value = " << value << endl;
}
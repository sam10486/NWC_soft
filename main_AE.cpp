#include <iostream>
#include <vector>
#include <algorithm>
#include <iterator>
#include <random>
#include <bitset>
#include "NWC_math.h"
#include "BitOperate.h"
#include "math.h"
#include "NWC.h"

using namespace std;

int main(){
    long long degree_N = 16;
    long long radix_r = 2;
    long long num_stage_p;
    long long bit_width_s;
    long long relocation_group_g;


    num_stage_p = log(degree_N)/log(radix_r);
    relocation_group_g = degree_N/pow(radix_r,2);
    bit_width_s = log(radix_r)/log(2);

    for(long long t=0; t<num_stage_p; t++){
        for(long long i=0; i<relocation_group_g; i++){
            for(long long j=0; j<radix_r; j++){
                long long BC = j*g + 
            }
        }
    }
}
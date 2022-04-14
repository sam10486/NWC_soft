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



long long BitOperate::BitReserve(long long DataToReverse, long long BitLength){
    long long result = 0;
    for(long long i = 0; i < BitLength; i++){
        if((DataToReverse >> i) & 1){
            result |= 1 << (BitLength - 1 -i);
        }
    }
    return result;
}

vector<long long> BitOperate::DecToBin(long long data, long long bit_width){
    vector<long long> BinVec(bit_width);
    for(long long int j=0; j<bit_width; j++){
        BinVec.at(j) = (data >> j) & 1;
    }
    return BinVec;
}

void BitOperate::IntToVec(long long integer, long long N, long long r, vector<long long> &bit_array){
    long long bit_width = (long long)ceil(log2(N));
    bit_array.resize(bit_width);
    for(long long j=0; j < bit_width; j++){
        bit_array[j] = (integer >> j) & 1;
        //cout << bit_array[j] << endl;
    }
    
}

long long BitOperate::VecToInt(vector<long long> &bit_array, long long N, long long r){
    long long bit_width = (long long)ceil(log2(N));
    long long integer = 0;
    for(long long j=0; j < bit_width; j++){
        integer += bit_array[j] << j;
        //cout << "bit_array[" << j << "] = " << bit_array[j] << endl;
    }
    return integer;
}

long long BitOperate::RR(long long BC, long long shift_bit, long long N, long long r){
    long long RR_out = 0;
    long long bit_width = (long long)ceil(log2(N/r));
    vector<long long> bit_array(bit_width);
    BitOperate DecToBin;
    bit_array = DecToBin.DecToBin(BC, bit_width);

    rotate(bit_array.begin(), bit_array.begin()+shift_bit, bit_array.end());
    for(long long j=0; j < bit_width; j++){
        RR_out += bit_array[j] << j;
        //cout << RR_out << endl ;
    }
    return RR_out;
}
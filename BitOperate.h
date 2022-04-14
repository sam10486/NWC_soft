#ifndef _BitOperate_
#define _BitOperate_

#include <iostream>
#include <vector>
using namespace std; 

class BitOperate{
public:
    long long BitReserve(long long DataToReverse, long long BitLength);
    vector<long long> DecToBin(long long data, long long bit_width);
    void IntToVec(long long integer, long long N, long long r, vector<long long> &bit_array);
    long long VecToInt(vector<long long> &bit_array, long long N, long long r);
    long long RR(long long BC, long long shift_bit, long long N, long long r);
};

#endif
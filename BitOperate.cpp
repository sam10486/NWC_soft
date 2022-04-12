#include <iostream>
#include <bitset>
#include "math.h"
#include "BitOperate.h"
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

/*
    long long n = 8;
    BitOperate rev;
    for (long long i = 0; i < n; i++){
        ans = rev.BitReserve(i, log2(n));
        cout << "i = " << i << " reverse index = " << ans << endl;
    }
*/
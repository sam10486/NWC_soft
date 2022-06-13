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
    long long radix_mem = 2;
    long long modular = 12289;
    long long phi = find_phi(n, modular);
    BitOperate BR;

    int BN = 2;
    int MA = (n/BN)/radix_mem;
    vector<vector<vector<long long> > > memory;

    memory.resize(BN);
    for(int i=0; i<BN; i++){
        memory[i].resize(MA);
        for(int j=0; j<MA; j++){    
            memory[i][j].resize(radix_mem);
        }
    }

    /*for(int i=0; i<BN; i++){
        for(int j=0; j<MA; j++){
            for(int k=0; k<radix_mem; k++){
                cout << "memory[" << i << "][" << j << "][" << k << "] = " << memory[i][j][k] << endl;
            }
        }
    }*/

    for(int i=0; i<n; i++){
        data_in[i] = i;
    }

    mem_init(memory, data_in, n, radix_mem);
    cout << "--------------" << endl;
    for(int i=0; i<BN; i++){
        for(int j=0; j<MA; j++){
            for(int k=0; k<radix_mem; k++){
                cout << "memory[" << i << "][" << j << "][" << k << "] = " << memory[i][j][k] << endl;
            }
        }
    }
    cout << "-----------------" << endl;
    //mem_AE(memory, data_in, n, radix_mem);


    
    /*mixed_radix_NWC(data_out, data_in, n, 2, 2, phi, modular);

    for(int i=0; i<n; i++){
        int idx = BR.BitReserve(i, log2(n));
        cout << data_out[idx] << " ";
    }*/
    cout << endl;
}
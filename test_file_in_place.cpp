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
    long long n = 32768;
    long long input_idx[n];
    long long data_in[n];
    long long data_out[n];
    long long radix = 16;
    long long modular = 65537;
    long long phi = find_phi(n, modular);
    cout << "phi = " << phi << endl;
    BitOperate BR;

    int BN = radix;
    int MA = (n/BN);

    vector<vector<long long> > memory;

    ofstream ofs;
    ofs.open("/home/ldap-users/siang/Desktop/NWC_software/check_in_place/memory_ans.txt");

    memory.resize(BN);
    for(int i=0; i<BN; i++){
        memory[i].resize(MA);
    }

    for(int i=0; i<BN; i++){
        for(int j=0; j<MA; j++){
            cout << "memory[" << i << "][" << j << "] = " << memory[i][j] << endl;
        }
    }

    for(int i=0; i<n; i++){
        input_idx[i] = i;
        data_in[i] = i;
        cout << "data_in[" << i << "] = " << data_in[i] << endl;
        mem_init_in_place(memory, data_in[i], input_idx[i], n, radix, modular);
        cout << "------------" << endl;
    }

    for(int i=0; i<BN; i++){
        cout << "------------" << endl;
        for(int j=0; j<MA; j++){
            cout << "memory[" << i << "][" << j << "] = " << memory[i][j] << endl;
        }
    }
    cout << "-----------------" << endl;

    long long log_radix = log2(radix);
    mixed_radix_NWC_in_place(data_out, data_in, n, log_radix, 3, phi, modular, memory);

    cout << "-----------------" << endl;
    for(int i=0; i<BN; i++){
        cout << "------------" << endl;
        for(int j=0; j<MA; j++){
            cout << "memory[" << i << "][" << j << "] = " << memory[i][j] << endl;
        }
    }


    if(!ofs.is_open()){
        cout << "failed to open file.\n" << endl;
    }else {
        for(int i=0; i<BN; i++){
            for(int j=0; j<MA; j++){
                 ofs << memory[i][j] << endl;
            }
        }
        ofs.close();
    }

}
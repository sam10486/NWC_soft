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
    long long radix_2 = 2;
    long long radix_16 = 16;
    long long m = 512;
    long long modular = 12289;

    long long input_vector[m];
    long long input_stage0[2][256];

    long long stage0_data_in[2][16][16];
    long long stage0_data_out[2][16][16];
    long long stage0_data_out_tmp[2][16][16];

    long long stage1_data_in_radix16[2][16][radix_16];
    long long stage1_data_out_radix16[2][16][radix_16];
    long long stage1_data_out_tmp[2][256];
   

    long long phi_512 = find_phi(m, modular);
    long long phi_256 = ExpMod(phi_512, 2, modular);
    long long phi_16 = ExpMod(phi_512, 32, modular);
    long long phi_2 = ExpMod(phi_512, 256, modular);

    long long prou_512 = ExpMod(phi_512, 2, modular);
    long long prou_256 = ExpMod(phi_512, 2, modular);
    long long prou_16 = ExpMod(prou_512, 32, modular);
    long long prou_2 = ExpMod(prou_16, 256, modular);

    for(int i=0; i<m; i++){
        input_vector[i] = i;
        cout << "input_vector[" << i  << "] "<< " = " << input_vector[i] << endl;
    }

    cout << "----------stage0_data_in-----------" << endl;
    for(int i=0; i<2; i++){
        for(int j=0; j<256; j++){
            input_stage0[i][j] = input_vector[i+2*j];
            cout << "input_stage0[" << i << "][" << j << "]" << " = " << input_stage0[i][j] << endl;
        }
    }

    cout << "----------stage0_radix16 input------------" << endl;
    for(int k=0; k<2; k++){
        for(int i=0; i<16; i++){
            for(int j=0; j<16; j++){
                stage0_data_in[k][i][j] = input_stage0[k][i+j*16];
                cout << "stage0_data_in[" << k << "][" << i << "][" << j << "]" << " = " << stage0_data_in[k][i][j] << endl;
            }
        }
    }

    cout << "----------stage0_radix 16 fft------------" << endl;
    for(int k=0; k<2; k++){
        for(int i=0; i<16; i++){
            for(int j=0; j<16; j++){
                NWC_forward_DIT(stage0_data_out[k][i],stage0_data_in[k][i], radix_16, phi_16, modular);
                cout << "stage0_data_out[" << k << "][" << i << "][" << j << "]" << " = " << stage0_data_out[k][i][j] << endl;
            }
        }
    }

    cout << "----------stage0_radix16 fft twiidle-----------" << endl;
    for(int k=0; k<2; k++){
        for(int i=0; i<16; i++){
            for(int j=0; j<16; j++){
                BitOperate BU;
                long long index_rev = BU.BitReserve(j, 4);
                long long exp_prou256 = ExpMod(prou_256, i*index_rev, modular);
                long long exp_phi256 = ExpMod(phi_256, index_rev, modular);
                long long twiddle = MulMod(exp_phi256, exp_prou256, modular);
                stage0_data_out_tmp[k][i][j] = MulMod(stage0_data_out[k][i][j] , twiddle, modular);
                cout << "stage0_data_out_tmp[" << k << "][" << i << "][" << j << "]" << " = " << stage0_data_out_tmp[k][i][j] << endl;
            }
        }
    }
    
    cout << "----------stage1_data_in-----------" << endl;
    for(int k=0; k<2; k++){
        for(int i=0; i<16; i++){
            for(int j=0; j<16; j++){
                stage1_data_in_radix16[k][i][j] = stage0_data_out_tmp[k][j][i];
                cout << "stage1_data_in_radix16[" << k << "][" << i << "][" << j << "]" << " = " << stage1_data_in_radix16[k][i][j] << endl;
            }
        }
    }

    cout << "----------stage1_radix16 fft-----------" << endl;
    for(int k=0; k<2; k++){
        for(int i=0; i<16; i++){
            for(int j=0; j<16; j++){
                NWC_forward_DIT(stage1_data_out_radix16[k][i], stage1_data_in_radix16[k][i], radix_16, phi_16, modular);
                cout << "stage1_data_out_radix16[" << k << "][" << i << "][" << j << "]" << " = " << stage1_data_out_radix16[k][i][j] << endl;
            }
        }
    }



}
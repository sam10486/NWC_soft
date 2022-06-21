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
    ifstream  ifs1, ifs2;
    ifs1.open("/home/ldap-users/siang/Desktop/NWC_software/check_AE/algo_asn.txt");
    ifs2.open("/home/ldap-users/siang/Desktop/NWC_software/check_AE/memory_ans.txt");
    
    vector<long long > algo_asn_array;
    vector<long long > memory_ans_array;
    

    if(!ifs1.is_open() || !ifs2.is_open()){
        cout << "failed to open file.\n" << endl;
    }else {
        long long algo_asn;
        long long memory_ans;
        while(ifs1 >> algo_asn){
            cout << algo_asn << endl;
            algo_asn_array.push_back(algo_asn);
        }
        cout << "-----------------" << endl;
        while(ifs2 >> memory_ans){
            cout << memory_ans << endl;
            memory_ans_array.push_back(memory_ans);
        }
    }
    ifs1.close();
    ifs2.close();
    int correct=0;
    int j;

    for(int i=0; i<memory_ans_array.size(); i++){
        long long mem_ans = memory_ans_array[i];
        //cout << "mem_ans = " << mem_ans << endl;
        for(j=0; j<algo_asn_array.size(); j++){
            long long algo_ans = algo_asn_array[j];
            if(mem_ans == algo_ans){
                correct++;
                break;
            }
        }
        if(correct != 0){
            cout << "correct = " << correct << endl;
            cout << "memory_ans[" << i << "] = " << memory_ans_array[i] << " === algo_asn_array[" << j << "] = " << algo_asn_array[j] << ", PASS!" << endl;
            correct = 0;
        }else{
            cout << "memory_ans[" << i << "] False XX " << endl;
            correct = 0;
        }
    }
}

#include <iostream>
#include "NWC_math.h"
#include "math.h"
#include "assert.h"
#include "BitOperate.h"
#include <vector>


using namespace std;


// q is prime
long long AddMod(long long in_1, long long in_2, long long q) {
    long long tmp;
    tmp = in_1 + in_2;
    tmp = tmp % q;
    return tmp;
}

long long SubMod(long long in_1, long long in_2, long long q) {
    long long tmp;
    tmp = in_1 - in_2;
    if (tmp < 0){
        tmp = tmp + q;
    }
    tmp = tmp % q;
    return tmp;
}

long long XGCD(long long a, long long q){
    long long x0 = 1;
    long long y0 = 0;
    long long x1 = 0;
    long long y1 = 1;
    long long u = a;
    long long v = q;
    while (u != 0){
        long long quotient = v / u ;
        long long r = v % u;
        v = u;
        u = r;
        long long x2 = x0 - (quotient * x1);
        long long y2 = y0 - (quotient * y1);
        x0 = x1;
        x1 = x2;
        y0 = y1;
        y1 = y2;
    }
    if(y0 < 0)
        y0 = y0 + q;
    
    return y0;
}

long long InvMod(long long a, long long q){
    long long tmp;
    tmp = XGCD(a, q);
    return tmp;
}

long long DivMod(long long a,long long b,long long q){
    long long Ib  = InvMod(b,q);  //inverse b in finite filed P  {GF(P)}
    long long tmp = a * Ib;
    tmp = tmp % q;
    return tmp;
}

long long MulMod(long long a,long long b,long long q){
    if(a < 0 || b < 0 ){
        cout << "a or b is negative number, return 0" << endl;
        return 0;
    }else{
        long long tmp = a * b;
        tmp = tmp % q;
        return tmp;
    }
}

long long ExpMod(long long a, long long exp, long long q){
    long long tmp = a;
    for(long long i=1; i < exp; i++){
        tmp = tmp * a;
        tmp = tmp % q;
    }
    if(exp == 0)
        tmp = pow(a, 0);

    return tmp; 
}

// p = 1 mod 2n
long long find_phi(long long n, long long modular){
    //output = primitive root of unity
    long long i,j;
    long long phi;
    long long phi_temp;
    long long double_n = n*2;
    for(j=2;j<modular;j++){
        if((modular % j) == 0){
            printf("---------------modular is no prime----------------\n");
            return 0;
        }
    }

    for(phi=2;phi<modular;phi++){
        phi_temp = 1;
        for(i=1;i<double_n;i++){
            phi_temp *= phi;
            phi_temp %= modular;
            if(phi_temp == 1){
                break;
            }
        }
        
        phi_temp *= phi;
        phi_temp %= modular;
        if(phi_temp == 1){
            break;
        }
    }
    
    if(phi == modular){
        return 0;
    }
    else{
        return phi;
    }
}

bool isPrime(long long n){
    if(n==1)
        return 0;
    long long i=2;
    for(; i*i<=n; i++){
        if(n%i==0){
            return 0;
        }
    }
    return 1;
}

long long find_prou(long long n, long long modular){
    long long phi = find_phi(n, modular);
    assert(phi != 0);
    long long prou ;
    prou = MulMod(phi, phi, modular);
    return prou;
}

//****************index reverse phi array*****************
vector<long long> phi_array(long long n, long long modular){
    vector<long long> phi_array(n);
    long long phi = find_phi(n, modular);
    cout << "phi = " << phi << endl;
    BitOperate rev;
    long long rev_index;
    /*cout << "n = " << n << endl;
    for(long longi=0; i<n; i++)
            cout << "phi_array = " << phi_array[i] << endl;*/
    for(long long i=0; i<n; i++){
        rev_index = rev.BitReserve(i, log2(n));
        //cout << "rev_index = " << rev_index << endl;     
        phi_array.at(rev_index) = ExpMod(phi, i, modular);
        
    }
    /*for(long longi=0; i<n;i++)
        cout << "phi_array = " << phi_array[i] << endl;*/
    return phi_array;
}
//**************index reverse phi inv array*********************
vector<long long> phi_array_inv(long long n, long long modular){
    vector<long long> phi_array_tmp = phi_array(n, modular);
    
    vector<long long> phi_array_inv(n);
    for(long long i=0; i<n; i++){
        phi_array_inv.at(i) =  InvMod(phi_array_tmp.at(i), modular);
        //cout << "phi_array_inv.at " << i << " = "<< phi_array_inv.at(i) << endl;
    }
        
    return phi_array_inv;
}
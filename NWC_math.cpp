
#include <iostream>
#include "NWC_math.h"
#include "math.h"
#include "assert.h"
#include "BitOperate.h"
#include <vector>
#include <NTL/ZZ.h>


using namespace std;
using namespace NTL;


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
    //cout << "phi = " << phi << endl;
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

long long barrett_reduction(long long a, long long b, long long modular){
    long long n = ceil(log2(modular));
    cout << "n = " << n << endl;
    long long mu = floor(pow(2,2*n)/modular);
    cout << "mu = " << mu << endl;
    long long z = a * b;
    cout << "z = " << z << endl;
    long long tmp1 = floor(z >> (long long)log2(pow(2,n)));
    long long tmp2 = mu >> (long long)log2(pow(2,n));
    cout << "tmp1 = " << tmp1 << endl;
    cout << "tmp2 = " << tmp2 << endl;
    long long r = z - (tmp1 * tmp2) * modular;;
    cout << "r = " << r << endl;
    while(r >= modular){
        r = r - modular;
    }
    return r;
}

long long modular_mul(long long x, long long y, long long modular){
    long long modular_BW = ceil(log2(modular)); //precompute
    long long tmp_y = floor((y*pow(2,modular_BW))/modular); //precompute
    long long z = (x * y);
    if(z > pow(2, modular_BW))
        z = z - pow(2, modular_BW);
    long long t = floor((x*tmp_y)>>modular_BW);
    long long z_e = t *  modular;
    if(z_e > pow(2, modular_BW))
        z_e = z_e - pow(2, modular_BW);
    z = z - z_e;
    cout << "z = " << z << endl;
    if(z >= modular)
        z = z - modular;
    
    return z;
}

long long precompute_value(long long modular, long long bit_width, long long alpha){
    long long tmp = pow(2, (bit_width + alpha));
    long long result = floor(tmp/modular);
    return result;
}

long long find_prime(long long m, long long powerof2){
	bool flag = 0;
	bool powerof2_flag = 0;
	bool prime_flag = 0;	
	long long i = 0;
	long long tmp ;
	long long init = m * pow(2,powerof2);
	while(flag == 0){
		i++ ;
		tmp = init * i;
		powerof2_flag = 0;
		prime_flag = 0;
		prime_flag = isPrime(tmp+1);
		if(prime_flag == 1){		
			flag = 1;
		}			
	}
	return tmp+1 ;	
}
//------------ZZ ----------------------

ZZ find_n_rou(ZZ base, long long m, ZZ modular) // a^(p-1) = 1 (mod p)  ---> base^(modular-1) = 1 (mod modular)
{
	//cout << " m = " << m << " modular = "<< modular << endl; 
	assert(( modular % m ) == 1);
	ZZ i;
	ZZ n_rou;
	i = (modular-1)/m ;   // base^(modular - 1) = base^( n * i ) = (base^i)^n = 1 (mod modular)
	PowerMod(n_rou, base, i, modular);
	//cout << " n_rou = " << n_rou << endl;
	return n_rou;
}

bool check_prou(ZZ n_rou, long long m, ZZ modular){ //check if n_rou^1, n_rou^2,...,n_rou^(m-1) is not equal 1;
	bool is_prou = true;
	ZZ tmp;
	for(int i = 1; i < m; i++){
		PowerMod(tmp, n_rou, i, modular);
		if(tmp == 1){
			is_prou = false;
			break;
		}
	}
	return is_prou;
}

ZZ find_phi(long long m, ZZ modular)
{   
	bool is_prou = false;
	ZZ i = (ZZ)2 ;
	ZZ n_rou;
	ZZ prou;
    long long degree = 2 * m;
	while(is_prou == false)
	{
		//cout << " 1 " << endl;
		n_rou = find_n_rou(i, degree, modular);
		//cout << " 2 " << endl;
		is_prou = check_prou(n_rou, degree, modular);
		//cout << " 3 " << endl;
		i = i + 1;		
	}
	//cout << " base " << i-1 << endl;
	prou = n_rou;
	return prou;
}




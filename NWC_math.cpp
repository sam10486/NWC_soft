
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

long long prou_power(long long data_in, long long power, long long modular){
    //output = data_in^power % modular
    long long ans;
    
    ans = 1;
    
    if(power >= 2)
    {
        if((power % 2) == 1)
        {
            ans = prou_power(data_in, (power - 1)/2, modular);
            ans *= ans;
            ans %= modular;
            ans *= data_in;
            ans %= modular;
        }
        else
        {
            ans = prou_power(data_in, power/2, modular);
            ans *= ans;
            ans %= modular;
        }
    }
    else if(power == 1)
    {
        ans = data_in;
    }
    
    return ans;
}

long long FFT(long long *DFT_data, long long *data_in, long long n, long long prou, long long modular){ //primitive root of unity in n-point FFT
	long long DFT_data_tmp_1[n];
	long long DFT_data_tmp_2[n];
    long long two_to_i, ind_j;
    long long i, j, k;
	
	long long check_n;// check if m | modular - 1
	check_n = (modular-1) % n ;
	assert(check_n == 0) ;    
	
    for(j=0;j<n;j++)
    {
        DFT_data_tmp_1[j] = data_in[j];
    }
        
    for(i=0;i<log2(n);i++)
    {
        two_to_i=1<<i;
        for(k=0;k<two_to_i;k++)
        {
            for(j=0;j<((n/two_to_i)/2);j++)
            {
                ind_j = j + k * (n/two_to_i);
                //BU2 up output
                DFT_data_tmp_2[ind_j] = DFT_data_tmp_1[ind_j] + DFT_data_tmp_1[ind_j + ((n/two_to_i)/2)];
                DFT_data_tmp_2[ind_j] %= modular;
                //BU2 down output
                DFT_data_tmp_2[ind_j + ((n/two_to_i)/2)] = DFT_data_tmp_1[ind_j] - DFT_data_tmp_1[ind_j + ((n/two_to_i)/2)];
                DFT_data_tmp_2[ind_j + ((n/two_to_i)/2)] *= prou_power(prou, j * two_to_i, modular);
                DFT_data_tmp_2[ind_j + ((n/two_to_i)/2)] %= modular;
            }
        }
        for(j=0;j<n;j++)
        {
            DFT_data_tmp_1[j] = DFT_data_tmp_2[j];
            DFT_data_tmp_1[j] %= modular;
        }
    } 	
    
    //output index
    for(i=0;i<n;i++)
    {
        ind_j = 0;
        for(k=0;k<log2(n);k++)
        {
           if(((i >> k) & (long long)1) == (long long)1)
           {
               ind_j |= (1 << (int)(log2(n) - k - 1));
           }
        }

        DFT_data[ind_j] = DFT_data_tmp_1[i]; //deal with negative
        if(DFT_data[ind_j] < 0)
        {
        	DFT_data[ind_j] += modular;
        }
    }
    
	return 0;
}
long long FFT_no_bit_reverse(long long *DFT_data, long long *data_in, long long n, long long prou, long long modular){ //primitive root of unity in n-point FFT
	long long DFT_data_tmp_1[n];
	long long DFT_data_tmp_2[n];
    long long two_to_i, ind_j;
    long long i, j, k;
	
	long long check_n;// check if m | modular - 1
	check_n = (modular-1) % n ;
	assert(check_n == 0) ;    
	
    for(j=0;j<n;j++)
    {
        DFT_data_tmp_1[j] = data_in[j];
    }
        
    for(i=0;i<log2(n);i++)
    {
        two_to_i=1<<i;
        for(k=0;k<two_to_i;k++)
        {
            for(j=0;j<((n/two_to_i)/2);j++)
            {
                ind_j = j + k * (n/two_to_i);
                //BU2 up output
                //cout << "ind_j = " << ind_j << endl << "ind_j + ((n/two_to_i)/2) = " << ind_j + ((n/two_to_i)/2)<< endl;
                DFT_data_tmp_2[ind_j] = DFT_data_tmp_1[ind_j] + DFT_data_tmp_1[ind_j + ((n/two_to_i)/2)];
                DFT_data_tmp_2[ind_j] %= modular;
                //BU2 down output
                DFT_data_tmp_2[ind_j + ((n/two_to_i)/2)] = DFT_data_tmp_1[ind_j] - DFT_data_tmp_1[ind_j + ((n/two_to_i)/2)];
                DFT_data_tmp_2[ind_j + ((n/two_to_i)/2)] *= prou_power(prou, j * two_to_i, modular);
                DFT_data_tmp_2[ind_j + ((n/two_to_i)/2)] %= modular;
            }
        }
        for(j=0;j<n;j++)
        {
            DFT_data_tmp_1[j] = DFT_data_tmp_2[j];
            DFT_data_tmp_1[j] %= modular;
        }
    } 	
    
    //output index
    for(i=0;i<n;i++)
    {
        DFT_data[i] = DFT_data_tmp_1[i]; //deal with negative
        if(DFT_data[i] < 0)
        {
        	DFT_data[i] += modular;
        }
    }
    
	return 0;
}
long long DFT(long long *DFT_data, long long *data_in, long long m, long long prou, long long modular){ //primitive root of unity in m-point DFT
	long long DFT_data_tmp[m];
    long long i, j, prou_tmp;
	
	long long check_m;// check if m | modular - 1
	check_m = (modular-1) % m ;
	assert(check_m == 0) ;
	
	//cout << "DFT_in = " << endl;	
    for(i=0;i<m;i++)
    {
        DFT_data_tmp[i] = 0;
		//cout << data_in[i] << endl ;		
    }
	//cout << endl; 

 
    for(i=0;i<m;i++)
    {
        prou_tmp = prou_power(prou, i, modular);
        for(j=m-1;j>0;j--)
        {
            DFT_data_tmp[i] += data_in[j];
            DFT_data_tmp[i] *= prou_tmp;
            DFT_data_tmp[i] %= modular;
        }
        DFT_data_tmp[i] = AddMod( DFT_data_tmp[i], data_in[0], modular);
    } 	
    
	//cout << "DFT_out = " << endl;
    for(i=0;i<m;i++)
    {
        DFT_data[i] = DFT_data_tmp[i];
		//cout << DFT_data[i] << endl ;
    }
	//cout << endl;
    
	return 0;
}

vector<long long> NWC(vector<long long> a, long long degree, long long phi, long long modular){
    long long t = degree;
    BitOperate Bitrev;
    for(long long m = 1; m < degree; m = m << 1){
        t = t / 2;
        for(long long i = 0; i < m; i++){
            long long j1 = 2*i*t;
            long long j2 = j1 + t - 1;
            long long bit_num = ceil(log2(degree));
            long long index = Bitrev.BitReserve((m+i), bit_num);
            long long S = ExpMod(phi, index, modular);
            for(long long j = j1; j <= j2; j++){
                long long U = a.at(j);
                long long V = MulMod(a.at(j+t), S, modular);
                a.at(j) = AddMod(U, V, modular);
                a.at(j+t) = SubMod(U, V, modular);
            }
        }
    }
    return a;
}

vector<long long> INWC(vector<long long> a, long long degree, long long phi, long long modular){
    long long t = 1;
    long long bit_num = ceil(log2(degree));
    BitOperate Bitrev;
    for(long long m=degree; m>1; m = m >> 1){
        long long j1 = 0;
        long long h = m / 2;
        for(long long i = 0; i < h; i++){
            long long j2 = j1 + t - 1;
            long long index = Bitrev.BitReserve((h+i), bit_num);
            long long S = ExpMod(phi, index, modular);
            for(long long j = j1; j <= j2; j++){
                long long U = a.at(j);
                long long V = a.at(j+t);
                a.at(j) = AddMod(U, V, modular);
                a.at(j+t) = MulMod(SubMod(U,V,modular), S, modular);
            }
            j1 = j1 + 2*t;
        }
        t = 2*t;
    }
    long long degree_inv = InvMod(degree, modular);
    for(long long j = 0; j < degree; j++){
        a.at(j) = MulMod(a.at(j), degree_inv, modular);
    }
    return a;
}

long long NWC_FFT_no_bit_reverse(long long *DFT_data, long long *data_in, long long n, long long prou, long long phi, long long modular){ //primitive root of unity in n-point FFT
	long long DFT_data_tmp_1[n];
	long long DFT_data_tmp_2[n];
    long long two_to_i, ind_j;
    long long i, j, k;
    long long phi_exp;
    long long prou_exp;
    long long twiddle;

    long long check_n;// check if m | modular - 1
	check_n = (modular-1) % (2*n);
    assert(check_n == 0) ; 

    for(j=0;j<n;j++){
        DFT_data_tmp_1[j] = data_in[j];
    }

    for(i=0;i<log2(n);i++){
        two_to_i=(n/2)>>i;
        for(k=0;k<two_to_i;k++){
            for(j=0;j<((n/two_to_i)/2);j++){
                ind_j = j + k * (n/two_to_i);
                phi_exp = ExpMod(phi, two_to_i, modular);
                prou_exp = ExpMod(prou, j*two_to_i, modular);
                twiddle = MulMod(phi_exp, prou_exp, modular);
                //cout << "i = " << i << "   two_to_i = " << two_to_i << endl;
                //cout << "phi_exp = " << phi_exp << "   " << "prou_exp = " << prou_exp << "   " << "twiddle = " << twiddle << endl;
                //cout << "ind_j = " << ind_j << "   " << "ind_j + ((n/two_to_i)/2) = " << ind_j + ((n/two_to_i)/2) << endl;
                //cout << "------------------------------------------" << endl;
                DFT_data_tmp_1[ind_j + ((n/two_to_i)/2)] = MulMod(DFT_data_tmp_1[ind_j + ((n/two_to_i)/2)], twiddle, modular);
                //BU2 up output
                DFT_data_tmp_2[ind_j] = AddMod(DFT_data_tmp_1[ind_j], DFT_data_tmp_1[ind_j + ((n/two_to_i)/2)], modular);
                //BU2 down output
                DFT_data_tmp_2[ind_j + ((n/two_to_i)/2)] = SubMod(DFT_data_tmp_1[ind_j], DFT_data_tmp_1[ind_j + ((n/two_to_i)/2)], modular);
                //cout << "phi_exp = " << phi_exp << "   " << "prou_exp = " << prou_exp << endl;
                //cout << "twiddle = " << twiddle << endl;
                //cout << "DFT_data_tmp_1[" << ind_j + ((n/two_to_i)/2) << "] = " <<  DFT_data_tmp_1[ind_j + ((n/two_to_i)/2)] << endl;
                //cout << "DFT_data_tmp_2[" << ind_j << "] = " << DFT_data_tmp_2[ind_j] << endl;
                //cout << "DFT_data_tmp_2[" << ind_j + ((n/two_to_i)/2) << "] = " << DFT_data_tmp_2[ind_j + ((n/two_to_i)/2)] << endl;
            }
        }
        //cout << "1231321231313" << endl;
        for(j=0;j<n;j++){
            DFT_data_tmp_1[j] = DFT_data_tmp_2[j];
            DFT_data_tmp_1[j] %= modular;
        }
    }
    //output index
    for(i=0;i<n;i++)
    {
        DFT_data[i] = DFT_data_tmp_1[i]; //deal with negative
        if(DFT_data[i] < 0)
        {
        	DFT_data[i] += modular;
        }
    }	

    return 0;
}


long long NWC_FFT_no_bit_reverse_DIF_ver(long long *DFT_data, long long *data_in, long long n, long long prou, long long phi, long long modular){ //primitive root of unity in n-point FFT
	long long DFT_data_tmp_1[n];
	long long DFT_data_tmp_2[n];
    long long two_to_i, ind_j;
    long long i, j, k;
    long long phi_exp;
    long long prou_exp;
    long long twiddle;

    long long check_n;// check if m | modular - 1
	check_n = (modular-1) % (2*n);
    assert(check_n == 0) ; 

    for(j=0;j<n;j++){
        DFT_data_tmp_1[j] = data_in[j];
    }

    for(i=0;i<log2(n);i++){
        two_to_i=(n/2)>>i;
        for(k=0;k<two_to_i;k++){
            for(j=0;j<((n/two_to_i)/2);j++){
                ind_j = j + k * (n/two_to_i);
                phi_exp = ExpMod(phi, two_to_i, modular);
                prou_exp = ExpMod(prou, j*two_to_i, modular);
                twiddle = MulMod(phi_exp, prou_exp, modular);
                //cout << "i = " << i << "   two_to_i = " << two_to_i << endl;
                //cout << "phi_exp = " << phi_exp << "   " << "prou_exp = " << prou_exp << "   " << "twiddle = " << twiddle << endl;
                //cout << "ind_j = " << ind_j << "   " << "ind_j + ((n/two_to_i)/2) = " << ind_j + ((n/two_to_i)/2) << endl;
                //cout << "------------------------------------------" << endl;
                DFT_data_tmp_1[ind_j + ((n/two_to_i)/2)] = MulMod(DFT_data_tmp_1[ind_j + ((n/two_to_i)/2)], twiddle, modular);
                //BU2 up output
                DFT_data_tmp_2[ind_j] = AddMod(DFT_data_tmp_1[ind_j], DFT_data_tmp_1[ind_j + ((n/two_to_i)/2)], modular);
                //BU2 down output
                DFT_data_tmp_2[ind_j + ((n/two_to_i)/2)] = SubMod(DFT_data_tmp_1[ind_j], DFT_data_tmp_1[ind_j + ((n/two_to_i)/2)], modular);
                //cout << "phi_exp = " << phi_exp << "   " << "prou_exp = " << prou_exp << endl;
                //cout << "twiddle = " << twiddle << endl;
                //cout << "DFT_data_tmp_1[" << ind_j + ((n/two_to_i)/2) << "] = " <<  DFT_data_tmp_1[ind_j + ((n/two_to_i)/2)] << endl;
                //cout << "DFT_data_tmp_2[" << ind_j << "] = " << DFT_data_tmp_2[ind_j] << endl;
                //cout << "DFT_data_tmp_2[" << ind_j + ((n/two_to_i)/2) << "] = " << DFT_data_tmp_2[ind_j + ((n/two_to_i)/2)] << endl;
            }
        }
        //cout << "1231321231313" << endl;
        for(j=0;j<n;j++){
            DFT_data_tmp_1[j] = DFT_data_tmp_2[j];
            DFT_data_tmp_1[j] %= modular;
        }
    } 

    //output index
    for(i=0;i<n;i++)
    {
        DFT_data[i] = DFT_data_tmp_1[i]; //deal with negative
        if(DFT_data[i] < 0)
        {
        	DFT_data[i] += modular;
        }
    }	

    return 0;
}


long long NWC_forward_DIT(long long *NWC_ans, long long *a, long long degree, long long phi, long long modular){
    BitOperate BR;
    // vector s a bit reverse
    long long bit_width = log2(degree);
    for(long long i=0; i<degree; i++){
        long long idx_BR = BR.BitReserve(i, bit_width);
        NWC_ans[i] = a[idx_BR];
        //cout << "input = " << NWC_ans[i] << endl;
    }

    for(long long s=1; s<=log2(degree); s++){
        //cout << "stage = " << s << endl;
        long long m = ExpMod(2, s, modular);
        for(long long j=0; j<=(m/2)-1; j++){
            long long twiddle_dg = ((2*j+1) * degree)/m;
            long long twiddle = ExpMod(phi, twiddle_dg, modular);
            //cout << " twiddle = " <<  twiddle << endl;
            for(long long k=0; k<=(degree/m)-1; k++){
                long long u = NWC_ans[k*m+j];
                long long t = MulMod(twiddle, NWC_ans[k*m+j+(m/2)], modular);
                //cout << " NWC_ans[k*m+j+(m/2)] = " <<  NWC_ans[k*m+j+(m/2)] << endl;
                NWC_ans[k*m+j] = AddMod(u, t, modular);
                NWC_ans[k*m+j+(m/2)] = SubMod(u, t, modular);
                //cout << " u = " <<  u << endl;
                //cout << " t = " <<  t << endl;
                //cout << " NWC_ans up = " <<  NWC_ans[k*m+j] << endl;
                //cout << " NWC_ans down = " <<  NWC_ans[k*m+j+(m/2)] << endl;
            }
        }
    }
}

long long power2_NTT(   long long *NTT_data, long long *data_in, long long n, 
                        long long *twiddle_array, long long modular){
    
    long long W;
    /*for(int i=0;i<n;i++)
        //cout << data_in[i] << " ";
        //cout << twiddle_array[i] << " ";
    cout << endl;*/
    //cout << "n = " << n << endl;
    for(int h=1; h<n; h = 2 * h){
        for(int j=0; j<h; j++){
            W = twiddle_array[h+j];
            //cout << "h = " << h << endl;
            //cout << "j = " << j << endl;
            //cout << "W[" << h+j << "]= " << W << endl;
            for(int i=(j*n)/h; i< ( (2*j+1)*n ) / (2*h) ; i++){
                long long tmp = MulMod(W, data_in[i + n/(2*h)], modular);
                //cout << "data_in[" << i + n/(2*h) << "] = " << data_in[i + n/(2*h)] << endl;
                //cout << "tmp = " << tmp << endl;
                data_in[i + n/(2*h)] = SubMod(data_in[i], tmp, modular);
                data_in[i] = AddMod(data_in[i], tmp, modular);
                //cout << "data_in[" << i + n/(2*h) << "] = " << data_in[i + n/(2*h)] << endl;
                //cout << "data_in[" << i << "] = " << data_in[i] << endl;
            }
        }
    }
    for(int i=0; i<n; i++){
        NTT_data[i] = data_in[i];
    }
}

long long mixed_radix_NWC(  long long *NWC_data, long long *NWC_data_in, 
                            long long n, long long radix_k1, long long radix_k2, long long phi, 
                            long long modular){

    long long k = ( log2(n) - radix_k2) / radix_k1;
    cout << "k = " << k << endl;
    long long parameter_check = radix_k1 * k + radix_k2;
    assert(parameter_check == log2(n));

    int element_num = pow(2, radix_k1);
    long long W[element_num - 1] = {0};
    long long tmp[element_num] = {0} ;
    long long tmp_array[element_num] = {0};
    BitOperate BR;
    //k radix 2^k1 NTTs
    for(long long l=0; l<k; l++){
        for(long long m=1; m<pow(2, radix_k1); m++){
            long long m_bar = BR.BitReserve( (pow(2, radix_k1*l)-1) * pow(2, floor(log2(m))) + m, log2(n) ); 
            //cout << "(pow(2, radix_k1*l)-1) * pow(2, floor(log2(m))) + m = " << (pow(2, radix_k1*l)-1) * pow(2, floor(log2(m))) + m << endl;
            //cout << "m_bar = " << m_bar << endl;
            W[m] = ExpMod(phi, m_bar, modular);
            //cout << "W[" << m << "] = " << W[m] << endl;
        }
        cout << "------------" << endl;
        for(int j=0; j<pow(2, radix_k1*l); j++){
            long long j_bar = BR.BitReserve(j, radix_k1*l);
            for(int i=0; i<pow(2, log2(n)- radix_k1*(l+1)); i++){
                for(int m=0; m<pow(2, radix_k1); m++){
                    long long index = j_bar * pow(2, log2(n) - radix_k1*l) + m * pow(2, log2(n) - radix_k1*(l+1)) + i;
                    tmp[m] = NWC_data_in[index];
                    cout << "BU[" << i << "] = " << index << endl;
                    //cout << "tmp[" << m << "] = " << tmp[m] << endl;
                }
                cout << "i = " << i << endl;
                /*for(int i=0; i<pow(2, radix_k1);i++)
                    cout << W[i] << " ";
                cout << endl;*/
                power2_NTT(tmp_array, tmp, pow(2, radix_k1), W, modular); 

                for(int m=0; m<pow(2, radix_k1); m++){
                    long long index = j_bar * pow(2, log2(n) - radix_k1*l) + m * pow(2, log2(n) - radix_k1*(l+1)) + i;
                    NWC_data_in[index] = tmp_array[m];
                    //cout << "tmp_array[" << m << "] = " << tmp_array[m] << endl;
                    //cout << "NWC_data_in[" << index << "] = " << NWC_data_in[index] << endl;
                }
            }

            for(int m=1; m<pow(2, radix_k1); m++){
                long long Wc_degree = log2(n) - radix_k1 * l - floor(log2(m));
                if(Wc_degree != log2(n))
                    Wc_degree = pow(2, Wc_degree);
                else 
                    Wc_degree = 0;
                //cout << "Wc_degree = " << Wc_degree << endl;
                long long Wc = ExpMod(phi, Wc_degree, modular);
                W[m] = MulMod(W[m], Wc, modular);
                //cout << "W[" << m << "] = " << W[m] << endl;
            }
        }
    }
    for(int i=0; i<n; i++)
        cout << NWC_data_in[i] << " ";
    cout << endl;

    cout << "------------------" << endl;
    int element_num_k2 = pow(2, radix_k2);
    long long W_k2[element_num_k2 - 1] = {0} ;
    long long tmp_k2[element_num_k2] = {0} ;
    long long tmp_array_k2[element_num_k2] = {0};
    // radix_2^k2 NTT
    for(int m=1; m<pow(2, radix_k2); m++){
        long long idx = (pow(2, radix_k1*k) - 1) * pow(2, floor(log2(m))) + m;
        //cout << "idx = " << idx << endl;
        long long m_bar = BR.BitReserve(idx, log2(n));
        //cout << "m_bar = " << m_bar << endl;
        W_k2[m] = ExpMod(phi, m_bar, modular);
        //cout << "W_k2[" << m << "] = " <<  W_k2[m] << endl;
    }
    for(int j=0; j<pow(2, radix_k1 * k); j++){
        long long j_bar = BR.BitReserve(j, radix_k1 * k);
        for(int m=0; m<pow(2, radix_k2); m++){
            int idx = j_bar * pow(2, radix_k2) + m;
            tmp_k2[m] = NWC_data_in[idx];
        }

        for(int i=0; i<pow(2, radix_k2); i++){
            //cout << "tmp_k2[" << i << "] = " << tmp_k2[i] << " ";
        }
        //cout << endl;
        power2_NTT(tmp_array_k2, tmp_k2, pow(2,radix_k2), W_k2, modular);
        for(int m=0; m<pow(2, radix_k2); m++){
            int idx = j_bar * pow(2, radix_k2) + m;
            NWC_data_in[idx] = tmp_array_k2[m];
            //cout << "NWC_data_in[" << idx << "] = " << NWC_data_in[idx] << endl;
        }
        for(int m=1; m<pow(2, radix_k2); m++){
            long long Wc_degree_k2 = radix_k2 - floor(log2(m));
            Wc_degree_k2 = pow(2, Wc_degree_k2);
            long long Wc_k2 = ExpMod(phi, Wc_degree_k2, modular);
            W_k2[m] = MulMod(W_k2[m], Wc_k2, modular);
        }
    }

    for(int i=0; i<n; i++){
        NWC_data[i] = NWC_data_in[i];
        //cout << "NWC_data[" << i << "] = " << NWC_data[i] << endl;
    }
    cout << endl;
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

bool isPrime(ZZ n){
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

ZZ find_prime(ZZ m, long long powerof2){
	bool flag = 0;
	bool powerof2_flag = 0;
	bool prime_flag = 0;	
	ZZ i ;
	ZZ tmp ;
	ZZ init;
	i = 0;	
	init = m * power((ZZ)2,powerof2);
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

ZZ precompute_value(ZZ modular, long long bit_width, long long alpha){
    ZZ tmp, result;
    tmp = power2_ZZ((bit_width + alpha));
    div(result, tmp, modular);
    return result;
}


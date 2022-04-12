#ifndef NWC_math
#define NWC_math

#include <iostream>
#include <vector>
using namespace std;

long long AddMod(long long, long long, long long);
long long SubMod(long long in_1, long long in_2, long long q);
long long XGCD(long long a, long long q);
long long InvMod(long long a, long long q);
long long DivMod(long long a,long long b,long long q);
long long MulMod(long long a,long long b,long long q);
long long ExpMod(long long a, long long exp, long long q);
long long find_phi(long long n, long long modular); //// p = 1 mod 2n
bool isPrime(long long n);
long long find_prou(long long n, long long modular);
vector<long long> phi_array(long long n, long long modular);
vector<long long> phi_array_inv(long long n, long long modular);
#endif
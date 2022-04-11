#ifndef NWC_math
#define NWC_math

#include <iostream>
#include <vector>
using namespace std;

int AddMod(int, int, int);
int SubMod(int in_1, int in_2, int q);
int XGCD(int a, int q);
int InvMod(int a, int q);
int DivMod(int a,int b,int q);
int MulMod(int a,int b,int q);
int ExpMod(int a, int exp, int q);
long long find_phi(long long n, long long modular); //// p = 1 mod 2n
bool isPrime(long long n);
long long find_prou(long long n, long long modular);
vector<long long> phi_array(long long n, long long modular);
#endif
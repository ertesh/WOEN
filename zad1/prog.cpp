#include <iostream>
#include <vector>
#include <cmath>
#include <cassert>
#include <cstdlib>
#include <ctime>

using namespace std;

typedef long long LL;
typedef vector<LL> (*ffun)(const vector<LL>&, double);

const int runs = 2; //20
const int iterations = 10; //100000
const int stepDistance = 100;
const int m = 50; // 50
const int L = 60;
const int minX = -32768;
const int maxX = 32768;


double ackley(const vector<double> &x) {
    const double a = 20;
    const double b = 0.2;
    assert(m == x.size());
    double sum1 = 0;
    double sum2 = 0;
    for (int i = 0; i < m; i++) {
        sum1 += x[i] * x[i];
        sum2 += cos(2 * M_PI * x[i]);
    }
    sum1 /= m;
    sum2 /= m;
    double ret = -a * exp(-b * sqrt(sum1)) - exp(sum2) + a + exp(1);
    cout << "ackley " << ret << endl;
    for (int i = 0; i < m; i++) cout << x[i] << " ";
    cout << "\n";
    return ret;
}

vector<LL> randomStart() {
    vector<LL> ret(m);
    for (int i = 0; i < m; i++) {
        LL p = 0;
        for (int j = 0; j < L; j++) {
            int q = rand() % 2;
            p = 2 * p + q;
        }
        ret[i] = p;
    }
    return ret;
}

LL decodeGray(LL j) {
    LL i = j;
    while (j > 1) {
        j /= 2;
        i ^= j;
    }
    return i;
}

vector<LL> mutateAdd(const vector<LL> &x, double p) {
    int prog = p * RAND_MAX;
    vector<LL> ret(m);
    for (int i = 0; i < m; i++) {
        LL p = 0;
        for (int j = 0; j < L; j++) {
            int q = rand();
            p = 2 * p + (q < prog);
        }
        LL mod = 1LL << L;
        if (rand() % 2) {
            ret[i] = (ret[i] + p) % mod;
        }
        else {
            ret[i] = (ret[i] + mod - p) % mod;
        }
//        cout << a << " " << b << " " << aa << " " << bb << endl;
    }
    return ret;
}

vector<LL> mutateBern(const vector<LL> &x, double p) {
    int prog = p * RAND_MAX;
    vector<LL> ret(m);
    for (int i = 0; i < m; i++) {
        LL p = 0;
        for (int j = 0; j < L; j++) {
            int q = rand();
            p = 2 * p + (q < prog);
        }
        ret[i] = x[i] ^ p;
        LL a = decodeGray(x[i]);
        LL b = decodeGray(ret[i]);
        double e = (double) (maxX - minX) / (1LL << L);
        double aa = minX + e * a;
        double bb = minX + e * b;
        cout << a << " " << b << " " << aa << " " << bb << endl;
    }
    return ret;
}

vector<double> decode(const vector<LL> &x, bool isGray) {
    double e = (double) (maxX - minX) / (1LL << L);
    cout << e << " " << (1LL << L) << endl;
    vector<double> ret(x.size(), 0);
    for (int i = 0; i < m; i++) {
        LL j = i;
        if (isGray) j = decodeGray(x[i]);
        ret[i] = minX + e * j;
    }
    return ret;
}

int testAlgorithm(double p, ffun mutate, bool isGray) {
    int size = iterations / stepDistance + 1;
    vector<double> E(size, 0);
    vector<double> Var(size, 0);
    vector<double> Min(size, 1e9);
    vector<double> Max(size, 0);

    for (int i = 0; i < runs; i++) {
        cout << "Run " << i << endl;
        vector<LL> x = randomStart();
        double val = ackley(decode(x, isGray));

        for (int j = 0; j <= iterations; j++) {
            if (j % stepDistance == 0) {
                int k = j / 100;
                E[k] += val;
                Var[k] += val * val;
                Min[k] = min(Min[k], val);
                Max[k] = max(Max[k], val);
            }
            vector<LL> y = (*mutate)(x, p);
            double nextVal = ackley(decode(y, isGray));
            if (nextVal < val) {
                val = nextVal;
                x.clear();
                x = y;
            }
        }
    }
    for (int i = 0; i < size; i++) {
        E[i] /= runs;
        Var[i] = sqrt((Var[i] - runs * E[i] * E[i]) / (runs - 1));
        cout << i << " " << E[i] << " " << Var[i] << " " << Min[i] << " " << Max[i] << "\n";
    }
}

int main() {
    srand(12345);
//    testAlgorithm(0.005, mutateBern, true);
    testAlgorithm(0.005, mutateBern, false);
//    testAlgorithm(0.005, mutateAdd, true);
//    testAlgorithm(0.005, mutateAdd, false);
}


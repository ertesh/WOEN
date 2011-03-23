#include <iostream>
#include <vector>
#include <cmath>
#include <cassert>
#include <cstdlib>
#include <ctime>

using namespace std;

typedef long long LL;
typedef double LD;
typedef vector<LL> (*ffun)(const vector<LL>&, LD);

const int runs = 2; //20
const int iterations = 20000; //100000
const int stepDistance = 2000;
const int m = 50; // 50
const int L = 60;
const LD minX = -32.768;
const LD maxX = 32.768;

const LD e = (maxX - minX) / (1LL << L);


LD ackley(const vector<LD> &x) {
    const LD a = 20;
    const LD b = 0.2;
    assert(m == (int) x.size());
    LD sum1 = 0;
    LD sum2 = 0;
    for (int i = 0; i < m; i++) {
        sum1 += x[i] * x[i];
        sum2 += cos(2 * M_PI * x[i]);
    }
    sum1 /= m;
    sum2 /= m;
    LD ret = -a * exp(-b * sqrt(sum1)) - exp(sum2) + a + exp(1);
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

vector<LL> mutateAdd(const vector<LL> &x, LD prob) {
    int prog = prob * RAND_MAX;
    vector<LL> ret(m);
    for (int i = 0; i < m; i++) {
        LL p = 0;
        for (int j = 0; j < L; j++) {
            int q = rand();
            p = 2 * p + (q < prog);
        }
        LL mod = 1LL << L;
        if (rand() % 2) {
            ret[i] = (x[i] + p) % mod;
        }
        else {
            ret[i] = (x[i] + mod - p) % mod;
        }
    }
    return ret;
}

vector<LL> mutateBern(const vector<LL> &x, LD prob) {
    int prog = prob * RAND_MAX;
    vector<LL> ret(m);
    for (int i = 0; i < m; i++) {
        LL p = 0;
        for (int j = 0; j < L; j++) {
            int q = rand();
            p = 2 * p + (q < prog);
        }
        ret[i] = x[i] ^ p;
    }
    return ret;
}

vector<LD> decode(const vector<LL> &x, bool isGray) {
    vector<LD> ret(x.size(), 0);
    for (int i = 0; i < m; i++) {
        LL j = x[i];
        if (isGray) j = decodeGray(x[i]);
        ret[i] = minX + e * j;
    }
    return ret;
}

void testAlgorithm(LD p, ffun mutate, bool isGray) {
    int size = iterations / stepDistance + 1;
    vector<LD> E(size, 0);
    vector<LD> Var(size, 0);
    vector<LD> Min(size, 1e9);
    vector<LD> Max(size, 0);

    for (int i = 0; i < runs; i++) {
        cout << "Run " << i << endl;
        vector<LL> x = randomStart();
        LD val = ackley(decode(x, isGray));

        for (int j = 0; j <= iterations; j++) {
            if (j % stepDistance == 0) {
                int k = j / stepDistance;
                E[k] += val;
                Var[k] += val * val;
                Min[k] = min(Min[k], val);
                Max[k] = max(Max[k], val);
            }
            vector<LL> y = (*mutate)(x, p);
            LD nextVal = ackley(decode(y, isGray));
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
    testAlgorithm(0.001, mutateBern, true);
    testAlgorithm(0.001, mutateBern, false);
    testAlgorithm(0.001, mutateAdd, true);
    testAlgorithm(0.001, mutateAdd, false);
}


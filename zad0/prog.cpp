#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <cmath>

using namespace std;

const int K = 100; // Rodzaje kuponow = 100

int main() {
    int N, M;
    int p = scanf("%d %d", &N, &M);
    vector<double> E(N + 1, 0);
    vector<double> Var(N + 1, 0);
    vector<int> Mi(N + 1, K);
    Mi[0] = 0;
    vector<int> Ma(N + 1, 0);
    srand(12345);
    for (int te = 0; te < M; te++) {
        vector<bool> v(K, 0);
        int counter = 0;
        for (int i = 1; i <= N; i++) {
            int k = rand() % K;
            if (!v[k]) {
                v[k] = 1;
                counter++;
            }
            E[i] += counter;
            Var[i] += counter * counter;
            Mi[i] = min(Mi[i], counter);
            Ma[i] = max(Ma[i], counter);
        }
    }
    for (int i = 0; i <= N; i++) {
        E[i] /= M;
        Var[i] -= M * E[i] * E[i];
        Var[i] /= M - 1; // now it's variation
        Var[i] = sqrt(Var[i]); // now it's standard deviation
        printf("%d %.5f %.5f %d %d\n", i, E[i], Var[i], Ma[i], Mi[i]);
    }
}


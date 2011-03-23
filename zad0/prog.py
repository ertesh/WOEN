#!/usr/bin/env python

import random
import Gnuplot

K = 10 # Rodzaje kuponow = 100
N = 5 * K # Ilosc losowan
M = 1000 # Liczba testow

E = [0 for i in range(N + 1)]
Var = [0 for i in range(N + 1)]
Mi = [K for i in range(N + 1)]
Ma = [0 for i in range(N + 1)]

def draw(*tab):
    gp = Gnuplot.Gnuplot(debug = 1)
    for el in tab:
        data = Gnuplot.Data(range(N + 1), el)
        data.set_option(='a')
        gp.replot(data)

    gp.hardcopy(filename="a.png", terminal="png")

random.seed(12345)

for test in range(M):
    tab = [0 for i in range(K)]
    counter = 0
    for i in range(1, N + 1):
        k = random.randint(0, K - 1)
        if tab[k] == 0:
            tab[k] = 1
            counter += 1
        E[i] += counter
        Var[i] += counter * counter
        Mi[i] = min(Mi[i], counter)
        Ma[i] = max(Ma[i], counter)

for i in range(N + 1):
    E[i] /= float(M)
    Var[i] -= M * E[i] * E[i]

draw(E, Mi, Ma)


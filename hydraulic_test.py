import start
import hydraulic as hc
import Data as dt
from math import pi
import saor_data as sdt

T = start.start_calc(dt.Gpb, 350 + 273.15)
G = dt.Gpb
p = 0
i = 0
# Подъемная ветвь
dx = dt.h_az / dt.n_az
p += hc.dp_az_pg_in_and_out(G, T[i], 0.5)
for j in range(dt.n_az):
    p += hc.dp_friction_az(G, T[i])
    if j % 5 == 0:
        p += hc.dp_az_spacer(G, T[i])
    i += 1
p += hc.dp_az_pg_in_and_out(G, T[i], 0.5)

dx = dt.h_1 / dt.n_1
for _ in range(dt.n_1):
    p += hc.dp_friction(G, T[i], dx, dt.f_1)
    i += 1
    
p += hc.dp_expansion(G, T[i], dt.f_1, dt.f_2)

dx = dt.h_2 / dt.n_2
for _ in range(dt.n_2):
    p += hc.dp_friction(G, T[i], dx, dt.f_2)
    i += 1
    
p += hc.dp_bend_90(G / 4, T[i], dt.f_2)

dx = dt.l_3 / dt.n_3
for _ in range(dt.n_3):
    p += hc.dp_friction(G / 4, T[i], dx, dt.f_3)
    i += 1
p += hc.dp_bend_90(G / 4, T[i], dt.f_3)
p += hc.dp_az_pg_in_and_out(G, T[i], 0.5)
dx = dt.h_pg / dt.n_pg
for _ in range(dt.n_pg):
    p += hc.dp_bundle_pg(G, T[i], dx)
    i += 1  
p += hc.dp_az_pg_in_and_out(G, T[i], 0.5)

p += hc.dp_bend_90(G / 4, T[i], dt.f_4)

dx = dt.h_4 / dt.n_4
for _ in range(dt.n_4):
    p += hc.dp_friction(G / 4, T[i], dx, dt.f_4)
    i += 1 
dx = dt.l_5 / dt.n_5

p += hc.dp_bend_90(G / 4, T[i], dt.f_5)

for _ in range(dt.n_5):
    p += hc.dp_friction(G / 4, T[i], dx, dt.f_5)
    i += 1

p += hc.dp_bend_90(G / 12, T[i], dt.f_5)
# САОР — нисходящая ветвь
dx = dt.h_saor / dt.n_saor
for _ in range(dt.n_saor):
    p += hc.dp_friction_saor(G / 12, T[i])
    i += 1
p += hc.dp_expansion(G / 12, T[i], pi * (sdt.d5 ** 2 - sdt.d4 ** 2) / 4, dt.f_6)
# Опускной участок
dx = dt.h_6 / dt.n_6
for _ in range(dt.n_6):
    p += hc.dp_friction(G / 12, T[i], dx, dt.f_6)
    i += 1
p += hc.dp_bend_90(G / 4, T[i], dt.f_6)
# Горизонталь до АЗ
dx = dt.l_7 / dt.n_7
for _ in range(dt.n_7):
    p += hc.dp_friction(G / 4, T[i], dx, dt.f_7)
    i += 1
p += hc.dp_bend_90(G / 4, T[i], dt.f_7)
print(p)
import start
import hydraulic as hc
import Data as dt

T = start.start_calc(dt.Gpb, 350 + 273.15)
G = dt.Gpb
p = 0
i = 0
# Подъемная ветвь
dx = dt.h_az / dt.n_az
for _ in range(dt.n_az):
    p += hc.dp_friction(G, T[i], dx, dt.f_az, dt.dg_az)
    i += 1

dx = dt.h_1 / dt.n_1
for _ in range(dt.n_1):
    p += hc.dp_friction(G, T[i], dx, dt.f_1, dt.dg_az)
    i += 1

dx = dt.h_2 / dt.n_2
for _ in range(dt.n_2):
    p += hc.dp_friction(G, T[i], dx, dt.f_2)
    i += 1

dx = dt.l_3 / dt.n_3
for _ in range(dt.n_3):
    p += hc.dp_friction(G, T[i], dx, dt.f_3)
    i += 1

dx = dt.h_pg / dt.n_pg
for _ in range(dt.n_pg):
    p += hc.dp_friction(G, T[i], dx, dt.f_pg, dt.dg_pg)
    i += 1
    
dx = dt.h_4 / dt.n_4
for _ in range(dt.n_4):
    p += hc.dp_friction(G, T[i], dx, dt.f_4)
    i += 1
    
dx = dt.l_5 / dt.n_5
for _ in range(dt.n_5):
    p += hc.dp_friction(G, T[i], dx, dt.f_5)
    i += 1
# САОР — нисходящая ветвь
dx = dt.h_saor / dt.n_saor
for _ in range(dt.n_saor):
    p += hc.dp_friction_saor(G / 12, T[i])
    i += 1
# Опускной участок
dx = dt.h_6 / dt.n_6
for _ in range(dt.n_6):
    p += hc.dp_friction(G / 12, T[i], dx, dt.f_6)
    i += 1
# Горизонталь до АЗ
dx = dt.l_7 / dt.n_7
for _ in range(dt.n_7):
    p += hc.dp_friction(G / 12, T[i], dx, dt.f_7)
    i += 1

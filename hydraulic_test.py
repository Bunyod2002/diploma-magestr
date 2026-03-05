# ruff: noqa: F403, F405
from Data import *
from hydraulic import *
from functions import *
import start
T = start.T
i = 0
p = 0
for j in range(n_az):
    p += dp_friction(T[i], h_az / n_az, f_az, dg_az)
    i += 1
for j in range(n_1):
    p += dp_friction(T[i], h_1 / n_1, f_1)
    i += 1
for j in range(n_2):
    p += dp_friction(T[i], h_2 / n_2, f_2)
    i += 1
for j in range(n_pg):
    p += dp_friction(T[i], h_pg / n_pg, f_pg, dg_pg)
    i += 1
for j in range(n_3):
    p += dp_friction(T[i], h_3 / n_3, f_3, dg_az)
    i += 1
for j in range(n_4):
    p += dp_friction(T[i], h_4 / n_4, f_4)
    i += 1
for j in range(n_5):
    p += dp_friction(T[i], l_5 / n_5, f_5)
    i += 1
for j in range(n_6):
    p += dp_friction(T[i], h_6 / n_6, f_6)
    i += 1


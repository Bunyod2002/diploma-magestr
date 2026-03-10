
# ruff: noqa: F403, F405
from start import T
from Data import *
from math import exp
from functions import *  

h = 2.25
T0 = T
dtime = 600
time = 0.05
n = 0

# Цикл естественной циркуляции
# 1.Активная зона  2.Область до тягового участка  3. Тяговый участок до отметки естественной циркуляции
# 4. Горизонтальный участок до САОР  5.Теплообменник Фильда  6. Опускной участок 7. Горизонтальный участок до АЗ 
T1 = [T0[-1]] + [0] * N
def part_x(n, f, h, G):
    global i, T1
    dx = h / n
    dt_dx = dt / dx
    for j in range(n):
        t_i = T0[i]  # T_i-1_k
        t_i_1 = T0[i-1]  # T_i_k
        r = ro(t_i_1 - 273.15)
        t_k_1 = round(t_i + dt_dx * (G* cp_Pb*(t_i_1 - t_i)) / (cp_Pb * r * f), 3)
        T1[i] = t_k_1
        i += 1
        
while time < dtime:
    # Активная зона
    if time < 120:
        G = Gpb * exp(-time/60)
    Q_veg = 0.065 * Q * (time ** (-0.2) - (time + 2592000) ** (-0.2))
    i = 1
    dx = h_az / n_az
    dt_dx = dt / dx
    h = 2.25
    z = 0
    for j in range(n_az):
        t_i = T0[i]  # T_i-1_k
        t_i_1 = T0[i-1]  # T_i_k
        r = ro(t_i_1 - 273.15)
        t_k_1 = round(t_i + dt_dx * (G* cp_Pb*(t_i_1 - t_i) + ql(Q_veg, z) * dx) / (cp_Pb * r * f_az), 3)
        T1[i] = t_k_1
        h += dx
        i += 1
        z += dx
    # Область до тягового участка
    part_x(n_1, f_1, h_1, G)
    # Тяговый участок
    part_x(n_2, f_2, h_2, G)
    # Парогенератор
    h = 7.75
    dx = h_pg / n_pg
    dt_dx = dt / dx

    for j in range(n_pg):
        t_i = T0[i]  # T_i-1_k
        t_i_1 = T1[i - 1]  # T_i_k
        r = ro(t_i_1 - 273.15)
        alfa = alphaPb(r, f_pg, dg_pg, G)
        t_k_1 = round(t_i + dt_dx * (G * cp_Pb * (t_i_1 - t_i) + alfa * s_pg * (T_pg - t_i)) / (r * f_pg * cp_Pb), 3)
        T1[i] = t_k_1
        h -= dx
        i += 1

    
    time += dt
    T0 = T1[:]
    T1 = [T1[-1]] + [0] * N









# ruff: noqa: F403, F405
from Data import *
from math import exp
from functions import * 
from saor import *

h = 2.25
T0 = [0] * (N + 1)
dtime = 1000
time = 0.05
n = 0

# Цикл естественной циркуляции
# 1.Активная зона  2.Область до тягового участка  3. Тяговый участок до отметки естественной циркуляции
# 4. Горизонтальный участок до САОР  5.Теплообменник Фильда  6. Опускной участок 7. Горизонтальный участок до АЗ 
T1 = [623,15] + [0] * N
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
    G = Gpb
    Q_veg = 0.065 * Q * (time ** (-0.2) - (time + 2592000) ** (-0.2))
    i = 1
    dx = h_az / n_az
    dt_dx = dt / dx
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
    # Горизонтальный участок до САОР
    part_x(n_3, f_3, l_3, G)
    # CАОР
    t_saor = saor_calc(T1[-1])
    for j in range(1, len(t_saor)):
        T1[i] = t_saor[j]
        i += 1
    # Опускной участок
    part_x(n_4, f_4, h_4, G)
    # Горизонтальный участок до АЗ
    part_x(n_5, f_5, l_5, G)

    
    time += dt
    T0 = T1[:]
    T1 = [T1[-1]] + [0] * N
print(T0)







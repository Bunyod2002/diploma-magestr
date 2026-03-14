
# ruff: noqa: F403, F405
import Data as dt
from math import exp
import functions as fc
import saor
import start
T0 = start.T[:dt.n_1 + dt.n_az + dt.n_2] + start.T[dt.n_1 + dt.n_az:dt.n_1 + dt.n_az + dt.n_2] * 2 + [0] * 51
h = 2.25
dtime = 500
time = 0.05
n = 0

# Цикл естественной циркуляции
# 1.Активная зона  2.Область до тягового участка  3. Тяговый участок до отметки естественной циркуляции
# 4. Горизонтальный участок до САОР  5.Теплообменник Фильда  6. Опускной участок 7. Горизонтальный участок до АЗ 
T1 = [623,15] + [0] * dt.N
def part_x(n, f, h, G):
    global i, T1
    dx = h / n
    dt_dx = dt.dt / dx
    for j in range(n):
        t_i = T0[i]  # T_i-1_k
        t_i_1 = T0[i-1]  # T_i_k
        r = fc.ro(t_i_1 - 273.15)
        t_k_1 = round(t_i + dt_dx * (G* dt.cp_Pb*(t_i_1 - t_i)) / (dt.cp_Pb * r * f), 3)
        T1[i] = t_k_1
        i += 1
        
while time < dtime:
    # Активная зона
    G = dt.Gpb
    Q_veg = 0.065 * dt.Q * (time ** (-0.2) - (time + 2592000) ** (-0.2))
    i = 1
    dx = dt.h_az / dt.n_az
    dt_dx = dt.dt / dx
    z = 0
    for j in range(dt.n_az):
        t_i = T0[i]  # T_i-1_k
        t_i_1 = T0[i-1]  # T_i_k
        r = fc.ro(t_i_1 - 273.15)
        t_k_1 = round(t_i + dt_dx * (G* dt.cp_Pb*(t_i_1 - t_i) + fc.ql(Q_veg, z) * dx) / (dt.cp_Pb * r * dt.f_az), 3)
        T1[i] = t_k_1
        h += dx
        i += 1
        z += dx
    # Область до тягового участка
    part_x(dt.n_1, dt.f_1, dt.h_1, G)
    # Тяговый участок
    part_x(dt.n_2, dt.f_2, dt.h_2, G)
    # Горизонтальный участок до САОР
    part_x(dt.n_3, dt.f_3, dt.l_3, G)
    # CАОР
    t_saor = saor.saor_calc(T1[-1], G/12)
    for j in range(1, len(t_saor)):
        T1[i] = t_saor[j]
        i += 1
    # Опускной участок
    part_x(dt.n_4, dt.f_4, dt.h_4, G)
    # Горизонтальный участок до АЗ
    part_x(dt.n_5, dt.f_5, dt.l_5, G)

    
    time += dt.dt
    T0 = T1[:]
    T1 = [T1[-1]] + [0] * dt.N
print(T0)







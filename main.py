# ruff: noqa: F403, F405
import Data as dt
from math import exp
import functions as fc
import saor
import start
import hydraulic as hc

T0 = start.start_calc(dt.Gpb, 773.15)
h = 2.25
dtime = 86400
time = 120
n = 0
# Цикл естественной циркуляции
# 1.Активная зона  2.Область до тягового участка  3. Тяговый участок до отметки естественной циркуляции
# 4. Горизонтальный участок до САОР  5.Теплообменник Фильда  6. Опускной участок 7. Горизонтальный участок до АЗ 
T1 = [T0[-1]] + [0] * dt.N
def part(n, f, h, G):
    global i, T1
    dx = h / n
    dt_dx = dt.dt / dx
    for j in range(n): 
        t_i = T0[i]  # T_i-1_k
        t_i_1 = T0[i-1]  # T_i_k
        r = fc.ro(t_i_1 - 273.15)
        t_k_1 = round(t_i + dt_dx * (G* dt.cp_Pb*(t_i_1 - t_i)) / (dt.cp_Pb * r * f), 3)
        T1[i] = t_k_1
        fc.kurrent(f, t_k_1, G, dx, dt.dt)
        i += 1
          
G = dt.Gpb    
while time < dtime:
    # Активная зона
    T0_old = T0.copy()
    T1_old = T1.copy()
    Q_veg = dt.Q
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
        fc.kurrent(dt.f_az, t_k_1, G, dx, dt.dt)
        h += dx
        i += 1
        z += dx
    # Область до тягового участка
    part(dt.n_1, dt.f_1, dt.h_1, G)
    # Тяговый участок
    part(dt.n_2, dt.f_2, dt.h_2, G)
    # Горизрнтальный участок до ПГ
    part(dt.n_3, dt.f_3, dt.l_3, G)
    # Парогенератор
    for j in range(dt.n_pg):
        t_i = T0[i]  # T_i-1_k
        t_i_1 = T1[i - 1]  # T_i_k
        r = fc.ro(t_i_1 - 273.15)
        alfa = fc.alphaPb(r, fc.f_pg, fc.dg_pg, G)
        t_k_1 = round(t_i + dt_dx * (G * dt.cp_Pb * (t_i_1 - t_i) + alfa * dt.s_pg * (dt.T_pg - t_i)) / (r * dt.f_pg * dt.cp_Pb), 3)
        T1[i] = t_k_1
        h -= dx
        i += 1
    # Подъемный участок через насос
    part(dt.n_4, dt.f_4, dt.h_4, G)
    # Горизонтальный участок до САОР
    part(dt.n_5, dt.f_5, dt.l_5, G)
    # CАОР
    t_saor = saor.saor_calc(T1[i-1], G/12)
    for j in range(len(t_saor)):
        T1[i] = t_saor[j]
        i += 1
    # Опускной участок
    part(dt.n_6, dt.f_6, dt.h_6, G)
    # Горизонтальный участок до АЗ
    part(dt.n_7, dt.f_7, dt.l_7, G)
    
    time += dt.dt
    T0 = T1[:]
    T1 = [T1[-1]] + [0] * dt.N








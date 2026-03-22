# ruff: noqa: F403, F405
import Data as dt
from math import exp
import functions as fc
import saor
import start
import hydraulic as hc

#from gui import create_default_temp_plot
#plotter = create_default_temp_plot()
T0 = start.start_calc(dt.Gpb, 773.15)
h = 2.25
dtime = 100
time = 60
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
step = 0  
G = dt.Gpb      
while time < dtime:
    # Активная зона
    T0_old = T0.copy()
    T1_old = T1.copy()
    Q_veg =  0.065 * dt.Q * (time ** (-0.2) - (time + 2592000) ** (-0.2))
    i = 1
    dx = dt.h_az / dt.n_az
    dt_dx = dt.dt / dx
    z = 0
    #plotter.push_point("AZ_in", time, T1[i-1])
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
    #plotter.push_point("AZ_out", time, T1[i-1])
    # Область до тягового участка
    part(dt.n_1, dt.f_1, dt.h_1, G)
    # Тяговый участок
    part(dt.n_2, dt.f_2, dt.h_2, G)
    # Горизонтальный участок до САОР
    part(dt.n_3, dt.f_3, dt.l_3, G)
    # CАОР
    t_saor = saor.saor_calc(T1[i-1], G/12)
    #plotter.push_point("PG_in", time, T1[i-1])
    for j in range(len(t_saor)):
        T1[i] = t_saor[j]
        i += 1
    #plotter.push_point("PG_out", time, T1[i-1])
    # Опускной участок
    part(dt.n_4, dt.f_4, dt.h_4, G)
    # Горизонтальный участок до АЗ
    part(dt.n_5, dt.f_5, dt.l_5, G)
    p = hc.p_full(G, T1)
    if abs(p[0] - p[1]) < 5:
        time += dt.dt
        T0 = T1[:]
        T1 = [T1[-1]] + [0] * dt.N
    else:
        G = G * (p[1] / p[0])
        T0 = T0_old
        T1 = T1_old
        print(G)
    T0 = T1[:]
    T1 = [T1[-1]] + [0] * dt.N
    #if step % 10 == 0:
        #plotter.redraw()
    #step += 1






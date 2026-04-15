# ruff: noqa: F403, F405
import Data as dt
from math import exp
import functions as fc
import saor
import start
import hydraulic as hc
from gui import create_default_temp_plot

T0 = start.start_calc(dt.Gpb, 350 + 273.15)
plotter = plotter = create_default_temp_plot()
plotter2 = create_default_temp_plot()
p = hc.pressure_balance(dt.Gpb, T0)
p_nasos = p[0] + p[1]
t_draw = [0, 0, 0, 0]
 
h = 2.25
dtime = 86400
time = 0
n = 0
# Цикл естественной циркуляции
# 1.Активная зона  2.Область до тягового участка  3. Тяговый участок до отметки естественной циркуляции
# 4. Горизонтальный участок до САОР  5.Теплообменник Фильда  6. Опускной участок 7. Горизонтальный участок до АЗ 
T1 = [T0[-1]] + [0] * dt.N
def part(n, f, h, G, dtt):
    global i, T1
    dx = h / n
    dt_dx = dtt / dx
    for j in range(n): 
        t_i = T0[i]  # T_i-1_k
        t_i_1 = T0[i-1]  # T_i_k
        r = fc.ro(t_i_1 - 273.15)
        t_k_1 = round(t_i + dt_dx * (G* dt.cp_Pb*(t_i_1 - t_i)) / (dt.cp_Pb * r * f), 3)
        T1[i] = t_k_1
        fc.kurrent(f, t_k_1, G, dx, dtt)
        i += 1          
G = dt.Gpb
Q_veg = dt.Q  
dtt = dt.dt 
step = 0 
while time < dtime:
    # Активная зона
    T0_old = T0.copy()
    T1_old = T1.copy()
    if time > 90:
        p_nasos = p_nasos * exp(-(time-90) / 60)
    if time > 100:
        Q_veg = fc.residual_power(time - 100)
    i = 1
    dx = dt.h_az / dt.n_az
    dt_dx = dtt / dx
    z = 0
    t_draw[0] = T1[i-1]
    for j in range(dt.n_az):
        t_i = T0[i]  # T_i-1_k
        t_i_1 = T0[i-1]  # T_i_k
        r = fc.ro(t_i_1 - 273.15)
        t_k_1 = round(t_i + dt_dx * (G* dt.cp_Pb*(t_i_1 - t_i) + fc.ql(Q_veg, z) * dx) / (dt.cp_Pb * r * dt.f_az), 3)
        T1[i] = t_k_1
        fc.kurrent(dt.f_az, t_k_1, G, dx, dtt)
        h += dx
        i += 1
        z += dx
    t_draw[1] = T1[i-1]
    # Область до тягового участка
    part(dt.n_1, dt.f_1, dt.h_1, G, dtt)
    # Тяговый участок
    part(dt.n_2, dt.f_2, dt.h_2, G, dtt)
    # Горизрнтальный участок до ПГ
    part(dt.n_3, dt.f_3, dt.l_3, G / 4, dtt)
    # Парогенератор
    dx = dt.h_pg / dt.n_pg
    dt_dx = dtt / dx
    if time > 90:
        part(dt.n_pg, dt.f_pg, dt.h_pg, G, dtt)
    else:
        for j in range(dt.n_pg):
            t_i = T0[i]  # T_i-1_k
            t_i_1 = T1[i - 1]  # T_i_k
            r = fc.ro(t_i_1 - 273.15)
            alfa = fc.alphaPb(r, dt.f_pg, dt.dg_pg, G)
            t_k_1 = round(t_i + dt_dx * (G * dt.cp_Pb * (t_i_1 - t_i) + alfa * dt.s_pg * (dt.T_pg - t_i)) / (r * dt.f_pg * dt.cp_Pb), 3)
            T1[i] = t_k_1
            h -= dx
            i += 1
    # Подъемный участок через насос
    part(dt.n_4, dt.f_4, dt.h_4, G / 4, dtt)
    # Горизонтальный участок до САОР
    part(dt.n_5, dt.f_5, dt.l_5, G / 4, dtt)
    # CАОР
    t_draw[2] = T1[i-1]
    t_saor = saor.saor_calc(T1[i-1], G/12)
    for j in range(len(t_saor)):
        T1[i] = t_saor[j]
        i += 1
    t_draw[3] = T1[i-1]
    # Опускной участок
    part(dt.n_6, dt.f_6, dt.h_6, G / 12, dtt)
    # Горизонтальный участок до АЗ
    part(dt.n_7, dt.f_7, dt.l_7, G / 4, dtt)
    
    p = hc.pressure_balance(G, T1)
    if abs(p_nasos - p[1] - p[0]) < 10:
        time += dtt
        T0 = T1[:]
        T1 = [T1[-1]] + [0] * dt.N
        plotter.push_point("AZ_in", time, t_draw[0])
        plotter.push_point("AZ_out", time, t_draw[1])
        plotter.push_point("PG_in", time, t_draw[2])
        plotter.push_point("PG_out", time, t_draw[3])
        plotter2.push_point("Flow_rate", time, G + 273.15)
        step += 1
    else:
        G = G * (abs(p_nasos - p[1])) / p[0]
        '''if time > 5:
            dtt = fc.new_dt(dt.f_az, T0_old[30], G, dt.h_az / dt.n_az)'''
        T0 = T0_old
        T1 = T1_old
        print(G)
    if step % 10 == 0:
        plotter.redraw()
        plotter2.redraw()
    







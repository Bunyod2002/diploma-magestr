# ruff: noqa: F403, F405
import Data as dt
import functions as fc
from math import pi, exp
import saor
import start
import hydraulic as hc
from gui import create_default_temp_plot
datas = start.start_calc(21000)
T0 = datas[0]

p = hc.pressure_balance(dt.Gpb, T0)
p_0 = p[0] - p[1]
p_nasos = p_0
din = fc.din_count()
#t_draw = [0, 0, 0, 0]
plotter = create_default_temp_plot()
 
h = 2.25
dtime = 10000
time = 20
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
        #fc.kurrent(f, t_k_1, G, dx, dtt)
        i += 1          
G = 21000
Q_veg = dt.Q  
dtt = dt.dt 
step = 0
t_f = datas[1]
while time < dtime:
    i = 1
    dx = dt.h_az / dt.n_az
    dt_dx = dtt / dx
    z = 0

    p_nasos = p_0 * exp(-(time) / 60)

    Q_veg = fc.residual_power(time - 10)
    plotter.push_point("AZ_in", time, T1[i-1])
    for j in range(dt.n_az):
        tc = t_f[j]
        t_w = tc[-1]
        if j == 5:
            plotter.push_point("Flow_rate", time, t_w)
        t_i = T0[i]  # T_i-1_k
        t_i_1 = T0[i-1]  # T_i_k
        r = fc.ro(t_i_1 - 273.15)
        alfa = fc.alphaPb(r, dt.f_az, dt.dg_az, G)
        t_k_1 = round(t_i + dt_dx * (G * dt.cp_Pb*(t_i_1 - t_i) + alfa * dt.s_tvel * (t_w - t_i)) / (dt.cp_Pb * r * dt.f_az), 3)
        T1[i] = t_k_1
        a_coef = []
        b_coef = []
        dr1 = (dt.r2 - dt.r1) / (dt.m - 1)
        dr2 = (dt.r4 - dt.r3) / (dt.n - 1)
        sr1 = dt.at * dtt / (dr1 ** 2)
        sr2 = dt.ac * dtt / (dr2 ** 2)
        qv = fc.ql(Q_veg, z) / (pi * (dt.r2 ** 2 - dt.r1**2))
        for k in range(9):
            if k == 0:
                ri14 = dt.r1 + dr1 / 4
                ri12 = dt.r1 + dr1 / 2
                a = 1 / (1 + (ri14 / (2 * ri12 * sr1)))
                b = tc[k] * ri14 / (2 * ri12 * sr1) + qv * dr1 ** 2 * ri14 / (2* ri12 * dt.clt)
                a_coef.append(a)
                b_coef.append(b)
            elif 0 < k < 4:
                ri = dt.r1 + dr1 * k
                ri1 = ri - dr1 / 2
                ri2 = ri + dr1 / 2
                a = 1 / (ri * (2 + 1 / sr1) / ri2 - ri1 * a_coef[k-1] / ri2)
                b = ri * tc[k] / (ri2 * sr1) + qv * dr1 ** 2 * ri / (ri2 * dt.clt) + ri1 * a_coef[k-1] * b_coef[k-1] / ri2
                a_coef.append(a)
                b_coef.append(b) 
            elif k == 4:
                alfgap1 = dt.clg / dt.delta
                alfgap2 = dt.sigma * dt.eps1 * (tc[k] + tc[k+1]) * (tc[k] ** 2 + tc[k + 1] ** 2)
                alfgap = alfgap1 + alfgap2
                ri = dt.r2 
                ri14 = ri - dr1 / 4
                ri12 = ri - dr1 / 2
                a = 1 / (1 + dt.clt * ri14 / (dr1 * alfgap * sr1 * 2 * ri) + dt.clt * ri12 * (1- a_coef[k-1]) / (dr1 * alfgap * ri))
                b = dt.clt * ri14 * tc[k] / (dr1 * alfgap * sr1 * 2 * ri) + qv * dr1 * ri14 / (alfgap * 2 * ri) + dt.clt * ri12 * a_coef[k-1] * b_coef[k-1] / (dr1 * alfgap * ri)
                a_coef.append(a)
                b_coef.append(b) 
            elif k == 5:
                ri = dt.r3
                ri14 = ri + dr2 / 4
                ri12 = ri + dr2 / 2
                a = 1 / (1 + (ri14 / (2 * sr2 * ri12)) + alfgap * dr2 * ri * (1 - a_coef[k-1]) / (dt.clc * ri12))
                b = ri14 * tc[k] / (sr2 * 2 * ri12) + alfgap * dr2 * ri * a_coef[k-1] * b_coef[k-1] / (dt.clc * ri12)
                a_coef.append(a)
                b_coef.append(b) 
            else:
                ri = dt.r3 + dr2 * (k - 5)
                ri1 = ri - dr2 / 2
                ri2 = ri + dr2 / 2
                a = 1 / (ri * (2 + 1 / sr2) / ri2 - ri1 * a_coef[k-1] / ri2)
                b = ri * tc[k] / (ri2 * sr2) + ri1 * a_coef[k-1] * b_coef[k-1] / ri2
                a_coef.append(a)
                b_coef.append(b)
            
        for k in range(len(tc) - 1, -1, -1):
            if k == len(tc) - 1:
                ri = dt.r4
                ri14 = ri - dr2 / 4
                ri12 = ri - dr2 / 2
                tc[k] = (tc[k] + alfa * t_i * sr2 * dr2 * 2 * ri / (dt.clc * ri14) + sr2 * 2 * ri12 * a_coef[k-1] * b_coef[k-1] / ri14) / (1 + sr2 * 2 * ri12 * (1 - a_coef[k-1]) / ri14 + alfa * sr2 * dr2 * 2 * ri / (dt.clc * ri14))
            else:
                tc[k] = a_coef[k] * (b_coef[k] + tc[k + 1])
        t_f[j] = tc   
        #fc.kurrent(dt.f_az, t_k_1, G, dx, dtt)
        h += dx
        i += 1
        z += dx
    plotter.push_point("AZ_out", time, T1[i-1])
    # Область до тягового участка
    part(dt.n_1, dt.f_1, dt.h_1, G, dtt)
    # Тяговый участок
    part(dt.n_2, dt.f_2, dt.h_2, G, dtt)
    # Горизрнтальный участок до ПГ
    part(dt.n_3, dt.f_3, dt.l_3, G / 4, dtt)
    # Парогенератор
    dx = dt.h_pg / dt.n_pg
    dt_dx = dtt / dx
    plotter.push_point("PG_in", time, T1[i-1])
    if time > 10:
        part(dt.n_pg, dt.f_pg, dt.h_pg, G / 4, dtt)
    else:
        for j in range(dt.n_pg):
            t_i = T0[i]  # T_i-1_k
            t_i_1 = T0[i - 1]  # T_i_k
            r = fc.ro(t_i_1 - 273.15)
            alfa = fc.alphaPb(r, dt.f_pg, dt.dg_pg, G / 4)
            t_k_1 = round(t_i + dt_dx * (G / 4 * dt.cp_Pb * (t_i_1 - t_i) + alfa * dt.s_pg * (dt.T_pg - t_i)) / (r * dt.f_pg * dt.cp_Pb), 3)
            T1[i] = t_k_1
            h -= dx
            i += 1
    
    # Подъемный участок через насос
    part(dt.n_4, dt.f_4, dt.h_4, G / 4, dtt)
    # Горизонтальный участок до САОР
    part(dt.n_5, dt.f_5, dt.l_5, G / 4, dtt)
    # CАОР
    t_saor = saor.saor_calc(T1[i-1], G/12)
    for j in range(len(t_saor)):
        T1[i] = t_saor[j]
        i += 1
    plotter.push_point("PG_out", time, T1[i-1])
    # Опускной участок
    part(dt.n_6, dt.f_6, dt.h_6, G / 12, dtt)
    # Горизонтальный участок до АЗ
    part(dt.n_7, dt.f_7, dt.l_7, G / 4, dtt)
    if step % 10 == 0:
        plotter.redraw()
    p = hc.pressure_balance(G, T1)
    time += dtt
    G += dtt * (p_nasos + p[1] - p[0]) / din
    dtt = fc.new_dt(dt.f_az, T1[30], G, dt.h_az / dt.n_az)
    T0 = T1[:]
    T1 = [T1[-1]] + [0] * dt.N
    step += 1
    print(G, time, p_nasos)
plotter.hold()








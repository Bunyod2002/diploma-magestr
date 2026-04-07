import Data as dt
from math import exp
import functions as fc
import saor

def start_calc(G, T):
    T0 = [0] * (dt.N + 1)
    dtime = 600
    time = 0
    n = 0
    # Цикл естественной циркуляции
    # 1.Активная зона  2.Область до тягового участка  3. Тяговый участок до отметки естественной циркуляции
    # 4. Горизонтальный участок до САОР  5.Теплообменник Фильда  6. Опускной участок 7. Горизонтальный участок до АЗ 
    T1 = [T] + [0] * dt.N
    h = 2.25
    def part_x(n, f, h, G, i):
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
        return i
            
    while time < dtime:
        # Активная зона
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
        i = part_x(dt.n_1, dt.f_1, dt.h_1, G, i)
        # Тяговый участок
        i = part_x(dt.n_2, dt.f_2, dt.h_2, G, i)
        # Горизрнтальный участок до ПГ
        i = part_x(dt.n_3, dt.f_3, dt.l_3, G, i)
        # Парогенератор
        dx = dt.h_pg / dt.n_pg
        dt_dx = dt.dt / dx
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
        i = part_x(dt.n_4, dt.f_4, dt.h_4, G, i)
        # Горизонтальный участок до САОР
        i = part_x(dt.n_5, dt.f_5, dt.l_5, G, i)
        # CАОР
        i = part_x(dt.n_saor, dt.d_test, dt.h_saor, G, i)
        '''t_saor = saor.saor_calc(T1[i-1], G/12)
        for j in range(len(t_saor)):
            T1[i] = t_saor[j]
            i += 1'''
        # Опускной участок
        i = part_x(dt.n_6, dt.f_6, dt.h_6, G, i)
        # Горизонтальный участок до АЗ
        i = part_x(dt.n_7, dt.f_7, dt.l_7, G, i)
    
        
        time += dt.dt
        T0 = T1[:]
        T1 = [T] + [0] * dt.N

    return T0

a = start_calc(23761, 350+273.15)        

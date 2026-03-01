from Data import *
from functions import * 

h = 2.25
T_k = [0] * (N + 1)
dtime = 250
time = 0
n = 0
T_k_1 = [T0] + [0] * N
def part(n, f, h):
    global T_k_1, i
    dx = h / n
    dt_dx = dt / dx
    for j in range(n):
        t_i = T_k[i]  # T_i-1_k
        t_i_1 = T_k[i-1]  # T_i_k
        r = ro(t_i_1 - 273.15)
        t_k_1 = round(t_i + dt_dx * (Gpb* cp_Pb*(t_i_1 - t_i)) / (cp_Pb * r * f), 3)
        T_k_1[i] = t_k_1
        i += 1

while time < dtime:
    # Активная зона
    i = 1
    dx = h_az / n_az
    dt_dx = dt / dx
    h = 2.25
    z = 0
    for j in range(n_az):
        t_i = T_k[i]  # T_i-1_k
        t_i_1 = T_k[i-1]  # T_i_k
        r = ro(t_i_1 - 273.15)
        t_k_1 = round(t_i + dt_dx * (Gpb* cp_Pb*(t_i_1 - t_i) + ql(Q, z) * dx) / (cp_Pb * r * f_az), 3)
        T_k_1[i] = t_k_1
        h += dx
        i += 1
        z += dx
    # Область до тягового участка
    part(n_1, f_1, h_1)
    # Тяговый участок
    part(n_2, f_2, h_2)
    # Парогенератор
    h = 7.75
    dx = h_pg / n_pg
    dt_dx = dt / dx
    for j in range(n_pg):
        t_i = T_k[i]  # T_i-1_k
        t_i_1 = T_k[i - 1]  # T_i_k
        r = ro(t_i_1 - 273.15)
        alfa = alphaPb(r, f_pg, dg_pg)
        t_k_1 = round(t_i + dt_dx * (Gpb * cp_Pb * (t_i_1 - t_i) + alfa * s_pg * (T_pg - t_i)) / (r * f_pg * cp_Pb), 3)
        T_k_1[i] = t_k_1
        h -= dx
        i += 1
    # Вертикальный участок с ГЦН
    part(n_3, f_3, h_3)
    # Опускной участок
    part(n_4, f_4, h_4)
    # Горизонтальный участок
    part(n_5, f_5, l_5)
    # Подъемный участок до активной зоны
    part(n_6, f_6, h_6)
    
    time += dt
    T_k = T_k_1
    T_k_1 = [T0] + [0] * N
    
T = T_k
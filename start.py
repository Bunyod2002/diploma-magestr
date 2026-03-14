import start_data as st
import functions as fc

h = 2.25
T_k = [0] * (st.N + 1)
dtime = 250
time = 0
n = 0
T_k_1 = [st.T0] + [0] * st.N
def part(n, f, h):
    global T_k_1, i
    dx = h / n
    dt_dx = st.dt / dx
    for j in range(n):
        t_i = T_k[i]  # T_i-1_k
        t_i_1 = T_k[i-1]  # T_i_k
        r = fc.ro(t_i_1 - 273.15)
        t_k_1 = round(t_i + dt_dx * (st.Gpb* st.cp_Pb*(t_i_1 - t_i)) / (st.cp_Pb * r * f), 3)
        T_k_1[i] = t_k_1
        i += 1

while time < dtime:
    # Активная зона
    i = 1
    dx = st.h_az / st.n_az
    dt_dx = st.dt / dx
    h = 2.25
    z = 0
    for j in range(st.n_az):
        t_i = T_k[i]  # T_i-1_k
        t_i_1 = T_k[i-1]  # T_i_k
        r = fc.ro(t_i_1 - 273.15)
        t_k_1 = round(t_i + dt_dx * (st.Gpb* st.cp_Pb*(t_i_1 - t_i) + fc.ql(st.Q, z) * dx) / (st.cp_Pb * r * st.f_az), 3)
        T_k_1[i] = t_k_1
        h += dx
        i += 1
        z += dx
    # Область до тягового участка
    part(st.n_1, st.f_1, st.h_1)
    # Тяговый участок
    part(st.n_2, st.f_2, st.h_2)
    # Парогенератор
    h = 7.75
    dx = st.h_pg / st.n_pg
    dt_dx = st.dt / dx
    for j in range(st.n_pg):
        t_i = T_k[i]  # T_i-1_k
        t_i_1 = T_k[i - 1]  # T_i_k
        r = fc.ro(t_i_1 - 273.15)
        alfa = fc.alphaPb(r, st.f_pg, st.dg_pg, st.Gpb)
        t_k_1 = round(t_i + dt_dx * (st.Gpb * st.cp_Pb * (t_i_1 - t_i) + alfa * st.s_pg * (st.T_pg - t_i)) / (r * st.f_pg * st.cp_Pb), 3)
        T_k_1[i] = t_k_1
        h -= dx
        i += 1
    # Вертикальный участок с ГЦН
    part(st.n_3, st.f_3, st.h_3)
    # Опускной участок
    part(st.n_4, st.f_4, st.h_4)
    # Горизонтальный участок
    part(st.n_5, st.f_5, st.l_5)
    # Подъемный участок до активной зоны
    part(st.n_6, st.f_6, st.h_6)
    
    time += st.dt
    T_k = T_k_1
    T_k_1 = [st.T0] + [0] * st.N
    
T = T_k

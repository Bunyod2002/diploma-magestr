# TODO 1 переписать программу start, сделав его стационарным расчетом всего контура с заданной температурой входа
# таким образом получить распределения температур по всему контуру в первом приближении
import Data as dt
import functions as fc
import saor

G = 500
T = [350 + 273.5]
t_cur = 350 + 273.5
Q = 0.065 * dt.Q
z = 0
dx = dt.h_az / dt.n_az
for j in range(dt.n_az):
    t_temp = t_cur + fc.ql(Q, z) * dx / (G * dt.cp_Pb)
    T.append(t_temp)
    t_cur = t_temp
    z += dx

T += [T[-1]] * 20
    
t_saor = saor.saor_calc(T[-1], G / 12)
T += t_saor
T += [T[-1]] * 15
print(T)
    
        


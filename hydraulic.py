import functions as fc
from Data import g
import Data as dt
import saor_data as sdt
from math import log10, pi
def dp_friction(G, T1, x, f, d = 0):
    if d == 0:
        d = (4 * f / pi) ** 0.5
    r = fc.ro(T1 - 273.15)
    w = fc.vel(f, T1, G)
    v = (43.8 - 7.57 * 0.01 * (T1 - 273.15) + 0.467 * 0.0001 * (T1 - 273.15)**2) * 10 **(-8)  # кинематическая вязкость
    Re = w * d / v
    if Re < 2300 and Re > 0:
        e_fric = 64 / Re
    elif Re > 0:
        e_fric = (1.82*log10(Re) - 1.64)**(-2) # коэффициент потери на трение
    else:
        e_fric = 0
    dp_fric = e_fric*x*r*w*w*0.5/d
    return dp_fric

def dp_polez(T, x, direction=None):
    r = fc.ro(T - 273.15)
    if direction == 'up':
        dp = r * g * x
    elif direction == 'down':
        dp = r * g * x
    else:
        dp = 0
    return dp

def dp_friction_saor(G, T1):
    x = dt.h_saor / dt.n_saor
    d_g = sdt.d5 - sdt.d4
    f = pi * (sdt.d5 ** 2 - sdt.d4 ** 2) / 4
    r = fc.ro(T1 - 273.15)
    w = fc.vel(G, f, T1)
    v = (43.8 - 7.57 * 0.01 * (T1 - 273.15) + 0.467 * 0.0001 * (T1 - 273.15)**2) * 10 **(-8)  # кинематическая вязкость
    Re = w * d_g / v
    if Re < 2300 and Re > 0:
        e_fric = 64 / Re
    elif Re > 0:
        e_fric = (1.82*log10(Re) - 1.64)**(-2) # коэффициент потери на трение
    else:
        e_fric = 0
    dp_fric = e_fric*x*r*w*w*0.5/d_g
    return dp_fric

def p_full(T):
    i = 0
    p1 = 0
    p2 = 0
    for j in range(dt.n_az):
        p1 += dp_friction(dt.Gpb, T[i], dt.h_az / dt.n_az, dt.f_az, dt.dg_az)
        p1 += dp_polez(T[i], dt.h_az / dt.n_az, 'up')
        i += 1
    for j in range(dt.n_1):
        p1 += dp_friction(dt.Gpb, T[i], dt.h_1 / dt.n_1, dt.f_1)
        p1 += dp_polez(T[i], dt.h_1 / dt.n_1, 'up')
        i += 1
    for j in range(dt.n_2):
        p1 += dp_friction(dt.Gpb, T[i], dt.h_2 / dt.n_2, dt.f_2)
        p1 += dp_polez(T[i], dt.h_2 / dt.n_2, 'up')
        i += 1
    for j in range(dt.n_3):
        p1 += dp_friction(dt.Gpb, T[i], dt.l_3 / dt.n_3, dt.f_3, dt.dg_az)
        i += 1
        # Здесь должен быть участок САОР
    for j in range(dt.n_saor):
        p1 += dp_friction_saor(dt.Gpb / 12, T[i])
        p2 += dp_polez(T[i], dt.h_saor / dt.n_saor, 'down')
        i += 1
    for j in range(dt.n_4):
        p1 += dp_friction(dt.Gpb, T[i], dt.h_4 / dt.n_4, dt.f_4)
        p2 += dp_polez(T[i], dt.h_4 / dt.n_4, 'down')
        i += 1
    for j in range(dt.n_5):
        p1 += dp_friction(dt.Gpb, T[i], dt.l_5 / dt.n_5, dt.f_5)
        i += 1
    return (p1, p2)
    
    
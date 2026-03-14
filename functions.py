# ruff: noqa: F403, F405
import Data as dt
from math import pi, e
from math import log10, cos

def mass(ro, x, f):  # функция считает массу свинца
    return round(ro * x * f, 3)

def alphaPb(ro, f, d, G):  # альфа свинца в зависимости от плотности
    vPb = 4 * G / (ro * f)  # скорость свинца
    PePb = vPb * d / dt.aPb  # число Пекле
    Nu = 5 + 0.025 * PePb ** 0.8  # число Нуссельта
    alfaPb = Nu * dt.lambdaPb / d
    return round(alfaPb, 3)

def ro(t: float):   # плотность свинца в зависимости от температуры в градусах
    return round(1000*(11.05 - 12.49 * t * 0.0001), 3)

def ql(Q, z):
    kz = 1.5                  # максимум = 1.5·среднего
    ql_avg = Q / dt.h_az         # средняя линейная мощность
    a = kz - 1.0              # a=0.5
    x = (0.5*dt.h_az - z) / (0.5*dt.h_az)
    Z = 1.0 + a * cos(pi * x) # среднее(Z)=1 → ∫ ql = Q
    return ql_avg * Z         # Вт/м

def vel(f, T, G):
    return G / (ro(T - 273.15) * f)


def dp(T1, x, f, direction, d = 0):
    if d == 0:
        d = (4 * f / pi) ** 0.5
    r = ro(T1 - 273.15)
    w = vel(f, T1)
    v = (43.8 - 7.57 * 0.01 * (T1 - 273.15) + 0.467 * 0.0001 * (T1 - 273.15)**2) * 10 **(-8)  # кинематическая вязкость
    Re = w * d / v
    if Re < 2300 and Re > 0:
        e_fric = 64 / Re
    elif Re > 0:
        e_fric = (1.82*log10(Re) - 1.64)**(-2) # коэффициент потери на трение
    else:
        e_fric = 0
    dp_fric = e_fric*x*r*w*w*0.5/d
    if direction == 'up':
        dp_polez = r * g * x
    elif direction == 'down':
        dp_polez = -r * g * x
    else:
        dp_polez = 0
    return dp_fric + dp_polez


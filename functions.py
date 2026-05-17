# ruff: noqa: F403, F405
import Data as dt
from math import pi, e
from math import log10, cos

def mass(ro, x, f):  # функция считает массу свинца
    return round(ro * x * f, 3)

def alphaPb(ro, f, d, G):  # альфа свинца в зависимости от плотности
    """На вход подается плотность, проходное сечение, гидравлический диаметр, расход"""
    vPb = 4 * G / (ro * f)  # скорость свинца
    PePb = vPb * d / dt.aPb  # число Пекле
    Nu = 5 + 0.025 * PePb ** 0.8  # число Нуссельта
    alfaPb = Nu * dt.lambdaPb / d
    return alfaPb

def kurrent(f, T, G, x, dt): # условие куррента
    v = vel(f, T, G) # считаем скорость
    if not x / v > dt:
        raise ValueError('слишком большой шаг по времени')

def new_dt(f, T, G, x):
    v = vel(f, T, G)
    return 0.9 * x / v

def ro(t: float):   # плотность свинца в зависимости от температуры в градусах
    return round(1000*(11.05 - 12.49 * t * 0.0001), 3)

def ql(Q, z):
    kz = 1.5                  # максимум = 1.5·среднего
    ql_avg = Q / dt.h_az         # средняя линейная мощность
    a = kz - 1.0              # a=0.5
    x = (0.5*dt.h_az - z) / (0.5*dt.h_az)
    Z = 1.0 + a * cos(pi * x) # среднее(Z)=1 → ∫ ql = Q
    return ql_avg * Z        # Вт/м

def vel(f, T, G): #скорость свинца, м^2, K, кг/c
    return G / (ro(T - 273.15) * f)

def residual_power(time): # остаточное тепловыделение в зависимости от времени
    return 0.065 * dt.Q * (time ** (-0.2) - (time + 2592000.0) ** (-0.2))

def din_count():
    din = dt.h_az / dt.f_az + dt.h_1 / dt.f_1 + dt.h_2 / dt.f_2 + dt.l_3 / dt.f_3 + dt.h_pg / dt.f_pg + dt.h_4 / dt.f_4 + dt.l_5 / dt.f_5 + dt.h_saor / dt.f_saor + dt.h_6 / dt.f_6 + dt.l_7 / dt.f_7
    return din


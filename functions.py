# ruff: noqa: F403, F405
from Data import *  
from math import pi, e
from math import log10, cos

def mass(ro, x, f):  # функция считает массу свинца
    return round(ro * x * f, 3)

def alphaPb(ro, f, d, G=Gpb):  # альфа свинца в зависимости от плотности
    vPb = 4 * G / (ro * f)  # скорость свинца
    PePb = vPb * d / aPb  # число Пекле
    Nu = 5 + 0.025 * PePb ** 0.8  # число Нуссельта
    alfaPb = Nu * lambdaPb / d
    return round(alfaPb, 3)

def ro(t: float):   # плотность свинца в зависимости от температуры в градусах
    return round(1000*(11.05 - 12.49 * t * 0.0001), 3)

def ql(Q, z):
    kz = 1.5                  # максимум = 1.5·среднего
    ql_avg = Q / h_az         # средняя линейная мощность
    a = kz - 1.0              # a=0.5
    x = (0.5*h_az - z) / (0.5*h_az)
    Z = 1.0 + a * cos(pi * x) # среднее(Z)=1 → ∫ ql = Q
    return ql_avg * Z         # Вт/м

def vel(G, f, T):
    return G / (ro(T - 273.15) * f)



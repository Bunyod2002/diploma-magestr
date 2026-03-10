# ruff: noqa: F403, F405
from functions import *
from math import log10
from Data import g

def dp_friction(G, T1, x, f, d = 0):
    if d == 0:
        d = (4 * f / pi) ** 0.5
    r = ro(T1 - 273.15)
    w = vel(G, f, T1)
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
    r = ro(T - 273.15)
    if direction == 'up':
        dp = r * g * x
    elif direction == 'down':
        dp = - r *g * x
    else:
        dp = 0
    return dp


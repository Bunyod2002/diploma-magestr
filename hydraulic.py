import functions as fc
from Data import g
import Data as dt
import saor_data as sdt
from math import log10, pi

def friction_koef(re):
    """Определяет коэффициент сопротивления на трение, на вход принимается число Рейнольдса"""
    if re <= 0.0:
        return 0.0
    if re < 2300.0:
        return 64.0 / re
    return (1.82 * log10(re) - 1.64) ** (-2)

def dp_friction(G, T, x, f, d = None):
    """Определяет потери на трение на определеноом участке
    В качестве входных данных нужен расход, температура
    длина участка, площадь участка и гидравлический диаметр"""
    if not d:
        d = (4 * f / pi) ** 0.5
    r = fc.ro(T - 273.15)
    w = fc.vel(f, T, G)
    v = (43.8 - 7.57 * 0.01 * (T - 273.15) + 0.467 * 0.0001 * (T - 273.15)**2) * 10 **(-8)  # кинематическая вязкость свинца
    Re = w * d / v # Число Рейнольдса
    e_fric = friction_koef(Re)
    dp_fric = e_fric * (x / d) * r * w * w * 0.5
    return dp_fric

def dp_gravity(T, x, direction=None):
    """Возвращает гравитационный напор, со знаком + если движение вверх
    со знаком - если движение вниз"""
    r = fc.ro(T - 273.15)
    if direction == 'up':
        dp = r * g * x
    elif direction == 'down':
        dp = r * g * x
    else:
        dp = 0
    return dp

def dp_friction_saor(G, T1):
    """Отдельная функция, которая считает сопротивление на трении
    для участка САОР"""
    x = dt.h_saor / dt.n_saor
    d_g = sdt.d5 - sdt.d4
    area = pi * (sdt.d5 ** 2 - sdt.d4 ** 2) / 4
    return dp_friction(G, T1, x, area, d_g)


def pressure_balance(G: float, T: list[float]):
    """Возвращает (потери, движущий напор, невязка)."""
    i = 0
    losses = 0.0
    driving = 0

    # Подъемная ветвь
    dx = dt.h_az / dt.n_az
    for _ in range(dt.n_az):
        losses += dp_friction(G, T[i], dx, dt.f_az, dt.dg_az)
        losses += dp_gravity(T[i], dx, 'up')
        i += 1

    dx = dt.h_1 / dt.n_1
    for _ in range(dt.n_1):
        losses += dp_friction(G, T[i], dx, dt.f_1)
        losses += dp_gravity(T[i], dx, 'up')
        i += 1

    dx = dt.h_2 / dt.n_2
    for _ in range(dt.n_2):
        losses += dp_friction(G, T[i], dx, dt.f_2)
        losses += dp_gravity(T[i], dx, 'up')
        i += 1

    # Горизонталь до САОР
    dx = dt.l_3 / dt.n_3
    for _ in range(dt.n_3):
        losses += dp_friction(G, T[i], dx, dt.f_3, dt.dg_az)
        i += 1

    # САОР — нисходящая ветвь
    dx = dt.h_saor / dt.n_saor
    for _ in range(dt.n_saor):
        losses += dp_friction_saor(G / 12, T[i])
        driving += dp_gravity(T[i], dx, 'down')
        i += 1

    # Опускной участок
    dx = dt.h_4 / dt.n_4
    for _ in range(dt.n_4):
        losses += dp_friction(G, T[i], dx, dt.f_4)
        driving += dp_gravity(T[i], dx, 'down')
        i += 1

    # Горизонталь до АЗ
    dx = dt.l_5 / dt.n_5
    for _ in range(dt.n_5):
        losses += dp_friction(G, T[i], dx, dt.f_5)
        i += 1

    return losses, driving

    
    
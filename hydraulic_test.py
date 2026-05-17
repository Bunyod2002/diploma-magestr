import functions as fc
from Data import g
import Data as dt
import saor_data as sdt
from math import log10, pi
import start
T = start.start_calc(23761)
# Гидравлические сопротивления активной зоны
def dp_friction_az(G, T):
    """Определяет потери на трение в активной зоне"""
    x = dt.h_az / dt.n_az
    d_g = dt.dg_az
    area = dt.f_az
    return dp_friction(G, T, x, area, d_g)

def dp_az_pg_in_and_out(G, T, koef):
    """Определяет потери на местные сопротивления на входе и выходе из активной зоны"""
    r = fc.ro(T - 273.15)
    w = fc.vel(dt.f_az, T, G)
    return koef * r * w * w * 0.5

def dp_az_spacer(G, T):
    """Определяет потери в АЗ для дистационирующих решеток"""
    rho = fc.ro(T - 273.15)
    w = G / (rho * dt.f_az)
    v = (43.8 - 7.57 * 0.01 * (T - 273.15) + 0.467 * 0.0001 * (T - 273.15)**2) * 10 **(-8)  # кинематическая вязкость свинца
    Re = w * dt.dg_az / v # Число Рейнольдса
    eps_spacer = 0.2
    cv = min(
        3.5 + 73.14 / (Re ** 0.264) + 2.79e10 / (Re ** 2.79),
        2.0 / (eps_spacer ** 2),
    )
    x = cv * (eps_spacer ** 2)
    return x * 0.5 * rho * w * w

# Гидравлические сопротивления для змеевика в парогенераторе
def dp_bundle_pg(G, T, x):
    """Определяет потери при оптекании змеевика в парогенераторе"""
    r = fc.ro(T - 273.15)
    w = fc.vel(dt.f_pg, T, G)
    v = (43.8 - 7.57 * 0.01 * (T - 273.15) + 0.467 * 0.0001 * (T - 273.15)**2) * 10 **(-8)  # кинематическая вязкость свинца
    Re = w * dt.dg_pg / v # Число Рейнольдса
    eu = 0.55 + 1800.0 / Re + 2.0e6 / (Re ** 2)
    return eu * (x / dt.dg_pg) * r * w * w * 0.5


# Гидравлические сопротивления для обычных труб
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

def dp_gravity(T, x):
    """Возвращает гравитационный напор, со знаком + если движение вверх
    со знаком - если движение вниз"""
    r = fc.ro(T - 273.15)
    return r * g * x

# Гидравлические сопротивления для САОР
def dp_friction_saor(G, T1):
    """Отдельная функция, которая считает сопротивление на трении
    для участка САОР"""
    x = dt.h_saor / dt.n_saor
    d_g = sdt.d5 - sdt.d4
    area = pi * (sdt.d5 ** 2 - sdt.d4 ** 2) / 4
    return dp_friction(G, T1, x, area, d_g)

# Местные сопротивления
def dp_expansion(G, T, A1, A2):
    """Внезапное расширение A1 -> A2 (A2 > A1)"""
    rho = fc.ro(T - 273.15)
    w = G / (rho * A1)
    beta = A1 / A2
    xi = (1 - beta) ** 2
    return xi * 0.5 * rho * w**2

def dp_contraction(G, T, A1, A2):
    """Внезапное сужениеA1 -> A2 (A2 < A1)"""
    rho = fc.ro(T - 273.15)
    w = G / (rho * A2)
    beta = A2 / A1
    xi = 0.5 * (1 / beta - 1)
    return xi * 0.5 * rho * w**2

def dp_bend_90(G, T, A, xi=0.5):
    """Поворот на 90 градусов"""
    rho = fc.ro(T - 273.15)
    w = G / (rho * A)
    return xi * 0.5 * rho * w**2

def dp_bend_180(G, T, A, xi=1.5):
    """Разворот на 180 градусов"""
    rho = fc.ro(T - 273.15)
    w = G / (rho * A)
    return xi * 0.5 * rho * w**2

# Общие гидравлическое сопротивление для всего контура
def pressure_balance(G: float, T: list[float]):
    """Возвращает (потери, движущий напор, невязка)."""
    p_loses = 0
    p_driving = 0
    i = 0
    # Подъемная ветвь
    dx = dt.h_az / dt.n_az
    p_loses += dp_az_pg_in_and_out(G, T[i], 0.5)
    for j in range(dt.n_az):
        p_loses += dp_friction_az(G, T[i])
        #p_loses += dp_gravity(T[i], dx)
        i += 1


    dx = dt.h_1 / dt.n_1
    for _ in range(dt.n_1):
        p_loses += dp_friction(G, T[i], dx, dt.f_1)
        #p_loses += dp_gravity(T[i], dx)
        i += 1
        

    dx = dt.h_2 / dt.n_2
    for _ in range(dt.n_2):
        p_loses += dp_friction(G, T[i], dx, dt.f_2)
        #p_loses += dp_gravity(T[i], dx)
        i += 1
        
    p_loses += dp_bend_90(G / 4, T[i], dt.f_2)

    dx = dt.l_3 / dt.n_3
    for _ in range(dt.n_3):
        p_loses += dp_friction(G / 4, T[i], dx, dt.f_3)
        i += 1
    dx = dt.h_pg / dt.n_pg
    for _ in range(dt.n_pg):
        p_loses += dp_bundle_pg(G / 4, T[i], dx)
        #p_driving += dp_gravity(T[i], dx)
        i += 1  

    dx = dt.h_4 / dt.n_4
    for _ in range(dt.n_4):
        p_loses += dp_friction(G / 4, T[i], dx, dt.f_4)
        #p_loses += dp_gravity(T[i], dx)
        i += 1 
    dx = dt.l_5 / dt.n_5

    for _ in range(dt.n_5):
        p_loses += dp_friction(G / 4, T[i], dx, dt.f_5)
        i += 1

    # САОР — нисходящая ветвь
    dx = dt.h_saor / dt.n_saor
    for _ in range(dt.n_saor):
        p_loses += dp_friction_saor(G / 12, T[i])
        #p_driving += dp_gravity(T[i], dx)
        i += 1
    p_loses += dp_expansion(G / 12, T[i], pi * (sdt.d5 ** 2 - sdt.d4 ** 2) / 4, dt.f_6)
    # Опускной участок
    dx = dt.h_6 / dt.n_6
    for _ in range(dt.n_6):
        #p_driving += dp_gravity(T[i], dx)
        p_loses += dp_friction(G / 4, T[i], dx, dt.f_6)
        i += 1
    # Горизонталь до АЗ
    dx = dt.l_7 / dt.n_7
    for _ in range(dt.n_7):
        p_loses += dp_friction(G / 4, T[i], dx, dt.f_7)
        i += 1

    return p_loses, p_driving


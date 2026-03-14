# ruff: noqa: F403, F405
from saor_data import *
from saor_func import *
def saor_calc(t_pb):
    T2 = 380.0 + 273.15
    T20 = 380.0 + 273.15
    while True:
        lst_t3 = []
        h1 = 7.0
        h2 = 5.0
        h = 12.0
        T1 = 39.0 + 273.15  # Температура входа воздуха
        T3 = t_pb + 273.15 # Температура входа свинца
        TAr = 600.0 + 273.15
        Tw2 = T3 - 122
        Tw1 = Tw2 - 220

        while h2 > 0:
            alpha21 = alfa21(T2)
            alpha22 = alfa22(T2)
            alpha1 = alfa1(T1)
            Tw1 = newton(Tw1, Tw2, T1, T2, TAr, alpha1, alphaAr, alpha21, alpha22)[0]
            Tw2 = newton(Tw1, Tw2, T1, T2, TAr, alpha1, alphaAr, alpha21, alpha22)[1]
            T1 = t1(T1, Tw2, Tw1, alpha21, T2, X2)
            T2 = t2(T2, alpha21, Tw1, alpha22, Tw2, X2)
            h -= X2
            h2 -= X2

        Tw2 = T3 - 122
        Tw1 = Tw2 - 220

        while True:
            alpha21 = alfa21(T2)
            alpha22 = alfa22(T2)
            alpha1 = alfa1(T1)
            Tw1 = newton(Tw1, Tw2, T1, T2, T3, alpha1, alphaPb, alpha21, alpha22)[0]
            Tw2 = newton(Tw1, Tw2, T1, T2, T3, alpha1, alphaPb, alpha21, alpha22)[1]
            lst_t3.append(T3-273.15)
            T1 = t1(T1, Tw2, Tw1, alpha21, T2, X1)
            T2 = t2(T2, alpha21, Tw1, alpha22, Tw2, X1)
            T3 = t3(T3, T2, alpha22, Tw2, Tw1, X1)
            h -= X1
            h1 -= X1
            if h1 < 0:
                Tw1 = newton(Tw1, Tw2, T1, T2, T3, alpha1, alphaPb, alpha21, alpha22)[0]
                Tw2 = newton(Tw1, Tw2, T1, T2, T3, alpha1, alphaPb, alpha21, alpha22)[1]
                lst_t3.append(T3 - 273.15)
                break
        if abs(T2-T1) > 0.5:
            T20 -= 0.5*(T2-T1)
            T2 = T20
            continue
        break
    return lst_t3

# ruff: noqa: F403, F405
from saor_data import *
from math import pi, log10
A12 = 1/(1/e1 + d2*(1/e2 - 1)/d3)

def alpha3():
    vPb = 4*Gpb/(RoPb*pi*(d5**2-d4**2)) #скорость свинца
    PePb = vPb*(d5-d4)/aPb #число Пекле
    Nu = 5 + 0.025*PePb**0.8 #число Нуссельта
    alfa = Nu*lambdaPb/(d5-d4)
    return alfa
alphaPb = alpha3() #коэффициент теплопередачи свинца

def Kt(alpha1):
    return 0.5*(d2-d1)/lamw1 + 1/alpha1

def B(alpha3):
    B = 0.5*(d4-d3)/lamw2 + 1/alpha3
    return B

def t1(T1, Tw2, Tw1, alpha21, T2, X):
    Q1 = C0*A12*(Tw2**4-Tw1**4) - alpha21*(Tw1-T2)
    T11 = T1 + Q1*pi*d1*X/(Gair*Cp_down)
    return T11

def t2(T2, alpha21, Tw1, alpha22, Tw2, X):
    T22 = T2 - ((alpha21*(Tw1-T2)*d2+alpha22*(Tw2-T2)*d3)*pi*X)/(Gair*Cp_up)
    return T22
def t3(T3, T2, alpha22, Tw2, Tw1, X):
    Qs = alpha22*(Tw2-T2)+C0*A12*(Tw2**4-Tw1**4)
    T33 = T3 - Qs*pi*d4*X/(Gpb*Cp_Pb)
    return T33

def alfa21(T2):
    s = pi*(d3**2-d2**2)*0.25
    ro = p/(R*T2)
    w = Gair/(ro*s)
    v = (4.031 + 0.0551*T1 - 2.2*10**(-5)*T1**2 + 5.43*10**(-9)*T1**3)*10**(-6)/ro
    Re = w*(d3-d2)/v
    Nu = 0.021*Re**0.8*Pr_air**0.43
    lam_air = (-6.364 + 0.137 * T2 - 1.13 * 10 ** (-4) * T2 ** 2 + 6.363 * 10 ** (-8) * T2 ** 3 - 1.146 * 10 ** (
        -11) * T2 ** 4) * 10 ** (-3)
    n = 0.16*Pr_air**(-0.15)
    QTOQ = 1
    TET = 22*(0.27*(d2/d3)**2 - 1)*Re**(-0.87)*Pr_air**(-1.08)
    Nu1 = Nu*(1- 0.45/(Pr_air + 2.4))*(d3/d2)**n
    Nu2 = Nu1/(1 + TET*Nu1*QTOQ)
    alfa = Nu2*lam_air/(d3-d2)
    return alfa
def alfa22(T2):
    s = pi*(d3**2-d2**2)*0.25
    ro = p/(R*T2)
    w = Gair/(ro*s)
    v = (4.031 + 0.0551*T1 - 2.2*10**(-5)*T1**2 + 5.43*10**(-9)*T1**3)*10**(-6)/ro
    Re = w*(d3-d2)/v
    Nu = 0.021*Re**0.8*Pr_air**0.43
    lam_air = (-6.364 + 0.137*T2 - 1.13*10**(-4)*T2**2 + 6.363*10**(-8)*T2**3 - 1.146*10**(-11)*T2**4)*10**(-3)
    QTOQ = 1
    TET = 22 * (0.27 * (d2 / d3) ** 2 - 1) * Re ** (-0.87) * Pr_air ** (-1.08)
    TET2 = TET*d2/d3
    Nu1 = Nu * (1 - 0.45 / (Pr_air + 2.4)) * (d2 / d3) ** 0.6
    Nu2 = Nu1 / (1 + TET2 * Nu1/QTOQ)
    alfa = Nu2 * lam_air / (d3 - d2)
    return alfa

def alfa1(T1):
    s = pi*d1**2*0.25
    ro = p/(R*T1)
    w = Gair/(ro*s)
    v = (4.031 + 0.0551*T1 - 2.2*10**(-5)*T1**2 + 5.43*10**(-9)*T1**3)*10**(-6)/ro
    Re = w*d1/v
    Nu = 0.021*Re**0.8*Pr_air**0.43
    lam_air = (-6.364 + 0.137 * T1 - 1.13 * 10 ** (-4) * T1 ** 2 + 6.363 * 10 ** (-8) * T1 ** 3 - 1.146 * 10 ** (
        -11) * T1 ** 4) * 10 ** (-3)
    alfa = Nu*lam_air/d1
    return alfa

def newton(Tw1, Tw2, T1, T2, T3, alpha1, alpha3, alpha21, alpha22):
    while True:
        dF1dTw1 = -4 * Kt(alpha1) * C0 * A12 * Tw1 ** 3 - Kt(alpha1) * alpha21 - 1
        dF1dTw2 = 4 * Kt(alpha1) * C0 * A12 * Tw2 ** 3
        dF2dTw1 = 4 * B(alpha3) * C0 * A12 * Tw1 ** 3
        dF2dTw2 = -4 * B(alpha3) * C0 * A12 * Tw2 ** 3 - B(alpha3) * alpha22 - 1
        F1 = T1 + Kt(alpha1) * C0 * A12 * (Tw2 ** 4 - Tw1 ** 4) - Kt(alpha1) * alpha21 * (Tw1 - T2) - Tw1
        F2 = T3 - B(alpha3) * alpha22 * (Tw2 - T2) - B(alpha3) * C0 * A12 * (Tw2 ** 4 - Tw1 ** 4) - Tw2
        det0 = dF1dTw1*dF2dTw2 -  dF1dTw2*dF2dTw1
        det1 = (-1)*F1*dF2dTw2 + F2*dF1dTw2
        det2 = (-1)*F2*dF1dTw1 + F1*dF2dTw1
        dTw1 = det1/det0
        dTw2 = det2/det0
        if max(abs(dTw1),abs(dTw2)) <= 0.5:
            lst = [Tw1, Tw2]
            break
        Tw1 += dTw1
        Tw2 += dTw2
    return lst
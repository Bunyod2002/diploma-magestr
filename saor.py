# ruff: noqa: F403, F405
import saor_data as sdt
import saor_func as sfc

def saor_calc(t_pb, G):
    T2 = 380.0 + 273.15
    T20 = 380.0 + 273.15
    while True:
        lst_t3 = []
        h1 = 7.0
        h2 = 5.0
        h = 12.0
        T1 = 39.0 + 273.15  # Температура входа воздуха
        T3 = t_pb 
        TAr = 600.0 + 273.15
        Tw2 = T3 - 122
        Tw1 = Tw2 - 220

        while h2 > 0:
            alpha21 = sfc.alfa21(T2)
            alpha22 = sfc.alfa22(T2)
            alpha1 = sfc.alfa1(T1)
            Tw1 = sfc.newton(Tw1, Tw2, T1, T2, TAr, alpha1, sdt.alphaAr, alpha21, alpha22)[0]
            Tw2 = sfc.newton(Tw1, Tw2, T1, T2, TAr, alpha1, sdt.alphaAr, alpha21, alpha22)[1]
            T1 = sfc.t1(T1, Tw2, Tw1, alpha21, T2, sdt.X2)
            T2 = sfc.t2(T2, alpha21, Tw1, alpha22, Tw2, sdt.X2)
            h -= sdt.X2
            h2 -= sdt.X2

        Tw2 = T3 - 122
        Tw1 = Tw2 - 220

        while True:
            alpha21 = sfc.alfa21(T2)
            alpha22 = sfc.alfa22(T2)
            alpha1 = sfc.alfa1(T1)
            Tw1 = sfc.newton(Tw1, Tw2, T1, T2, T3, alpha1, sfc.alphaPb(G), alpha21, alpha22)[0]
            Tw2 = sfc.newton(Tw1, Tw2, T1, T2, T3, alpha1, sfc.alphaPb(G), alpha21, alpha22)[1]
            lst_t3.append(T3)
            T1 = sfc.t1(T1, Tw2, Tw1, alpha21, T2, sdt.X1)
            T2 = sfc.t2(T2, alpha21, Tw1, alpha22, Tw2, sdt.X1)
            T3 = sfc.t3(T3, T2, alpha22, Tw2, Tw1, sdt.X1, G)
            h -= sdt.X1
            h1 -= sdt.X1
            if h1 < 0:
                Tw1 = sfc.newton(Tw1, Tw2, T1, T2, T3, alpha1, sfc.alphaPb(G), alpha21, alpha22)[0]
                Tw2 = sfc.newton(Tw1, Tw2, T1, T2, T3, alpha1, sfc.alphaPb(G), alpha21, alpha22)[1]
                lst_t3.append(T3)
                break
        if abs(T2-T1) > 0.5:
            T20 -= 0.5*(T2-T1)
            T2 = T20
            continue
        break
    return lst_t3


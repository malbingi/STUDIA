from math import pow, exp

Q = 238_000.0
R = 8.3144
D = 30.0
_KIRCHOFF = 45000.0
_ITERACJE = 1_000
_BURGERS = 0.25*pow(10, -9)
_B2 = _BURGERS*_BURGERS
_TAU = (_KIRCHOFF * _B2 * 10.0**6)/2

def ivm_calculations(vector_a, Temp, _EPSILON_DOT, time):
    Temp += 273.15
    _KROK_CZASOWY = time/_ITERACJE

    RO_0 = 10_000.0
    historia_ro = [RO_0]
    _Z = _EPSILON_DOT * exp(Q/(R*Temp))
    _SWOBODNA_DROGA_DYSLOKACJI = vector_a[0]/pow(_Z, vector_a[12])
    _B_L = _BURGERS * _SWOBODNA_DROGA_DYSLOKACJI
    A1 = 1/(_B_L)
    A2 = vector_a[1] * pow(_EPSILON_DOT, (-vector_a[8])) * exp(-vector_a[2]/(R*Temp))
    A3 = vector_a[3] * (_TAU/D) * exp(-vector_a[4]/(R*Temp))
    RO_CR = -vector_a[10] + vector_a[11] * pow(_Z, vector_a[9])

    T_CR = -1
    index_of_a3 = 0

    for i in range(_ITERACJE):
        if historia_ro[-1] > RO_CR and historia_ro[-1] >= 0.:
            temp = 1.0
            if T_CR == -1:
                T_CR = i
            index_of_a3 = i - T_CR
            try:
                last_step = A3*temp*pow(historia_ro[-1], vector_a[7])*historia_ro[index_of_a3]
            except:
                print(A3, temp, historia_ro[-1], vector_a[7], historia_ro[index_of_a3] ,index_of_a3)
        else:
            temp = 0.0
            T_CR = -1
            index_of_a3 = 0
            last_step = 0.0
        delta_ro = A1*_EPSILON_DOT - A2 * _EPSILON_DOT*historia_ro[-1] - last_step
        RO_0 += delta_ro*_KROK_CZASOWY
        
        historia_ro.append(RO_0)
    boolean = penalty(vector_a)
    return RO_0, boolean

def penalty(vector):
    ret = False
    for i in range(13):
        if limit_vector(vector, i) == True:
            ret = True
    return ret

def limit_vector(vector, i):
    ret = False
    if i == 0:
        if vector[i] < 0.05*pow(10,-3):
            ret = True
        if vector[i] > 0.15*pow(10, -3):
            ret = True

    if i == 1:
        if vector[i] < 15_000.:
            ret = True
        if vector[i] > 22_000.:
            ret = True

    if i == 2:
        if vector[i] < 50.0*pow(10, 3):
            ret = True
        if vector[i] > 100.0*pow(10, 3):
            ret = True

    if i == 3:
        if vector[i] < 0.01 *3 * pow(10, 10):
            ret = True
        if vector[i] > 0.09 *3 * pow(10, 10):
            ret = True

    if i == 4:
        if vector[i] < 100.0*pow(10, 3):
            ret = True
        if vector[i] > 150.0*pow(10, 3):
            ret = True

    if i == 5:
        if vector[i] < 1.5:
            ret = True
        if vector[i] > 2.5:
            ret = True

    if i == 6:
        if vector[i] != 0.0:
            ret = True

    if i == 7:
        if vector[i] < 0.2:
            ret = True
        if vector[i] > 0.8:
            ret = True

    if i == 8:
        if vector[i] < 0.05:
            ret = True
        if vector[i] > 0.25:
            ret = True

    if i == 9:
        if vector[i] < 0.1:
            ret = True
        if vector[i] > 0.9:
            ret = True

    if i == 10:
        if vector[i] != 0:
            ret = True

    if i == 11:
        if vector[i] < 0.00001*pow(10,13):
            ret = True
        if vector[i] > 0.00009*pow(10, 13):
            ret = True
            
    if i == 12:
        if vector[i] < 0.01:
            ret = True
        if vector[i] > 0.09:
            ret = True
    return ret

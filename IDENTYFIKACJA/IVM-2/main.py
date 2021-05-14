from math import log, pow, exp, sqrt
from saveToFile import saveToFile as save

RO_0 = 10_000.0

historia_ro = [RO_0]
napr_uplast = 0.0

vector_a = [2.1 * pow(10, -3),
            176.0,
            19.5 * pow(10, 3),
            0.000148 * 3 * pow(10, 10),
            151.00 * pow(10, 3),
            0.973,
            5.77,
            1.0,
            0.0,
            0.262,
            0.0 * pow(10, 13),
            0.000605 * pow(10, 13),
            0.167]
Q = 238_000.0
R = 8.3144
Temp = 675.0 + 273.15
D = 30.0
_EPSILON_DOT = 1.0
_TIME = 1.0
_KIRCHOFF = 45000.0
_ITERACJE = 100
_KROK_CZASOWY = _TIME / _ITERACJE
_Z = _EPSILON_DOT * exp(Q/(R*Temp))
_BURGERS = 0.25*pow(10, -9)
_SWOBODNA_DROGA_DYSLOKACJI = vector_a[0]/pow(_Z, vector_a[12])
_B2 = _BURGERS*_BURGERS
_TAU = (_KIRCHOFF * _B2 * 10.0**6)/2
_B_L = _BURGERS * _SWOBODNA_DROGA_DYSLOKACJI

A1 = 1/(_B_L)
A2 = vector_a[1] * pow(_EPSILON_DOT, (-vector_a[8])) * exp(-vector_a[2]/(R*Temp))
A3 = vector_a[3] * (_TAU/D) * exp(-vector_a[4]/(R*Temp))
RO_CR = -vector_a[10] + vector_a[11] * pow(_Z, vector_a[9])
T_CR = -1
index_of_a3 = 0

save(RO_0, vector_a[6] + vector_a[5]*_KIRCHOFF*_BURGERS*sqrt(RO_0))

for i in range(_ITERACJE):
    if historia_ro[-1] >= RO_CR:
        temp = 1.0
        if T_CR == -1:
            T_CR = i
        index_of_a3 = i - T_CR
    else: 
        temp = 0.0
        T_CR = -1
        index_of_a3 = 0
    
    delta_ro = A1*_EPSILON_DOT - A2*_EPSILON_DOT*historia_ro[-1] - A3*temp*pow(historia_ro[-1], vector_a[7])*historia_ro[index_of_a3]
    
    RO_0 += delta_ro*_KROK_CZASOWY
    
    historia_ro.append(RO_0)

    napr_uplast = vector_a[6] + vector_a[5]*_KIRCHOFF*_BURGERS*sqrt(RO_0)
    
    save(RO_0, napr_uplast)
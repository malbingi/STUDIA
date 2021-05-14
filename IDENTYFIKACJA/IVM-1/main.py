from math import log, pow
from saveToFile import saveToFile as save
# WARTOŚCI A1, A2, A3, RO_CT I A8 PODMIENIAM DLA KAŻDEGO PRZYPADKU
A1 = 3.5
A2 = 5.0
A3 = 1
epsilon_dot = 1.0
ro_0 = 0.0
ro_cr = 0.4
a8 = 0.0
number_of_steps = 1000
historia_ro = [0.0]

t_cr = (1/A2) * log((ro_0-(A1/A2))/(ro_cr-(A1/A2)))
t_cr_step = int(t_cr*number_of_steps)

save(ro_0)

for i in range(number_of_steps):
    
    if i < t_cr_step:
        temp = 0.0
        index = 0
    else:
        temp = 1.0
        index = i-t_cr_step
    delta = A1*epsilon_dot - A2*historia_ro[-1]*epsilon_dot - A3*pow(historia_ro[-1], a8)*temp*historia_ro[index]
    ro_0 += delta*(1/number_of_steps)
    historia_ro.append(ro_0)
    save(ro_0)

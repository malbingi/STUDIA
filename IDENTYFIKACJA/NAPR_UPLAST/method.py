from func import square_error as c # funkcja optymizowana
from saveToFile import saveToFile as save

def hooke_jeeves(x0, epsilon, step, alpha):
    it = 0
    xB = x0
    a = b = float('inf')
    while c(x0) > epsilon or it < 1_000: #
        # warunek stopu
        if it == 0:
            a = c(x0)
        else:
            a = b
            b = c(x0)
        if a == b:
            break
        it += 1
        
        xB = x0
        x0 = try_procedure(xB, step)
        if c(x0) < c(xB):
            while c(x0) >= c(xB):
                xB_temp = xB
                xB = x0
                x0 = 2*xB-xB_temp
                x0 = try_procedure(x0, step)
        else:
            step = alpha*step
        save(c(x0),0)
    save(xB,0)
    return xB

def try_procedure(x,step):
    size = len(x)
    for i in range(size):
        temp = x
        
        # ograniczenia współczynników
        if temp[i] + step*temp[i] < dict_constraints[i][1] and temp[i] + step*temp[i] > dict_constraints[i][0]:
            temp[i] = temp[i] + step*temp[i]
        else:
            continue
        
        # po korekcji współczynników
        if c(temp) < c(x):
            x = temp
        else:
            temp[i] = temp[i] - 2.0*step*temp[i]
            if c(temp) < c(x):
                x = temp
    return x

dict_constraints = {
    0: [1.0, 1000.0],
    1: [0.0, 1.0],
    2: [0.0, 1.0],
    3: [1.0, 10_000.0],
    4: [0.0, 1.0],
    5: [1.0, 90_000.0],
    6: [0.0, 1.0]
}
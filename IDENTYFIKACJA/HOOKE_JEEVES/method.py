from func import call as c
from saveToFile import saveToFile as save

def hooke_jeeves(x0, epsilon, s, alpha, function):
    it = 0
    xB = x0
    a = b = float('inf')
    while s > epsilon or it < 1_000:
        if it == 0:
            a = c(x0, function)
        else:
            a = b
            b = c(x0, function)
        save(function, c(x0, function))
        if a == b:
            break
        it += 1
        xB = x0
        x0 = try_procedure(xB, s, function)
        if c(x0, function) < c(xB, function):
            while c(x0, function) >= c(xB, function):
                xB_temp = xB
                xB = x0
                x0 = 2*xB-xB_temp
                x0 = try_procedure(x0, s)
        else:
            s = alpha*s
    return xB

def try_procedure(x,s, function):
    size = len(x)
    for i in range(size):
        temp = x
        temp_val = list(temp)
        temp_val[i] = temp_val[i]+s
        temp = tuple(temp_val)
        if c(temp, function) < c(x, function):
            x = temp
        else:
            temp_val = list(temp)
            temp_val[i] = temp_val[i] - 2.0*s
            temp = tuple(temp_val)
            if c(temp, function) < c(x, function):
                x = temp
    return x
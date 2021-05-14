import method
from func import sigma_p as sigma
from saveToFile import saveToFile as save

dict_T_epsilon_dot = {
    0: [850, 0.1],
    1: [850, 1],
    2: [850, 10],
    3: [1000, 0.1],
    4: [1000, 1],
    5: [1000, 10],
    6: [1150, 0.1],
    7: [1150, 1],
    8: [1150, 10]
}

def main():
    x0 = [499.5, 0.5, 0.5, 4499.5, 0.5, 44999.5, 0.5]
    epsilon = 0.01
    step = 0.4
    alpha = 0.5
    val = method.hooke_jeeves(x0, epsilon, step, alpha)
    print(val)
    return val
    

if __name__ == '__main__':
    vector_a = main()
    for i in range(9):
        T, epsilon_dot = dict_T_epsilon_dot[i]
        for j in range(21):
           save(sigma(vector_a, 0.05*(j+1), epsilon_dot, T), i+1)

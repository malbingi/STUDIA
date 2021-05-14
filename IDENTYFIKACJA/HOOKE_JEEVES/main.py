import method

def main():
    x0 = (20.5,30.2)
    epsilon = 0.01
    step = 0.4
    alpha = 0.5
    val = method.hooke_jeeves(x0, epsilon, step, alpha, 0)
    print(val)
    val1 = method.hooke_jeeves(x0, epsilon, step, alpha, 1)
    print(val1)
    x0 = (10.35, 10.35)
    val2 = method.hooke_jeeves(x0, epsilon, step, alpha, 2)
    print(val2)

if __name__ == '__main__':
    main()

from math import sqrt as sq

_input = [(0,1.5), (1,2.8), (-1,-0.5), (3,6)]
_epsilon = 0.001
_K = (sq(5)-1)/2

def least_squares(a):
    val = 0.0
    for i in range(4):
        val = val + (_input[i][1] - model(_input[i][0],a) )**2
    return val
    
def model(x,a):
    return 1 + a*x

def main():
    _a = -10
    _b = 10
    _p = _b - _K*(_b-_a)
    _k = _a + _K*(_b-_a)
    _it = 0
    while (_b-_a) > _epsilon or _it < 100:
        _it += 1
        if least_squares(_p) < least_squares(_k):
            _b = _k
            _k = _p
            _p = _b - _K*(_b-_a)          
        else:
            _a = _p
            _p = _k
            _k = _a + _K*(_b-_a)
    return (_a+_b)/2

if __name__ == '__main__':
    print(main())        
        
        
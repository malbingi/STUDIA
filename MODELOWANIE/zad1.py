import numpy as np
from saveToFile import dumpStateVTR
import math
from numba import jit

_LX = 128
_LY = 128
_F = 9

_CS_SQ = 1.0/3.0 # pr 

_wC = 4.0/9.0 
_wS = 1.0/9.0
_wL = 1.0/36.0

_CX = (0, 0,  0, -1, 1, 1, -1, -1,  1)
_CY = (0, 1, -1,  0, 0, 1,  1, -1, -1)

_W = (_wC, _wS, _wS, _wS, _wS, _wL, _wL, _wL, _wL)

_Rho = np.zeros((_LX, _LY), dtype=np.float)
_Ux = np.zeros((_LX, _LY), dtype=np.float)
_Uy = np.zeros((_LX, _LY), dtype=np.float)
_map = np.zeros((_LX, _LY), dtype=np.float)

_CELLS = np.zeros((_LX, _LY, _F), dtype=np.float)
_temp_CELLS = np.zeros((_LX, _LY, _F), dtype=np.float)
_PHI_1 = np.zeros((_LX,_LY), dtype=np.float)

_ULBI = 0.0
_VLBI= 0.0

_Re = 100.0
_Ma = 0.1
_rhoiLB = 1.0

_rho = 0.0

_ULB = 0.1
_VLB = 0.0
_PHI_LBI = 1.0

_rho_LB_IN = 1.0
_ux_LB_IN = _Ma*math.sqrt(_CS_SQ)
_uy_LB_IN = 0.0

_NU_LB = _ux_LB_IN * _LY / _Re
_TAU = _NU_LB/_CS_SQ + 0.5

@jit(nopython=True)
def fEq(i, rho, ux, vy):
    cu = _CX[i] * ux + _CY[i] * vy
    u2 = ux*ux + vy*vy
    return rho * _W[i] * (1.0 + cu / _CS_SQ + 0.5*cu*cu/(_CS_SQ*_CS_SQ) - 0.5*u2/_CS_SQ)

@jit()    
def collideAndStream(c, tc):
    _FSTAR = np.zeros(9, dtype=np.float)
    for y in range(_LY):
        for x in range(_LX):
            # 0. BC
            
            #if x == _LX-1:
            #    for i in range(_F):
            #        c[x,y,i] = fEq(i, _Rho[x-1,y], _ux_LB_IN,_uy_LB_IN)
            #if x == 0:
            #    for i in range(_F):
            #        c[x,y,i] = c[x+1,y,i]
            _rho = _ux = _uy = 0.0
            if not _map[x,y] == 1:
            
                if y == _LY-1:
                    for i in range(_F):
                        c[x,y,i] = fEq(i, _rho_LB_IN, _ux_LB_IN, _uy_LB_IN)
            #if x == _LX-1:
            #    for i in range(_F):
            #        c[x,y,i] = c[x-1,y,i]
            
            # wartości makroskopowe
            
                for i in range(_F):
                    _rho += c[x,y,i]
                    _ux += c[x,y,i]*_CX[i]
                    _uy += c[x,y,i]*_CY[i]
                #print("rho ux uy", _rho, _ux, _uy)
                _ux /= _rho
                _uy /= _rho
                        
                _Rho[x,y] = _rho
                _Ux[x,y] = _ux
                _Uy[x,y] = _uy
            
            #_PHI_1[x,y] = _PHI
            # kolizja 
            
                for k in range(_F):
                    #(_TAU)
                    _FSTAR[k] = (1.0 - 1.0/_TAU)*c[x,y,k] + 1.0/_TAU*fEq(k,_rho, _ux, _uy)
                    #print(_FSTAR)
            else:
                _FSTAR[0] = c[x,y,0]
                _FSTAR[4] = c[x,y,3]
                _FSTAR[3] = c[x,y,4]
                _FSTAR[2] = c[x,y,1]
                _FSTAR[1] = c[x,y,2]
                _FSTAR[5] = c[x,y,7]
                _FSTAR[6] = c[x,y,8]
                _FSTAR[7] = c[x,y,5]
                _FSTAR[8] = c[x,y,6]
            
            #for i in range(_F):
            #    _FSTAR[i] = (1.0 - 1.0/_TAU)*c[x,y,i] + 1.0/_TAU*fEq(i, _rho, _ULB, _VLB)
            # streaming 
            _XP = 0 if x == _LX-1 else x+1 
            _YP = 0 if y == _LY-1 else y+1 
            _XM = _LX-1 if x == 0 else x-1 
            _YM = _LY-1 if y == 0 else y-1 
            
            #print(_FSTAR)
            tc[x,y,0] = _FSTAR[0]
            tc[x,_YP,1] = _FSTAR[1]
            tc[x,_YM,2] = _FSTAR[2]
            tc[_XM,y,3] = _FSTAR[3]
            tc[_XP,y,4] = _FSTAR[4]
            tc[_XP,_YP,5] = _FSTAR[5]
            tc[_XM,_YP,6] = _FSTAR[6]
            tc[_XM,_YM,7] = _FSTAR[7]
            tc[_XP,_YM,8] = _FSTAR[8]
            #_CX = (0, 0,  0, -1, 1, 1, -1, -1,  1)
            #_CY = (0, 1, -1,  0, 0, 1,  1, -1, -1)
    return tc

def setInitialConditions(rhoi, uxi, vyi):
    _DX = 1.0 / _LY
    for y in range(_LY):
        for x in range(_LX):
            #_Rho[x,y] = rhoi
            #_Ux[x,y] = uxi
            #_Uy[x,y] = vyi
            
           
            #if x == _LX/2:
            #    _PHI = phiH
            #else:
            #    _PHI = phiL
            
            #_PHI_1[x,y] = _PHI
            #_XP = _DX * ( x + 0.5)/2
            #phi = 0.5 * (1 + math.sin(2 * math.pi * _XP))
            for i in range(_F):
                _CELLS[x,y,i] = fEq(i, rhoi, uxi, vyi)
                #print(x,y,i,_CELLS[x,y,i])
            if x == 0 or x == _LX-1 or y == 0 :
                _map[x,y] = 1
                #_Rho[x,y] = 0
                #_Ux[x,y] = _Uy[x,y] = 0.0
    #for i in range(_F):
    #    _CELLS[_LX//2, _LY//2, 1] = fEq(i, 1.0, ui, vi)

def main():
    iter = 0
    ITERMAX = 100_001
    setInitialConditions(_PHI_LBI, _ULBI, _VLBI)
    dumpStateVTR(_LX, _LY, _F, _PHI_1, f"./0405/state_prev.vtk", _Ux, _Uy, _map, _Rho)
    
    fi0 = fi1 = 0.0
    
    while iter < ITERMAX:
        for i in range(_LX):
            for j in range(_LY):
                fi1 = _Ux[i,j]**2 + _Uy[i,j]**2
        #print(math.fabs(fi1))
        
        #if iter >= 100_001 and math.fabs(fi1-fi0) < 1e-3:
        #    break
        
        fi0 = fi1
        fi1 = 0.0
        """
        if iter%1000 == 0:
            print(iter)
        """
        if iter%2 == 0:
            collideAndStream(_CELLS, _temp_CELLS)
        else:
            collideAndStream(_temp_CELLS, _CELLS)
        if iter%100 == 0:
            dumpStateVTR(_LX, _LY, _F, _PHI_1, f"./0405/state{iter}.vtk", _Ux, _Uy, _map, _Rho)
        iter += 1
    
    return

if __name__ == "__main__":
    main()
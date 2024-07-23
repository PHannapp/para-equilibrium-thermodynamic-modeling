import numpy as np

def F11(p):
    return -1.8692e-09*p
def F12(p):
    return -2.3753e-10*p
def F13(p):
    return -2.50627e-11*p
def F14(p):
    return -3.4483e-08*p
def F15(p):
    return -1.2469e-08*p
def EXPO1(p):
    return np.exp(F11(p))
def EXPO2(p):
    return np.exp(F12(p))
def EXPO3(p):
    return np.exp(F13(p))
def EXPO4(p):
    return np.exp(F14(p))
def EXPO5(p):
    return np.exp(F15(p))
def FF1(p):
    return 5.35e+08*EXPO1(p)
def FF2(p):
    return 4.21e+09*EXPO2(p)
def FF3(p):
    return 3.99e+10*EXPO3(p)
def FF4(p):
    return 29000000*EXPO4(p)
def FF5(p):
    return 80200000*EXPO5(p)

# reference -------------------------------------------------------------------------------
def GHSERH2(T):
    return -9522.97393+78.5273873*T-31.35707*T*np.log(T)+0.0027589925*T**2-0.000000746390667*T**3+56582.3*T**(-1)
def GHSERNI(T):
    return -5179.159+117.854*T-22.096*T*np.log(T)-0.0048407*T**2
def GHSERCE(T):
    return -7160.519+84.23022*T-22.3664*T*np.log(T)-0.0067103*T**2-0.000000320773*T**3-18117*T**(-1)
def GHSERLA(T):
    return -7968.403 + 120.284604*T - 26.34*T*np.log(T) - 0.001295165*T**2

# Endmembers -------------------------------------------------------------------------------
def GCE_NI_V_V(T):
    return -197460 + (24.3084) * T + GHSERCE(T) + 5*GHSERNI(T)
def GCE_NI_H_V(T):
    return -128606.71962211627 + (24.3084 + 61) * T + GHSERCE(T) + 5*GHSERNI(T)+0.5*GHSERH2(T)
def GCE_NI_V_H(T):
    return -230655.58499331266 + (24.3084 + 366) * T + GHSERCE(T) + 5*GHSERNI(T)+3*GHSERH2(T)
def GCE_NI_H_H(T):
    return -258147.0997064103 + (24.3084 + 427) * T + GHSERCE(T) + 5*GHSERNI(T) + 3.5*GHSERH2(T)

def GLA_NI_V_V(T):
    return -154674 + 30.6215 * T + GHSERLA(T) + 5*GHSERNI(T)
def GLA_NI_H_V(T):
    return -153101.98965182144 + (30.6215 + 61) * T + GHSERLA(T) + 5*GHSERNI(T) + 0.5*GHSERH2(T)
def GLA_NI_V_H(T):
    return  -216752.56601723537 + (30.6215 + 366) * T + GHSERLA(T) + 5*GHSERNI(T) + 3*GHSERH2(T)
def GLA_NI_H_H(T):
    return -268858.9951716817 + (30.6215 + 427) * T + GHSERLA(T) + 5*GHSERNI(T) + 3.5*GHSERH2(T)

# Interaction parameters -------------------------------------------------------------------
def LCE_NI_HV_H_0(T):
    return -22772.324506943554
def LCE_NI_H_HV_0(T):
    return -141445.99185666817
def LCE_NI_V_HV_1(T):
    return 3565.7555723479286
def LLA_NI_H_HV_1(T):
    return -22942.766915992015

def LCELA_NI_H_H_0(T):
    return -38562
def LCELA_NI_H_V_0(T):
    return 74202
def LCELA_NI_H_V_1(T):
    return -61141

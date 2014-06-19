import numpy as np

def EOM_inertial_2bodyconstant_thrust(t, X, Thrust, Mdot, mu):
    r = np.linalg.norm(X[0:3])
    r3 = r**3
    dX = np.zeros(7)
    g0 = 9.80665

    dX[0] = X[3]
    dX[1] = X[4]
    dX[2] = X[5]
    dX[3] = -mu*X[0]/r3 + Thrust[0]/X[6]
    dX[4] = -mu*X[1]/r3 + Thrust[1]/X[6]
    dX[5] = -mu*X[2]/r3 + Thrust[2]/X[6]
    dX[6] = -Mdot

    return dX

def EOM_inertial_2body(t, X, mu):
    r = np.linalg.norm(X[0:3])
    r3 = r**3
    dX = np.zeros(6)

    dX[0] = X[3]
    dX[1] = X[4]
    dX[2] = X[5]
    dX[3] = -mu*X[0]/r3
    dX[4] = -mu*X[1]/r3
    dX[5] = -mu*X[2]/r3

    return dX
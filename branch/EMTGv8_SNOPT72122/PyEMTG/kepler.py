#Kepler solver for PyEMTG
#Ryne Beeson 7-9-2014

#  Orbital Elements -> r, v Arrays
def coe2rv(oe, mu):
    #  import statements
    from math  import cos, sin, sqrt
    from numpy import matrix, array, zeros
    
    #  if 'a' is set to zero, then body is at rest in the
    #+ current frame of reference.
    #+ return a zero vector for r and v
    if oe[0] == 0.0: return zeros(3,1), zeros(3,1)
    
    a  = oe[0]
    e  = oe[1]
    i  = oe[2]
    Om = oe[3]
    om = oe[4]
    f  = oe[5]
    
    p  = a*(1 - e*e)
    r  = p/(1 + e*cos(f))
    rv = matrix([r*cos(f), r*sin(f),   0])
    vv = matrix([-sin(f),  e + cos(f), 0])
    vv = sqrt(mu/p)*vv
    
    c0 = cos(Om); s0 = sin(Om)
    co = cos(om); so = sin(om)
    ci = cos(i);  si = sin(i)
    
    R  = matrix([[c0*co - s0*so*ci, -c0*so - s0*co*ci,  s0*si],
                 [s0*co + c0*so*ci, -s0*so + c0*co*ci, -c0*si],
                 [so*si,             co*si,             ci]])
                 
    ri = array(R*rv.T); ri = ri.reshape(3)
    vi = array(R*vv.T); vi = vi.reshape(3)
    
    return ri, vi


#  laguerre_conway solver
def laguerre_conway(e, M):
    #  import statements
    from math  import sqrt, sin, cos
    from numpy import abs
    #  generate an initial guess for eccentric anomaly (E)
    En = (M*(1 - sin(M + e)) + (M + e)*sin(M)) / (1 + sin(M) - sin(M + e))
    n  = 4
    #  for-loop
    for i in range(1,20):
        f      =  M - En + e*sin(En)
        fdash  = -1 + e*cos(En)
        fddash = -e*sin(En)
        g = sqrt(((n - 1)**2) * (fdash**2) - n*(n - 1)*f*fddash)
        if fdash > 0:
            En1 = En - (n*f/(fdash + g))
        else:
            En1 = En - (n*f/(fdash - g))
        #  calculate error
        error = abs(En1 - En)
        En = En1
        if error <= 1E-4:
            break
    
    #  return the eccentric anomaly
    return En


#  solve kepler's equation using laguerre_conway
def kepler(sma, ecc, inc, RAAN, AOP, MA, reference_time, epoch, mu):
    #  input delta_time in (secs)
    #  import statements
    from math import pi, sqrt, tan, atan
    #  delta_time
    delta_time = (epoch - reference_time)*86400.0
    n = 1.0 / sqrt((sma**3)/mu)
    MA = (MA + n*delta_time) % (2.0*pi)
    #  update E
    E   = laguerre_conway(ecc, MA) % (2.0*pi)
    #  calculate f
    f = 2.0*atan(tan(E/2.0) * sqrt((1 + ecc) / (1 - ecc))) % (2.0*pi)
    #  coe2rv
    r, v = coe2rv([sma, ecc, inc, RAAN, AOP, f], mu)
    #  return r and v

    return r, v


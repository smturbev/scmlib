''' util.py
    created by Sami Turbeville
    on 9/21/23
'''

# universal constants
R    = 287.058
G    = 9.80665

def omega2w(omega, p, t):
    ''' converts vertical velocity from pressure coords
        to m/s.
        
        input:
         - omega : array in units of Pa/s
         - p     : pressure array in Pa
         - t     : temperature array in K
        output:
         - w     : vertical velocity array in units of m/s
    '''
    rho = p / (R*t)
    w = -omega/(rho*G)
    return w

def w2omega(w, p, t):
    ''' converts vertical velocity from m/s
        to pressure coords (Pa/s).
        
        input:
         - w     : array in units of m/s
         - p     : pressure array in Pa
         - t     : temperature array in K
        output:
         - omega : vertical velocity array in units of Pa/s
    '''
    rho = p / (R*t)
    omega = -w*rho*G
    return omega
    
    
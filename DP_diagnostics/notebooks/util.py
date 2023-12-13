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
    
def calc_rice(qi, ni):
    rho = 920 # kg/m3
    qi = qi # kg/kg
    ni =  ni # 1/kg
    r_ice = np.where((ni>1e-5),(3*qi/(4*np.pi*rho*ni))**(1/3),0)*1e6 # um
    return r_ice

def calc_ni(numice, qv, p, t):
    ni = numice * calc_rho(qv,p,t)  # 1/kg * kg/m3 = 1/m3
    return ni

def calc_rho(qv, p, t):
    """calculates density of air in kg/m3 for given input"""
    R =  287 # (Gas constant of air) J/(kg*K)
    Tv = (1 + 0.61*qv) * t # K
    rho = p / (R*Tv) # kg/m3
    return rho

def calc_rhice(ds, varQ="Q", varT="T", z_units="hPa"):
    """ input: xarray with variables Q and T
        output: xarray of rh wrt ice
    """
    e_si = np.exp(9.550426 - 5723.265/ds[varT] \
                  + 3.53068*np.log(ds[varT]) - 0.00728332*ds[varT])
    if z_units=="hPa":
        z = ds.lev*100
    else:
        z=ds.lev
    w_si = (0.622 * e_si) / (z - e_si)
    w_i  = ds[varQ] / (1 - ds[varQ])
    rh_ice = w_i/w_si * 100
    return rh_ice
    
############################################
####      Unit Conversion Functions     ####
############################################
# These are obviously unneccessary if we just use metric, but I'll leave them here for convenience later.

### LENGTH/DISTANCE
def mm2in(dim_mm):              # Convert from mm to in for some reason.
    """Return inches from input millimeter value"""
    dim_in = dim_mm / 25.4
    return dim_in

def in2mm(dim_in):              # In case we need to convert from in to mm for some reason.
    """Return millimeters from input inch value"""
    dim_mm = dim_in * 25.4
    return dim_mm

### STRESS/PRESSURE/MODULUS
ksi2mpa = 6.89475729        # constant to convert imperial to metric units

def psi2MPa(E_psi):        # convert modulus in psi to MPa
    """Return megapascals (MPa) from input psi value"""
    E_MPa = E_psi * 1e-3 * ksi2mpa     # convert from psi to ksi, then ksi to MPa
    return E_MPa

def Mbar2MPa(Mbar):
    """Convert pressure from Mbar to MPa"""
    MPa = Mbar*100000
    return MPa

def MPa2Mbar(MPa):
    """Convert pressure from MPa to Mbar"""
    Mbar = MPa/100000
    return Mbar

def kbar2MPa(kbar):
    """Convert pressure from Mbar to MPa"""
    MPa = kbar*100
    return MPa

def MPa2kbar(MPa):
    """Convert pressure from MPa to Mbar"""
    kbar = MPa/100
    return kbar

### VELOCITY
def in_s2mm_us(c_in_s):         # convert wave speeds in in/s to mm/microsecond (km/s).
    """Return millimeter/microsecond from input inches/second"""
    c_mm_us = c_in_s * 25.4 * 1e-6      # convert from in to mm and /s to /us
    return c_mm_us

def m_s2mm_us(c_m_s):
    """Return millimeter/microsecond from input meters/second"""
    c_mm_us = c_m_s * 1e-3         # (* 1e3 * 1e-6) to convert from m to mm and /s to /us
    return c_mm_us

def mm_us2m_s(c_mm_us):
    """Return millimeter/microsecond from input meters/second"""
    c_m_s = c_mm_us * 1000         # (*1e-3 * 1e6) to convert from mm to m and /us to /s
    return c_m_s

### TEMPERATURE
def temp_C2K(temp_C):
    """Convert temperature from celcius (C) to absolute temperature in Kelvin (K)"""
    temp_K = temp_C + 273.15
    return temp_K

def temp_F2C(temp_F):
    """Convert temperature from Fahrenheit (F) to Celcius (C)"""
    temp_C = (temp_F -32)*(5/9)
    return temp_C

def temp_C2F(temp_C):
    """Convert temperature from Fahrenheit (F) to Celcius (C)"""
    temp_F = (temp_C*(9/5))+32
    return temp_F

### TIME
def time_hr2s(time_hr):
    """Convert time from hours to seconds"""
    time_s = time_hr * 3600
    return time_s

def time_min2s(time_min):
    """Convert time from hours to seconds"""
    time_s = time_min * 60
    return time_s

if __name__ == '__main__':
    print('unit conversion tool currently has no standalone functionality')
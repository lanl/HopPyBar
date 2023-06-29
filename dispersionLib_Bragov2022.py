# dispersionLib.py
# Dispersion Correction via Python
# @author: Sasha
 # (received via personal correspondence from A.M. Bragov)
# Created on Mon Feb 07 15:42:59 2011
# See Bragov et al., "Dispersion correction in split-Hopkinson pressure bar: theoretical and experimental analysis." Continuum Mechanics and Thermodynamics 34 (2022) 895-907, (https://doi.org/10.1007/s00161-019-00776-0) for more info on dispersion correction. 
# Many comments added by B.M. Morrow (for clarity), 2020.

import mpmath as mm
import numpy as np
import json

class barWithDispersion:
    """Dispersion Correction Class for Kolsky Bar.
    Methodology: Perform dispersion shift of a pulse to obtain an explicit dependence of the wave number ξ, AKA xi/ksi, on frequency omega. 
    This method is based on the direct solution of the Pochhammer-Chree frequency equation. Functional dependence ksi(omega) is determined by solving the frequency (equation 1 in paper).
    The algorithm for solving this equation is as follows: 
        1. The interval [0, omegamax] is selected to solve the given equation
        2. The step delta omega is set
        3. Value ksi_n, satisfying Eq. (1), is numerically determined for frequency omega_n = omega_n-1 + delta omega
        4. The linear extrapolation of the last two solutions is chosen as the initial approximation for ksi_n:
          (ksi_n)^0 = ksi_(n-1) + (ksi_(n-1) - ksi_(n-2))/(omega_(n-1) - omega_(n-2))*(omega_n - omega_(n-1))
        5. If the solution with given accuracy could not be found, then step (delta omega) is to be divided into half and steps 3 and 4 are to be repeated.

    Usage: barWithDispersion(E = 194310000000, rho = 8082., nu = 0.31025, r0 = (d_bar_i_equiv/2)/1000)

    Input Definitions:
        r - bar radius (m)
        rho - bar density (kg/m^3)
        E - Young's modulus (Pa)
        nu - Poisson's ratio
        cb = 'bar' velocity (m/s)"""

    def __init__(self, r0 = 5e-3, rho = 7800., E = 200e9, nu = 0.29):
        self.r0 = r0                      # bar radius (m)
        self.rho = rho                    # density (kg/m^3)
        self.E = E                        # Young's Modulus (Pa)
        self.nu = nu                      # Poisson ratio
        self.l = E*nu/(1+nu)/(1-2*nu)     # lambda (Lame's first parameter)
        self.mu = E/2./(1+nu)             # shear modulus (mu or G) (Lame's second parameter)
        self.cb = np.sqrt(E/rho)          # 'bar' velocity (m/s); units reduce from Pa/(kg/m^3) to m/s.
        self.cs = np.sqrt(self.mu/rho)    # shear velocity
        self.cl = np.sqrt((self.l+2*self.mu)/rho)   # velocity of longitudinal waves (P-waves) (m/s)
        self.cr = (0.862+1.14*nu)/(1+nu)*self.cs      # velocity of Raleigh waves (dispersionless)

    def f(self, ksi, omega):    # function relateing harmonic phase velocity xi and frequency omega
        # Remember: ksi is lowercase greek letter xi ("ksi"), not a unit of pressure
        alpha_sq = self.rho*omega**2/(self.l + 2*self.mu) - ksi**2  # Alpha squared in Pochhammer-Chree (P-C) equation
        beta_sq = self.rho*omega**2/self.mu - ksi**2  # Beta squared in P-C eq
        alpha = mm.sqrt(alpha_sq)
        beta = mm.sqrt(beta_sq)
        # The equation in Bragov 2022 is misprinted. The equations for alpha_sq and beta_sq differ between Bragov2022 and Bacon1998 by sign of ksi**2.
        # Here we use the intent of Bragov et al. 
        freq = (
            ((2 * alpha/self.r0) * (beta_sq + ksi**2) * mm.besselj(1, alpha * self.r0) * mm.besselj(1, beta * self.r0)) - 
            ((beta_sq - ksi**2)**2 * mm.besselj(0, alpha * self.r0) * mm.besselj(1, beta * self.r0)) -
            (4 * ksi**2 * alpha * beta * mm.besselj(1, alpha * self.r0) * mm.besselj(0, beta*self.r0))
            )
        return abs(freq)

    def solve(self, omegamax = 5e6, domega0 = lambda x: 10e3):
        # this function numerically solves the Pochamer-Cree equation in the frequency range up to omegamax with domega increments
        omega = [0]  # [0,1000.,2000.]
        ksi = [0]  # [0,omega[1]/cb, omega[2]/cb]
        domega = domega0(0)
        
        while omega[-1] <= omegamax and domega > 1e-6 * domega0(omega[-1]):
            xx = omega[-1] + domega
            def q(x): return self.f(x, xx)/xx**2
            try:
                if len(omega) >= 3:
                    y0 = ((xx - omega[-2])/(omega[-1] - omega[-2]) * (ksi[-1] - ksi[-2]) + ksi[-2])
                else:
                    y0 = 1.0
                y = mm.findroot(q, y0, solver='mnewton')
                omega.append(omega[-1] + domega)
                ksi.append(y.real)
                #print("Dispersion Correction; Converged:", omega[-1])  # we only need to know when there's an error (below)
                domega = domega0(omega[-1])
            except ValueError:
                print(f'Dispersion Correction Not Converged: {omega[-1] + domega} {domega}')
                domega *= 0.5
        self.ksi_omega = lambda x: np.interp(x, omega, ksi)
        self.omega = np.array(omega)
        self.ksi = np.array(ksi)
        self.omegamax = max(self.omega)

# These two properties are unused currently.
    # @property
    # def cp_c0(self):
    #     rez = []
    #     for i in range(len(self.omega)):
    #         if self.ksi[i] == 0:
    #             rez.append(1)
    #         else:
    #             rez.append(self.omega[i]/self.ksi[i]/self.cb)
    #     return rez

    # @property
    # def r_l(self):
    #     return self.r0*self.ksi/2/np.pi       # should this be 2*pi?

    def Ksi(self, o):
        if o <= self.omegamax:
            return self.ksi_omega(o)
        else:
            return self.ksi[-1]+(o-self.omega[-1])/self.cr
#            return self.ksi[-1]+(omega-omega[-1])/self.cr

    def save(self, fname):
        rez = {}
        rez['r0'] = self.r0                # bar radius (m, I think)
        rez['rho'] = self.rho              # density (g/cm^3)
        rez['E'] = self.E                  # Young's Modulus
        rez['nu'] = self.nu                # Poisson Ratio
        rez['omegamax'] = self.omegamax    # maximum frequency?
        rez['omega'] = list(map(float, self.omega))    # list of frequencies solved?
        rez['ksi'] = list(map(float, self.ksi))        # list of ksi values corresponding to those omega values?
        json.dump(rez, open(fname, 'w'))

    def load(self, fname='', data=None):
        if data == None:
            rez = json.load(open(fname, 'r'))
        else:
            rez = json.loads(data)
        self.r0 = rez['r0']                # bar radius? (m?)
        self.rho = rez['rho']              # density (g/cm^3)
        self.E = rez['E']                  # Young's Modulus
        self.nu = rez['nu']                # Poisson Ratio
        self.omegamax = rez['omegamax']    # Maximum frequecy
        self.omega = np.array(rez['omega'])    # list of frequencies solved?
        self.ksi = np.array(rez['ksi'])        # list of ksi values corresponding to those omega values?
        self.l = self.E*self.nu/(1+self.nu)/(1-2*self.nu)      # lambda (Lame's first parameter)
        self.mu = self.E/2./(1+self.nu)        # shear modulus (mu or G) (Lame's second parameter)
        self.cb = np.sqrt(self.E/self.rho)     # 'bar' velocity
        self.cs = np.sqrt(self.mu/self.rho)    # shear velocity (I'm assuming the "t" is for "transverse", because this is defo c_shear)
        self.cl = np.sqrt((self.l+2*self.mu)/self.rho)         # velocity of P-wave (longitudinal)
        self.cr = (0.862+1.14*self.nu)/(1+self.nu)*self.cs     # velocity of Raleigh waves (dispersionless)
        self.ksi_omega = lambda x: np.interp(x, self.omega, self.ksi)   # interpolating between discrete values of ksi and omega?

    def disp_shift_wave(self, t, y, dz, backshift = False):
        # shifts the signal y(t) by the distance dz, taking into account the dispersion effect
        fy = np.fft.rfft(y)
        n = len(fy)
        dz = -dz
        omega = 2 * np.pi * np.arange(n) / max(t)
        fy2 = []
        
        for i in range(n):
            fy2.append(np.exp(complex(0, self.Ksi(omega[i]) * dz)) * fy[i])
        
        rez = np.fft.irfft(fy2).tolist()
        if backshift:
            dt = dz/self.cb
            rez = self.simple_shift_wave(t, rez, dt)
        
        rez = rez + [rez[-1]] * (len(y) - len(rez))
        return rez

    def simple_shift_wave(self, t, y, tshift):
        dt = t[-1] - t[0]
        tshift = tshift - np.floor(tshift/dt) * dt
        dt = t[1] - t[0]
        n = int(np.round(tshift/dt))
        if isinstance(y, np.ndarray):
            y = y.tolist()
        return y[-n:]+y[:-n]

if __name__ == '__main__':
    import matplotlib.pylab as plt
    c = barWithDispersion(E=200e9, rho=7850, nu=0.29, r0=10e-3)
    c.solve()
    c.save('bar_dispersion.txt')
#    c1=barWithDispersion()
#    c1.load('bar1.txt')
#    plt.plot(c.omega, c.ksi)
#    plt.show()
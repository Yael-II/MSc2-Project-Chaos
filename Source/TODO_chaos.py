import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad

import potentials as pot
import integrator as itg
import initial_conditions as init

DEFAULT_N_iter = 30000
DEFAULT_N_part = 100
DEFAULT_h = 0.01
DEFAULT_c = 1.7

E_all = np.array([1/100, 1/12, 1/10, 1/8, 1/6])



def compute_coordinates(E: float,
                        N_iter: int = DEFAULT_N_iter,
                        N_part: int = DEFAULT_N_part,
                        h: float = DEFAULT_h) -> tuple:
    W_part = init.n_energy_part(pot.hh_potential, N_part, E)

    # Perform integration
    t_part, coord_part = itg.rk4(0, W_part, h, N_iter, pot.hh_evolution)

    # Extract positions and velocities
    x_part = coord_part[:, 0, 0]
    y_part = coord_part[:, 0, 1]
    u_part = coord_part[:, 1, 0]
    v_part = coord_part[:, 1, 1]

    return t_part, x_part, y_part, u_part, v_part

def x(t: float,
      x_part: np.ndarray,
      h: float = DEFAULT_h):
    return x_part[t//h]

def theta(t: float, 
          x,
          x_part: np.ndarray,
          c: float = DEFAULT_c, 
          h: float = DEFAULT_h,
          ds: float = DEFAULT_h):
    """theta(t) = c*t + int_0^t phi(x(s)) ds"""
    I = quad(lambda s: x(s, x_part), 0, t)
    return c*t + I

def p(t: float,
      x,
      theta,
      x_part: np.ndarray,
      h: float = DEFAULT_h,
      ds: float = DEFAULT_h):
    I = quad(lambda s: x(s, x_part)*cos(theta(s, x, x_part)), 0, t)
    return I

def M(t,
      p,
      x_part: np.ndarray,
      h: float = DEFAULT_h,
      dtau: float = DEFAULT_h):
    I = quad(lambda tau: (p(t+tau, x, theta, x_part) 
                          - p(tau, x, theta, x_part)))**2, 0, np.max(t))
    return I

if __name__ == "__main__":
    for i in range(len(E_all)):
        E = E_all[i]
        t_part, x_part, y_part, v_part, u_part = compute_coordinates(E)
        I = M(t_part, x_part)


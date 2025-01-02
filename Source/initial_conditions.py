import numpy as np
import matplotlib.pyplot as plt
POS_MIN = -1
POS_MAX = +1
VEL_MIN = -1
VEL_MAX = +1
N_PART = 1

def mesh_grid(N: int = N_PART,
              xmin: float = POS_MIN, 
              xmax: float = POS_MAX, 
              ymin: float = POS_MIN, 
              ymax: float = POS_MAX) -> np.ndarray:
    """Generates a set of regularly sampled particles with no velocity
    :param N: number of particles
    :parma xmin: Minimum value for position x
    :param xmax: Maximum value for position x
    :param ymin: Minimum value for position y
    :param ymax: Maximum value for position y
    :returns W: Phase-space vector
    """
    X = np.linspace(xmin, xmax, N)
    Y = np.linspace(ymin, ymax, N)
    X,Y = np.meshgrid(X,Y, indexing="ij")
    return np.array([X,Y])

def one_part(x0: float = 0, 
             y0: float = 0, 
             u0: float= 0, 
             v0: float = 0) -> np.ndarray:
    """Generates a particle at position (x0, y0) 
    with velocities (u0,v0)
    :param x0: Initial position x
    :param y0: Initial position y
    :param u0: Initial velocity u
    :param v0: Initial velocity v
    """
    X = x0
    Y = y0
    U = u0
    V = v0
    return np.array([[X,Y], [U,V]])

def n_energy_part(potential,
                  N: int = N_PART,
                  E: float = 0,
                  xmin: float = -1,
                  xmax: float = +1,
                  ymin: float = -0.5,
                  ymax: float = +1):
    X = []
    Y = []
    POT = []
    U = []
    V = []
    while len(X) < N:
        x = np.random.random()*(xmax-xmin)+xmin
        y = np.random.random()*(ymax-ymin)+ymin
        w = np.array([x, y])
        pot = potential(w, position_only=True)[0]
        if pot <= E:
            X.append(x)
            Y.append(y)
            POT.append(pot)
    X = np.array(X)
    Y = np.array(Y)
    POT = np.array(POT)
    U = np.zeros_like(X)
    V = np.zeros_like(Y)
    C = np.sqrt(2 * (E - POT))
    THETA = np.random.random(N)*2*np.pi
    U = C*np.cos(THETA)
    V = C*np.sin(THETA)
    return np.array([[X, Y], [U, V]])

def n_energy_2part(potential,
                   N: int = N_PART,
                   E: float = 0,
                   sep: float = 1e-7,
                   xmin: float = -1,
                   xmax: float = +1,
                   ymin: float = -0.5,
                   ymax: float = +1):
    W_1 = n_energy_part(potential, N, E)
    W_2 = np.zeros_like(W_1)
    alpha = np.random.uniform(0, 2*np.pi, N)
    W_2[0, 0] = W_1[0, 0] + sep*np.cos(alpha)
    W_2[0, 1] = W_1[0, 1] + sep*np.sin(alpha)
    W_2[1, 0] = W_1[1, 0]
    W_2[1, 1] = W_1[1, 1]
    return (W_1, W_2)






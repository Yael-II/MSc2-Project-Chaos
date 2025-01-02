import numpy as np

def kinetic(W: np.ndarray) -> np.ndarray:
    """Computes the kinetic energy.
    :param W: Phase-space vectors
    :returns T: Kinetic energy
    """
    U = W[1,0]
    V = W[1,1]
    # If U or V is not an array (or a list), but rather a scalar, then we
    # create a list of one element so that it can work either way
    if np.ndim(U) == 0: U = np.array([U])
    if np.ndim(V) == 0: V = np.array([V])
    
    return (U**2 + V**2)/2

def total(W: np.ndarray,
          potential,
          kinetic = kinetic) -> np.ndarray:
    return potential(W) + kinetic(W)

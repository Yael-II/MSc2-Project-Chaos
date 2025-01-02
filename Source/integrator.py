import numpy as np

def euler(x0, y0, h, n, func): # FIXME cannot be used with vectors
    """
    Euler method adapted for state vector [[x, y], [u, v]]
    :param x0: initial time value
    :param y0: initial state vector [[x, y], [u, v]]
    :param h: step size (time step)
    :param n: number of steps
    :param func: RHS of differential equation
    :returns: x array, solution array
    """
    x_values = np.zeros(n)
    y_values = np.zeros((n, 2, 2))  #  to accommodate the state vector

    for i in range(n):
        dydt = func(x0, y0)
        y0 = y0 + h * dydt
        x0 = x0 + h

        x_values[i] = x0
        y_values[i, :, :] = y0

    return x_values, y_values

# Updated RK2 integrator
def rk2(x0, y0, h, n, func): # FIXME cannot be used with vectors
    """
    RK2 method adapted for state vector [[x, y], [u, v]]
    :param x0: initial time value
    :param y0: initial state vector [[x, y], [u, v]]
    :param h: step size (time step)
    :param n: number of steps
    :param func: RHS of differential equation
    :returns: x array, solution array
    """
    x_values = np.zeros(n)
    y_values = np.zeros((n, 2, 2))  #  to accommodate the state vector

    for i in range(n):
        k1 = func(x0, y0)
        k2 = func(x0 + h / 2., y0 + h / 2. * k1)

        y0 = y0 + h * (k1 / 2. + k2 / 2.)
        x0 = x0 + h

        x_values[i] = x0
        y_values[i, :, :] = y0

    return x_values, y_values

def rk4(t0: float, 
        W0: np.ndarray, 
        h: float, 
        n: int, 
        func):
    """
    RK4 method adapted for state vector [[x, y], [u, v]]
    :param x0: initial time value
    :param y0: initial state vector [[x, y], [u, v]]
    :param h: step size (time step)
    :param n: number of steps
    :param func: RHS of differential equation
    :returns: x array, solution array
    """
    time = np.zeros(n)
    W = np.zeros((n,) + np.shape(W0))
    # to accommodate the state vector
    t = t0
    w = W0
    for i in range(n):
        k1 = func(t, w)
        k2 = func(t + h / 2., w + h / 2. * k1)
        k3 = func(t + h / 2., w + h / 2. * k2)
        k4 = func(t + h, w + h * k3)

        w = w + h * (k1 / 6. + k2 / 3. + k3 / 3. + k4 / 6.)
        t = t + h

        time[i] = t
        W[i] = w

    return t, W


def integrator_type(x0, y0, h, n, func, int_type):
    return int_type(x0, y0, h, n, func)

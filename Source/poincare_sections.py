import numpy as np
def pcs_find(pos_x, pos_y, vel_x, vel_y):
    if np.ndim(pos_x) == 1: 
        pos_x = np.array([pos_x])
        pos_y = np.array([pos_y])
        vel_x = np.array([vel_x])
        vel_y = np.array([vel_y])
    pcs_pos_y = []
    pcs_vel_y = []
    for j in range(len(pos_x[0])): # for each particle
        i = 0
        x = pos_x[:,j]
        y = pos_y[:,j]
        u = vel_x[:,j]
        v = vel_y[:,j]
        while i < len(x) - 1:
            if x[i] * x[i+1] < 0:
                y0 = y[i] \
                        + (y[i+1] - y[i])/(x[i+1] - x[i]) \
                        * (0 - x[i])
                v0 = v[i] \
                        + (v[i+1] - v[i])/(x[i+1] - x[i])\
                        * (0 - x[i])
                pcs_pos_y.append(y0)
                pcs_vel_y.append(v0)
            i += 1
    return pcs_pos_y, pcs_vel_y
def pcs_find_legacy(pos_x, pos_y, vel_x, vel_y):
    """Depreciated legacy function that should not be used
    """
    if np.ndim(pos_x) == 1: 
        pos_x = np.array([pos_x])
        pos_y = np.array([pos_y])
        vel_x = np.array([vel_x])
        vel_y = np.array([vel_y])
    pcs_pos_y = []
    pcs_vel_y = []
    for j in range(len(pos_x)): # for each particle
        i = 0
        x = pos_x[j]
        y = pos_y[j]
        u = vel_x[j]
        v = vel_y[j]
        while i < len(x) - 1:
            if x[i] * x[i+1] < 0:
                y0 = y[i] \
                        + (y[i+1] - y[i])/(x[i+1] - x[i]) \
                        * (0 - x[i])
                v0 = v[i] \
                        + (v[i+1] - v[i])/(x[i+1] - x[i])\
                        * (0 - x[i])
                pcs_pos_y.append(y0)
                pcs_vel_y.append(v0)
            i += 1
    return pcs_pos_y, pcs_vel_y

import matplotlib.pyplot as plt
import numpy as np
import P2L3


def plot(resid, parms, fdv):
    '''
    A function to plot the result of TestP2L1.py

    The customizability of this function is fairly limited without editing.

    L2Norm is defined as the norm of the vector residuals array.
    '''

    EV, IV, vols, Mzev, Mziv = P2L3.P2L3(parms, fdv)
    tot = IV[0:2]*vols[0] + EV[0:2]*vols[1] + EV[2:4]*vols[2]

    residNorm = np.linalg.norm(resid)

    fig, ax = plt.subplots(2, figsize = (12,6))

    # Top plot
    ax[0].title.set_text(
        fdv['Name'] + ': L2Norm = ' + str(residNorm)
        )
    ax[0].plot(fdv['taxis'], Mziv[0]*vols[0], color = 'green', linestyle = 'solid', label = 'Pyriv')
    ax[0].plot(fdv['taxis'], Mzev[0]*vols[1], color = 'green', linestyle = 'dashdot', label = 'Pyree')
    ax[0].plot(fdv['taxis'], Mzev[2]*vols[1], color = 'green', linestyle = 'dotted', label = 'Pyric')
    ax[0].plot(fdv['taxis'], Mzev[1]*vols[2], color = 'blue', linestyle = 'solid', label = 'Lacee')
    ax[0].plot(fdv['taxis'], Mzev[3]*vols[2], color = 'blue', linestyle = 'dotted', label = 'Lacic')

    ax[0].grid()
    ax[0].legend()
    ax[0].set_ylabel('Mz')
    ax[0].set_xlabel('Time (s)')

    # Bottom plot
    ax[1].title.set_text(
        'kpl = ' + str(parms[0]) + ', kve = ' + str(fdv['knowns']['kve']) + ', vb = ' + str(fdv['knowns']['vb'])
        )
    ax[1].plot(fdv['taxis'], tot[0], color = 'green', label = 'Pyr')
    ax[1].plot(fdv['taxis'], tot[1], color = 'blue', label = 'Lac')
    ax[1].scatter(fdv['taxis'], fdv['data'][0], color = 'green', marker = 'x')
    ax[1].scatter(fdv['taxis'], fdv['data'][1], color = 'blue', marker = 'x')

    ax[1].grid()
    ax[1].legend()
    ax[1].set_ylabel('Mxy')

    plt.subplots_adjust(hspace = .5)
    plt.show()

    return fig
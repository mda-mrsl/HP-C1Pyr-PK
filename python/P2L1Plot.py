import matplotlib.pyplot as plt
import numpy as np

import P2L1

def plot(resid, parms, fdv):
    '''
    A function to plot the result of TestP2L1.py

    The customizability of this function is fairly limited without editing.

    L2Norm is defined as the norm of the vector residuals array.
    '''

    # Just to keep the indices straight
    pyr = 0
    lac = 1

    Mxy, Mtot = P2L1.P2L1(parms, fdv)

    resn = np.linalg.norm(resid)

    fig, ax = plt.subplots(2, figsize = (12, 6))
    fig.suptitle('Vertically Stacked plot')

    ax[0].title.set_text(fdv['Name'] + ': L2Norm = ' + str(resn))
    ax[0].plot(fdv['taxis'], Mtot[pyr,:], color = 'green', label = 'Pyr')
    ax[0].plot(fdv['taxis'], Mtot[lac,:], color = 'blue', label = 'Lac')

    ax[1].title.set_text('kpl = '+ str( parms[0] ))
    ax[1].plot(fdv['taxis'], Mxy[pyr,:], color = 'green', label = 'Pyr')
    ax[1].plot(fdv['taxis'], Mxy[lac,:], color = 'blue', label = 'Lac')
    ax[1].scatter(fdv['taxis'], fdv['data'][pyr], color = 'green', marker = 'x')
    ax[1].scatter(fdv['taxis'], fdv['data'][lac], color = 'blue', marker = 'x')

    ax[0].legend(fontsize = 'x-large')
    ax[0].grid()
    ax[0].set_ylabel('Mz = Mtot')

    ax[1].set_ylabel('Mxy')
    ax[1].grid()
    ax[1].set_xlabel('Time (s)', fontsize = 'x-large')

    plt.subplots_adjust(hspace = .5)
    plt.show() 

    return fig
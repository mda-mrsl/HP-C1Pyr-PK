import matplotlib.pyplot as plt
import numpy as np
import P2L2

def plot(resid, parms, fdv):
    '''
    A function to plot the result of TestP2L1.py

    The customizability of this function is fairly limited without editing.

    L2Norm is defined as the norm of the vector residuals array.
    '''
    pyr = 0 
    lac = 1

    EV , IV, vb, Mzev, Mziv = P2L2.P2L2(parms, fdv)
    tot = IV*vb + (1-vb)*EV
    Mz = (Mziv*vb) + Mzev*(1-vb)

    residNorm = np.linalg.norm(resid)

    fig, ax = plt.subplots(2, figsize = (12, 6))
    fig.suptitle('Vertically stacked plot')

    ax[0].title.set_text(fdv['Name'] + ': L2 Norm = ' + str(residNorm))
    ax[0].plot(fdv['taxis'], Mz[pyr], color = 'green', label = 'Pyr')
    ax[0].plot(fdv['taxis'], Mz[lac], color = 'blue', label = 'Lac')

    ax[1].title.set_text('kpl = ' + str(parms[0])
        + ', kve = ' + str( fdv['knowns']['kve'] )
        + ', vb = ' + str( vb ))
    ax[1].plot(fdv['taxis'], tot[pyr], color = 'green', label = 'Pyr')
    ax[1].plot(fdv['taxis'], tot[lac], color = 'blue', label = 'Lac')
    ax[1].scatter(fdv['taxis'], fdv['data'][pyr], color = 'green', marker = 'x')
    ax[1].scatter(fdv['taxis'], fdv['data'][lac], color = 'blue', marker = 'x')

    ax[0].legend(fontsize = 'x-large')
    ax[0].grid()
    ax[0].set_ylabel('Mz')
    
    ax[1].set_xlabel('Time (s)')
    ax[1].grid()
    ax[1].set_ylabel('Mxy')

    plt.subplots_adjust(hspace = .5)
    plt.show()

    return fig
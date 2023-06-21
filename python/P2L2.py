import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats

import helperFunctions as hf

def P2L2(vars, fdv):
    '''
    input:
        vars:               Value for parameter to be fit as array of size (2,)
        fdv:                Dictionary with the following keys
            fdv['fitvarNames']      Names of variables to be fit as string
            fdv['knowns']           Dictionary object containing known values
                                        (see TestP2L1.py for default)
            fdv['ntp']              Number of time points as int
            fdv['NSeg']             Number of segments per time point as int
                                        Note: for version 0.0.1, NSeg = 1
            fdv['NFlips']           Total number of excitations = ntp * NSeg as int
            fdv['FlipAngle']        Numpy array of excitation angles in degrees size = (2, NFlips)
            fdv['TR']               Single valued numpy array of repetition time size = (1, NFlips)
            fdv['data']             Observed data in the format [ [pyr], [lac] ] with size = (2, ntp)
            fdv['taxis']            Time axis used for the data with size = (1, NFlips)
            fdv['UseVIF']           Bool type. True if VIFs have been measured
            fdv['VIFP']             **Mz** VIF for pyruvate with shape = (1, NFlips)
            fdv['VIFL']             **Mz** VIF for lactate with shape = (1, NFlips), usually zeros
            fdv['Name']             Title for the output. Only affects plot title on TestP2L1.py
            fdv['verbose']          Bool type. True plots results of P2L1 directly

    output:
        Mxyev:              Extravascular observed with shape = (2, NTP) as [ [pyr], [lac] ]
        Mxyiv:              Intravascular observed with shape = (2, NTP) as [ [pyr], [lac] ]    
        vb:                 vascular blood vol fraction
        Mzev:               Extravascular longitudinal with shape = (2, NTP) as [ [pyr], [lac] ]
        Mziv:               Intravascular longitudinal with shape = (2, NTP) as [ [pyr], [lac] ]

    Note: Mzev, Mziv reflects the z-directed magnetization before the excitation pulse,
        while Mxy reflects transverse magnetization due to pulse

    Total observed signal is [Pyr, Lac] (1-vb)*Mxyev + vb*Mxyiv
    '''

    eps = 1e-6

    # Unpack fit variables into better naming scheme
    if len(vars) != 2:
        print(
            'P2L2 requires two parameters to be fit. Please enter vars as list with shape (1, 2)'
        )
        return 1

    i = 0
    for item in fdv['fitvarNames']:
        fdv[item] = vars[i]
        i += 1

    kvedve = fdv['knowns']['kve'] / (1 - fdv['knowns']['vb'])

    # Assuming that the VIF is given
    if fdv['UseVIF']:
        ff = np.array( [ fdv['VIFScale']*fdv['VIFP'], fdv['VIFScale']*fdv['VIFL'] ] )
    else:
        ff = np.array( [ fdv['VIFScale']*stats.gamma(
            fdv['taxis'], fdv['knowns']['Gam1'], fdv['knowns']['Gam2']), np.zeros(1, fdv['FlipAngle']) ] )

    # IV magnetization at each segment
    MivSeg = ff
    MxyivSeg = ff *  np.sin( np.radians(fdv['FlipAngle']) )

    # Initialize return variables
    Mxyev = np.zeros( (2,fdv['ntp']) )
    Mzev = np.zeros( (2,fdv['ntp']) )
    vb = fdv['knowns']['vb']
    Mxyiv = np.zeros( (2,fdv['ntp']) )
    Mziv = np.zeros( (2,fdv['ntp']) )

    #Check for trivial cases
    if (vb < eps):
        print('\nNo signal acquired\n')
        return (Mxyev, Mxyiv, vb, Mzev, Mziv)
    elif vb > (1-eps):
        print('\nAll vessel signal\n')
        # If segmented scquisition with NFlips > NTP, where VIF has
        #   been measured or interpolated to each excitation
        #   combine segments: ** THIS ASSUMES EQUAL WEIGHTING **
        for i in range(0, fdv['ntp']):
            Mxyiv[:,i] = np.mean( MxyivSeg[:, (i-1)*fdv['Nseg']+1 : i*fdv['NSeg']], 2 )
            Mziv[:, i] = np.mean(MivSeg[:, (i-1)*fdv['NSeg']+1 : i*fdv['NSeg']], 2)
        return 

    # Calculate non trivial result:
    # [Pev, Lev] = [Pev0, Lev0]*exp(At)+kpl*integral(0,t:exp(A*(t-T))*Piv(T), Liv(T)*dT)
    # A is a matrix with shape (2,2)
    a11 = - ( kvedve + fdv['kpl'] + (1/fdv['knowns']['T1Pyr']) )
    a12 = fdv['knowns']['klp']
    a21 = fdv['kpl']
    a22 = - ( kvedve + fdv['knowns']['klp'] + (1/fdv['knowns']['T1Lac']) )
    A = [ [a11, a12], [a21, a22] ]
    
    dD, P = np.linalg.eig(A)
    # dD should be array of eigenvalues with shape (1,2)
    # P is a matrix containing the eigenvectors of A

    ### Calculate signal evolution in TR ###
    MzevSegIC = np.array( [fdv['knowns']['P0'], fdv['knowns']['L0']] ) # Set initial conditions
    for i in range(fdv['ntp']):
        MxyevSeg = np.zeros( (2, fdv['NSeg']) )
        MzevSeg = np.zeros( (2, fdv['NSeg']) )
        for j in range(fdv['NSeg']):
            iSeg = (i-1)*(fdv['NSeg']+j)
            TR = fdv['TR'][0][iSeg]

            # First account for signal already present and its evolution
            #  Longitudinal magnetization available at the start of each segment
            MzevSeg[:,j] = MzevSegIC
            # Signal observed at each excitation, at start of segment
            MxyevSeg[:,j] = np.multiply(
                MzevSegIC, np.sin( np.radians(fdv['FlipAngle'][:,iSeg]) )
            )

            # At the end of this TR, Mz evolves to contribute to IC for next
            if iSeg<fdv['NFlips']: # Don't calculate after last datapoint
                MzevSeg1 = np.exp(dD*TR).flatten() * hf.mrdivide(P, np.multiply(
                        MzevSegIC, np.cos( np.radians(fdv['FlipAngle'][:, iSeg]) )
                    ).T ).flatten()


            # Now calculate new spins flowing into the system
            #  Assume piecewise linear VIF and diagonalize
            dff1 = hf.mrdivide(P, ff[:,iSeg]) # Diag VIF @ start of TR
            dff2 = hf.mrdivide(P, ff[:,iSeg+1]) # Diag VIF @ end of TR            
            
            # Get slope and y-intercept for diagonalized pyruvate forcing function
            b = dff1
            m = (dff2-dff1) / TR

            # At end of this TR inflowing spins will cause
            MevSeg2a = - (b.flatten() / dD) * (1-np.exp(dD*TR))
            MevSeg2b = np.multiply(m.flatten(),
                (-np.divide(TR, dD)- ((1/dD)/dD) ) + (np.exp(dD*TR) * ((1/dD)/dD)) )
            
            # Total signal at end of TR is combination of inflowing and already present signals
            MzevSegIC = np.matmul(P, (MzevSeg1 + kvedve*(MevSeg2a+MevSeg2b)))



        if fdv['NSeg'] > 1:
            # In version 0.0.1, fdv['NSeg']>1 is not supported
            print('Script does not work for segmented acquisition yet')
        else:
            Mxyev[:,i-1] = np.reshape(MxyevSeg, (2,) )
            Mzev[:,i-1] = np.reshape(MzevSeg, (2,) )

            Mxyiv[:,i-1] = MxyivSeg[:,i-1]
            Mziv[:,i-1] = MivSeg[:,i-1]
    ### END OF CALCULATION LOOP ###

    if fdv['verbose']:
        fig, ax = plt.subplots(figsize = (12, 6))

        ax.plot(fdv['taxis'], Mxyev[0,:], color = 'green')
        ax.plot(fdv['taxis'], Mxyev[1,:], color = 'blue')

        ax.grid()

        plt.show()

    return (Mxyev, Mxyiv, vb, Mzev, Mziv)
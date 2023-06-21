import numpy as np 
import matplotlib.pyplot as plt
import scipy.stats as stats

import helperFunctions as hf

def P2L3(vars, fdv):
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
        Mxyev:              Extravascular observed with shape = (4, NTP) as [ [pyr], [lac] ]
        Mxyiv:              Intravascular observed with shape = (4, NTP) as [ [pyr], [lac] ]    
        vols:               Vascular blood vol fractions [vb, ve, vc]
        Mzev:               Extravascular longitudinal with shape = (4, NTP) as [ [pyr], [lac] ]
        Mziv:               Intravascular longitudinal with shape = (4, NTP) as [ [pyr], [lac] ]

    Note: Mz reflects the z-directed magnetization before the excitation pulse,
        while Mxy reflects transverse magnetization due to pulse

    Total observed signal is [Pyr, Lac] = IV[0:2]*vols[0] + EV[0:2]*vols[1] + EV[2:4]*vols[2]
    '''

    eps = 1e-6

    # Check if user entered the right number of vars
    if len(vars) != 2:
        print('Please enter vars as a list or array with shape (2, 1)')
        return 1
    
    # Unpack variables into better naming scheme
    i = 0
    for item in fdv['fitvarNames']:
        fdv[item] = vars[i]
        i += 1
    
    # Check if Kecl is given
    if fdv['knowns'].get('kecl') == None:
        fdv['knowns']['kecl'] = fdv['knowns']['kecp']

    if fdv.get('FA') != None:
        fdv['FlipAngle'] = fdv['FA'] * np.ones( np.shape(fdv['FlipAngle']) )

    vb = fdv['knowns']['vb']
    ve = (1-vb)*fdv['knowns']['vef']
    vc = 1-vb-ve
    vols = np.array( [vb, ve, vc] )
    R1Pyr = 1 / fdv['knowns']['T1Pyr']
    R1Lac = 1 / fdv['knowns']['T1Lac']
    kvedve = fdv['knowns']['kve'] / ve
    kecpdve = fdv['knowns']['kecp'] / ve
    kecpdvc = fdv['knowns']['kecp'] / vc
    kecldve = fdv['knowns']['kecl'] / ve
    kecldvc = fdv['knowns']['kecl'] / vc

    # Assuming that the VIF is given
    if fdv['UseVIF']:
        ff = np.array(
            [ fdv['VIFScale']*fdv['VIFP'].flatten(), 
             fdv['VIFScale']*fdv['VIFL'].flatten(), 
             np.zeros( (1, fdv['NFlips']) ).flatten(),
             np.zeros( (1, fdv['NFlips']) ).flatten() 
            ]
        )
    else:
        ff = np.array(
            [ fdv['VIFScale']*stats.gamma.pdf(a=fdv['knowns']['Gam1'], scale=fdv['known']['Gam2'], x = fdv['taxis']).flatten(),
            np.zeros( (1, fdv['NFlips']) ).flatten(),
            np.zeros( (1, fdv['NFlips']) ).flatten(),
            np.zeros( (1, fdv['NFlips']) ).flatten()
            ]
        )

    # IV magnetization at each segment
    MzivSeg = ff
    MxyivSeg = ff * np.sin( np.radians(fdv['FlipAngle']) )

    # Initialize return variables
    Mxyev = np.zeros( (4, fdv['ntp']) )
    Mxyiv = np.zeros( (4, fdv['ntp']) )
    Mzev = np.zeros( (4, fdv['ntp']) )
    Mziv = np.zeros( (4, fdv['ntp']) )
    

    # Calculate non-trivial result
    #  Diff Eq: y'(t) = A y(t) + ff(t)
    #  Soln: y(t) = exp(A*t) * integral(0,t: exp(-A*T)*ff(T) dT)
    #  A is 4x4 matrix
    a11 = - (kvedve + kecpdve + R1Pyr)
    a13 = kecpdve
    a22 = - (kvedve + kecldve + R1Lac)
    a24 = kecldve
    a31 = kecpdvc
    a33 = - (kecpdvc + fdv['kpl'] + R1Pyr)
    a34 = fdv['knowns']['klp']
    a42 = kecldvc
    a43 = fdv['kpl']
    a44 = - (kecldvc + fdv['knowns']['klp'] + R1Lac)

    A = [ [a11, 0, a13, 0],
        [0, a22, 0, a24],
        [a31, 0, a33, a34],
        [0, a42, a43, a44] ]

    # Diagonalize to permit matrix integral
    dD, P = np.linalg.eig(A)
    # dD is (4,1) array containing eigenvalues
    # P is 2D array containing matrix of eigenvecters

    ### Calculate signal evolution in TR ###
    MzevSegIC = np.array(
        [fdv['knowns']['Pe0'], fdv['knowns']['Le0'], fdv['knowns']['Pi0'], fdv['knowns']['Li0']]
    ) # Set initial conditions
    for i in range( fdv['ntp'] ):
        MxyevSeg = np.zeros( (4, fdv['NSeg']) )
        MzevSeg = np.zeros( (4, fdv['NSeg']) )
        for j in range( fdv['NSeg'] ):
            iSeg = (i-1) * fdv['NSeg'] + j

            TR = fdv['TR'][0,iSeg]

            # First account for EV signal already present and its evolution
            #  Longitudinal magnetization available at the start of each segment
            MzevSeg[:, j] = MzevSegIC
            # Signal observed at each excitation, at start of segment
            MxyevSeg[:,j] = MzevSegIC * np.sin( np.radians(fdv['FlipAngle'][:, iSeg]) )

            #At the end of this TR, Mz evolves to contribute to IC for the next segment
            if iSeg < fdv['NFlips']: # Don't calculate after last datapoint

                MzevSeg1 = np.exp(dD*TR) * hf.mrdivide(
                    P, 
                    MzevSegIC * np.cos( np.radians(fdv['FlipAngle'][:, iSeg]) )
                ).flatten()

                # Now calculate the new spins flowing into the system
                #  Assume piecewise linear VIF. Diagonalize:
                dff1 = hf.mrdivide(P, ff[:, iSeg])
                dff2 = hf.mrdivide(P, ff[:,iSeg+1])

                # Get the slope and y-intercept for diagonalized forcing function
                b = dff1.flatten()
                m = ( (dff2-dff1)/TR ).flatten()

                # At the end of the segment inflowing spins will cause
                MzevSeg2a = (-b/dD) * ( 1-np.exp(dD*TR) )
                MzevSeg2b = m * ( (-TR/dD-((1/dD)/dD)) + np.exp(dD*TR)*((1/dD)/dD) )

                # Total signal at the end of TR equals the IC for the next TR
                MzevSegIC = np.matmul( P,(MzevSeg1 + kvedve*(MzevSeg2a+MzevSeg2b)) )

        if fdv['NSeg'] > 1:
            # In version 0.0.1, fdv['NSeg']>1 is not supported
            print('Segmented acquisition not set up yet')
        else:
            Mxyev[:,i-1] = np.reshape( MxyevSeg, (4,) )
            Mzev[:, i-1] = np.reshape( MzevSeg, (4,) )

            Mxyiv[:, i] = MxyivSeg[:, i]
            Mziv[:, i] = MzivSeg[:, i]
    ### END OF CALCULATION LOOP ###

    if fdv['verbose']:

        fig, ax = plt.subplots(figsize = (12,6))

        ax.plot(fdv['taxis'], Mxyev[0], color = 'green', linestyle = 'solid', label = 'Pyrev')
        ax.plot(fdv['taxis'], Mxyev[2], color = 'green', linestyle = 'dotted', label = 'Lacev')
        ax.plot(fdv['taxis'], Mxyev[1], color = 'blue', linestyle = 'solid', label = 'Pyric')
        ax.plot(fdv['taxis'], Mxyev[3], color = 'blue', linestyle = 'dotted', label = 'Lacic')

        ax.set_xlabel('Time (s)')
        ax.grid()
        ax.legend()

        plt.show()

    return (Mxyev, Mxyiv, vols, Mzev, Mziv)
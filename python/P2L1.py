import numpy as np
import matplotlib.pyplot as plt

def P2L1(vars, fdv):

    '''
    input:
        vars:               Value for parameter to be fit as array of size (1,)
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
            fdv['Name']             Title for the output. Only affects plot title on TestP2L1.py
            fdv['verbose']          Bool type. True plots results of P2L1 directly

    output:
        Mxy:                Mxy observed [ [Pyr], [Lac] ] with shape = (2, ntp)
        Mz:                 Mz [ [Pyr], [Lac] ] with shape = (2, ntp)

    Note: Mz is the magnetization before the excitation pulse,
        while Mxy is the transverse magnetization due to the pulse

    Total observed signal is [ [Mxy], [Mz] ]
    '''

    # Unpack fit variables into better naming scheme
    i = 0
    for item in fdv['fitvarNames']:
        fdv[item] = vars[i]
        i += 1
        
    # Data from each excitation
    PxySeg = fdv['data'][0]
    PzSeg = np.divide( PxySeg, np.sin( np.radians(fdv['FlipAngle'][0]) ))
    
    # Initialize return variables 
    Mxy = np.zeros( (2, fdv['ntp']) )
    Mz = np.zeros( (2, fdv['ntp']) )

    # Set up for the following equation
    #  Lac = Lac0 * exp(At) + kpl * inte(0,t: exp(A*(t-T)) *[Pyr(T)] *dT)
    A = -(fdv['knowns']['klp'] + (1/fdv['knowns']['T1Lac']))

    #### Caluclate signal evolution in each TR ###
    LzSegIC = fdv['knowns']['L0']  # Initial condition before acquisition starts
    for i in range(0, fdv['ntp']):
        LxySeg = np.zeros( (1, fdv['NSeg']) )
        LzSeg = np.zeros( (1, fdv['NSeg']) )
    
        for j in range(0, fdv['NSeg']):
            iSeg = (i)*(fdv['NSeg'])+j
            TR = fdv['TR'][0][iSeg]
    
            # First account for signal already in the slice and its evolution

            # Longitudinal M available at start of each cycle
            LzSeg[j] = LzSegIC
            # Signal observed at start of each cycle
            LxySeg[j] = np.multiply(LzSegIC, np.sin( np.radians(fdv['FlipAngle'][1, iSeg]) ))

            # Evolution of this cycle becomes the IC of the next
            if iSeg < fdv['NFlips']: # No evolution after last datapoint
                LzSeg1 = np.exp(A*TR)*LzSegIC*np.cos(
                    np.radians(fdv['FlipAngle'][1, iSeg] ))

                # Now account for new spins flowing nto the system during TR

                # Obtain parameters for linear pyruvate forcing function
                b = PzSeg[iSeg-1]
                m = (PzSeg[iSeg]-PzSeg[iSeg-1]) / TR 

                # Contribution from inflowing spins during this TR
                LzSeg2 = (((m/A + b)* np.exp(A*TR)) - (m*(TR+1/A)) - b)/A

                #Total signal at the end of TR is sum of inflowing and already present
                LzSegIC = LzSeg1 + fdv['kpl'] * LzSeg2
                
        
        if fdv['NSeg']>1:
            # In version 0.0.1, fdv['NSeg']>1 will never be evaluated
            print('Segmented acquisition not available in version 0.0.1')
            #Pyr
            Mxy[0, i-1] = np.mean(
                PxySeg[range(fdv['NSeg']*(i-1)+ 1, i*fdv['NSeg'])], 2)
            Mz[0, i-1] = np.mean( 
                PzSeg[range(fdv['NSeg']*(i-1)+ 1, i*fdv['NSeg'])], 2)

            #Lac
            Mxy[1:i-1] = np.mean(LxySeg)
            Mz[1:i-1] = np.mean(LzSeg)
        else:
            Mxy[:,i-1] = [PxySeg[i-1], LxySeg[0][0]]
            Mz[:,i-1] = [PzSeg[i-1], LzSeg[0][0]]
    ### END OF CALCULATION LOOP ###

    if fdv['verbose']:
        fig, ax = plt.subplots(figsize = (12, 6))

        ax.plot(fdv['taxis'], Mxy[0], color = 'green', label = 'Pyr')
        ax.plot(fdv['taxis'], Mxy[1], color = 'blue', label = 'Lac')

        ax.grid()
        ax.legend(fontsize = 'x-large')

        ax.set_xlabel('Time (s)', fontsize = 'x-large')

        plt.show()

    return (Mxy, Mz)
import numpy as np
import scipy.stats as stats
import scipy.optimize as opt

import P2L2
import P2L2Err
import P2L2Plot


def TestP2L2(parms=None, keys=None, values=None, knowns=None, knownVals=None):
    '''
    Script to generate synthetic data and verify that P2L2.py works as expected

    Optional inputs:
        parms                   Changes the default parameters for the calculation
        keys                    Keys for which the values of the dictionary fdv should be changed
        values                  Corresponding values to the keys of the dictionary to change
        knowns                  Names of known variables to be changed from default
        knownVals               Corresponding known variables to be changed from default

    Output:                     pyplot figure
    '''

    ### Generate synthetic data ###
    # Setting default parameters for single shot typical MDA
    #  For description of fdv keys call help(P2L2)
    fdv = {}

    fdv['fitvarNames'] = ['kpl', 'VIFScale']

    fdv['knowns'] = {
        'kve': .02,
        'vb': .05,
        'T1Pyr': 43,
        'T1Lac': 33,
        'klp': 0,
        'Gam1': 2.8,
        'Gam2': 4.5,
        'tdel': 10,
        'P0': 0,
        'L0': 0
    }
    # Check for user entered parameters
    if knowns != None:
        i = 0
        for item in knowns:
            fdv['knowns'][item] = knownVals[i]
            i+=1

    # Independent parameters
    fdv['Name'] = 'HP 2-Compartment Model Demo'
    fdv['ntp'] = 64    
    fdv['NSeg'] = 1   
    fdv['TR'] = 2  
    fdv['verbose'] = False
    # Check for user entered parameters
    if keys != None:
        i = 0
        for item in keys:
            fdv[item] = values[i]
            i+=1

    # Describe acquisition scheme 
    fdv['NFlips'] = fdv['ntp'] * fdv['NSeg'] # Total number of excitations

    # Describe temporal sampling scheme
    fdv['TR'] = np.ones( (1, fdv['NFlips']) )*fdv['TR']
    fdv['taxis'] = np.cumsum(fdv['TR']) - fdv['TR'][0]

    # Describe excitation scheme
    fdv['FlipAngle'] = 20 * np.ones( (2, fdv['NFlips']) )

    # Describe vascular input function
    #  Test data is generated synthetically from gamma distribution
    if fdv.get('UseVIF') == None:
        fdv['UseVIF'] = True
    if fdv.get('VIFP') == None:
        fdv['VIFP'] = stats.gamma.pdf(a=2.8, scale=4.5, x=fdv['taxis']-6)
    if fdv.get('VIFL') == None:
        fdv['VIFL'] = np.zeros( (1, fdv['NFlips']) )[0]


    # Placeholder data to be fit
    fdv['data'] = np.ones( (2, fdv['NFlips']) )

    # Generate synthetic data
    if parms == None:
        parms = np.array( [0.1, 1000] )
    else:
        parms = np.array( parms )
    EV, IV, vb, Mzev, Mziv = P2L2.P2L2(parms, fdv) 
    fdv['data'] = IV*vb + EV*(1-vb)

    ### Create a fit for the generated data and try to recover parameters above ###
    LB = [0, 0]
    UB = [1, np.Inf]
    guess = (parms*10) * np.random.rand(2)

    # Optimize parameters by solving least squares problem for residuals
    fit = opt.least_squares(
        fun = P2L2Err.P2L2ErrOuter(fdv),
        x0 = guess,
        bounds = (LB, UB)
    )

    popt = fit['x']
    resid = fit['fun']

    ### Plot the fits ###
    fig = P2L2Plot.plot(resid, popt, fdv)

    return fig


if __name__ == '__main__':
    TestP2L2()


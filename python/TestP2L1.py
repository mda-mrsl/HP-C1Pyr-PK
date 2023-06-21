import numpy as np
import scipy.stats as stats
import scipy.optimize as opt

import P2L1
import P2L1Err
import P2L1Plot


def TestP2L1(parms=None, keys=None, values=None, knowns=None, knownVals=None):
    '''
    Script to generate synthetic data and verify that P2L1.py works as expected
    
        Optional inputs:
        parms                   Changes the default parameters for the calculation
        keys                    Keys for which the values of the dictionary fdv should be changed
        values                  Corresponding values to the keys of the dictionary to change
        knowns                  Names of known variables to be changed from default
        knownVals               Corresponding known variables to be changed from default

    Output:                     pyplot figure
    '''

    ### Generate synthetic data ###
    # Setting default parameters for single shot, typical MDA
    #  For description of fdv keys call help(P2L1) or see docstring for P2L1.py
    fdv = {}

    fdv['kpl'] = None
    fdv['fitvarNames'] = ['kpl']

    fdv['knowns'] = {
        'T1Lac': 33, 
        'klp': 0, 
        'L0': 0
        }
    # Check for user entered parameters
    if knowns != None:
        i = 0
        for item in knowns:
            fdv['knowns'][item] = knownVals[i]
            i+=1

    # Independent parameters
    fdv['Name'] = 'Test P2L1'
    fdv['ntp'] = 60
    fdv['NSeg'] = 1
    fdv['TR'] = 2
    fdv['FA'] = 20
    fdv['verbose'] = False
    # Check for user entered parameters
    if keys != None:
        i = 0
        for item in keys:
            fdv[item] = values[i]
            i+=1

    # Describe acquisition scheme
    fdv['NFlips'] = (fdv['ntp']*fdv['NSeg'])

    #Describe temporal sampling scheme
    fdv['TR'] = fdv['TR'] * np.ones( (1, fdv['NFlips']) )
    fdv['taxis'] = np.cumsum(fdv['TR']) - fdv['TR'][0]

    #Describe excitation scheme
    fdv['FlipAngle'] = fdv['FA']*np.ones( (2,fdv['NFlips']) )

    # To test, we generate synthetic data
    randomArray = stats.gamma.pdf(a=2.8, scale=4.5, x=fdv['taxis']-10)
    oscillatorArray = np.sin( np.radians( fdv['FlipAngle'][0] ))
    fdv['data'] = np.zeros( (2, fdv['NFlips']) )
    fdv['data'][0] = np.multiply( randomArray, oscillatorArray )


    # Check for user entered parameters
    if parms == None:
        parms = [.05] # This is Kpl only
    else:
        np.array(parms)
    Mxy, Mz = P2L1.P2L1(parms, fdv)
    fdv['data'] = Mxy 

    ### Create a fit for the generated data ###
    guess = [.01]
    LB = 0
    UB = np.Inf

    # Optimize parameters by solving least squares problem for residuals
    fit = opt.least_squares(
        fun = P2L1Err.P2L1ErrOuter(fdv),
        x0 = guess,
        bounds = (LB, UB)
    )

    popt = fit['x']
    resid = fit['fun']

    ### Plot the fits ###
    fig = P2L1Plot.plot(resid, popt, fdv)

    return fig

if __name__ == '__main__':
    TestP2L1()
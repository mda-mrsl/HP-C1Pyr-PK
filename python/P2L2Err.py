import numpy as np

import P2L2

def P2L2ErrOuter(fdv):
    '''
    Function to calculate the residuals for use with scipy least_squares() func

    Input:
        fdv:        For information about the keys in this dictionary, call help(P2L2)

    Output:
        P2L1Inner:   A function of only vars that plays nice with least_squares
    '''

    def P2L2Err2(vars):
        '''
        Function of one variable to calculate the residuals for use with scipy least_squares() func

        Input:
            vars:       A two element array of test parameters

        Output:
            resid:       A 1D array containing residuals 
        '''

        EV, IV, vb, Mzev, Mziv = P2L2.P2L2(vars, fdv)
        tot = IV*vb + (1-vb)*EV

        resid = fdv['data'] - tot

        resid[0] /= np.max(fdv['data'][0])
        resid[1] /= np.max(fdv['data'][1])

        resid = resid.reshape( (1, np.prod(np.shape(resid))) ).flatten()

        return resid 

    return P2L2Err2


def P2L2Err(fdv, vars):
    '''
    Function of two variables to calculate the residuals to display or for further processing

    Input:
        fdv:        For information about the keys in this dictionary, call help(P2L2)
        pamrs:      A two element array of test parameters

    Output:
        resid:      A 1D array containing residuals
    '''

    inner = P2L2ErrOuter(fdv)
    resid = inner(vars)
    
    return resid

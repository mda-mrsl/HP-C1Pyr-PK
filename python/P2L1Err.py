import numpy as np
import P2L1

def P2L1ErrOuter(fdv):
    '''
    Function to calculate the residuals for use with scipy least_squares() func

    Input:
        fdv:        For information about the keys in this dictionary, call help(P2L1)

    Output:
        P2L1Inner:   A function of only parms that plays nice with least_squares
    '''

    def P2L1Inner(parms):
        '''
        Function of one variable to calculate the residuals for use with scipy least_squares() func

        Input:
            parms:       A single element array of test parameters

        Output:
            resid:       A 1D array containing residuals   
        '''

        Mxy, Mz = P2L1.P2L1(parms, fdv)
        resid = fdv['data'][1] - Mxy[1]

        resid[0] /= np.max( fdv['data'][0] )

        return resid.flatten()
    
    return P2L1Inner


def P2L1Err(parms, fdv):
    '''
    Function of two variables to calculate the residuals to display or for further processing

    Input:
        fdv:        For information about the keys in this dictionary, call help(P2L1)
        pamrs:      A signle element array of test parameters

    Output:
        resid:      A 1D array containing residuals
    '''


    inner = P2L1ErrOuter(fdv)
    resid = inner(parms)

    return resid
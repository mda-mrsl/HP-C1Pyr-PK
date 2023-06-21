import numpy as np
import P2L3

def P2L3ErrOuter(fdv):
    '''
    Function to calculate the residuals for use with scipy least_squares() func

    Input:
        fdv:        For information about the keys in this dictionary, call help(P2L3)

    Output:
        P2L1Inner:   A function of only vars that plays nice with least_squares
    '''

    def P2L3ErrInner(parms):
            '''
            Function of one variable to calculate the residuals for use with scipy least_squares() func

            Input:
                vars:       A two element array of test parameters

            Output:
                resid:       A 1D array containing residuals 
            '''

            EV, IV, vols, Mzev, Mziv = P2L3.P2L3(parms, fdv)
            tot = IV[0:2]*vols[0] + EV[0:2]*vols[1] + EV[2:4]*vols[2]

            resid = fdv['data'] - tot
            # Normalize to max of measured data. May need to reweight
            resid[0] /= np.max( fdv['data'][0] )
            resid[1] /= np.max( fdv['data'][1] )

            resid = resid.reshape(( np.prod(np.shape(resid)), ))
            
            return resid

    return P2L3ErrInner


def P2L3Err(parms, fdv):
    '''
    Function of two variables to calculate the residuals to display or for further processing       
    Input:
        fdv:        For information about the keys in this dictionary, call help(P2L3)
        pamrs:      A two element array of test parameters      
    Output:
        resid:      A 1D array containing residuals
    '''
      
    inner = P2L3ErrOuter(fdv)
    resid = inner(parms)

    return resid
      
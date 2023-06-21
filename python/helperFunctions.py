import numpy as np
import pandas as pd

def writeToDataframe(dataColumns, names = ['Time', 'Pyr', 'Lac'], path = 'pythonTest.txt'):
    '''
    Takes a list of columns and a list of headers and puts
     them into a pandas Dataframe to be written to a text file for
     later testing. This is mostly for debugging purposes.

    By default the DataFrame is written to a file called pythonTest.txt in the directory from which the script is run
    Inputs:
        dataColumns = List of lists or arrays. Each element should be one column
        names = List of strings or key for dataframes
        path = String containing location of directory to save the DataFrame

    Returns: 
        pd.DataFrame() object
    '''

    if len( dataColumns ) != len( names ):
        print('Error: missing column names.')
        return
    else:
        data = pd.DataFrame()
        for i in range( len( names ) ):
            data[names[i]] = dataColumns[i]

    data.to_csv(path, sep = '\t')

    return data

def mldivide(B, A):
    '''
    Designed to mimic the forward slash notation in MATLAB
    For a MATLAB statement of the form B/A = x or an equation of the form xA = B
    '''
    shapeA = np.shape(A)
    shapeB = np.shape(B)

    if len(shapeA) == 1:
        A = np.reshape(A, (1, shapeA[0]))

    if len(shapeB) == 1:
        B = np.reshape(B, (shapeB[0], 1))

    # lstsq returns a bunch of other information so the [0] index 
    #  selects out only the answer to the equation
    almostX = np.linalg.lstsq(A.T, B.T, rcond = None)[0]

    x = almostX.T

    return x

def mrdivide(B, A):

    '''
    Designed to mimic the backslash notation in MATLAB
    For a MATLAB statement of the form B\A = x or an equation of the form Ax = B
    '''

    shapeA = np.shape(A)
    shapeB = np.shape(B)

    if len(shapeA) == 1:
        A = np.reshape(A, (shapeA[0], 1))
    if len(shapeB) == 1:
        B = np.reshape(B, (1, shapeB[0]))
        
    # lstsq returns a bunch of other information so the [0] index
    #  extracts out only the answer to the equation
    x = np.linalg.lstsq(B, A, rcond = None)[0]

    
    return x

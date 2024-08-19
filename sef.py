import numpy as np

def add_slack_variable(A, c, irow, identity_val=1):
    newA = np.array([np.append(A[i], [(identity_val if i == irow else 0)]) for i in range(len(A))])
    newc = np.append(c, [0])
    return newA, newc

def to_equality_form(A, c, constraint_types):
    """
    Converts the given LP to Equality Form.
    """
    for row in range(len(A)):
        if constraint_types[row] == "=":
            continue
        identity_val = -1 if constraint_types[row] == ">=" else 1
        A, c = add_slack_variable(A, c, row, identity_val)
    return A, c

def sef(A, c, constraint_types):
    """
    Converts the given LP to Standard Equality Form.
    """
    

import numpy as np

def add_slack_variable(A, c, irow, identity_val=1):
    A[irow] += [identity_val]
    c += [0]
    for other_row in range(len(A)):
        if other_row != irow:
            other_row += [0]

def to_equality_form(A, c, constraint_types):
    """
    Converts the given LP to Equality Form.
    """
    for row in range(len(A)):
        if constraint_types[row] == ">=":
            add_slack_variable(A, c, row, identity_val=-1)
        elif constraint_types[row] == "<=":
            add_slack_variable(A, c, row, identity_val=1)

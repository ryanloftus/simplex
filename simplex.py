import argparse
import json
import numpy as np

FEASIBLE_BOUNDED = 0
UNBOUNDED = 1
INFEASIBLE = 2

def read_lp_file(filename: str):
    with open(filename) as f:
        return json.loads(f.read())
    
def nonnegative(v):
    for x in v:
        if x < 0:
            return False
    return True

def nonpositive(v):
    for x in v:
        if x > 0:
            return False
    return True

def certificate_of_optimality(A, c, B):
    return np.matmul(np.transpose(A[:][B]), c[B])

def canonical_form(A, b, c, z, B):
    AB = [[A[i][j] for j in range(len(A[i])) if j in B] for i in range(len(A))]
    AB_inv = np.linalg.inv(AB)
    A_tilde = np.matmul(AB_inv, A)
    b_tilde = np.matmul(AB_inv, b)
    z_tilde = np.dot(c[B], b_tilde) + z
    c_tilde = c - np.matmul(np.transpose(A_tilde), c[B])
    return A_tilde, b_tilde, c_tilde, z_tilde
    
def phase_one_simplex(A, b, c, z, B):
    """
    Phase 1: find a basis with non-negative objective value or,
    return the certificate of optimality if there is no such basis
    """

    canonical_form(A, b, c, z, B)

    return (INFEASIBLE, [], [], 0, B)

def phase_two_simplex(A, b, c, z, B):
    # convert LP to canonical form w.r.t. B
    A, b, c, z = canonical_form(A, b, c, z, B)
    
    # find the basic feasible solution
    x = [0] * len(c)
    for i in range(len(B)):
        x[B[i]] = b[i]

    # check for optimal solution
    if nonpositive(c):
        return {"outcome": FEASIBLE_BOUNDED, "optimal_bfs": x, "optimal_basis": B, "obj_val": z}
    
    # find which k enters the basis with Bland's Rule
    k = 0
    for i in range(len(c)):
        if c[i] > 0:
            k = i
            break
    if nonpositive(A[:][k]):
        return UNBOUNDED, x # TODO: certificate d

    # find which i leaves the basis with Bland's Rule
    t = min([b[i] / A[i][k] for i in range(len(b))])
    for i in range(len(b)):
        if b[i] / A[i][k] == t:
            B[i] = k
            break

    # iterate again with new basis
    return phase_two_simplex(A, b, c, z, B)

def auxillary_lp(A, b):
    dA = []
    for r in range(len(A)):
        identity_row = [0] * len(A)
        identity_row[r] = 1
        dA.append([a for a in A[r]] + identity_row)
    dc = [0] * len(A[0]) + [-1] * len(A)
    db = [bi for bi in b]
    dz = 0
    B = [len(A[0]) + i for i in range(len(A))]
    return dA, db, dc, dz, B

def two_phase_simplex(lp):
    c, z, A, b = lp['c'], lp['z'], lp['A'], lp['b']

    # First phase - find a basis or a certificate of infeasibility
    auxA, auxb, auxc, auxz, auxB = auxillary_lp(A, b)
    B, obj_val = phase_one_simplex(auxA, auxb, auxc, auxz, auxB)
    if obj_val < 0:
        print("Certificate of infeasibility:", certificate_of_optimality(auxA, auxc, B))

    # Second phase - find an optimal solution and certificate of optimality or a certificate of unboundedness
    res = phase_two_simplex(A, b, c, z, B)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--input_file',
        help='file that defines a linear program (LP) in Standard Equality Form (SEF)',
        required=True,
    )
    args = parser.parse_args()

    LP = read_lp_file(args.input_file)
    two_phase_simplex(LP) 

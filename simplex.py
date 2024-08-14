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
    return np.matmul(np.transpose(A[:,B]), c[B])

def choose_i_to_enter_basis(c):
    """
    Select k that should enter the basis in the next iteration by Bland's Rule
    """
    for i in range(len(c)):
        if c[i] > 0:
            return i
    return None

def choose_i_to_leave_basis(col, b):
    """
    Select k that should leave the basis in the next iteration by Bland's Rule
    """
    t = min([b[i] / col[i] for i in range(len(b)) if col[i] > 0])
    for i in range(len(b)):
        if col[i] > 0 and b[i] / col[i] == t:
            return i

def canonical_form(A, b, c, z, B):
    AB = A[:,B]
    AB_inv = np.linalg.inv(AB)
    A_tilde = np.matmul(AB_inv, A)
    b_tilde = np.matmul(AB_inv, b)
    z_tilde = np.dot(c[B], b_tilde) + z
    c_tilde = c - np.matmul(np.transpose(A_tilde), c[B])
    bfs = np.zeros_like(c)
    bfs[B] = b_tilde
    return A_tilde, b_tilde, c_tilde, z_tilde, bfs

def simplex(A, b, c, z, B):
    while True:
        A, b, c, z, bfs = canonical_form(A, b, c, z, B)

        # check for optimal solution
        if nonpositive(c):
            return {"outcome": FEASIBLE_BOUNDED, "optimal_bfs": bfs, "optimal_basis": B, "obj_val": z}
        
        # update the basis and check for unboundedness
        k = choose_i_to_enter_basis(c)
        if nonpositive(A[:,k]):
            return {"outcome": UNBOUNDED, "certificate_x": bfs, "certificate_d": np.zeros_like(bfs)} # TODO: certificate d
        B[choose_i_to_leave_basis(A[:,k], b)] = k

def auxillary_lp(A, b):
    dA = []
    for r in range(len(A)):
        identity_row = [0] * len(A)
        identity_row[r] = 1
        dA.append([(a if b[r] >= 0 else -a) for a in A[r]] + identity_row)
    dc = np.array([0] * len(A[0]) + [-1] * len(A))
    db = np.abs(np.copy(b))
    dz = 0
    B = np.array([len(A[0]) + i for i in range(len(A))])
    return np.array(dA), db, dc, dz, B

def two_phase_simplex(lp):
    c, z, A, b = np.array(lp['c']), lp['z'], np.array(lp['A']), np.array(lp['b'])

    # First phase - find a basis or a certificate of infeasibility
    auxA, auxb, auxc, auxz, auxB = auxillary_lp(A, b)
    res = simplex(auxA, auxb, auxc, auxz, auxB)
    if res["obj_val"] < 0:
        res["outcome"] = INFEASIBLE
        res["certificate"] = certificate_of_optimality(auxA, auxc, res["optimal_basis"])
        return res
    B = res["optimal_basis"]

    # Second phase - find an optimal solution and certificate of optimality or a certificate of unboundedness
    res = simplex(A, b, c, z, B)
    if res["outcome"] == FEASIBLE_BOUNDED:
        res["certificate"] = certificate_of_optimality(A, c, res["optimal_basis"])
    return res

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--input_file',
        help='file that defines a linear program (LP) in Standard Equality Form (SEF)',
        required=True,
    )
    args = parser.parse_args()

    LP = read_lp_file(args.input_file)
    res = two_phase_simplex(LP)
    if res["outcome"] == INFEASIBLE:
        print("Certificate of infeasibility:", res["certificate"])
    elif res["outcome"] == UNBOUNDED:
        print("Certificate of unboundedness:")
        print("x = {x}".format(x=res["certificate_x"]))
        print("d = {d}".format(d=res["certificate_d"]))
    elif res["outcome"] == FEASIBLE_BOUNDED:
        print("Optimal solution:", res["optimal_bfs"])
        print("Optimal basis:", res["optimal_basis"])
        print("Optimal value:", res["obj_val"])
        print("Certificate of optimality:", res["certificate"])

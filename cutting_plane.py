import argparse
import simplex
import json
import numpy as np
import sef

def is_integral(x):
    for n in x:
        if np.ceil(n) != np.floor(n):
            return False
    return True

def find_cutting_plane(A, b):
    for row in range(len(A)):
        if np.floor(b[row]) != np.ceil(b[row]):
            new_row = np.array(([np.floor(a) for a in A[row]]))
            A = np.vstack([A, new_row])
            b = np.append(b, np.floor(b[row]))
            return A, b, "<="
    return None

def solve_ip(ip):
    c, z, A, b = np.array(ip['c']), ip['z'], np.array(ip['A']), np.array(ip['b'])
    num_int = len(c)
    while True:
        # First phase - find a basis or a certificate of infeasibility
        auxA, auxb, auxc, auxz, auxB = simplex.auxillary_lp(A, b)
        res = simplex.simplex(auxA, auxb, auxc, auxz, auxB)
        if res["obj_val"] < 0:
            return simplex.INFEASIBLE, None
        B = res["optimal_basis"]

        # Second phase - find an optimal solution and certificate of optimality or a certificate of unboundedness or a cutting plane
        res = simplex.simplex(A, b, c, z, B)
        if res["outcome"] == simplex.UNBOUNDED:
            return simplex.UNBOUNDED, None
        elif is_integral(res["optimal_bfs"][:num_int]):
            return simplex.FEASIBLE_BOUNDED, res["optimal_bfs"][:num_int]
        else:
            A, b, c, z, x = simplex.canonical_form(A, b, c, z, B)
            A, b, constraint_type = find_cutting_plane(A, b)
            constraint_types = ["="] * (len(A)-1) + [constraint_type]
            A, c = sef.to_equality_form(A, c, constraint_types)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--input_file',
        help='file that defines an integer program (IP) in Standard Equality Form (SEF)',
        required=True,
    )
    args = parser.parse_args()

    IP = simplex.read_lp_file(args.input_file)
    res = solve_ip(IP)
    print(res)

# simplex

Implementation of the simplex algorithm. Solves Linear Programs given as a system $Ax = b$ of linear constraints with objective function $c^Tx$ assuming $x_i \geq 0$. Constraint types can be specified in the input so the system can have inequalities as well as equalities. Input examples are shown in the json files in this repo.

This implementation uses the 2-Phase Simplex with Bland's Rule so termination is guarenteed. If the LP is infeasible, the script prints a certificate of infeasiblity. If the LP is unbounded, the script prints a certificate of unboundedness. If the LP has an optimal solution, the script prints an optimal solution, optimal basis, optimal value, and a certificate of optimality.

## Solving Linear Programs (LPs)

To run the script, use:

```
$ python3 simplex.py --input_file feasible_bounded_LP.json
```

## Solving Integer Programs (IPs)

To run the script, use:

```
$ python3 cutting_plane.py --input_file IP.json
```

# simplex

Implementation of the simplex algorithm. Solves Linear Programs given in Standard Equality Form. Input examples are shown in the json files in this repo.

This implementation uses the 2-Phase Simplex with Bland's Rule so termination is guarenteed. If the LP is infeasible, the script prints a certificate of infeasiblity. If the LP is unbounded, the script prints a certificate of unboundedness. If the LP has an optimal solution, the script prints an optimal solution, optimal basis, optimal value, and a certificate of optimality.

## Running the Script

To run the script, use:

```
$ python3 simplex.py --input_file feasible_bounded_LP.json
```

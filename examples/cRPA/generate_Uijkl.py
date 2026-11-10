"""Generate a five-orbital static Kanamori tensor <ij|W|kl>, in eV."""
from pathlib import Path

import numpy as np

norb = 5
U = 4.0
J = 1.0

tensor = np.zeros((norb,) * 4)
for m in range(norb):
    tensor[m, m, m, m] = U
    for n in range(norb):
        if m != n:
            tensor[m, n, m, n] = U - 2 * J
            tensor[m, n, n, m] = J
            tensor[m, m, n, n] = J

indices = np.indices(tensor.shape).reshape(4, -1).T
np.savetxt(Path("crpa/Uijkl_impurity0.dat"),
           np.column_stack((indices, tensor.ravel())),
           fmt=["%d"] * 4 + ["%.12g"], header="i j k l U", comments="")

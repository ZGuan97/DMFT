"""Average a static spatial-orbital tensor <ij|W|kl> into Kanamori U and J.

First line: i j k l U. Subsequent rows: i j k l value, with optional # comments.
U is the mean of U[m,m,m,m]; J is the mean of U[m,n,n,m] for m != n.
"""

from pathlib import Path

import numpy as np


tensor_file = Path("crpa/Uijkl_impurity0.dat")
output_file = Path("crpa/UJ.dat")
index_base = 0
energy_to_eV = 1.0  # Use 27.211386245 for input values in Hartree.

data = np.loadtxt(tensor_file, skiprows=1, ndmin=2)
indices = data[:, :4].astype(int) - index_base
norb = indices.max() + 1
tensor = np.zeros((norb, norb, norb, norb))
tensor[tuple(indices.T)] = data[:, 4] * energy_to_eV

m = np.arange(norb)
U = tensor[m, m, m, m].mean()
m, n = np.indices((norb, norb))
J = tensor[m, n, n, m][m != n].mean() if norb > 1 else 0.0

np.savetxt(output_file, [[0, U, J]], fmt=["%d", "%.12g", "%.12g"],
           header="impurity_id U J (eV)")

print(f"U = {U:.12g} eV, J = {J:.12g} eV -> {output_file}")

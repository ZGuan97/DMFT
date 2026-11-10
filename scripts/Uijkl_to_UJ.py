"""Average a static spatial-orbital tensor <ij|W|kl> into Kanamori U and J.

First line is a header. Rows: i j k l real [imag], with optional # comments.
U is the mean of U[m,m,m,m]; J is the mean of U[m,n,n,m] for m != n.
"""

from pathlib import Path

import numpy as np


tensor_file = Path("crpa/Uijkl_impurity0.dat")
output_file = Path("crpa/UJ.dat")
index_base = 0
energy_to_eV = 1.0  # Use 27.211386245 for input values in Hartree.
imag_tolerance = 1e-10  # eV, applied to the averaged U and J.

data = np.loadtxt(tensor_file, skiprows=1, ndmin=2)
indices = data[:, :4].astype(int) - index_base
norb = indices.max() + 1
values = data[:, 4].astype(complex)
if data.shape[1] == 6:
    values += 1j * data[:, 5]
tensor = np.zeros((norb, norb, norb, norb), dtype=complex)
tensor[tuple(indices.T)] = values * energy_to_eV

m = np.arange(norb)
U = tensor[m, m, m, m].mean()
m, n = np.indices((norb, norb))
J = tensor[m, n, n, m][m != n].mean() if norb > 1 else 0.0

if not (abs(U.imag) <= imag_tolerance and abs(J.imag) <= imag_tolerance):
    raise ValueError(f"U/J are not real: Im(U)={U.imag:.12g} eV, Im(J)={J.imag:.12g} eV")
U, J = U.real, J.real

np.savetxt(output_file, [[0, U, J]], fmt=["%d", "%.12g", "%.12g"],
           header="impurity_id U J (eV)")

print(f"U = {U:.12g} eV, J = {J:.12g} eV -> {output_file}")

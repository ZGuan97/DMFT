## Input files
Two kinds of input files, i.e., DFT and DMFT files, are required if one wants to do DFT+DMFT calculations. The DFT files include all input and output files of DFT, and they must be grouped together into the folder named `dft/` under the root directory where the DFT+DMFT calculations are carried out. The DMFT files include all parameters controlling the DFT+DMFT calculations.

### How to prepare DFT files
  -i) Create a folder named dft under the root working directory. 

  -ii) Do a scf DFT calculation in the `dft/` directory as introduced in the manual of FHI-aims and ABACUS. In the DFT calculations, several extra parameters are required in the DFT input files, e.g., the INPUT file in ABACUS and control.in file in FHI-aims. Thoese parameters are listed in the manual of FHI-aims and ABACUS with explanations. There will be a folder named outputs_to_DMFT, which contains all data need to be passed to DMFT by DFT, after the DFT calculation finished.

### How to prepare DMFT files
  - `DMFT.in:` This file controls how the DMFT calulation will be carried out. An example is presented as follows
    ```bash
    DFT_solver   AIMS
    calculation scf
    temperature   1160
    impurity_solver rutgers-cthyb
    magnetism     para
    max_charge_step 10
    max_DMFT_step     1
    max_DFT_step    10
    charge_mix_param  0.05
    delta_sigma 0.01
    delta_rho  1.0e-4
    MC_step      50000000
    energy_window -2.0 1.0
    local_symmetry 1
    ```
    All available parameters are explictly explained [here](list_of_parameters.md). Note that all parameters and their values are case-insensitive. The symbol `#` uses for annotating.


### cRPA interface files

The cRPA interface uses the following fixed paths relative to the DMFT calculation directory:

| File | Purpose |
| --- | --- |
| `crpa/space.json` | Correlated-space definition and orbital/impurity mappings |
| `crpa/projectors.dat` | Bloch-to-local projector matrices |
| `crpa/UJ.dat` | Scalar U and J for a single inequivalent impurity, in eV |
| `crpa/Uijkl_impurityN.dat` | Interaction tensor for inequivalent impurity `N`, numbered from zero |

These paths are defined in `src/para/crpa.h` and `src/para/crpa.cpp`.
The reader uses the independently bundled nlohmann/json library; it does not depend
on the ComCTQMC JSON implementation.
When `crpa true` is set in `DMFT.in`, `space.json` and `projectors.dat` are read during
initialization and projection setup respectively. DMFT uses the supplied projector instead of constructing NAO/CWF
projectors, both for SCF and for spectra. The existing `KS_bands_window()` selects
the bands from `projection_window`, or reads the saved window for restart and spectra.
Integration weights come from the existing DFT `k_weight.dat`.
The projector construction settings do not alter the imported matrices.
The spectral band window uses the existing DOS-window handling.
`crpa/UJ.dat` supplies scalar U and J for a single inequivalent impurity. They replace
the values read from `correlated_atoms.info` for that impurity and its equivalent atoms.
The existing interaction construction and double counting use these parameters.

The JSON structure is:

```json
{
  "window_bands": [[8, 9, 10, 11, 12]],
  "correlated_atoms": [
    {
      "atom_id": 0,
      "impurity_id": 0,
      "orbitals": ["dxy", "dyz", "dz2", "dxz", "dx2-y2"]
    }
  ]
}
```

`atom_id` is the zero-based atom index in the DFT structure.
`impurity_id` is the zero-based inequivalent impurity index.
`orbitals` contains the ordered orbital labels.
`window_bands` is required: one array per spin channel, listing zero-based KS
band indices in projector-row order. These indices define Q_W, the projector
onto the selected KS window. The same list applies at every k point. DMFT
uses these indices directly for eigenvalues, band mappings and window electron
count, including restarts; `projection_window` does not select bands in cRPA mode
and may be omitted. Existing inputs must add this field. Each spin list must be
nonempty, ascending, consecutive and within the DFT band range; the number of
lists must match DFT spins. These are input conventions; the reader does not validate them. K-dependent/disconnected windows are not supported.
`dos_window` still selects the spectra bands and must contain the imported window.

The order of `correlated_atoms` defines the correlated-atom indices. Each `orbitals`
array defines the local orbital order; concatenating these arrays defines the global
projector column order. Atom, impurity and orbital counts are derived from the contents.

The existing `correlated_atoms.info` and `symmetry.dat` are still read for shell,
interaction and magnetic metadata. Counts, global orbital offsets and impurity mappings used in the
calculation are taken from `space.json`. Local symmetry averaging is disabled for
the imported orbitals, so their supplied order and basis are retained.

### Projector file format

`crpa/projectors.dat` is whitespace-separated numeric text, without keywords or
comments. Its layout, in reading order, is:

The file contains only `real imag` pairs, ordered by spin, k point, band, then
global correlated orbital (orbital index varies fastest).
For each spin channel there are `nk * nbands[spin] * norb` pairs.

The number of spin channels and `nk` come from DFT, `nbands[spin]` and matrix row
order come from `space.json.window_bands`, and `norb` comes from `space.json`.
There are no dimensions or band indices in this numeric file; indices are supplied by `space.json.window_bands`. The in-memory `projectors`
structure stores only the matrices. Nonmagnetic and collinear up/down channels
are supported; SOC spinors are not supported.

Each pair represents the dimensionless complex matrix element
`tilde P[n,a](k,spin) = <psi_nk,spin | w_ak,spin>`; no conjugation or orthonormalization is applied on input. Columns must already be orthonormal.
The global orbital order is obtained by concatenating the orbital lists in `space.json`.
K-point coordinates and weights are not stored in this file; k-point indices follow
the existing DFT ordering, and integration uses the DFT weights.
Each matrix uses the project's row-major Eigen `Matrix` type and is accessed as
`crpa_projectors.P[spin][k](band, orb)`.
Only fixed-DFT calculations (`max_charge_step 1`) are supported.

For each atom block DMFT computes `P^dagger diag(epsilon) P` in Hartree and uses its
diagonal entries through the existing local-level interface. This integration retains
the current solver's diagonal local-level and local-self-energy approximation; it does
not add a general off-diagonal impurity model. Atomic overlap and eigenvector files
are not read by the imported-projector path.

### U and J file format

`crpa/UJ.dat` contains one data row with three columns: `impurity_id U J`.
The impurity index is 0, and U and J are in eV. Blank lines and comments beginning
with `#` are skipped. DMFT converts the values to Hartree when reading.
Use `interaction_type density_density` or `interaction_type kanamori` for this scalar
interface. Tensor conversion and multiple inequivalent impurities are not implemented.

## Output files
After the DFT+DMFT calculations finished, there will be following output files

  -`DMFT_running.log`: This is the log file of the DFT+DMFT calculations. It contains information about the DFT+DMFT calculations.

  -`dmft/`:This folder contains all outputs given by impurity solvers.

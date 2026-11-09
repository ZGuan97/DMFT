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
| `crpa/Uijkl_impurityN.dat` | Interaction tensor for inequivalent impurity `N`, numbered from zero |

These paths are defined in `src/para/crpa.h` and `src/para/crpa.cpp`.
The reader uses the independently bundled nlohmann/json library; it does not depend
on the ComCTQMC JSON implementation.
When `crpa true` is set in `DMFT.in`, `space.json` and `projectors.dat` are read and
stored during initialization. Interaction import, and using the external space and
projectors in place of the existing DMFT orbital construction, are not yet implemented.

The JSON structure is:

```json
{
  "correlated_atoms": [
    {
      "atom_id": 0,
      "impurity_id": 0,
      "orbitals": ["dxy", "dyz", "dz2", "dxz", "dx2-y2"]
    }
  ]
}
```

`atom_id` is the zero-based atom index in the DFT structure, not the index among
correlated atoms. Atom IDs must be distinct and nonnegative; checking their correspondence
to the DFT structure is deferred until the space is used in the calculation.
`impurity_id` values must be consecutive starting at zero; atoms sharing an impurity
must have identical ordered orbital labels. Labels are stored as supplied, without
interpreting them as atomic basis indices or applying rotations.
The order of `correlated_atoms` defines the correlated-atom indices. Each `orbitals`
array defines the local orbital order; concatenating these arrays defines the global
projector column order. Atom, impurity and orbital counts are derived from the contents.
The reader directly extracts these fields and derives the counts; it does not perform
additional index-range, duplicate-atom or equivalent-orbital consistency checks.

### Projector file format

`crpa/projectors.dat` is whitespace-separated numeric text, without keywords or
comments. Its layout, in reading order, is:

1. Three integers: `nspin nk norb`. `nspin` is 1 (nonmagnetic) or 2 (collinear
   up/down); SOC spinors are not supported. `norb` equals the total orbital count
   from `space.json`.
2. `nk` rows: `kx ky kz weight`, in the same k-point order as the DFT data.
   Coordinates are fractional in reciprocal lattice vectors including 2*pi;
   weights are dimensionless and normalized to sum to one, without a spin factor.
3. For each spin channel, in order:
   - One integer `nbands`, then `nbands` zero-based DFT band indices in matrix row order.
     The same band list applies to every k point in that spin channel.
   - `nk * nbands * norb` pairs `real imag`, ordered by k point, then band,
     then global correlated orbital (orbital index varies fastest).

Each pair represents the dimensionless complex matrix element
`P[n,a](k,spin) = <psi_nk,spin | w_ak,spin>`; no conjugation is applied on input.
The global orbital order is obtained by concatenating the orbital lists in `space.json`.
Each matrix uses the project's row-major Eigen `Matrix` type and is accessed as
`crpa_projectors.P[spin][k](band, orb)`.
Matrices are expected to satisfy `P^dagger P = I`; the reader neither orthogonalizes
nor rotates them. It retains the supplied k-point weights and band indices without
matching them against the DFT data at this stage. The same DFT wavefunctions and
their phase/band conventions must be used when these projectors are applied.

## Output files
After the DFT+DMFT calculations finished, there will be following output files

  -`DMFT_running.log`: This is the log file of the DFT+DMFT calculations. It contains information about the DFT+DMFT calculations.

  -`dmft/`:This folder contains all outputs given by impurity solvers.

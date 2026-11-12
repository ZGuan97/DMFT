#pragma once

#include <string>
#include <vector>
#include "../matrix_types.h"

namespace DMFT
{
  namespace crpa
  {
    // Paths relative to the DMFT calculation directory.
    extern const char space_file[];
    extern const char projector_file[];
    extern const char UJ_file[];

    // Single inequivalent impurity; file values are eV, returned values Hartree.
    void read_UJ(double& U, double& J);

    // Zero-based inequivalent impurity index.
    std::string interaction_file(int ineq);

    struct correlated_atom
    {
      int atom_id;
      int impurity_id;
      std::vector<std::string> orbitals;
    };

    struct corr_space
    {
      // Array order defines corr_id; concatenated orbital lists define global_orb.
      std::vector<correlated_atom> correlated_atoms;
      int nimpurities = 0;
      int norb = 0;
      // Zero-based KS indices, ordered by spin then projector row; defines Q_W.
      std::vector<std::vector<int>> window_bands;

      void read();
    };

    struct projectors
    {
      // Stores tilde P[spin][k](band, orb) = <psi_nk|w_ak>, dimensionless.
      std::vector<std::vector<DFT_plus_DMFT::Matrix>> P;

      void read(const std::vector<int>& nbands, int nk, int norb);
    };
  }
}

#include "crpa.h"

#include <algorithm>
#include <fstream>
#include <stdexcept>
#include <utility>

#include <nlohmann/json.hpp>

namespace DMFT
{
  namespace crpa
  {
    const char space_file[] = "crpa/space.json";
    const char projector_file[] = "crpa/projectors.dat";

    std::string interaction_file(int ineq)
    {
      return "crpa/Uijkl_impurity" + std::to_string(ineq) + ".dat";
    }

    void corr_space::read()
    {
      std::ifstream input(space_file);
      if(!input) throw std::runtime_error(std::string("cannot open ") + space_file);
      const auto document = nlohmann::json::parse(input);
      const auto atoms = document.at("correlated_atoms").get<nlohmann::json::array_t>();

      corr_space result;
      for(const auto& entry : atoms)
      {
        correlated_atom atom;
        atom.atom_id = entry.at("atom_id").get<int>();
        atom.impurity_id = entry.at("impurity_id").get<int>();
        atom.orbitals = entry.at("orbitals").get<std::vector<std::string>>();
        result.norb += atom.orbitals.size();
        result.nimpurities = std::max(result.nimpurities, atom.impurity_id + 1);
        result.correlated_atoms.push_back(atom);
      }
      *this = std::move(result);
    }

    void projectors::read()
    {
      std::ifstream input(projector_file);
      if(!input) throw std::runtime_error(std::string("cannot open ") + projector_file);
      projectors result;
      input >> result.nspin >> result.nk >> result.norb;

      result.kpoints.resize(result.nk);
      result.weights.resize(result.nk);
      for(int ik=0; ik<result.nk; ik++)
        input >> result.kpoints[ik][0] >> result.kpoints[ik][1]
              >> result.kpoints[ik][2] >> result.weights[ik];

      result.bands.resize(result.nspin);
      result.P.resize(result.nspin);
      for(int is=0; is<result.nspin; is++)
      {
        int nbands = 0;
        input >> nbands;
        result.bands[is].resize(nbands);
        for(int& band : result.bands[is]) input >> band;
        result.P[is].resize(result.nk);
        for(int ik=0; ik<result.nk; ik++)
        {
          result.P[is][ik].resize(nbands, result.norb);
          for(int iband=0; iband<nbands; iband++)
            for(int iorb=0; iorb<result.norb; iorb++)
            {
              double real = 0.0, imag = 0.0;
              input >> real >> imag;
              result.P[is][ik](iband, iorb) = std::complex<double>(real, imag);
            }
        }
      }
      *this = std::move(result);
    }
  }
}

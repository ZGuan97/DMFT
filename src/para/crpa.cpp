#include "crpa.h"
#include "../constants.h"

#include <algorithm>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <utility>

#include <nlohmann/json.hpp>

namespace DMFT
{
  namespace crpa
  {
    const char space_file[] = "crpa/space.json";
    const char projector_file[] = "crpa/projectors.dat";
    const char UJ_file[] = "crpa/UJ.dat";

    void read_UJ(double& U, double& J)
    {
      std::ifstream input(UJ_file);
      if(!input) throw std::runtime_error(std::string("cannot open ") + UJ_file);
      std::string line;
      int impurity = 0;
      while(std::getline(input, line))
      {
        line = line.substr(0, line.find('#'));
        std::istringstream row(line);
        if(row >> impurity >> U >> J) break;
      }
      U /= GLC::Hartree_to_eV;
      J /= GLC::Hartree_to_eV;
    }

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

    void projectors::read(const std::vector<int>& nbands, int nk, int norb)
    {
      std::ifstream input(projector_file);
      if(!input) throw std::runtime_error(std::string("cannot open ") + projector_file);
      projectors result;
      const int nspin = nbands.size();
      result.P.resize(nspin);
      for(int is=0; is<nspin; is++)
      {
        result.P[is].resize(nk);
        for(int ik=0; ik<nk; ik++)
        {
          result.P[is][ik].resize(nbands[is], norb);
          for(int iband=0; iband<nbands[is]; iband++)
            for(int iorb=0; iorb<norb; iorb++)
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

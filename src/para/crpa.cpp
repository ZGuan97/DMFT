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
  }
}

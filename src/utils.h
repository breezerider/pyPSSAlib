#include <memory>
#include <iostream>

#include <boost/config.hpp>
#include <boost/smart_ptr/make_unique.hpp>

#include <pssalib/datamodel/detail/Model.h>
#include <pssalib/datamodel/detail/Species.h>
#include <pssalib/datamodel/detail/Reaction.h>
#include <pssalib/datamodel/detail/SpeciesReference.h>

#ifndef PYPSSA_UTILS_H_
#define PYPSSA_UTILS_H_

#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

// inspired by https://stackoverflow.com/questions/30670909/boostpython-extract-c-value-from-numpy-ndarray-return-value
#define PY_ERRMSG(msg) \
  { PyErr_SetString(PyExc_TypeError, (boost::format("PY_ERRMSG(%1%:%2%): %3%") % (__FILE__) % (__LINE__) % (msg)).str().c_str()); };

/*
 * Prints reaction network
 */
std::ostream & printReactionNetwork(const pssalib::datamodel::detail::Model & model, std::ostream & os)
{
  os << "'" << model.toString() << "' model:\n\n";

  os << "Compartment volume: " << model.getCompartmentVolume() << "\n\n";

  os << "Initial species population:\n";

  for(UINTEGER si = 0; si < model.getSpeciesCount(); ++si)
  {
    const pssalib::datamodel::detail::Species * s = model.getSpecies(si);

    os << s->toString() << " = " << s->getInitialAmount() << "\t";
  }

  os << "\n\nReaction network:\n";

  for(UINTEGER ri = 0; ri < model.getReactionsCount(); ++ri)
  {
    const pssalib::datamodel::detail::Reaction * r = model.getReaction(ri);

    os << r->toString() << ": ";

    for(UINTEGER sri = 0; sri < r->getSpeciesReferencesCount(); ++sri)
    {
      const pssalib::datamodel::detail::SpeciesReference * sr = r->getSpeciesReferenceAt(sri);

      if(sri == r->getReactantsCount()) {
        if(r->isReversible()) {
          os << " <-" << r->getReverseRate() << "-" << r->getForwardRate() << "-> ";
        } else {
          os << " -" << r->getForwardRate() << "-> ";
        }
      } else if(sri > 0) {
        os << " + ";
      }

      os << sr->toString();
    }
    os << "\n";
  }

  return os;
}

/*
 * Compute the reaction rate for mass-action kinetics
 */
std::unique_ptr<double[]>
computeMassactioRates(
  const pssalib::datamodel::detail::Model & model,
  const unsigned int * paruPopulation
)
{
  REAL Vinv = model.getCompartmentVolume();
  if(Vinv <= 0.0)
    Vinv = 1.0;
  else
    Vinv = 1.0 / Vinv;

  std::unique_ptr<double[]> rates = boost::make_unique<double[]>(model.getReactionsCount());

  // std::cerr << "rates = [";

  for(UINTEGER ri = 0; ri < model.getReactionsCount(); ++ri)
  {
    const pssalib::datamodel::detail::Reaction * r = model.getReaction(ri);

    rates[ri] = r->getForwardRate(); // fill with reaction rates

    // compute mass-action rates
    for(UINTEGER sri = 0; sri < r->getReactantsCount(); ++sri)
    {
      const pssalib::datamodel::detail::SpeciesReference * sr = r->getSpeciesReferenceAt(sri);

      if(!sr->isReservoir() && (sri < r->getReactantsCount()))
        rates[ri] *= std::pow(REAL(paruPopulation[sr->getIndex()]) * Vinv, sr->getStoichiometryAbs());
    }

    // std::cerr << rates[ri] << ", ";;
  }
  // std::cerr << "]" << std::endl;

  return rates;
}

/*
 * Compute the value of mass-action system of ODEs at a given point
 */
void
computeODEs(
  const pssalib::datamodel::detail::Model & model,
  const unsigned int * paruPopulation,
  double * pardODEs
)
{
  std::unique_ptr<double[]> rates = computeMassactioRates(model, paruPopulation);

  // printReactionNetwork(model, std::cerr);
  //
  // std::cerr << "before odes = [";
  // for(UINTEGER si = 0; si < model.getSpeciesCount(); ++si)
  //   std::cerr << pardODEs[si] << ", ";
  // std::cerr << "]" << std::endl;

  for(UINTEGER si = 0; si < model.getSpeciesCount(); ++si)
    pardODEs[si] = 0.0;

  for(UINTEGER ri = 0; ri < model.getReactionsCount(); ++ri)
  {
    const pssalib::datamodel::detail::Reaction * r = model.getReaction(ri);

    for(UINTEGER sri = 0; sri < r->getSpeciesReferencesCount(); ++sri)
    {
      const pssalib::datamodel::detail::SpeciesReference * sr = r->getSpeciesReferenceAt(sri);

      // fill the ODEs vector
      if(!sr->isReservoir())
      {
        if(sri < r->getReactantsCount())
          pardODEs[sr->getIndex()] -= sr->getStoichiometryAbs() * rates[ri];
        else
          pardODEs[sr->getIndex()] += sr->getStoichiometryAbs() * rates[ri];
      }
    }
  }

  // std::cerr << "after odes = [";
  // for(UINTEGER si = 0; si < model.getSpeciesCount(); ++si)
  //   std::cerr << pardODEs[si] << ", ";
  // std::cerr << "]" << std::endl;
}

/*
 * Compute the Jacobian of mass-action system of ODEs at a given point
 */
void
computeJacobian(
  pssalib::datamodel::detail::Model & model,
  const unsigned int * paruPopulation,
  double * pardJ
)
{
    double Vinv = model.getCompartmentVolume();
    if(Vinv <= 0.0)
        Vinv = 1.0;
    else
        Vinv = 1.0 / Vinv;

    size_t lenJ = model.getSpeciesCount() * model.getSpeciesCount();
    for(UINTEGER si = 0; si < lenJ; ++si)
        pardJ[si] = 0.0;

    for(UINTEGER ri = 0; ri < model.getReactionsCount(); ++ri)
    {
        const pssalib::datamodel::detail::Reaction * r = model.getReaction(ri);

        for(UINTEGER sri = 0; sri < r->getReactantsCount(); ++sri)
        {
            const pssalib::datamodel::detail::SpeciesReference * srI = r->getSpeciesReferenceAt(sri);

            double jri = r->getForwardRate();

            if(srI->isReservoir()) break;

            for(UINTEGER srj = 0; srj < r->getReactantsCount(); ++srj)
            {
                const pssalib::datamodel::detail::SpeciesReference * srJ = r->getSpeciesReferenceAt(srj);

                if(srI->getIndex() != srJ->getIndex())
                    jri *= std::pow(REAL(paruPopulation[srJ->getIndex()]) * Vinv, srJ->getStoichiometryAbs());
                else if(srJ->getStoichiometryAbs() > 1)
                    jri *= double(srJ->getStoichiometryAbs()) * std::pow(REAL(paruPopulation[srJ->getIndex()]) * Vinv, srJ->getStoichiometryAbs() - 1);
            }

            for(UINTEGER srj = 0; srj < r->getSpeciesReferencesCount(); ++srj)
            {
                const pssalib::datamodel::detail::SpeciesReference * srJ = r->getSpeciesReferenceAt(srj);

                if(srJ->isReservoir()) continue;

                if(srj < r->getReactantsCount())
                    pardJ[srI->getIndex() + srJ->getIndex() * model.getSpeciesCount()] -= jri * double(srJ->getStoichiometryAbs());
                else
                    pardJ[srI->getIndex() + srJ->getIndex() * model.getSpeciesCount()] += jri * double(srJ->getStoichiometryAbs());
            }
        }
    }
}

/*
 * Compute the Q-term for Lyapunov Equation
 */
void
computeQ(
  pssalib::datamodel::detail::Model & model,
  const unsigned int * paruPopulation,
  double * pardQ
)
{
    std::unique_ptr<double[]> rates = computeMassactioRates(model, paruPopulation);

    size_t lenQ = model.getSpeciesCount() * model.getSpeciesCount();
    for(UINTEGER si = 0; si < lenQ; ++si)
      pardQ[si] = 0.0;

    for(UINTEGER ri = 0; ri < model.getReactionsCount(); ++ri)
    {
        pssalib::datamodel::detail::Reaction * r = model.getReaction(ri);

        // Scale with compartment volume
        rates[ri] *= model.getCompartmentVolume();

        for(UINTEGER sri = 0; sri < r->getSpeciesReferencesCount(); ++sri)
        {
            const pssalib::datamodel::detail::SpeciesReference * srI = r->getSpeciesReferenceAt(sri);

            for(UINTEGER srj = 0; srj < r->getSpeciesReferencesCount(); ++srj)
            {
                const pssalib::datamodel::detail::SpeciesReference * srJ = r->getSpeciesReferenceAt(srj);

                // fill the Lyapunov eq matrix Q
                if(!srI->isReservoir()&&!srJ->isReservoir())
                {
                    double dq = rates[ri] * srI->getStoichiometryAbs() * srJ->getStoichiometryAbs();
                    if((sri < r->getReactantsCount())^(srj < r->getReactantsCount()))
                        pardQ[srI->getIndex() + srJ->getIndex() * model.getSpeciesCount()] -= dq;
                    else
                        pardQ[srI->getIndex() + srJ->getIndex() * model.getSpeciesCount()] += dq;
                }
            }
        }
    }
}

#endif

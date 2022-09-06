/*****************************************************************************\
* (c) Copyright 2019 CERN for the benefit of the LHCb Collaboration           *
*                                                                             *
* This software is distributed under the terms of the GNU General Public      *
* Licence version 3 (GPL Version 3), copied verbatim in the file "COPYING".   *
*                                                                             *
* In applying this licence, CERN does not waive the privileges and immunities *
* granted to it by virtue of its status as an Intergovernmental Organization  *
* or submit itself to any jurisdiction.                                       *
\*****************************************************************************/
#ifndef TupleToolDOCA_H
#define TupleToolDOCA_H

#include "DecayTreeTupleBase/TupleToolBase.h"
#include "Kernel/IParticleTupleTool.h"
#include "LoKi/ChildSelector.h"

struct IDVAlgorithm;
struct IDistanceCalculator;

namespace LHCb {
  class Particle;
}

/** @class TupleToolDOCA TupleToolDOCA.h
 *
 * @brief compute DOCAs between particles in the decay chain DecayTreeTuple
 *
 * Thanks to Vanya for the pointer to LoKi::Child::Selectors.
 *
 * a small example:
 * @code
 *
 * Bu2JpsiKTuple.Decay    = '[B- -> ^(J/psi(1S)-> ^mu+ ^mu-) ^K- ]CC'
 * Bu2JpsiKTuple.addBranches( {
 *     "muon_p"   : '[B- ->  (J/psi(1S)-> ^mu+  mu-)  K- ]CC'
 *    ,"muon_m"  : '[B- ->  (J/psi(1S)->  mu+ ^mu-)  K- ]CC'
 *    ,"kaon_m"  : '[B- ->  (J/psi(1S)->  mu+  mu-) ^K- ]CC'
 *    ,"Jpsi"    : '[B- -> ^(J/psi(1S)->  mu+  mu-)  K- ]CC'
 *    ,"Bu"      : '[B- ->  (J/psi(1S)->  mu+  mu-)  K- ]CC'
 *    } )
 *
 * budocas = Bu2JpsiKTuple.Bu.addTupleTool('TupleToolDOCA')
 * budocas.Name       = [ "Kmu_OS", "Kmu_SS", "mumu" ]
 * budocas.Location1  = [ "[B- ->  (J/psi(1S)-> mu+  mu-)  ^K-]CC",
 *                        "[B- ->  (J/psi(1S)-> mu+  mu-)  ^K-]CC",
 *                        "[B- ->  (J/psi(1S)-> mu+  ^mu-)  K-]CC"]
 * budocas.Location2  = [ "[B- ->  (J/psi(1S)-> ^mu+  mu-)  K-]CC",
 *                        "[B- ->  (J/psi(1S)-> mu+  ^mu-)  K-]CC",
 *                        "[B- ->  (J/psi(1S)-> ^mu+  mu-)  K-]CC"]
 *
 * @endcode
 * creates the branches: Bu_DOCA_Kmu_OS, Bu_DOCA_Kmu_SS, Bu_DOCA_mumu
 * and: Bu_DOCACHI2_Kmu_OS, Bu_DOCACHI2_Kmu_SS, Bu_DOCACHI2_mumu
 *
 *
 * @sa DecayTreeTuple
 *
 *  @author Paul Seyfert
 *  @date   2016-11-08
 */
class TupleToolDOCA : public TupleToolBase, virtual public IParticleTupleTool {
public:
  /// Standard constructor
  TupleToolDOCA( const std::string& type, const std::string& name, const IInterface* parent );

  virtual ~TupleToolDOCA(){}; ///< Destructor

  virtual StatusCode initialize() override;
  virtual StatusCode finalize() override;

  StatusCode fill( const LHCb::Particle*, const LHCb::Particle*, const std::string&, Tuples::Tuple& ) override;

private:
  IDVAlgorithm*              m_dva;
  const IDistanceCalculator* m_dist;

  /// Well, I wanted a std::map<std::string,std::pair<std::string,std::string>>
  /// here, but Gaudi cannot parse that. Now I'm vulnerable against name
  /// collisions in m_name and lenght mismatches between the three vectors

  std::vector<std::string> m_locations1;
  std::vector<std::string> m_locations2;
  std::vector<std::string> m_name;

  std::map<std::string, std::pair<LoKi::Child::Selector*, LoKi::Child::Selector*>> m_childSelectors;
};

#endif // TupleToolDOCA_H

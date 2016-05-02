// Dear emacs, this is -*- c++ -*-
#ifndef GenJetProps_H
#define GenJetProps_H


#include "UHHAnalysis/NtupleWriter/Objects/Jet.h"
#include "UHHAnalysis/NtupleWriter/Objects/TopJet.h"
#include "UHHAnalysis/NtupleWriter/Objects/GenJetWithParts.h"
#include "UHHAnalysis/NtupleWriter/Objects/GenParticle.h"
#include "UHHAnalysis/NtupleWriter/Objects/PFParticle.h"
#include "UHHAnalysis/NtupleWriter/interface/Njettiness.h"


/**
 *  @short for the calculation of jet properties
 *   most of the functions here use FastJet to calculate
 *   observables like N-subjettiness, mass-drop,
 *   sub-jets or jet grooming variables
 *   
 *   If you use this functionality, make sure that the 
 *   jet constituents are provided
 *
 *
 *  @author Roman Kogler, Tobias Lapsien
 */

class GenJetProps
{
 public:
  GenJetProps();
  GenJetProps(GenJetWithParts* jet);
  GenJetProps(GenJetWithParts* jet,  std::vector<GenParticle>* parts);
  ~GenJetProps();

  double GetNsubjettiness(int N, Njettiness::AxesMode mode, double beta, double R0, double Rcutoff=std::numeric_limits<double>::max());

  double GetPrunedNsubjettiness(int N, Njettiness::AxesMode mode, double beta, double R0, double Rcutoff=std::numeric_limits<double>::max());

  std::vector<fastjet::PseudoJet> GetFastJet(double R0, fastjet::JetAlgorithm jet_algo=fastjet::cambridge_algorithm);
  std::vector<fastjet::PseudoJet> GetFastJetcorrected(double R0, fastjet::JetAlgorithm jet_algo=fastjet::cambridge_algorithm);
  fastjet::PseudoJet GetPrunedJet(fastjet::PseudoJet injet);

  double GetQjetVolatility(int seed, double R0);

  double FindMean( std::vector< double > qjetmasses );
  double FindRMS( std::vector< double > qjetmasses );

  std::vector<fastjet::PseudoJet> GetJetConstituents();

  void set_pf_cands( std::vector<GenParticle>* ptr){m_pf_candidates = ptr;}

 private:

  std::vector<GenParticle>* m_pf_candidates;

  GenJetWithParts* m_jet;
  fastjet::ClusterSequence* m_JetFinder;
  fastjet::JetDefinition* m_JetDef ;

};


#endif // GenJetProps_H

#include "UHHAnalysis/NtupleWriter/interface/GenJetProps.h"

#include <vector>
#include <iostream>

#include "UHHAnalysis/NtupleWriter/interface/Nsubjettiness.h"
#include "UHHAnalysis/NtupleWriter/interface/QjetsPlugin.h"
#include "fastjet/tools/Pruner.hh"

#include "SFrameTools/include/EventCalc.h"
#include "fastjet/ClusterSequenceArea.hh"
#include "SFrameTools/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "SFrameTools/JetMETObjects/interface/JetCorrectorParameters.h"

#include "TString.h"
#include "TObject.h"
#include "TMath.h"

using namespace std;

GenJetProps::GenJetProps()
{
  // standard constructor
  m_jet = NULL;
  m_pf_candidates = NULL;
  m_JetFinder = NULL;
  m_JetDef = NULL;
}

GenJetProps::GenJetProps(GenJetWithParts* jet)
{
  // standard constructor
  m_jet = jet;
  m_pf_candidates = NULL;
  m_JetFinder = NULL;
  m_JetDef = NULL;
}

GenJetProps::GenJetProps(GenJetWithParts* jet, std::vector<GenParticle>* parts)
{
  // standard constructor
  m_jet = jet;
  m_pf_candidates = parts;
  m_JetFinder = NULL;
  m_JetDef = NULL;
}

GenJetProps::~GenJetProps()
{
  // destructor, clean up
  if (m_JetFinder) delete m_JetFinder;
  if (m_JetDef) delete m_JetDef;
  m_pf_candidates = NULL;
}

double GenJetProps::GetNsubjettiness(int N, Njettiness::AxesMode mode, double beta, double R0, double Rcutoff)
{
  /// get N-subjettiness
  /// the actual calculation is done by using spartyjet
  /// the second argument is the way the sub-jet axes are found
  /// from NsubjettinessPlugin.hh:
  /// kt_axes              : exclusive kT
  /// ca_axes              : exclusive CA
  /// min_axes             : minimize N-jettiness by iteration (100 passes default)
  /// onepass_{kt,ca}_axes : one-pass minimization seeded by kt or CA (pretty good)

  if (!m_jet){
    return 0.;
  }

  // first: create a fastjet-jet from the TopJet
  std::vector<fastjet::PseudoJet> FJparticles = GetJetConstituents();

  fastjet::PseudoJet fjet = join(FJparticles);
  
  fastjet::Nsubjettiness Nsubj(N, mode, beta, R0, Rcutoff);
  
  double sub=Nsubj.result(fjet);
  
  return sub;

}

double GenJetProps::GetPrunedNsubjettiness(int N, Njettiness::AxesMode mode, double beta, double R0, double Rcutoff)
{
  /// get N-subjettiness, same as above but for pruned jets
  /// pruning ist first called via GetPrunedJet (see below)

  if (!m_jet){
    return 0.;
  }

  std::vector<fastjet::PseudoJet> jets = GetFastJet(2.0);   // something large to make sure jet is inside radius
  if(jets.size()<1) return 0.;
  fastjet::PseudoJet pjet = GetPrunedJet(jets[0]);

  // first: create a fastjet-jet from the TopJet
  std::vector<fastjet::PseudoJet> FJparticles = pjet.constituents();
  fastjet::PseudoJet fjet = join(FJparticles);
  
  fastjet::Nsubjettiness Nsubj(N, mode, beta, R0, Rcutoff);
  
  double sub=Nsubj.result(fjet);
  
  return sub;

}

std::vector<fastjet::PseudoJet> GenJetProps::GetJetConstituents()
{
  // get the jet constituents as Fastjet pseudojets

  std::vector<fastjet::PseudoJet> FJparticles;
  if (m_pf_candidates==NULL){
    cout << "Error, can not get jet constituents without list of PF candidates"<< endl;
    return FJparticles;
  }

  std::vector<unsigned int> jet_parts_ind = m_jet->genparticles_indices();

  for(unsigned int ic=0; ic<jet_parts_ind.size(); ++ic){
    GenParticle p = m_pf_candidates->at(jet_parts_ind[ic]);
    fastjet::PseudoJet fj( p.pt()*cos(p.phi()),
			   p.pt()*sin(p.phi()),
			   p.pt()*sinh(p.eta()),
			   p.energy() );
    
    FJparticles.push_back( fj );   
  }
  return FJparticles;
}

vector<fastjet::PseudoJet> GenJetProps::GetFastJet(double R0, fastjet::JetAlgorithm jet_algo)
{
  /// get fastjet jet from the input particles of the object's jet
  /// could be more than one jet, if the jet definition is different 
  /// than the one being used for the initial jet or the radius is smaller

  std::vector<fastjet::PseudoJet> InputParticles = GetJetConstituents();

  if (m_JetDef){
    delete m_JetDef;
  }
  m_JetDef = new fastjet::JetDefinition(jet_algo, R0);  
  if (m_JetFinder){
    delete m_JetFinder;
  }
  m_JetFinder = new fastjet::ClusterSequence(InputParticles, *m_JetDef);

  vector<fastjet::PseudoJet> Jets = m_JetFinder->inclusive_jets(10.);

  // return if no jets were found
  if (Jets.size() == 0){
    return Jets; 
  }

  // sort jets by increasing pt
  vector<fastjet::PseudoJet> SortedJets = sorted_by_pt(Jets);

  return SortedJets;
  
}

vector<fastjet::PseudoJet> GenJetProps::GetFastJetcorrected(double R0, fastjet::JetAlgorithm jet_algo)
{
  /// get fastjet jet from the input particles of the object's jet
  /// could be more than one jet, if the jet definition is different 
  /// than the one being used for the initial jet or the radius is smaller

  std::vector<fastjet::PseudoJet> InputParticles = GetJetConstituents();

  if (m_JetDef){
    delete m_JetDef;
  }
  m_JetDef = new fastjet::JetDefinition(jet_algo, R0);  
  if (m_JetFinder){
    delete m_JetFinder;
  }
  m_JetFinder = new fastjet::ClusterSequence(InputParticles, *m_JetDef);


  vector<fastjet::PseudoJet> Jets = m_JetFinder->inclusive_jets(10.);
  double maxrap = 5.0;
   unsigned int n_repeat = 3; // default is 1
   double ghost_area = 0.01; // this is the default
   fastjet::GhostedAreaSpec area_spec(maxrap, n_repeat, ghost_area);
 
   fastjet::AreaDefinition area_def(fastjet::active_area, area_spec); 

   fastjet::ClusterSequenceArea clust_seq(InputParticles, *m_JetDef, area_def);
  // return if no jets were found
  if (Jets.size() == 0){
    return Jets; 
  }

  // sort jets by increasing pt
 
  double ptmin=20;
  vector<fastjet::PseudoJet> jets = sorted_by_pt(clust_seq.inclusive_jets(ptmin));
  // for(unsigned int i=0;i<jets.size();i++) std:cout<<jets[i].area()<<std::endl;
  vector<fastjet::PseudoJet> SortedJets = sorted_by_pt(jets);

  EventCalc* calc = EventCalc::Instance();
  BaseCycleContainer* bcc = calc->GetBaseCycleContainer();
  std::vector<JetCorrectorParameters> pars;
   pars.push_back(JetCorrectorParameters("/nfs/dust/cms/user/peiffer/CMSSW_5_3_10/src/START53_V27_L1FastJet_AK5PFchs.txt"));
   pars.push_back(JetCorrectorParameters("/nfs/dust/cms/user/peiffer/CMSSW_5_3_10/src/START53_V27_L2L3Residual_AK5PFchs.txt"));
   pars.push_back(JetCorrectorParameters("/nfs/dust/cms/user/peiffer/CMSSW_5_3_10/src/START53_V27_L3Absolute_AK5PFchs.txt"));
    FactorizedJetCorrector* corrector = new FactorizedJetCorrector(pars);

   for(unsigned int j=0;j<SortedJets.size();++j){
     fastjet::PseudoJet pseudojet=SortedJets[j];
     corrector->setJetPt(pseudojet.pt());
     corrector->setJetEta(pseudojet.eta());
     corrector->setJetE(pseudojet.E());
     corrector->setJetA(pseudojet.area());
     corrector->setRho(bcc->rho);
     corrector->setNPV(bcc->pvs->size());
     float correctionfactor = corrector->getCorrection();
     //  std::cout<<correctionfactor<<std::endl;
     //  std::cout<<SortedJets[j].pt()<<" "<<SortedJets[j].eta()<<" "<<SortedJets[j].pt()<<std::endl;
     SortedJets[j].operator*=(correctionfactor);
     //   std::cout<<"After correction: "<<SortedJets[j].pt()<<" "<<SortedJets[j].eta()<<" "<<SortedJets[j].pt()<<std::endl;
   }
  return SortedJets;
  
}


fastjet::PseudoJet GenJetProps::GetPrunedJet(fastjet::PseudoJet injet)
{
  // calculate the pruned jet

  double zcut = 0.1;   
  double rcut_factor = 0.5;   
  fastjet::Pruner pruner(fastjet::cambridge_algorithm, zcut, rcut_factor);   
  fastjet::PseudoJet pruned_jet = pruner(injet);

  return pruned_jet;
}



double GenJetProps::GetQjetVolatility(int seed, double R0)
{
    
    fastjet::JetDefinition jetDef(fastjet::cambridge_algorithm, R0);    
    int mQJetsN = 50;
    int mQJetsPreclustering = 30;
    
    std::vector< double > qjetmasses;
    
    std::vector<fastjet::PseudoJet> FJparticles = GetJetConstituents();

    fastjet::ClusterSequence thisClustering_basic(FJparticles, jetDef);
    std::vector<fastjet::PseudoJet> out_jets_basic = sorted_by_pt(thisClustering_basic.inclusive_jets(50.0));
    if (out_jets_basic.size()==0){ 
      return -1.; // no jets found, something wrong
    }

    int j = 0; // the hardest jet
    
    double zcut(0.1), dcut_fctr(0.5), exp_min(0.), exp_max(0.), rigidity(0.1);          
    
    vector<fastjet::PseudoJet> constits;
    unsigned int nqjetconstits = out_jets_basic.at(j).constituents().size();
    if (nqjetconstits < (unsigned int) mQJetsPreclustering) constits = out_jets_basic.at(j).constituents();
    else constits = out_jets_basic.at(j).associated_cluster_sequence()->exclusive_subjets_up_to(out_jets_basic.at(j),mQJetsPreclustering);
    for(unsigned int ii = 0 ; ii < (unsigned int) mQJetsN ; ii++){
        QjetsPlugin qjet_plugin(zcut, dcut_fctr, exp_min, exp_max, rigidity);
	qjet_plugin.SetRandSeed(seed+ii); // new feature in Qjets to set the random seed
        fastjet::JetDefinition qjet_def(&qjet_plugin);
        fastjet::ClusterSequence qjet_seq(constits, qjet_def);
        vector<fastjet::PseudoJet> inclusive_jets2 = sorted_by_pt(qjet_seq.inclusive_jets(50.0));
        
        if (inclusive_jets2.size()>0) { qjetmasses.push_back( inclusive_jets2[0].m() ); }
        
    }
    
        // find RMS of a vector
    double qjetsRMS = FindRMS( qjetmasses );
        // find mean of a vector
    double qjetsMean = FindMean( qjetmasses );
        // compute QJets volatility
    double qjetsVolatility = qjetsRMS/qjetsMean;
    return qjetsVolatility;
}
    // -------------------------------------------
    // -------------------------------------------



double GenJetProps::FindRMS( std::vector< double > qjetmasses )
{
  
  double total = 0.;
  double ctr = 0.;
  for (unsigned int i = 0; i < qjetmasses.size(); i++){
    total = total + qjetmasses[i];
    ctr++;
  }
  double mean =  total/ctr;
  
  double totalsquared = 0.;
  for (unsigned int i = 0; i < qjetmasses.size(); i++){
    totalsquared += (qjetmasses[i] - mean)*(qjetmasses[i] - mean) ;
  }
  double RMS = sqrt( totalsquared/ctr );
  return RMS;
}

double GenJetProps::FindMean( std::vector< double > qjetmasses )
{
  double total = 0.;
  double ctr = 0.;
  for (unsigned int i = 0; i < qjetmasses.size(); i++){
    total = total + qjetmasses[i];
    ctr++;
  }
  return total/ctr;
}

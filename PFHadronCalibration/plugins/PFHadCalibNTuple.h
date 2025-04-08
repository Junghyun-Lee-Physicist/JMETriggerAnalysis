#ifndef PFHadCalibNTuple_H
#define PFHadCalibNTuple_H

#include <string>
#include <vector>

#include "FWCore/Framework/interface/Frameworkfwd.h" // CMSSW Framework "forward" declaration
#include "FWCore/Framework/interface/one/EDAnalyzer.h" // edm::one::EDAnalyzer
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h" // reco::GenParticleCollection
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h" // reco::PFCandidateCollection
#include "DataFormats/ParticleFlowReco/interface/PFSimParticleFwd.h"// reco::PFSimParticleCollection

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include <TFile.h>
#include <TTree.h>

class PFHadCalibNTuple : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit PFHadCalibNTuple(const edm::ParameterSet&);
  ~PFHadCalibNTuple() override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void analyze(const edm::Event&, const edm::EventSetup&) override;

  // Module name
  std::string moduleLabel_;

  // TTree name and title
  std::string ttreeName_;
  std::string ttreeTitle_;

  // Flag to print detailed log
  bool printDetailLog_;

  // Input EDM collections
  edm::EDGetTokenT<reco::GenParticleCollection> genPartsToken_;
  edm::EDGetTokenT<reco::PFSimParticleCollection> pfSimPartsToken_;
  edm::EDGetTokenT<reco::PFCandidateCollection> recoPFCandsToken_;

  // Configured parameters
      // Status code of selected GEN particles (e.g., final-state particles for PFHC calibration)
      // Typical values: 1 (stable particle), 2 (unstable decay product), etc.
      int genParticleStatus_;

      // PDG ID of the selected GEN particles (e.g., quarks, leptons, hadrons)
      // PDG code uniquely identifies the particle type. 
      //   --> Examples: 11 (electron), 13 (muon), 211 (pion), 22 (photon)
      int genParticlePdgId_;

      // Minimum delta-R between isolated GEN particles and other GEN particles
      // delta-R represents the angular distance in eta-phi space; 
      // smaller delta-R indicates stronger interaction between particles.
      double genParticleIsoMinDeltaR_;

      // Parameters for charged hadrons
      double minPt_;
      double minTrackP_;  
      double minTrackPt_;
      double minCaloEnergy_;
      double maxECalEnergy_;
      std::vector<uint> minPixelHits_;
      std::vector<uint> minTrackerHits_;
      std::vector<double> maxEtaForMinTrkHitsCuts_;

      bool usePFBlockElements_;

  // Event counters
  std::vector<uint> globalCounter_;

  // File and TTree pointers
  TFile* file_;
  TTree* ttree_;

  // Ntuple variables
      // For charged hadrons
      float pfsim_true_energy_, pfsim_eta_, pfsim_phi_;
      float pfcan_trackRef_p_, pfcan_Eecal_, pfcan_Ehcal_;
      float pfcan_selected_trackRef_p_, pfcan_selected_Eecal_, pfcan_selected_Ehcal_;
      int charge_; 
      float ho_; //  HO energy
      float Ccorrecal_, Ccorrhcal_;  // Corrected ECAL/HCAL energy
      std::vector<float> correcal_, corrhcal_;  // Vectors for corrected energies

      // For neutral particles
      std::vector<int> pfcID_with_neutrals_;
      std::vector<float> dr_with_neutrals_;
      std::vector<float> Eecal_with_neutrals_;
      std::vector<float> Ehcal_with_neutrals_;
      std::vector<float> correcal_with_neutrals_; 
      std::vector<float> corrhcal_with_neutrals_; 
    
      // For event information
      edm::RunNumber_t run;
      edm::EventNumber_t evt;
      edm::LuminosityBlockNumber_t lumiBlock;
      edm::Timestamp time;
      size_t orun, oevt, olumiBlock, otime;

  void reset_variables();

};

#endif


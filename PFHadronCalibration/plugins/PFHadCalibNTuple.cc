
#include "JMETriggerAnalysis/PFHadronCalibration/plugins/PFHadCalibNTuple.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowReco/interface/PFSimParticle.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlock.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlockElementTrack.h"

#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

using namespace std;

PFHadCalibNTuple::PFHadCalibNTuple(const edm::ParameterSet& iConfig)
    : moduleLabel_(iConfig.getParameter<std::string>("@module_label")),
      ttreeName_(iConfig.getParameter<std::string>("TTreeName")),
      ttreeTitle_(iConfig.getParameter<std::string>("TTreeTitle")),
      printDetailLog_(iConfig.getParameter<bool>("printDetailLog")),
      genPartsToken_(consumes<reco::GenParticleCollection>(
          iConfig.getParameter<edm::InputTag>("genParticles"))),
      pfSimPartsToken_(consumes<reco::PFSimParticleCollection>(
          iConfig.getParameter<edm::InputTag>("pfSimParticles"))),
      recoPFCandsToken_(consumes<reco::PFCandidateCollection>(
          iConfig.getParameter<edm::InputTag>("recoPFCandidates"))),
      genParticleStatus_(iConfig.getParameter<int>("genParticleStatus")),
      genParticlePdgId_(iConfig.getParameter<int>("genParticlePdgId")),
      genParticleIsoMinDeltaR_(iConfig.getParameter<double>("genParticleIsoMinDeltaR")),
      minPt_(iConfig.getParameter<double>("minPt")),
      minTrackP_(iConfig.getParameter<double>("minTrackP")),
      minTrackPt_(iConfig.getParameter<double>("minTrackPt")),
      minCaloEnergy_(iConfig.getParameter<double>("minCaloEnergy")),
      maxECalEnergy_(iConfig.getParameter<double>("maxECalEnergy")),
      minPixelHits_(iConfig.getParameter<std::vector<uint>>("minPixelHits")),
      minTrackerHits_(iConfig.getParameter<std::vector<uint>>("minTrackerHits")),
      maxEtaForMinTrkHitsCuts_(
          iConfig.getParameter<std::vector<double>>("maxEtaForMinTrkHitsCuts")),
      usePFBlockElements_(iConfig.getParameter<bool>("usePFBlockElements")) {
  assert(minPixelHits_.size() == minTrackerHits_.size());
  assert(minPixelHits_.size() == maxEtaForMinTrkHitsCuts_.size());
  for (uint idx = 0; 1 + idx < maxEtaForMinTrkHitsCuts_.size(); ++idx) {
    assert(maxEtaForMinTrkHitsCuts_[idx] < maxEtaForMinTrkHitsCuts_[idx + 1]);
  }

  globalCounter_ = std::vector<uint>(17, 0);

  std::string outputfileName = iConfig.getParameter<std::string>("rootOutputFile");
  file_ = new TFile(outputfileName.c_str(), "RECREATE");
  ttree_ = new TTree(ttreeName_.c_str(), ttreeTitle_.c_str());

  ttree_->Branch("true", &pfsim_true_energy_, "true/F");
  ttree_->Branch("p", &pfcan_selected_trackRef_p_, "p/F");
  ttree_->Branch("eta", &pfsim_eta_, "eta/F");
  ttree_->Branch("phi", &pfsim_phi_, "phi/F");

  ttree_->Branch("ecal", &pfcan_selected_Eecal_, "ecal/F");
  ttree_->Branch("hcal", &pfcan_selected_Ehcal_, "hcal/F");
  ttree_->Branch("ho", &ho_, "ho/F");  

  ttree_->Branch("Ccorrecal", &Ccorrecal_, "Ccorrecal/F");  
  ttree_->Branch("Ccorrhcal", &Ccorrhcal_, "Ccorrhcal/F");  

  ttree_->Branch("charge", &charge_, "charge/I"); 

  ttree_->Branch("dr", &dr_with_neutrals_);
  ttree_->Branch("Eecal", &Eecal_with_neutrals_);
  ttree_->Branch("Ehcal", &Ehcal_with_neutrals_);
  ttree_->Branch("pfcID", &pfcID_with_neutrals_);

  ttree_->Branch("correcal", &correcal_with_neutrals_);  
  ttree_->Branch("corrhcal", &corrhcal_with_neutrals_);  
}

PFHadCalibNTuple::~PFHadCalibNTuple() {

  // Print results and save file
  cout << "-----------------------------------------------------------------------------------------------------" << endl;
  cout << "[ The result of making PFHC Ntuple ]" << endl;
  cout << moduleLabel_ << endl;
  cout << "-----------------------------------------------------------------------------------------------------" << endl;
  cout << "< Number of Events cut-flow >" << endl;
  cout << "Total number of events: " << globalCounter_[0] << endl;
  cout << "Number of events with 1 PFSimParticle: " << globalCounter_[1] << endl;
  cout << "Number of events with No daughter particles of first PFsimParticle( = Pion minus without decay ) and also has negative charge:  " << globalCounter_[2] << endl;
  cout << "  --> Case a) Number of event identified Neutral Particles by PF Algorithm: " << globalCounter_[3] << endl;
  cout << "                Number of events with PFSimParticle's eta > 1E-10: " << globalCounter_[4] << endl;
  cout << "  --> Case b) Number of event identified Charged Particles by PF Algorithm: " << globalCounter_[7] << endl;
  cout << "-----------------------------------------------------------------------------------------------------" << endl;
  cout << "< PF Candidates cut-flow >" << endl;
  cout << "    --> Case a) Identified Neutral Particles" << endl;
  cout << "          Total number of [case a] PF Candidates: " << globalCounter_[5] << endl;
  cout << "            - with dR between ECAL entrance Sim value & PF Candidates: " << globalCounter_[6] << endl;
  cout << endl;
  cout << "    --> Case b) Identified Charged Particles" << endl;
  cout << "          Total number of [case b] PF Candidates: " << globalCounter_[8] << endl;
  cout << "            - identify charged hadron: " << globalCounter_[9] << endl;
  cout << "            - with minimum pT [GeV] : " << globalCounter_[10] << endl;
  cout << "            - with mininum calo energy : " << globalCounter_[11] << endl;
  if (usePFBlockElements_)
    cout << "            - with only 1 track in the block: " << globalCounter_[12] << endl;
  cout << "            - with track-p > " << minTrackP_ << " GeV and track-pT > " << minTrackPt_
       << " GeV: " << globalCounter_[13] << endl;
  cout << "            - with min number of " << minPixelHits_[0] << " of pixel hits: " << globalCounter_[14] << endl;
  cout << "            - with min number of " << minTrackerHits_[0] << " pixel+strip hits( = tracker hit ): "
       << globalCounter_[15] << endl;
  cout << "Above 2 lines about [ pixel hits & track hits ], the values of minimum hits vary eta range" << endl;
  cout << "Please see the configuration file to check certain values" << endl;
  cout << "            - with E_ECAL < " << maxECalEnergy_ << " GeV: " << globalCounter_[16] << endl;
  cout << "-----------------------------------------------------------------------------------------------------" << endl;

  file_->cd();
  ttree_->Write();

  file_->Write();
  file_->Close();

}

void PFHadCalibNTuple::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  if (printDetailLog_) {
    cout << "----------------------------------------------------------" << endl;
    cout << "  [Event begin]  " << endl;
    cout << endl;
  }

  //------------------------------------------------------------------------
  // Declare collections
  //   [ edm::Hadnole --> Pointer ]
  edm::Handle<reco::GenParticleCollection> genPartsHandle;
  edm::Handle<reco::PFSimParticleCollection> pfSimPartsHandle;
  edm::Handle<reco::PFCandidateCollection> recoPFCandsHandle;

  iEvent.getByToken(genPartsToken_, genPartsHandle);
  iEvent.getByToken(pfSimPartsToken_, pfSimPartsHandle);
  iEvent.getByToken(recoPFCandsToken_, recoPFCandsHandle);

  const reco::GenParticleCollection& genParts = *genPartsHandle;
  const reco::PFSimParticleCollection& pfSimParts = *pfSimPartsHandle;
  const reco::PFCandidateCollection& recoPFCands = *recoPFCandsHandle;

  //------------------------------------------------------------------------

  // Get event information
  run = iEvent.id().run();
  evt = iEvent.id().event();
  lumiBlock = iEvent.id().luminosityBlock();
  time = iEvent.time();

  orun = (size_t)run;
  oevt = (size_t)evt;
  olumiBlock = (size_t)lumiBlock;
  otime = (size_t)((iEvent.time().value()) >> 32);

  if (printDetailLog_) {
    cout << "[analyze] Event info.." << endl;
    cout << "[analyze] Current run         : [" << orun << "]" << endl;
    cout << "[analyze] Current evt number  : [" << oevt << "]" << endl;
    cout << "[analyze] Current lumiBlock   : [" << olumiBlock << "]" << endl;
    cout << "[analyze] Current time        : [" << otime << "]" << endl;
    cout << endl;
  }

  reset_variables();

  ++globalCounter_[0];  // Total number of event

  if (pfSimParts.size() != 1) return;
  ++globalCounter_[1];  // Number of events with 1 PFsimParticle

  if (pfSimParts[0].daughterIds().size() != 0) return;
  if (pfSimParts[0].charge() >= 0) return;
  ++globalCounter_[2];  // Number of events with No daughter particles and negative charge

  bool is_current_PFCandidate_charged_particle = false;
  for (auto const& pfc : recoPFCands) {
    if (pfc.particleId() < 4) {
      // pfc id = 0 -> unknown, or dummy
      // pfc id = 1 -> charged hadron
      // pfc id = 2 -> electron
      // pfc id = 3 -> muon

      is_current_PFCandidate_charged_particle = true;
      break;
    }
  }

  //------------------------------------------------------------------------
  // Identified Neutral Particles by PF Algorithm
  if (!is_current_PFCandidate_charged_particle) {
    ++globalCounter_[3];  // case a) Number of event identified Neutral Particles by PF Algorithm

    typedef reco::PFTrajectoryPoint pfPathPoint;

    pfPathPoint::LayerType ecalEntrance = pfPathPoint::ECALEntrance;
    const pfPathPoint& tp_at_ecal = (pfSimParts[0]).extrapolatedPoint(ecalEntrance);

    pfsim_eta_ = tp_at_ecal.positionREP().Eta();
    pfsim_phi_ = tp_at_ecal.positionREP().Phi();
    pfsim_true_energy_ = std::sqrt(tp_at_ecal.momentum().Vect().Mag2());

    if (fabs(pfsim_eta_) < 1E-10) return;
    ++globalCounter_[4];  // Number of events with PFSimParticle's eta > 1E-10

    // Iteration over PF candidates
    for (auto const& pfc : recoPFCands) {
      ++globalCounter_[5];  // Total number of [Case a] PF Candidates

      double dR = reco::deltaR(pfsim_eta_, pfsim_phi_, pfc.eta(), pfc.phi());
      if (dR < 1.2) {
        ++globalCounter_[6];

        dr_with_neutrals_.push_back(dR);
        pfcID_with_neutrals_.push_back(pfc.particleId());
        Eecal_with_neutrals_.push_back(pfc.rawEcalEnergy());
        Ehcal_with_neutrals_.push_back(pfc.rawHcalEnergy());
        correcal_with_neutrals_.push_back(pfc.ecalEnergy());  
        corrhcal_with_neutrals_.push_back(pfc.hcalEnergy());  
      }

      // Modified to match energy calculation in offline code
      if (pfc.particleId() == 4 && dR < 0.2) pfcan_selected_Eecal_ += pfc.rawEcalEnergy();
      if (pfc.particleId() == 5 && dR < 0.4) pfcan_selected_Ehcal_ += pfc.rawHcalEnergy();
    }

    ttree_->Fill();
    return;
  }

  //------------------------------------------------------------------------
  // Identified Charged Particles by PF Algorithm
  ++globalCounter_[7];  // Number of event identified Charged Particles by PF Algorithm

  for (auto const& pfc : recoPFCands) {
    ++globalCounter_[8];  // Total number of [Case b] PFCandidates

    // Demand charged hadrons only
    if (pfc.particleId() != 1) continue;
    ++globalCounter_[9];  // - identify charged hadron

    // Charged hadron minimum pt
    if (pfc.pt() < minPt_) continue;
    ++globalCounter_[10];  // - with minimum pT

    pfcan_Eecal_ = pfc.rawEcalEnergy();
    pfcan_Ehcal_ = pfc.rawHcalEnergy();
    ho_ = pfc.rawHoEnergy();           
    Ccorrecal_ = pfc.ecalEnergy();     
    Ccorrhcal_ = pfc.hcalEnergy();     
    charge_ = pfc.charge();            

    if ((pfcan_Eecal_ + pfcan_Ehcal_) < minCaloEnergy_) continue;
    ++globalCounter_[11];  // - with mininum calo energy

    auto nTracks = 0u;
    auto const& theElements = pfc.elementsInBlocks();
    if (theElements.empty()) {
      if (not usePFBlockElements_)
        nTracks = 1;
    } else {
      auto const& elements = theElements[0].first->elements();
      for (unsigned iEle = 0; iEle < elements.size(); ++iEle) {
        if (elements[iEle].type() == reco::PFBlockElement::TRACK) {
          ++nTracks;
        }
      }
    }

    // Only one track in the block
    if (nTracks != 1) continue;
    ++globalCounter_[12];  // - with only 1 track in the block:

    auto trackRef = pfc.trackRef();
    auto const& hp = trackRef->hitPattern();
    uint const track_nValidPixelHits = hp.numberOfValidPixelHits();
    uint const track_nValidTrackerHits = trackRef->numberOfValidHits();

    pfcan_trackRef_p_ = trackRef->p();
    auto const pfcan_trackRef_pt = trackRef->pt();
    auto const pfcan_trackRef_eta = trackRef->eta();

    if (pfcan_trackRef_p_ < minTrackP_ || pfcan_trackRef_pt < minPt_) continue;
    ++globalCounter_[13];

    auto hasMinPixelHits = false;
    auto hasMinTrackerHits = false;
    for (uint ieta = 0; ieta < maxEtaForMinTrkHitsCuts_.size(); ++ieta) {
      auto const etaMin = ieta ? maxEtaForMinTrkHitsCuts_[ieta - 1] : 0.;
      auto const etaMax = maxEtaForMinTrkHitsCuts_[ieta];

      if (std::abs(pfcan_trackRef_eta) >= etaMin and std::abs(pfcan_trackRef_eta) < etaMax) {
        hasMinPixelHits = track_nValidPixelHits >= minPixelHits_[ieta];
        hasMinTrackerHits = track_nValidTrackerHits >= minTrackerHits_[ieta];
        break;
      }
    }

    if (not hasMinPixelHits) continue;
    ++globalCounter_[14];

    if (not hasMinTrackerHits) continue;
    ++globalCounter_[15];

    if (pfcan_Eecal_ > maxECalEnergy_) continue;
    ++globalCounter_[16];

    pfcan_selected_trackRef_p_ = pfcan_trackRef_p_;
    pfcan_selected_Eecal_ = pfcan_Eecal_;
    pfcan_selected_Ehcal_ = pfcan_Ehcal_;

    // Initialize neutral particle information
    dr_with_neutrals_.clear();
    pfcID_with_neutrals_.clear();
    Eecal_with_neutrals_.clear();
    Ehcal_with_neutrals_.clear();
    correcal_with_neutrals_.clear();  
    corrhcal_with_neutrals_.clear();  

    // Collect neutral particle information
    for (auto const& pfc_neutral : recoPFCands) {
      if (pfc_neutral.charge() != 0) continue;  // Select only neutral particles

      double dR = reco::deltaR(pfsim_eta_, pfsim_phi_, pfc_neutral.eta(), pfc_neutral.phi());
      if (dR < 1.2) {
        dr_with_neutrals_.push_back(dR);
        pfcID_with_neutrals_.push_back(pfc_neutral.particleId());
        Eecal_with_neutrals_.push_back(pfc_neutral.rawEcalEnergy());
        Ehcal_with_neutrals_.push_back(pfc_neutral.rawHcalEnergy());
        correcal_with_neutrals_.push_back(pfc_neutral.ecalEnergy());  
        corrhcal_with_neutrals_.push_back(pfc_neutral.hcalEnergy());  
      }
    }
  }

  typedef reco::PFTrajectoryPoint pfPathPoint;

  pfPathPoint::LayerType ecalEntrance = pfPathPoint::ECALEntrance;
  const pfPathPoint& tp_at_ecal = (pfSimParts[0]).extrapolatedPoint(ecalEntrance);

  pfsim_eta_ = tp_at_ecal.positionREP().Eta();
  pfsim_phi_ = tp_at_ecal.positionREP().Phi();
  pfsim_true_energy_ = std::sqrt(tp_at_ecal.momentum().Vect().Mag2());

  ttree_->Fill();

  if (printDetailLog_) {
    cout << "  [Event end]  " << endl;
    cout << "----------------------------------------------------------" << endl;
  }
}

void PFHadCalibNTuple::reset_variables() {
  pfsim_true_energy_ = -1234.0;
  pfsim_eta_ = -1234.0;
  pfsim_phi_ = -1234.0;

  pfcan_trackRef_p_ = 0.00;
  pfcan_Eecal_ = 0.00;
  pfcan_Ehcal_ = 0.00;

  pfcan_selected_trackRef_p_ = 0.00;
  pfcan_selected_Eecal_ = 0.00;
  pfcan_selected_Ehcal_ = 0.00;

  ho_ = 0.00;        
  Ccorrecal_ = 0.00; 
  Ccorrhcal_ = 0.00; 
  charge_ = 0;        

  pfcID_with_neutrals_.clear();
  dr_with_neutrals_.clear();
  Eecal_with_neutrals_.clear();
  Ehcal_with_neutrals_.clear();

  correcal_with_neutrals_.clear();  
  corrhcal_with_neutrals_.clear();  
}

void PFHadCalibNTuple::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(PFHadCalibNTuple);

// class to produce 2 pat::MuonCollections

// one matched to the selected triggers
// another fitered wrt the selected triggers

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Utilities/interface/StreamID.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "DataFormats/PatCandidates/interface/Muon.h"

#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "DataFormats/PatCandidates/interface/TriggerPath.h"
#include "DataFormats/PatCandidates/interface/TriggerEvent.h"
#include "DataFormats/PatCandidates/interface/TriggerAlgorithm.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"

#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

#include <TLorentzVector.h>
#include "helper.h"
#include "TVectorD.h"    // for fixing tracks
#include "TMatrixDSym.h" // for fixing tracks

using namespace std;

constexpr bool debug = false;

class TriMuonTriggerSelector : public edm::EDProducer {
  
public:
    
  explicit TriMuonTriggerSelector(const edm::ParameterSet &iConfig);
    
  ~TriMuonTriggerSelector() override {};
  
private:
  
  virtual void produce(edm::Event&, const edm::EventSetup&);

  reco::Track fix_track(const reco::Track *tk, double delta) const;  

  const edm::ESGetToken<MagneticField, IdealMagneticFieldRecord> bFieldToken_;  
  edm::EDGetTokenT<std::vector<pat::Muon>> muonSrc_;
  const edm::EDGetTokenT<reco::BeamSpot> beamSpotSrc_;
  const edm::EDGetTokenT<std::vector<reco::Vertex>> primaryVtxSrc_;
  edm::EDGetTokenT<edm::TriggerResults> triggerBits_;
  edm::EDGetTokenT<std::vector<pat::TriggerObjectStandAlone>> triggerObjects_;
  
  
  // for trigger match
  std::vector<std::string> HLTPaths_;
  const double drForTriggerMatch_;
  
  // Offline muons selection  
  const double ptMin_;           
  const double absEtaMax_;       
};


TriMuonTriggerSelector::TriMuonTriggerSelector(const edm::ParameterSet &iConfig):
  bFieldToken_(esConsumes<MagneticField, IdealMagneticFieldRecord>()),
  //
  muonSrc_( consumes<std::vector<pat::Muon>> ( iConfig.getParameter<edm::InputTag>( "muonCollection" ) ) ),
  //
  beamSpotSrc_(consumes<reco::BeamSpot>( iConfig.getParameter<edm::InputTag>( "beamSpot" ) ) ),
  primaryVtxSrc_(consumes<std::vector<reco::Vertex>>( iConfig.getParameter<edm::InputTag>( "primaryVtx" ) ) ),
  triggerBits_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("bits"))),
  triggerObjects_(consumes<std::vector<pat::TriggerObjectStandAlone>>(iConfig.getParameter<edm::InputTag>("objects"))),
  HLTPaths_(iConfig.getParameter<std::vector<std::string>>("HLTPaths")),
  drForTriggerMatch_(iConfig.getParameter<double>("drForTriggerMatch")),    
  //
  ptMin_(iConfig.getParameter<double>("ptMin")),
  absEtaMax_(iConfig.getParameter<double>("absEtaMax"))
{
  // produce the SelectedMuons collection (all muons passing the preselection)
  produces<pat::MuonCollection>("SelectedMuons");
  produces<TransientTrackCollection>("SelectedTransientMuons");  
  // produce the triggering muons collection
  produces<pat::MuonCollection>("trgMuons");
}

void TriMuonTriggerSelector::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  // Inputs
  const auto& bField = iSetup.getData(bFieldToken_);

  edm::Handle<reco::BeamSpot> beamSpotHandle;
  iEvent.getByToken(beamSpotSrc_, beamSpotHandle);
  const reco::BeamSpot& beamSpot = *beamSpotHandle;

  edm::Handle<std::vector<reco::Vertex>> primaryVtxHandle;
  iEvent.getByToken(primaryVtxSrc_, primaryVtxHandle);
  const reco::Vertex &PV = primaryVtxHandle->front();

  edm::Handle<edm::TriggerResults> triggerBits;
  iEvent.getByToken(triggerBits_, triggerBits);
  const edm::TriggerNames &trigNames = iEvent.triggerNames(*triggerBits); 

  edm::Handle<std::vector<pat::TriggerObjectStandAlone>> triggerObjects;
  iEvent.getByToken(triggerObjects_, triggerObjects);

  edm::Handle<std::vector<pat::Muon>> muons;
  iEvent.getByToken(muonSrc_, muons);

  // Outputs
  std::unique_ptr<pat::MuonCollection>      muons_out      ( new pat::MuonCollection );
  std::unique_ptr<TransientTrackCollection> trans_muons_out( new TransientTrackCollection );
  std::unique_ptr<pat::MuonCollection>      trgmuons_out   ( new pat::MuonCollection );
  
  std::vector<int> muonIsTrigger(muons->size(), 0);
  std::vector<std::vector<int>> fires;
  std::vector<std::vector<float>> matcher; 
  std::vector<std::vector<float>> DR;
  

  // Trigger debug
  if(debug) {
    const edm::TriggerNames &names = iEvent.triggerNames(*triggerBits);
    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << "-----------------------------------------------" << std::endl;
    std::cout << std::endl;
    std::cout << "\n == TRIGGER PATHS= " << std::endl;
    for (unsigned int i = 0, n = triggerBits->size(); i < n; ++i) {
      if (triggerBits->accept(i)) 
	std::cout << "Event = " << (iEvent.id()).event() << ", Trigger " << names.triggerName(i) 
		  << ": Pass = " << (triggerBits->accept(i)) 
		  << ", Was Run = " << (triggerBits->wasrun(i))
		  << std::endl;
    }
  }

  // Loop over reconstructed muons
  for(const pat::Muon &muon : *muons){
    
    if(debug) std::cout << "Slimmed muons: muon Pt = " << muon.pt() 
			<< " Eta = " << muon.eta() << " Phi = " << muon.phi()  <<endl;
    
    // These vectors have one entry per HLT path
    std::vector<int> frs(HLTPaths_.size(),0);              
    std::vector<float> temp_matched_to(HLTPaths_.size(),1000.);
    std::vector<float> temp_DR(HLTPaths_.size(),1000.);

    // Loop over trigger paths
    int ipath=-1;
    for (const std::string path: HLTPaths_){
      
      if(debug) std::cout << "ipath = " << ipath << ", path = " << path << std::endl;
      if(debug) std::cout << std::endl;
      ipath++;
      
      // Here we loop over trigger objects
      float minDr = 1000.;
      float minPt = 1000.;

      if(debug) std::cout << std::endl;
      if(debug) std::cout << "Now start loop over trigger objects" << std::endl;

      // Loop over trigger objects matching the reference path
      for (pat::TriggerObjectStandAlone obj : *triggerObjects) {
      
	if(debug) std::cout << "New object" << std::endl;

	// consider only objects which match the ref path    
	obj.unpackPathNames(trigNames);
	obj.unpackFilterLabels(iEvent, *triggerBits);
	std::vector<std::string> pathNamesAll  = obj.pathNames(false);
	bool isPathExist = false;
	for (unsigned h = 0, n = pathNamesAll.size(); h < n; ++h) {

	  // remove the path version
	  string pathNameStart;
	  pathNameStart = pathNamesAll[h].substr(0,pathNamesAll[h].find("_v")-0);
	  if(debug) std::cout << "In loop over trigger object: this is ipath = " << pathNamesAll[h] << ", I need path = " << path << std::endl;	  
	  if(pathNameStart==path) isPathExist = true;     
	}
	if(!isPathExist) continue;
      
	if(debug) std::cout << "Path found" << std::endl;	  

	// muObjNumberGM: hlt candidate firing the dimuon part of the trigger, requiring global muons
	// muObjNumberTM: hlt candidate firing the request for tracker muons
	int muObjNumberGM = -1;
	int muObjNumberTM = -1;
	int muObjNumberUniv = -1;
	for (unsigned hh = 0; hh < obj.filterLabels().size(); ++hh){	

	  if(debug) std::cout << "Event: Filter " << hh << " => " << obj.filterLabels()[hh] << ": pt =  " << obj.pt() << ", eta = " << obj.eta() << ", phi = " << obj.phi() << std::endl;
	  if(debug) std::cout << "" << std::endl;	  
	  
	  if (path=="HLT_Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15_Charge1" || path=="HLT_Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15") {
	    if(obj.filterLabels()[hh].find("hltTau3MuTriMuon1filter") != std::string::npos) {  
	      muObjNumberTM = hh;
	    }
	    if(obj.filterLabels()[hh].find("hltDiMuonForTau3MuDzFiltered0p3") != std::string::npos) {  
	      muObjNumberGM = hh;
	    }
	  }
	  if (path=="HLT_DoubleMu4_3_LowMass") {
	    if(obj.filterLabels()[hh].find("hltDoubleMu43LowMassL3Filtered") != std::string::npos) {  
	      muObjNumberUniv = hh;
	    }
	  }
	}
	if(debug && muObjNumberTM>=0) std::cout << "Loop over filters done, my TM filter found" << std::endl;
	if(debug && muObjNumberGM>=0) std::cout << "Loop over filters done, my GM filter found" << std::endl;
	if(debug && muObjNumberUniv>=0) std::cout << "Loop over filters done, my Univ filter found" << std::endl;
	if(debug && muObjNumberGM<0 && muObjNumberTM<0 && muObjNumberUniv<0) std::cout << "Loop over filters done, my filters NOT found" << std::endl;
      
	// here HLT obj vs reco muon candidates
	TVector3 muTV3, objTV3;
	muTV3.SetPtEtaPhi( muon.pt(), muon.eta(), muon.phi() );
	objTV3.SetPtEtaPhi( obj.pt(), obj.eta(), obj.phi() );
	Float_t deltaR = fabs(muTV3.DeltaR(objTV3));
    
	// here HLT-muon candidates
	if (muObjNumberTM>=0 || muObjNumberGM>=0 || muObjNumberUniv>=0) {  
	  if(debug) std::cout<< "DeltaR = " << deltaR << std::endl;
	  if(debug) std::cout << "This is a muon HLT candidate" << endl;
	  if(deltaR < drForTriggerMatch_){    
	    frs[ipath]=1;    
	    if (deltaR < minDr){
	      minDr = deltaR;
	      minPt = obj.pt(); 
	    }
	    if(debug) std::cout << "This object is matched with muon: minDr = " << minDr << std::endl;
	    if(debug) std::cout << "Offline: " << muon.pt() << " " << muon.eta() << " " << muon.phi() << std::endl;
	    if(debug) std::cout << "HLT: " << obj.pt() << " " << obj.eta() << " " << obj.phi() << std::endl;
	    if(debug) std::cout << muObjNumberTM << " " << muObjNumberGM << " " << muObjNumberUniv << std::endl;
	  }
	}
	
      } // Loop over trigger object

      // Here we store the minimum between reco muon and all its matched HLT objects for this HLT path
      temp_DR[ipath]=minDr;
      temp_matched_to[ipath]=minPt;
      
    } // Loop over HLT paths

    // One vector per muon. Each vector : one element per path, corresponding to the closest HLT object
    fires.push_back(frs);                 // This is used in order to see if a reco muon fired a Trigger (1) or not (0).
    matcher.push_back(temp_matched_to);   // PT of the HLT object matched to the reco muon
    DR.push_back(temp_DR);                // Minimum dR between online and offline

    if(debug) {
      std::cout << std::endl;
      std::cout << "Summary for this reco muon: " << std::endl;
      int size1 = frs.size();
      int size2 = temp_matched_to.size();
      int size3 = temp_DR.size();
      if (size1!=size2 || size1!=size3) 
	std::cout << "problem with size: " << size1 << " " << size2 << " " << size3 << std::endl;
      else {
	std::cout << "size ok: " << size1 << std::endl;	
	for (int jj=0; jj<size1; jj++) std::cout << "fired = " << frs[jj] 
						 << ", pT HLT obj = " << temp_matched_to[jj] 
						 << ", pT RECO obj = " << muon.pt()
						 << ", dR = " << temp_DR[jj] << std::endl;
      }					 
    } // debug
    
  } // Loop over reco muons

  if(debug) std::cout << std::endl;

  // Now check for different reco muons that are matched to the same HLTObject.
  for(unsigned int path=0; path<HLTPaths_.size(); path++){
    for(unsigned int iMuo=0; iMuo<muons->size(); iMuo++){
      for(unsigned int im=(iMuo+1); im<muons->size(); im++){
	if(matcher[iMuo][path]!=1000. && matcher[iMuo][path]==matcher[im][path]){ 
	  if(DR[iMuo][path]<DR[im][path]){     // Keep the one that has the minimum DR with the HLT object
	    fires[im][path]=0;
	    matcher[im][path]=1000.;
	    DR[im][path]=1000.;                       
	  }
	  else{
	    fires[iMuo][path]=0;
	    matcher[iMuo][path]=1000.;
	    DR[iMuo][path]=1000.;                       
	  }
	}              
      }
      if(matcher[iMuo][path]!=1000.) muonIsTrigger[iMuo]=1;
    }
  }

  // Create a collection with all trg muons
  for(const pat::Muon & muon : *muons){
    unsigned int iMuo(&muon -&(muons->at(0)));
    if(muonIsTrigger[iMuo]==1){
      pat::Muon recoTriggerMuonCand(muon);
      trgmuons_out->emplace_back(recoTriggerMuonCand);
    }
  }

  // Save all reco muons passing the selection
  for(const pat::Muon & muon : *muons){
    unsigned int iMuo(&muon - &(muons->at(0)) );
    
    if( muon.pt()<ptMin_ ) continue;
    if( fabs(muon.eta())>absEtaMax_ ) continue;

    //const reco::TransientTrack muonTT((*(muon.bestTrack())),&(*bFieldHandle));  
    const reco::TransientTrack muonTT( fix_track( &(*muon.bestTrack()), 1e-8 ), &bField);

    if(!muonTT.isValid()) continue; 
    
    muons_out->emplace_back(muon);

    int isPFcand = (int) muon.isPFMuon();
    muons_out->back().addUserInt("isPFcand", isPFcand);    
    int isGlobal = (int) muon.isGlobalMuon();
    muons_out->back().addUserInt("isGlobal", isGlobal);    
    int isLoose = (int) muon.isLooseMuon();
    muons_out->back().addUserInt("isLoose", isLoose);    
    int isMedium = (int) muon.isMediumMuon();
    muons_out->back().addUserInt("isMedium", isMedium);    
    muons_out->back().addUserInt("isSoft", muon.isSoftMuon(PV));    
    muons_out->back().addUserInt("isTight", muon.isTightMuon(PV));     
    int isTracker = (int) muon.isTrackerMuon();
    muons_out->back().addUserInt("isTracker", isTracker); 
    muons_out->back().addUserInt("charge", muon.charge());
    if( muon.innerTrack().isNull() ) 
      muons_out->back().addUserInt("trackQuality", 999);
    else
      muons_out->back().addUserInt("trackQuality", (muon.innerTrack()->quality(reco::TrackBase::highPurity)));
    muons_out->back().addUserFloat("dZpv", muon.muonBestTrack()->dz(PV.position()));
    muons_out->back().addUserFloat("err_dZpv", muon.muonBestTrack()->dzError());

    // dr cut (same quantity as in HLTMuonDimuonL3Filter, to emulate HLT)
    float mudr = fabs( (- (muon.vx()-beamSpot.x0()) * muon.py() + (muon.vy()-beamSpot.y0()) * muon.px() ) / muon.pt() );
    muons_out->back().addUserFloat("dr", mudr);  

    for(unsigned int i=0; i<HLTPaths_.size(); i++){
      muons_out->back().addUserInt(HLTPaths_[i],fires[iMuo][i]);  // fired HLT or not
      std::string namedr = HLTPaths_[i]+"_dr";
      muons_out->back().addUserFloat(namedr,DR[iMuo][i]);         // dR between HLT and offline, to be studied
    }
    trans_muons_out->emplace_back(muonTT);
  }

  iEvent.put(std::move(muons_out),       "SelectedMuons");
  iEvent.put(std::move(trans_muons_out), "SelectedTransientMuons");
  iEvent.put(std::move(trgmuons_out),    "trgMuons");
}

// O. Cerri's code to deal with not positive definite covariance matrices
// https://github.com/ocerri/BPH_RDntuplizer/blob/master/plugins/VtxUtils.cc
/* Check for a not positive definite covariance matrix. If the covariance matrix is not positive definite, we force it to be positive definite by
 * adding the minimum eigenvalue to the diagonal of the covariance matrix plus `delta`.
 * See https://nhigham.com/2020/12/22/what-is-a-modified-cholesky-factorization/ */

reco::Track TriMuonTriggerSelector::fix_track(const reco::Track *tk, double delta) const {

  unsigned int i, j;
  double min_eig = 1;

  // Get the original covariance matrix. 
  reco::TrackBase::CovarianceMatrix cov = tk->covariance();

  // Convert it from an SMatrix to a TMatrixD so we can get the eigenvalues. 
  TMatrixDSym new_cov(cov.kRows);
  for (i = 0; i < cov.kRows; i++) {
    for (j = 0; j < cov.kRows; j++) {
    // Need to check for nan or inf, because for some reason these
    // cause a segfault when calling Eigenvectors().
    //
    // No idea what to do here or why this happens. 
    if (std::isnan(cov(i,j)) || std::isinf(cov(i,j)))
	cov(i,j) = 1e-6;
      new_cov(i,j) = cov(i,j);
    }
  }

  // Get the eigenvalues. 
  TVectorD eig(cov.kRows);
  new_cov.EigenVectors(eig);
  for (i = 0; i < cov.kRows; i++)
    if (eig(i) < min_eig)
      min_eig = eig(i);

  // If the minimum eigenvalue is less than zero, then subtract it from the
  // diagonal and add `delta`. 
  if (min_eig < 0) {
    for (i = 0; i < cov.kRows; i++)
      cov(i,i) -= min_eig - delta;
  }

  return reco::Track(tk->chi2(), tk->ndof(), tk->referencePoint(), tk->momentum(), tk->charge(), cov, tk->algo(), (reco::TrackBase::TrackQuality) tk->qualityMask());
}

DEFINE_FWK_MODULE(TriMuonTriggerSelector);

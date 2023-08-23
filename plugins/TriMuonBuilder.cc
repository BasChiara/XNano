#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"
#include "../interface/IsolationComputer.h"

#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"


#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/PatternTools/interface/ClosestApproachInRPhi.h"

#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicVertex.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicState.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/MultiTrackKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/TwoTrackMassKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/KinematicConstrainedVertexFitter.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h" 

#include <vector>
#include <string>
#include <cmath>
#include "TLorentzVector.h"

#include "helper.h"
#include "diMuonResonances.h"
#include "KinVtxFitter.h"

class TriMuonBuilder : public edm::global::EDProducer<> {

public:

  typedef std::vector<pat::Muon> MuonCollection;
  typedef std::vector<reco::TransientTrack> TransientTrackCollection;
  typedef std::vector<pat::PackedCandidate> PackedCandidatesCollection;

  explicit TriMuonBuilder(const edm::ParameterSet &cfg):
    l1_selection_{cfg.getParameter<std::string>("lep1Selection")},
    l2_selection_{cfg.getParameter<std::string>("lep2Selection")},
    l3_selection_{cfg.getParameter<std::string>("lep3Selection")},
    pre_vtx_selection_{cfg.getParameter<std::string>("preVtxSelection")},
    post_vtx_selection_{cfg.getParameter<std::string>("postVtxSelection")},
    src_{consumes<MuonCollection>( cfg.getParameter<edm::InputTag>("src") )},
    ttracks_src_{consumes<TransientTrackCollection>( cfg.getParameter<edm::InputTag>("transientTracksSrc") )},
    pkdCand_src_{consumes<PackedCandidatesCollection>( cfg.getParameter<edm::InputTag>("packedCandidatesSrc") )},
    met_{consumes<pat::METCollection>( cfg.getParameter<edm::InputTag>("met") )},
    PuppiMet_{consumes<pat::METCollection>( cfg.getParameter<edm::InputTag>("PuppiMet") )},
    beamSpotSrc_{consumes<reco::BeamSpot>( cfg.getParameter<edm::InputTag>("beamSpot") )},
    isoRadius_{cfg.getParameter<double>("isoRadius")},
    dBetaCone_{cfg.getParameter<double>("dBetaCone")},
    dBetaValue_{cfg.getParameter<double>("dBetaValue")}
    {
      produces<pat::CompositeCandidateCollection>("SelectedTriMuons");
    }

  ~TriMuonBuilder() override {}
  
  void produce(edm::StreamID, edm::Event&, const edm::EventSetup&) const override;
  
  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions) {}
  
private:

  bool  vetoResonances(edm::Event&, const std::vector<size_t>, pat::CompositeCandidate*) const; 

  const StringCutObjectSelector<pat::Muon> l1_selection_;    
  const StringCutObjectSelector<pat::Muon> l2_selection_;    
  const StringCutObjectSelector<pat::Muon> l3_selection_;    
  const StringCutObjectSelector<pat::CompositeCandidate> pre_vtx_selection_; 
  const StringCutObjectSelector<pat::CompositeCandidate> post_vtx_selection_;
  const edm::EDGetTokenT<MuonCollection> src_;
  const edm::EDGetTokenT<TransientTrackCollection> ttracks_src_;
  const edm::EDGetTokenT<PackedCandidatesCollection> pkdCand_src_;
  const edm::EDGetTokenT<pat::METCollection> met_;
  const edm::EDGetTokenT<pat::METCollection> PuppiMet_;
  const edm::EDGetTokenT<reco::BeamSpot> beamSpotSrc_;
  const double isoRadius_;
  const double dBetaCone_;
  const double dBetaValue_;

  const bool debug = false;
};

void TriMuonBuilder::produce(edm::StreamID, edm::Event &evt, edm::EventSetup const &) const {

  float the_MUON_SIGMA = 0.0000001;
  
  // input
  edm::Handle<MuonCollection> muons;
  evt.getByToken(src_, muons);
  
  edm::Handle<TransientTrackCollection> ttracks;
  evt.getByToken(ttracks_src_, ttracks);

  edm::Handle<PackedCandidatesCollection> pkdPFcand_hdl;
  evt.getByToken(pkdCand_src_, pkdPFcand_hdl);
  //const PackedCandidatesCollection& pkdPFcand = *pkdPFcand_hdl;

  edm::Handle<pat::METCollection> Met;
  evt.getByToken(met_, Met);
  const pat::MET &met = Met->front();// get PF Type-1 corrected MET directly available from miniAOD

  edm::Handle<pat::METCollection> P_Met;
  evt.getByToken(PuppiMet_, P_Met);
  const pat::MET &PuppiMet = P_Met->front();

  edm::Handle<reco::BeamSpot> beamSpotHandle;
  evt.getByToken(beamSpotSrc_, beamSpotHandle);
  const reco::BeamSpot& beamSpot = *beamSpotHandle;

  // output
  std::unique_ptr<pat::CompositeCandidateCollection> ret_value(new pat::CompositeCandidateCollection());
  
  for(size_t l1_idx = 0; l1_idx < muons->size(); ++l1_idx) {
    edm::Ptr<pat::Muon> l1_ptr(muons, l1_idx);
    if(!l1_selection_(*l1_ptr)) continue;
    
    for(size_t l2_idx = l1_idx + 1; l2_idx < muons->size(); ++l2_idx) {
      edm::Ptr<pat::Muon> l2_ptr(muons, l2_idx);
      if(!l2_selection_(*l2_ptr)) continue;
      if (l1_idx==l2_idx) continue;  // Muons must be different

      for(size_t l3_idx = l2_idx + 1; l3_idx < muons->size(); ++l3_idx) {
        edm::Ptr<pat::Muon> l3_ptr(muons, l3_idx);
        if(!l3_selection_(*l3_ptr)) continue;
        if(l3_idx == l1_idx || l3_idx == l2_idx) continue; // Muons must be different

        // tau candidate : set the kinematics and its daughters
        pat::CompositeCandidate muon_triplet;
        muon_triplet.setP4(l1_ptr->p4() + l2_ptr->p4() + l3_ptr->p4());
        muon_triplet.setCharge(l1_ptr->charge() + l2_ptr->charge() + l3_ptr->charge());
        muon_triplet.addUserInt("charge", muon_triplet.charge());
        //muon_triplet.addUserFloat("lep_deltaR", reco::deltaR(*l1_ptr, *l2_ptr));
    
        // Put the lepton passing the corresponding selection
        muon_triplet.addUserInt("l1_idx", l1_idx );
        muon_triplet.addUserInt("l2_idx", l2_idx );
        muon_triplet.addUserInt("l3_idx", l3_idx );

        // Use UserCands as they should not use memory but keep the Ptr itself
        muon_triplet.addUserCand("l1", l1_ptr );
        muon_triplet.addUserCand("l2", l2_ptr );    
        muon_triplet.addUserCand("l3", l3_ptr );

        // 1st KinVtx fit
        if( !pre_vtx_selection_(muon_triplet) ) continue;
        if(debug) std::cout << "  muon_triplet charge " << muon_triplet.charge() << std::endl;
        KinVtxFitter fitter(
                {ttracks->at(l1_idx), ttracks->at(l2_idx), ttracks->at(l3_idx)},
                {l1_ptr->mass(), l2_ptr->mass(), l3_ptr->mass()},
                {LEP_SIGMA, LEP_SIGMA, LEP_SIGMA} //some small sigma for the particle mass
                );
        if ( !fitter.success() ) continue;
        // save intermediate quantities after 1st fit needed for selection and to be saved in the final ntuples
        muon_triplet.addUserFloat("vtx_prob", fitter.prob());
        muon_triplet.addUserFloat("vtx_chi2", fitter.chi2());
        muon_triplet.addUserFloat("vtx_Ndof", fitter.dof());
        KinematicState fitted_cand = fitter.fitted_candidate();
        muon_triplet.addUserFloat("fitted_wovc_mass", fitter.success() ? fitted_cand.mass() : -1);
        RefCountedKinematicVertex fitted_vtx = fitter.fitted_refvtx();
        muon_triplet.addUserInt("vtx_isValid", fitted_vtx->vertexIsValid());
        if( !post_vtx_selection_(muon_triplet) ) continue;
                
        // 2nd KinVtx fit with vertex costraint
        KinematicParticleFactoryFromTransientTrack factory;
        std::vector<RefCountedKinematicParticle> particles;
        particles.emplace_back(factory.particle( ttracks->at(l1_idx), l1_ptr->mass(), 0., 0., the_MUON_SIGMA));
        particles.emplace_back(factory.particle( ttracks->at(l2_idx), l2_ptr->mass(), 0., 0., the_MUON_SIGMA));
        particles.emplace_back(factory.particle( ttracks->at(l3_idx), l3_ptr->mass(), 0., 0., the_MUON_SIGMA));
        std::vector<KinematicState> kinstate_muons;
        kinstate_muons.push_back(particles.at(0)->currentState());
        kinstate_muons.push_back(particles.at(1)->currentState());
        kinstate_muons.push_back(particles.at(2)->currentState());

        MultiTrackKinematicConstraint * vtxCostraint = new VertexKinematicConstraint(); 
        // from : https://github.com/cms-sw/cmssw/blob/master/RecoVertex/KinematicFit/src/VertexKinematicConstraint.cc
        vtxCostraint->value(kinstate_muons, fitted_vtx->position());
        KinematicConstrainedVertexFitter vc_fitter;
        RefCountedKinematicTree vc_FitTree = vc_fitter.fit(particles, vtxCostraint);
        if (vc_FitTree->isEmpty() || !vc_FitTree->isValid() || !vc_FitTree->isConsistent()) continue;
        vc_FitTree->movePointerToTheTop(); 
        KinematicState refitted_cand = vc_FitTree->currentParticle()->currentState();
        TLorentzVector Tau_vc;
        Tau_vc.SetPtEtaPhiM(refitted_cand.globalMomentum().perp(), 
			                      refitted_cand.globalMomentum().eta(),
			                      refitted_cand.globalMomentum().phi(),
			                      refitted_cand.mass());
        // Lxy BS - 3mu-vtx
        // from : https://cmssdt.cern.ch/lxr/source/HLTrigger/btau/plugins/HLTDisplacedmumuFilter.cc
        GlobalPoint fittedVtxPoint(fitted_vtx->position().x(), fitted_vtx->position().y(), fitted_vtx->position().z());
        GlobalError fittedVtxError( fitted_vtx->error().cxx(), 
                                    fitted_vtx->error().cyx(), fitted_vtx->error().cyy(), 
                                    fitted_vtx->error().czx(), fitted_vtx->error().czy(), fitted_vtx->error().czz()); 
        //GlobalError fittedVtxError = 
        GlobalPoint displacementfitVtxBS( -1.*((beamSpot.x0() - fittedVtxPoint.x()) + (fittedVtxPoint.z() - beamSpot.z0()) * beamSpot.dxdz()), 
                                          -1.*((beamSpot.y0() - fittedVtxPoint.y()) + (fittedVtxPoint.z() - beamSpot.z0()) * beamSpot.dydz()),
                                          0.0 );
        float Lxy_3muVtxBS = displacementfitVtxBS.perp();
        float errLxy_3muVtxBS = std::sqrt(fittedVtxError.rerr(displacementfitVtxBS));
        float sigLxy_3muVtxBS = Lxy_3muVtxBS/errLxy_3muVtxBS;

        // CosAlpha2D( Lxy(BS-3muV) ; P_3mu)
        math::XYZVector TauCand_refittedP3(Tau_vc.Px(), Tau_vc.Py(), 0 );
        reco::Vertex::Point dist2D_3muVtxBS(displacementfitVtxBS.x(), displacementfitVtxBS.y(), 0.);
        float CosAlpha2D_LxyP3mu = dist2D_3muVtxBS.Dot(TauCand_refittedP3)/(dist2D_3muVtxBS.R()*TauCand_refittedP3.R());
      


       // TRANSVERSE MASS
       float Tau_mT = std::sqrt(2. * Tau_vc.Perp()* met.pt() * (1 - std::cos(Tau_vc.Phi()-met.phi())));
       float Tau_Puppi_mT = std::sqrt(2. * Tau_vc.Perp()* PuppiMet.pt() * (1 - std::cos(Tau_vc.Phi()-PuppiMet.phi())));

       // VETO di-muon resonances 
       bool isToVeto = vetoResonances(evt, {l1_idx,l2_idx,l3_idx}, &muon_triplet);
       muon_triplet.addUserInt("diMuVtxFit_toVeto", isToVeto);

      // Tau candidate ISOLATION
      // custom class ... to check carefully
      IsolationComputer isoComputer = IsolationComputer(pkdPFcand_hdl, isoRadius_, 0.2, dBetaCone_);
      isoComputer.addMuonsToVeto({l1_ptr, l2_ptr, l3_ptr});
      float ptChargedFromPV = isoComputer.pTcharged_iso(muon_triplet);
      float ptChargedFromPU = isoComputer.pTcharged_PU(muon_triplet);
      float ptPhotons = isoComputer.pTphoton(muon_triplet);
      // class initiated with outer beta cone radius (NOT WORKING!!)
      //heppy::IsolationComputer isoComputer = heppy::IsolationComputer(dBetaCone_);
      //isoComputer.setPackedCandidates(pkdPFcand, -1, 0.2, 9999, true); // std::vector<pat::PackedCandidate>, fromPV_thresh, dz_thresh, dxy_thresh, also_leptons
      //float ptChargedFromPV = isoComputer.chargedAbsIso(muon_triplet, isoRadius_, 0., 0.5);
      //float ptChargedFromPU = isoComputer.puAbsIso(muon_triplet, dBetaCone_, 0., 0.5);
      //float ptPhotons       = isoComputer.photonAbsIsoRaw(muon_triplet, dBetaCone_, 0., 0.5);
      float TauAbsIsolation = ptChargedFromPV + std::max(0., ptPhotons - dBetaValue_*ptChargedFromPU);

      // useful quantities for BDT
      // muons longitudinal distance
      float dz_mu12 = ( (l1_ptr->hasUserFloat("dZpv") && l2_ptr->hasUserFloat("dZpv")) ? fabs(l1_ptr->userFloat("dZpv") - l2_ptr->userFloat("dZpv")) : -1 );
      float dz_mu13 = ( (l1_ptr->hasUserFloat("dZpv") && l3_ptr->hasUserFloat("dZpv")) ? fabs(l1_ptr->userFloat("dZpv") - l3_ptr->userFloat("dZpv")) : -1 );
      float dz_mu23 = ( (l2_ptr->hasUserFloat("dZpv") && l3_ptr->hasUserFloat("dZpv")) ? fabs(l2_ptr->userFloat("dZpv") - l3_ptr->userFloat("dZpv")) : -1 );

        // 1st KIN FIT WITHOUT VTX COSTRAINT
        //   Tau infos after 1st fit
        TVector3 Tau_wovc(fitted_cand.globalMomentum().x(),
                          fitted_cand.globalMomentum().y(),
                          fitted_cand.globalMomentum().z());
        muon_triplet.addUserFloat("fitted_wovc_pt",  Tau_wovc.Pt());
        muon_triplet.addUserFloat("fitted_wovc_eta", Tau_wovc.Eta());
        muon_triplet.addUserFloat("fitted_wovc_phi", Tau_wovc.Phi());
        // Tau vertex after fit
        muon_triplet.addUserFloat("fitted_vtxX",  fitted_vtx->position().x());
        muon_triplet.addUserFloat("fitted_vtxY",  fitted_vtx->position().y());
        muon_triplet.addUserFloat("fitted_vtxZ",  fitted_vtx->position().z());
        muon_triplet.addUserFloat("fitted_vtxEx", fitted_vtx->error().cxx());
        muon_triplet.addUserFloat("fitted_vtxEy", fitted_vtx->error().cyy());
        muon_triplet.addUserFloat("fitted_vtxEz", fitted_vtx->error().czz());
        // 2nd KIN FIT WITH VTX COSTRAINT
        muon_triplet.addUserFloat("fitted_vc_pt"  , Tau_vc.Perp()); 
        muon_triplet.addUserFloat("fitted_vc_eta" , Tau_vc.Eta());
        muon_triplet.addUserFloat("fitted_vc_phi" , Tau_vc.Phi());
        muon_triplet.addUserFloat("fitted_vc_mass", refitted_cand.mass()); 

        // MET infos
        muon_triplet.addUserFloat("MET_pt", met.pt()); 
        muon_triplet.addUserFloat("mT", Tau_mT); 
        muon_triplet.addUserInt("MET_isPf", met.isPFMET()); 
        muon_triplet.addUserFloat("PuppiMET_pt", PuppiMet.pt()); 
        muon_triplet.addUserFloat("Puppi_mT", Tau_Puppi_mT); 
        muon_triplet.addUserInt("PuppiMET_isPf", PuppiMet.isPFMET()); 

        // ISOLATION info
        muon_triplet.addUserFloat("iso_ptChargedFromPV", ptChargedFromPV);
        muon_triplet.addUserFloat("iso_ptChargedFromPU", ptChargedFromPU);
        muon_triplet.addUserFloat("iso_ptPhotons", ptPhotons);
        muon_triplet.addUserFloat("absIsolation",TauAbsIsolation);

        // useful quantities for BDT
        muon_triplet.addUserFloat("dZmu12", dz_mu12); 
        muon_triplet.addUserFloat("dZmu13", dz_mu13);
        muon_triplet.addUserFloat("dZmu23", dz_mu23);
        muon_triplet.addUserFloat("sigLxy_3muVtxBS", sigLxy_3muVtxBS);
        muon_triplet.addUserFloat("Cos2D_LxyP3mu",CosAlpha2D_LxyP3mu),
      

        // save further quantities, to be saved in the final ntuples: muons before fit
        // Muons post fit are saved only after the very final B fit
        muon_triplet.addUserFloat("mu1_pt",  l1_ptr->pt());
        muon_triplet.addUserFloat("mu1_eta", l1_ptr->eta());
        muon_triplet.addUserFloat("mu1_phi", l1_ptr->phi());
        muon_triplet.addUserFloat("mu1_dr",  l1_ptr->userFloat("dr"));
        muon_triplet.addUserInt("mu1_charge" ,l1_ptr->charge());
        muon_triplet.addUserInt("mu1_trackQuality",  l1_ptr->userInt("trackQuality"));
        muon_triplet.addUserFloat("mu2_pt",  l2_ptr->pt());
        muon_triplet.addUserFloat("mu2_eta", l2_ptr->eta());
        muon_triplet.addUserFloat("mu2_phi", l2_ptr->phi());
        muon_triplet.addUserFloat("mu2_dr",  l2_ptr->userFloat("dr"));   
        muon_triplet.addUserInt("mu2_charge" ,l2_ptr->charge());
        muon_triplet.addUserInt("mu2_trackQuality",  l2_ptr->userInt("trackQuality"));
        muon_triplet.addUserFloat("mu3_pt",  l3_ptr->pt());
        muon_triplet.addUserFloat("mu3_eta", l3_ptr->eta());
        muon_triplet.addUserFloat("mu3_phi", l3_ptr->phi());
        muon_triplet.addUserFloat("mu3_dr",  l3_ptr->userFloat("dr"));   
        muon_triplet.addUserInt("mu3_charge" ,l3_ptr->charge());
        muon_triplet.addUserInt("mu3_trackQuality",  l3_ptr->userInt("trackQuality"));

        // save further quantities, to be saved in the final ntuples: fired paths
        muon_triplet.addUserInt("mu1_fired_Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15_Charge1", l1_ptr->userInt("HLT_Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15_Charge1"));
        muon_triplet.addUserInt("mu1_fired_Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15", l1_ptr->userInt("HLT_Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15"));
        muon_triplet.addUserInt("mu1_fired_DoubleMu4_3_LowMass", l1_ptr->userInt("HLT_DoubleMu4_3_LowMass"));


        muon_triplet.addUserInt("mu2_fired_Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15_Charge1", l2_ptr->userInt("HLT_Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15_Charge1"));
        muon_triplet.addUserInt("mu2_fired_Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15", l2_ptr->userInt("HLT_Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15"));
        muon_triplet.addUserInt("mu2_fired_DoubleMu4_3_LowMass", l2_ptr->userInt("HLT_DoubleMu4_3_LowMass"));

        muon_triplet.addUserInt("mu3_fired_Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15_Charge1", l3_ptr->userInt("HLT_Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15_Charge1"));
        muon_triplet.addUserInt("mu3_fired_Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15", l3_ptr->userInt("HLT_Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15"));
        muon_triplet.addUserInt("mu3_fired_DoubleMu4_3_LowMass", l3_ptr->userInt("HLT_DoubleMu4_3_LowMass"));
        
        // push in the event
        ret_value->push_back(muon_triplet);
      }
    }
  }

  evt.put(std::move(ret_value),  "SelectedTriMuons");
}


bool TriMuonBuilder::vetoResonances(edm::Event& iEvt, const std::vector<size_t> tauMu_idcs, pat::CompositeCandidate* tau_cand) const {

    bool isMatchingResonance = false; 
    bool debug = false;
    if (tauMu_idcs.size() != 3){
      std::cout << "ERROR in TriMuonBuilder::vetoResonances() : Tau-cand must be made of exactly 3 muons" << std::endl;
      return -1;
    }
    edm::Handle<MuonCollection> all_muons;
    iEvt.getByToken(src_, all_muons);
    edm::Handle<TransientTrackCollection> all_muTtracks;
    iEvt.getByToken(ttracks_src_, all_muTtracks);

    const float fitProb_min = 0.05;
    float best_prob = -1., best_mass = 0.;

    for(size_t mu_idx = 0; mu_idx < all_muons->size(); ++mu_idx) {
      // not in the tau candidate
      if ( std::find(tauMu_idcs.begin(), tauMu_idcs.end(), mu_idx) != tauMu_idcs.end()){
        if(debug) std::cout << "  Skip muon " << mu_idx << std::endl;
        continue;
      }
      edm::Ptr<pat::Muon> mu(all_muons, mu_idx);
      for(size_t Tmu_idx = 0; Tmu_idx < 3; ++Tmu_idx){
          edm::Ptr<pat::Muon> TauMu(all_muons, tauMu_idcs[Tmu_idx]);
          // opposite charge
          if((TauMu->charge() + mu->charge()) != 0) continue;
          if(debug) std::cout << " Veto resonance mu_tau" << tauMu_idcs[Tmu_idx] <<" + mu_" << mu_idx << std::endl;
          // fit
          KinVtxFitter fitter(
                  {all_muTtracks->at(tauMu_idcs[Tmu_idx]), all_muTtracks->at(mu_idx)},
                  {TauMu->mass(), mu->mass()},
                  {LEP_SIGMA, LEP_SIGMA} //some small sigma for the particle mass
                  );
          if ( !fitter.success() ) continue;
          if(debug) std::cout << " Fit-done : fit prob = "<< fitter.prob() << " di-muon mass = " << fitter.fitted_candidate().mass() << std::endl; 
          // require good di-muon vertex (prob > 5%)
          if(fitter.prob() < fitProb_min) continue;
          float fitMass = fitter.fitted_candidate().mass();
          if(fitter.prob() > best_prob){
            best_prob = fitter.prob();
            best_mass = fitMass;
          }
          // check if compatibility with di-muon resonances 
          if(debug) std::cout << " checking for matching resonance..." << std::endl;
          for(std::vector< std::pair<float, float> >::iterator reso = resonancesToVeto.begin(); reso != resonancesToVeto.end(); ++reso){
            if( fabs( (fitMass - reso->first)/reso->second)  < SIGMA_TO_EXCLUDE ){
               std::cout << " matching found " << fabs( (fitMass - reso->first)/reso->second) << std::endl; 
               isMatchingResonance = true;
            }
          }// loop on known di-muon resonances  
      }// loop on mu in Tau cand 
    }// loop on mu outside the triplet
    tau_cand->addUserFloat("diMuVtxFit_bestProb", best_prob);
    tau_cand->addUserFloat("diMuVtxFit_bestMass", best_mass);
    return isMatchingResonance;

}//vetoResonances({


#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(TriMuonBuilder);

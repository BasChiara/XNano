#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
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
#include "TLorentzVector.h"

#include "helper.h"
#include "KinVtxFitter.h"

class TriMuonBuilder : public edm::global::EDProducer<> {

public:

  typedef std::vector<pat::Muon> MuonCollection;
  typedef std::vector<reco::TransientTrack> TransientTrackCollection;

  explicit TriMuonBuilder(const edm::ParameterSet &cfg):
    l1_selection_{cfg.getParameter<std::string>("lep1Selection")},
    l2_selection_{cfg.getParameter<std::string>("lep2Selection")},
    l3_selection_{cfg.getParameter<std::string>("lep3Selection")},
    pre_vtx_selection_{cfg.getParameter<std::string>("preVtxSelection")},
    post_vtx_selection_{cfg.getParameter<std::string>("postVtxSelection")},
    src_{consumes<MuonCollection>( cfg.getParameter<edm::InputTag>("src") )},
    ttracks_src_{consumes<TransientTrackCollection>( cfg.getParameter<edm::InputTag>("transientTracksSrc") )},
    beamSpotSrc_{consumes<reco::BeamSpot>( cfg.getParameter<edm::InputTag>("beamSpot") )}
    {
      produces<pat::CompositeCandidateCollection>("SelectedTriMuons");
    }

  ~TriMuonBuilder() override {}
  
  void produce(edm::StreamID, edm::Event&, const edm::EventSetup&) const override;
  
  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions) {}
  
private:
  const StringCutObjectSelector<pat::Muon> l1_selection_;    
  const StringCutObjectSelector<pat::Muon> l2_selection_;    
  const StringCutObjectSelector<pat::Muon> l3_selection_;    
  const StringCutObjectSelector<pat::CompositeCandidate> pre_vtx_selection_; 
  const StringCutObjectSelector<pat::CompositeCandidate> post_vtx_selection_;
  const edm::EDGetTokenT<MuonCollection> src_;
  const edm::EDGetTokenT<TransientTrackCollection> ttracks_src_;
  const edm::EDGetTokenT<reco::BeamSpot> beamSpotSrc_;
  const bool debug = false;
};

void TriMuonBuilder::produce(edm::StreamID, edm::Event &evt, edm::EventSetup const &) const {

  float the_MUON_SIGMA = 0.0000001;
  
  // input
  edm::Handle<MuonCollection> muons;
  evt.getByToken(src_, muons);
  
  edm::Handle<TransientTrackCollection> ttracks;
  evt.getByToken(ttracks_src_, ttracks);

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
        KinematicState fitted_cand = fitter.fitted_candidate();
        muon_triplet.addUserFloat("fitted_wovc_mass", fitter.success() ? fitted_cand.mass() : -1);
        RefCountedKinematicVertex fitted_vtx = fitter.fitted_refvtx();
        muon_triplet.addUserInt("vtx_isValid", fitted_vtx->vertexIsValid());
        if( !post_vtx_selection_(muon_triplet) ) continue;
        //if(!fitted_vtx->vertexIsValid()) continue;
        
        // 2nd KinVtx fit with vertex costraint
        KinematicParticleFactoryFromTransientTrack factory;
        std::vector<RefCountedKinematicParticle> particles;
        particles.emplace_back(factory.particle( ttracks->at(l1_idx), l1_ptr->mass(), 0., 0., the_MUON_SIGMA));
        particles.emplace_back(factory.particle( ttracks->at(l2_idx), l2_ptr->mass(), 0., 0., the_MUON_SIGMA));
        particles.emplace_back(factory.particle( ttracks->at(l3_idx), l3_ptr->mass(), 0., 0., the_MUON_SIGMA));
        //std::vector<RefCountedKinematicParticle> fitted_muons = fitter.fitted_children();
        //kinstate_muons.push_back(fitter.fitted_daughter(0));
        //kinstate_muons.push_back(fitter.fitted_daughter(1));
        //kinstate_muons.push_back(fitter.fitted_daughter(2));
        std::vector<KinematicState> kinstate_muons;
        kinstate_muons.push_back(particles.at(0)->currentState());
        kinstate_muons.push_back(particles.at(1)->currentState());
        kinstate_muons.push_back(particles.at(2)->currentState());

        MultiTrackKinematicConstraint * vtxCostraint = new VertexKinematicConstraint(); // from : https://github.com/cms-sw/cmssw/blob/master/RecoVertex/KinematicFit/src/VertexKinematicConstraint.cc
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


        //// DCA between the two muons (used at HLT)
        //float DCA = 10.;
        //TrajectoryStateClosestToPoint mu1TS = (ttracks->at(l1_idx)).impactPointTSCP();
        //TrajectoryStateClosestToPoint mu2TS = (ttracks->at(l2_idx)).impactPointTSCP();
        //if (mu1TS.isValid() && mu2TS.isValid()) {
        //    ClosestApproachInRPhi cApp;
        //    cApp.calculate(mu1TS.theState(), mu2TS.theState());
        //    if (cApp.status()) DCA = cApp.distance();
        //}
        //muon_triplet.addUserFloat("DCA", DCA); 

        //// Lxy (used at HLT)  
        //// HLTrigger/btau/plugins/HLTDisplacedmumuFilter.cc
        //math::XYZVector pperp(l1_ptr->px() + l2_ptr->px(), l1_ptr->py() + l2_ptr->py(), 0.);
        //GlobalError fitted_vtx_err( fitted_vtx->error().cxx(), fitted_vtx->error().cyx(), fitted_vtx->error().cyy(), fitted_vtx->error().czx(), fitted_vtx->error().czy(), fitted_vtx->error().czz() );
        //GlobalPoint dispFromBS( -1*( (beamSpot.x0() - fitted_vtx->position().x()) + (fitted_vtx->position().z() - beamSpot.z0()) * beamSpot.dxdz()), -1*((beamSpot.y0() - fitted_vtx->position().y()) + (fitted_vtx->position().z() - beamSpot.z0()) * beamSpot.dydz()), 0);
        //float lxy = dispFromBS.perp();
        //float lxyerr = sqrt(fitted_vtx_err.rerr(dispFromBS));
        //float lxySign = lxy/lxyerr;
        //muon_triplet.addUserFloat("LxySign", lxySign); 

        //// CosAlpha (used at HLT)  
        //reco::Vertex::Point vperp(dispFromBS.x(),dispFromBS.y(),0.);
        //float cosAlpha = vperp.Dot(pperp)/(vperp.R()*pperp.R());
        //muon_triplet.addUserFloat("cosAlpha", cosAlpha); 

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

        //muon_triplet.addUserFloat("mu1_dr_Dimuon25_Jpsi",        l1_ptr->userFloat("HLT_Dimuon25_Jpsi_dr"));
        //muon_triplet.addUserFloat("mu2_dr_Dimuon25_Jpsi",        l2_ptr->userFloat("HLT_Dimuon25_Jpsi_dr"));
        
        // push in the event
        ret_value->push_back(muon_triplet);
      }
    }
  }

  evt.put(std::move(ret_value),  "SelectedTriMuons");
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(TriMuonBuilder);

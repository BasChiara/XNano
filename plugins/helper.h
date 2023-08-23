#ifndef PhysicsTools_XNano_helpers
#define PhysicsTools_XNano_helpers

#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/GeometryCommonDetAlgo/interface/Measurement1D.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/GeometryCommonDetAlgo/interface/GlobalError.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "DataFormats/GeometryVector/interface/PV3DBase.h"
#include "Math/LorentzVector.h"

#include "TrackingTools/GeomPropagators/interface/AnalyticalImpactPointExtrapolator.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicTree.h"
#include "RecoVertex/VertexPrimitives/interface/ConvertToFromReco.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "TVector3.h"

#include <vector>
#include <algorithm>
#include <limits>
#include <memory>

typedef std::vector<reco::TransientTrack> TransientTrackCollection;

// RK code
constexpr float MUON_MASS = 0.10565837;
constexpr float ELECTRON_MASS = 0.000511;
constexpr float K_MASS = 0.493677;
constexpr float PI_MASS = 0.139571;
constexpr float LEP_SIGMA = 0.0000001;
constexpr float K_SIGMA = 0.000016;
constexpr float PI_SIGMA = 0.000016;
// PDG
constexpr float KSHORT_MASS = 0.497611;
constexpr float X3872_MASS = 3.87165;
constexpr float RHO_MASS = 0.77526;
constexpr float JPSI_MASS = 3.096900;
constexpr float B0_MASS = 5.27965;
constexpr float W_MASS = 80.377;
constexpr float W_WIDTH = 2.085;
// S.P.'s code
constexpr float KSHORT_SIGMA = 0.000013;
// Guess
constexpr float RHO_SIGMA = 0.000016;
constexpr float JPSI_SIGMA = 0.000016;


inline std::pair<float, float> min_max_dr(const std::vector< edm::Ptr<reco::Candidate> > & cands) {
  float min_dr = std::numeric_limits<float>::max();
  float max_dr = 0.;
  for(size_t i = 0; i < cands.size(); ++i) {
    for(size_t j = i+1; j < cands.size(); ++j) {
      float dr = reco::deltaR(*cands.at(i), *cands.at(j));
      min_dr = std::min(min_dr, dr);
      max_dr = std::max(max_dr, dr);
    }
  }
  return std::make_pair(min_dr, max_dr);
}

template<typename FITTER, typename LORENTZ_VEC>
inline double cos_theta_2D(const FITTER& fitter, const reco::BeamSpot &bs, const LORENTZ_VEC& p4) {
  if(!fitter.success()) return -2;
  GlobalPoint point = fitter.fitted_vtx();
  auto bs_pos = bs.position(point.z());
  math::XYZVector delta(point.x() - bs_pos.x(), point.y() - bs_pos.y(), 0.);
  math::XYZVector pt(p4.px(), p4.py(), 0.);
  double den = (delta.R() * pt.R());
  return (den != 0.) ? delta.Dot(pt)/den : -2;
}

template<typename LORENTZ_VEC>
inline double cos_theta_2D(float vx, float vy, float vz, const reco::BeamSpot &bs, const LORENTZ_VEC& p4) {
  auto bs_pos = bs.position(vz);
  math::XYZVector delta(vx - bs_pos.x(), vy - bs_pos.y(), 0.);
  math::XYZVector pt(p4.Px(), p4.Py(), 0.);
  double den = (delta.R() * pt.R());
  return (den != 0.) ? delta.Dot(pt)/den : -2;
}

template<typename FITTER, typename LORENTZ_VEC>
inline double cos_theta_2D(const FITTER& fitter, const reco::Vertex &pv, const LORENTZ_VEC& p4) {
  if(!fitter.success()) return -2;
  GlobalPoint point = fitter.fitted_vtx();
  math::XYZVector delta(point.x() - pv.x(), point.y() - pv.y(), 0.);
  math::XYZVector pt(p4.px(), p4.py(), 0.);
  double den = (delta.R() * pt.R());
  return (den != 0.) ? delta.Dot(pt)/den : -2;
}

template<typename LORENTZ_VEC>
inline double cos_theta_2D(float vx, float vy, float vz, const reco::Vertex &pv, const LORENTZ_VEC& p4) {
  math::XYZVector delta(vx - pv.x(), vy - pv.y(), 0.);
  math::XYZVector pt(p4.Px(), p4.Py(), 0.);
  double den = (delta.R() * pt.R());
  return (den != 0.) ? delta.Dot(pt)/den : -2;
}

template<typename FITTER>
inline Measurement1D l_xy(const FITTER& fitter, const reco::BeamSpot &bs) {
  if(!fitter.success()) return {-2, -2};
  GlobalPoint point = fitter.fitted_vtx();
  GlobalError err = fitter.fitted_vtx_uncertainty();
  auto bs_pos = bs.position(point.z());
  GlobalPoint delta(point.x() - bs_pos.x(), point.y() - bs_pos.y(), 0.);  
  return {delta.perp(), sqrt(err.rerr(delta))};
}

template<typename FITTER>
inline Measurement1D l_xy(const FITTER& fitter, const reco::Vertex &pv) {
  if(!fitter.success()) return {-2, -2};
  GlobalPoint point = fitter.fitted_vtx();
  GlobalError err = fitter.fitted_vtx_uncertainty();
  GlobalPoint delta(point.x() - pv.x(), point.y() - pv.y(), 0.);  
  return {delta.perp(), sqrt(err.rerr(delta))};
}

inline float CosA(GlobalPoint & dist, ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> & Bp4)
{
    math::XYZVector vperp(dist.x(),dist.y(),0);
    math::XYZVector pperp(Bp4.Px(),Bp4.Py(),0); 
    return std::move(vperp.Dot(pperp)/(vperp.R()*pperp.R()));
}


inline std::pair<double,double> computeDCA(const reco::TransientTrack& trackTT,
					   const reco::BeamSpot& beamSpot)
{
  double DCABS    = -1.;
  double DCABSErr = -1.;

  TrajectoryStateClosestToPoint theDCAXBS = 
    trackTT.trajectoryStateClosestToPoint(GlobalPoint(beamSpot.position().x(),beamSpot.position().y(),beamSpot.position().z()));
  if (theDCAXBS.isValid()) {
    DCABS    = theDCAXBS.perigeeParameters().transverseImpactParameter();
    DCABSErr = theDCAXBS.perigeeError().transverseImpactParameterError();
  }

  return std::make_pair(DCABS,DCABSErr);
}

inline std::pair<bool, Measurement1D> absoluteImpactParameter(const TrajectoryStateOnSurface& tsos,
                                                              RefCountedKinematicVertex vertex,
                                                              VertexDistance& distanceComputer){
  if (!tsos.isValid()) {
      return std::pair<bool, Measurement1D>(false, Measurement1D(0., 0.));
  }
  GlobalPoint refPoint = tsos.globalPosition();
  GlobalError refPointErr = tsos.cartesianError().position();
  GlobalPoint vertexPosition = vertex->vertexState().position();
  GlobalError vertexPositionErr = RecoVertex::convertError(vertex->vertexState().error());
  return std::pair<bool, Measurement1D>(true,
                                        distanceComputer.distance(VertexState(vertexPosition, vertexPositionErr), VertexState(refPoint, refPointErr)));
}


inline std::pair<bool, Measurement1D> absoluteImpactParameter3D(const TrajectoryStateOnSurface& tsos,
                                                                RefCountedKinematicVertex vertex){
  VertexDistance3D dist;
  return absoluteImpactParameter(tsos, vertex, dist);
}


inline std::pair<bool, Measurement1D> absoluteTransverseImpactParameter(const TrajectoryStateOnSurface& tsos,
                                                                        RefCountedKinematicVertex vertex){
  VertexDistanceXY dist;
  return absoluteImpactParameter(tsos, vertex, dist);
}


inline std::pair<bool, Measurement1D> signedImpactParameter3D(const TrajectoryStateOnSurface& tsos,
                                                              RefCountedKinematicVertex vertex,
                                                              const reco::BeamSpot &bs, double pv_z){
  VertexDistance3D dist;

  std::pair<bool,Measurement1D> result = absoluteImpactParameter(tsos, vertex, dist);
  if (!result.first)
    return result;

  //Compute Sign
  auto bs_pos = bs.position(vertex->vertexState().position().z());
  GlobalPoint impactPoint = tsos.globalPosition();
  GlobalVector IPVec(impactPoint.x() - vertex->vertexState().position().x(),       
                     impactPoint.y() - vertex->vertexState().position().y(),        
                     impactPoint.z() - vertex->vertexState().position().z());

  GlobalVector direction(vertex->vertexState().position().x() - bs_pos.x(), 
                         vertex->vertexState().position().y() - bs_pos.y(), 
                         vertex->vertexState().position().z() - pv_z);

  double prod = IPVec.dot(direction);
  double sign = (prod >= 0) ? 1. : -1.;

  //Apply sign to the result
  return std::pair<bool, Measurement1D>(result.first, Measurement1D(sign * result.second.value(), result.second.error()));

}

 
inline std::pair<bool, Measurement1D> signedTransverseImpactParameter(const TrajectoryStateOnSurface& tsos,
                                                                      RefCountedKinematicVertex vertex,
                                                                      const reco::BeamSpot &bs){
  VertexDistanceXY dist;

  std::pair<bool,Measurement1D> result = absoluteImpactParameter(tsos, vertex, dist);
  if (!result.first)
    return result;

  //Compute Sign
  auto bs_pos = bs.position(vertex->vertexState().position().z());
  GlobalPoint impactPoint = tsos.globalPosition();
  GlobalVector IPVec(impactPoint.x() - vertex->vertexState().position().x(), impactPoint.y() - vertex->vertexState().position().y(), 0.);
  GlobalVector direction(vertex->vertexState().position().x() - bs_pos.x(), 
                         vertex->vertexState().position().y() - bs_pos.y(), 0);

  double prod = IPVec.dot(direction);
  double sign = (prod >= 0) ? 1. : -1.;

  //Apply sign to the result
  return std::pair<bool, Measurement1D>(result.first, Measurement1D(sign * result.second.value(), result.second.error()));

}

#endif

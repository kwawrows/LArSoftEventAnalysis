///////////////////////////////////////////////////////////////////////
// Class:       EventAnaTree
// Plugin Type: analyzer (art v3_05_01)
// File:        EventAnaTree_module.cc
// Generated Feb 14 04:39:22 2023 by Klaudia Wawrowska
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "art_root_io/RootIOPolicy.h"
#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"

#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Persistency/Common/FindManyP.h" // for associations
#include "canvas/Persistency/Common/PtrVector.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/WireGeo.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_vectors.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCFlux.h"
#include "lardataobj/MCBase/MCShower.h"
#include "lardataobj/MCBase/MCTrack.h"
#include "lardataobj/MCBase/MCStep.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "lardataobj/Simulation/AuxDetSimChannel.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RawData/raw.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/PhotonBackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/TrackHitMeta.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "larreco/SpacePointSolver/Solver.h"
#include "lardataobj/RecoBase/PointCharge.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RecoBase/EndPoint2D.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "lardataobj/RawData/RDTimeStamp.h"
#include "larreco/RecoAlg/TrackMomentumCalculator.h"
#include "lardataobj/AnalysisBase/CosmicTag.h"
#include "lardataobj/AnalysisBase/FlashMatch.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "lardataobj/AnalysisBase/BackTrackerMatchingData.h"
#include "larevt/SpaceChargeServices/SpaceChargeService.h"
#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"


//Additional framework
#include "TTree.h"
#include <string>
#include "TFile.h"
#include "TDirectory.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TProfile.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TMath.h"
#include "TGraphErrors.h"
#include "TMinuit.h"
#include "TString.h"
#include "TTimeStamp.h"
#include "TVectorD.h"
#include "TCanvas.h"
#include "TFrame.h"
#include "TLine.h"
#include "TAxis.h"
#include "TTimeStamp.h"
#include "TVector3.h"
#include <vector>
#include <fstream>
#include "TPaveStats.h"
#include <iostream>
#include <string>
#include "math.h"
#include "stdio.h"
#include <iterator>
#include "TGeoMatrix.h"

namespace dune {
  class EventAnaTree;
}


class dune::EventAnaTree : public art::EDAnalyzer {
public:
  explicit EventAnaTree(fhicl::ParameterSet const& p);

  EventAnaTree(EventAnaTree const&) = delete;
  EventAnaTree(EventAnaTree&&) = delete;
  EventAnaTree& operator=(EventAnaTree const&) = delete;
  EventAnaTree& operator=(EventAnaTree&&) = delete;

  // function initializations
  void beginJob() override;
  void endJob() override;
  void analyze(const art::Event& e) override;

private:

  void ResetVariables();

  TTree *fTree; 

  bool fSaveGenieInfo;

  std::string fTruthLabel; //which module produced simulation
  std::string fParticleLabel; //PFParticles 
  std::string fTrackLabel;
  std::string fShowerLabel; // 01.12
  std::string fHitLabel;
  std::string fGenieLabel;
  std::string fCaloLabel; 
  std::string fFlashT0Label;
  std::string fMCT0Label;
  std::string fOpHitModuleLabel;
  std::string fPandoraVertexLabel;
  std::string fSpacePointLabel; 

  // ===========================================================
  //                          EVENT VARIABLES 
  // ============================================================

  unsigned int fEvent; 
  unsigned int fRun;
  unsigned int fSubRun;


  //Genie / generator info
  unsigned int fgenie_no_primaries; 
  std::vector<int>    fgenie_primaries_pdg; 
  std::vector<float>  fgenie_Eng; 
  std::vector<float>  fgenie_Px;
  std::vector<float>  fgenie_Py;
  std::vector<float>  fgenie_Pz;
  std::vector<float>  fgenie_P;
  std::vector<int>    fgenie_status_code; 

  std::vector<int>    fnuPDG_truth; 
  std::vector<int>    fccnc_truth; // 0 = CC, 1 = NC
  std::vector<int>    fmode_truth; // 0 = QE/El, 1=RES, 2=DIS, 3= Coherent production
  std::vector<float>  fhitnuc_truth;//hit nucleon 
  std::vector<float>  fnu_vx_truth; // neutrino vertex x 
  std::vector<float>  fnu_vy_truth; // neutrino vertex y 
  std::vector<float>  fnu_vz_truth; // neutrino vertex z 
  std::vector<float>  fnu_dcosx_truth; 
  std::vector<float>  fnu_dcosy_truth; 
  std::vector<float>  fnu_dcosz_truth; 
  std::vector<float>  flep_mom_truth;
  std::vector<float>  flep_dcosx_truth;
  std::vector<float>  flep_dcosy_truth;
  std::vector<float>  flep_dcosz_truth;


  //Geant/truth info
  unsigned int fnPrimaries; // no. of primary particles 
  unsigned int fnGeantParticles; //total number of geant particles 
  std::vector<int>     fTrackId;
  std::vector<float>   fMother; 
  std::vector<float>   fEng;
  std::vector<float>   fEkin;
  std::vector<float>   fMass;
  std::vector<int>     fPdg;
  std::vector<float>   fP;
  std::vector<float>   fPx;
  std::vector<float>   fPy;
  std::vector<float>   fPz;
  std::vector<int>     fND; //# daughters 
  std::vector<double>  fstartX;
  std::vector<double>  fstartY;
  std::vector<double>  fstartZ;
  std::vector<double>  fendX;
  std::vector<double>  fendY;
  std::vector<double>  fendZ;

  //TPC info
  unsigned int fnHits; //number of hits 
  unsigned int fnColHits; // # of hits on collection plane 
  std::vector<int>     fhit_MCparticle; //pdg code of MC particle that produced this hit 
  std::vector<int>     fhit_tpc;
  std::vector<int>     fhit_channel;
  std::vector<float>   fhit_time;
  std::vector<float>   fhit_SADC;
  std::vector<float>   fhit_wire;
  std::vector<float>   fhit_charge;
  std::vector<int>     fhit_plane;
  std::vector<float>   fhit_width; //width RMS in time
  std::vector<int>     fhit_TOT; // hit width in ticks
  std::vector<float>   fhit_trueX;
  std::vector<float>   fhit_trueY;
  std::vector<float>   fhit_trueZ;
  std::vector<int>     fhit_clusterId; 

  //collection info
  std::vector<int>     fcolhit_tpc;
  std::vector<int>     fcolhit_channel;
  std::vector<float>   fcolhit_time;
  std::vector<float>   fcolhit_SADC;
  std::vector<float>   fcolhit_wire;
  std::vector<float>   fcolhit_charge;
  std::vector<float>   fcolhit_width; 
  std::vector<float>   fcolhit_trueX;
  std::vector<float>   fcolhit_trueY;
  std::vector<float>   fcolhit_trueZ;

  //Pandora info 
  unsigned int         fnPFParticles_pandora;
  unsigned int         fnPrimaries_pandora; 
  std::vector<int>     fNDPrimary_pandora; //number of reco primary daughters 
  std::vector<int>     fPrimaryID_pandora; 
  std::vector<int>     fPFPpdg_pandora;
  std::vector<int>     fPFPparent_pandora;
  std::vector<int>     fPFPID_pandora; 

 //Pandora track-like PFP
  unsigned int         fnRecoTracks_pandora; //# of reco tracks per event
  std::vector<int>     fTrackID_pandora;
  std::vector<float>   fTrackLength_pandora;
  std::vector<float>   fTrackStartx_pandora; 
  std::vector<float>   fTrackStarty_pandora; 
  std::vector<float>   fTrackStartz_pandora; 
  std::vector<float>   fTrackEndx_pandora;
  std::vector<float>   fTrackEndy_pandora;
  std::vector<float>   fTrackEndz_pandora;

  //Pandora PFP vertex
  unsigned int fnVertices_pandora;
  unsigned int fnPrimaryVertices_pandora;
  std::vector<float> fVrtxX_pandora; 
  std::vector<float> fVrtxY_pandora;
  std::vector<float> fVrtxZ_pandora;
  std::vector<float> fPrimaryVrtxX_pandora; 
  std::vector<float> fPrimaryVrtxY_pandora;
  std::vector<float> fPrimaryVrtxZ_pandora;
  std::vector<float> fPrimaryVrtxPDG_pandora;


  //Pandora Cluster information
  unsigned int fnClusters_pandora;
  std::vector<float>  fClusterID_pandora;
  std::vector<float>  fClusterView_pandora;
  std::vector<float>  fClusterNHits_pandora;
  std::vector<float>  fClusterIntegral_pandora;
  std::vector<float>  fClusterIntegralAverage_pandora;
  std::vector<float>  fClusterSummedADC_pandora;
  std::vector<float>  fClusterSummedADCaverage_pandora;
  std::vector<float>  fClusterStartAngle_pandora;
  std::vector<float>  fClusterEndAngle_pandora;
  std::vector<float>  fClusterStartOpeningAngle_pandora;
  std::vector<float>  fClusterEndOpeningAngle_pandora;

  //Reco hits associated with the above reco tracks 
  std::vector<int>         fnTrkHits_pandora;  
  std::vector<int>         fTrkHit_trkID; // ID of the track the hit belongs to
  std::vector<int>         fTrkHit_MCparticle;
  std::vector<int>         fTrkHit_PFParticle; 
  std::vector<float>       fTrkHit_tpc;
  std::vector<float>       fTrkHit_channel;
  std::vector<float>       fTrkHit_time;
  std::vector<float>       fTrkHit_SADC;
  std::vector<float>       fTrkHit_wire;
  std::vector<float>       fTrkHit_charge;
  std::vector<int>         fTrkHit_plane;
  std::vector<float>       fTrkHit_trueX;
  std::vector<float>       fTrkHit_trueY;
  std::vector<float>       fTrkHit_trueZ;
  std::vector<int>         fTrkHit_clusterID; 
  

  //Pandora shower-like PFP info
  unsigned int          fnRecoShowers_pandora;
  std::vector<int>      fShowerID_pandora;
  std::vector<double>   fShowerLength_pandora;
  std::vector<float>    fShowerStartx_pandora;
  std::vector<float>    fShowerStarty_pandora;
  std::vector<float>    fShowerStartz_pandora;

  //Reco hits associated with the above reco showers 
  std::vector<int>      fnShowerHits_pandora; 
  std::vector<int>      fShowerHit_ShowerID;
  std::vector<int>      fShowerHit_MCparticle;
  std::vector<int>      fShowerHit_PFParticle;
  std::vector<float>    fShowerHit_tpc;
  std::vector<float>    fShowerHit_channel;
  std::vector<float>    fShowerHit_time;
  std::vector<float>    fShowerHit_SADC;
  std::vector<float>    fShowerHit_wire;
  std::vector<float>    fShowerHit_charge;
  std::vector<int>      fShowerHit_plane;
  std::vector<float>    fShowerHit_trueX;
  std::vector<float>    fShowerHit_trueY;
  std::vector<float>    fShowerHit_trueZ;
  std::vector<int>      fShowerHit_clusterID; 

  //Pandora Calorimetry 
  std::vector<int>                     fCalo_TrackID;
  std::vector< std::vector<float> >    fCalo_TrackResidualRange;
  std::vector< std::vector<float> >    fCalo_TrackdEdx;
    
 
  art::ServiceHandle<geo::Geometry> geo;
  art::ServiceHandle<cheat::PhotonBackTrackerService> pbt_serv;
};

// Constructor with config. parameters 
dune::EventAnaTree::EventAnaTree(fhicl::ParameterSet const& p) :
  EDAnalyzer{p}//,
{
  //Module label configs
  fTruthLabel          = p.get<std::string>("TruthLabel");
  fHitLabel            = p.get<std::string>("HitLabel");
  fGenieLabel          = p.get<std::string>("GenieLabel");
  fParticleLabel       = p.get<std::string>("ParticleLabel");
  fTrackLabel          = p.get<std::string>("TrackLabel");
  fCaloLabel           = p.get<std::string>("CaloLabel");
  fFlashT0Label        = p.get<std::string>("FlashT0Label");
  //fMCT0Label         = p.get<std::string>("MCT0Label");
  fOpHitModuleLabel    = p.get<std::string>("OpHitModuleLabel");
  fShowerLabel         = p.get<std::string>("ShowerLabel");
  fPandoraVertexLabel  = p.get<std::string>("PandoraVertexLabel");
  fSpacePointLabel     = p.get<std::string>("SpacePointLabel");
  fSaveGenieInfo       = p.get<bool>("SaveGenieInfo", 0);
}


void dune::EventAnaTree::ResetVariables()
{

  fEvent = fRun = fSubRun = -1;
  fgenie_no_primaries = 0;
  fnPrimaries = fnGeantParticles = 0;
  fnHits = fnColHits = 0;
  fnPFParticles_pandora = fnPrimaries_pandora =  0; 
  fnVertices_pandora = fnPrimaryVertices_pandora = 0; 
  fnClusters_pandora = 0;
  fnRecoTracks_pandora = 0;
  fnRecoShowers_pandora = 0;  

    
  fgenie_primaries_pdg.clear();
  fgenie_Eng.clear();
  fgenie_Px.clear();
  fgenie_Py.clear();
  fgenie_Pz.clear();
  fgenie_P.clear();
  fgenie_status_code.clear();

  fnuPDG_truth.clear();
  fccnc_truth.clear();
  fmode_truth.clear();
  fhitnuc_truth.clear();
  fnu_vx_truth.clear();
  fnu_vy_truth.clear();
  fnu_vz_truth.clear();
  fnu_dcosx_truth.clear();
  fnu_dcosy_truth.clear();
  fnu_dcosz_truth.clear();
  flep_mom_truth.clear();
  flep_dcosx_truth.clear();
  flep_dcosy_truth.clear();
  flep_dcosz_truth.clear();


  fTrackId.clear();
  fMother.clear();
  fEng.clear();
  fEkin.clear();
  fMass.clear();
  fPdg.clear();
  fP.clear();
  fPx.clear();
  fPy.clear();
  fPz.clear();
  fND.clear();
  fstartX.clear();
  fstartY.clear();
  fstartZ.clear();
  fendX.clear();
  fendY.clear();
  fendZ.clear();


  fhit_MCparticle.clear();
  fhit_tpc.clear();
  fhit_channel.clear();
  fhit_time.clear();
  fhit_SADC.clear();
  fhit_wire.clear();
  fhit_charge.clear();
  fhit_plane.clear();
  fhit_width.clear();
  fhit_TOT.clear();
  fhit_trueX.clear();
  fhit_trueY.clear();
  fhit_trueZ.clear();
  fcolhit_tpc.clear();
  fcolhit_channel.clear();
  fcolhit_time.clear();
  fcolhit_SADC.clear();
  fcolhit_wire.clear();
  fcolhit_charge.clear();
  fcolhit_width.clear();
  fcolhit_trueX.clear();
  fcolhit_trueY.clear();
  fcolhit_trueZ.clear();
 
  
  fPFPpdg_pandora.clear();
  fPFPID_pandora.clear();
  fNDPrimary_pandora.clear();
  fPrimaryID_pandora.clear();
  fPFPparent_pandora.clear();
  
  fTrackID_pandora.clear();
  fTrackLength_pandora.clear();
  fTrackStartx_pandora.clear();
  fTrackStarty_pandora.clear();
  fTrackStartz_pandora.clear();
  fTrackEndx_pandora.clear();
  fTrackEndy_pandora.clear();
  fTrackEndz_pandora.clear();

  fVrtxX_pandora.clear();
  fVrtxY_pandora.clear();
  fVrtxZ_pandora.clear();
  fPrimaryVrtxX_pandora.clear();
  fPrimaryVrtxY_pandora.clear();
  fPrimaryVrtxZ_pandora.clear();
  fPrimaryVrtxPDG_pandora.clear();

  fClusterID_pandora.clear();
  fClusterView_pandora.clear();
  fClusterNHits_pandora.clear();
  fClusterIntegral_pandora.clear();
  fClusterIntegralAverage_pandora.clear();
  fClusterSummedADC_pandora.clear();
  fClusterSummedADCaverage_pandora.clear();
  fClusterStartAngle_pandora.clear();
  fClusterEndAngle_pandora.clear();
  fClusterStartOpeningAngle_pandora.clear();
  fClusterEndOpeningAngle_pandora.clear();

  fnTrkHits_pandora.clear();
  fTrkHit_trkID.clear();
  fTrkHit_MCparticle.clear();
  fTrkHit_PFParticle.clear();
  fTrkHit_tpc.clear();
  fTrkHit_channel.clear();
  fTrkHit_time.clear();
  fTrkHit_SADC.clear();
  fTrkHit_wire.clear();
  fTrkHit_charge.clear();
  fTrkHit_plane.clear();
  fTrkHit_trueX.clear();  
  fTrkHit_trueY.clear();  
  fTrkHit_trueZ.clear();  
  fTrkHit_clusterID.clear();
  
  
  fShowerID_pandora.clear();
  fShowerLength_pandora.clear();
  fShowerStartx_pandora.clear();
  fShowerStarty_pandora.clear();
  fShowerStartz_pandora.clear();

  fnShowerHits_pandora.clear();
  fShowerHit_ShowerID.clear();
  fShowerHit_MCparticle.clear();
  fShowerHit_PFParticle.clear();
  fShowerHit_tpc.clear();
  fShowerHit_channel.clear();
  fShowerHit_time.clear();
  fShowerHit_SADC.clear();
  fShowerHit_wire.clear();
  fShowerHit_charge.clear();
  fShowerHit_plane.clear();
  fShowerHit_trueX.clear();  
  fShowerHit_trueY.clear();  
  fShowerHit_trueZ.clear();  
  fShowerHit_clusterID.clear(); 

  fCalo_TrackID.clear();
  fCalo_TrackResidualRange.clear();
  fCalo_TrackdEdx.clear(); 
  
  art::ServiceHandle<geo::Geometry> geo;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///
///                                       ANALYZER
///
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void dune::EventAnaTree::analyze(art::Event const& e) 
{
  ResetVariables();
  
  fEvent   = e.id().event();
  fRun     = e.run();
  fSubRun  = e.subRun();

  bool isMC = !e.isRealData();
  auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(e);
  art::ServiceHandle<cheat::BackTrackerService> bt_serv; // MC cheat for mapping back to truth 

  if (isMC){

    if ( fSaveGenieInfo ) {

      art::Ptr<simb::MCTruth> mctruth;  

      // ----------------------- Genie:: neutrino truth info ---------------------------------------------
      art::ValidHandle<std::vector <simb::MCTruth> > mctruthListHandle = e.getValidHandle< std::vector<simb::MCTruth> >(fGenieLabel);
      std::vector< art::Ptr<simb::MCTruth> > mclist;

      if (mctruthListHandle.isValid()){
	art::fill_ptr_vector(mclist, mctruthListHandle);
      }

      if (!mclist.empty()){
	mctruth = mclist[0];
      }
      if (mctruth->NeutrinoSet() ){
	fgenie_no_primaries = mctruth->NParticles();
    
	for (size_t iPart = 0; iPart < fgenie_no_primaries; ++iPart){

	  const simb::MCParticle& part(mctruth->GetParticle(iPart));
	  fgenie_primaries_pdg.push_back( part.PdgCode() );
	  fgenie_Eng.push_back( part.E() );
	  fgenie_Px.push_back( part.Px() );
	  fgenie_Px.push_back( part.Py() );
	  fgenie_Px.push_back( part.Pz() );
	  fgenie_Px.push_back( part.P() );
	  fgenie_status_code.push_back( part.StatusCode() );
	  
	} // loop over genie particles 
      } // if (mctruth->NeutrinoSet()) 

      //Save neutrino interaction information 31.01.2022
      if (mclist.size() > 0) {
	int neutrino_i = 0;
	for (unsigned int iList = 0; iList < mclist.size(); ++iList){
	  if ( mclist[iList]->NeutrinoSet() ){
	    fnuPDG_truth.push_back( mclist[iList]->GetNeutrino().Nu().PdgCode() );
	    fccnc_truth.push_back( mclist[iList]->GetNeutrino().CCNC() );
	    fmode_truth.push_back( mclist[iList]->GetNeutrino().Mode() );
	    fhitnuc_truth.push_back( mclist[iList]->GetNeutrino().HitNuc() );
	    fnu_vx_truth.push_back( mclist[iList]->GetNeutrino().Nu().Vx() );
	    fnu_vy_truth.push_back( mclist[iList]->GetNeutrino().Nu().Vy() );
	    fnu_vz_truth.push_back( mclist[iList]->GetNeutrino().Nu().Vz() );

	    if ( mclist[iList]->GetNeutrino().Nu().P() ){
	      fnu_dcosx_truth.push_back( mclist[iList]->GetNeutrino().Nu().Px()/(mclist[iList]->GetNeutrino().Nu().P()) );
	      fnu_dcosy_truth.push_back( mclist[iList]->GetNeutrino().Nu().Py()/(mclist[iList]->GetNeutrino().Nu().P()) );
	      fnu_dcosz_truth.push_back( mclist[iList]->GetNeutrino().Nu().Pz()/(mclist[iList]->GetNeutrino().Nu().P()) );
		} //angular info for neutrino

	    flep_mom_truth.push_back( mclist[iList]->GetNeutrino().Lepton().P() );
	    if ( mclist[iList]->GetNeutrino().Lepton().P()){
	      flep_dcosx_truth.push_back(  mclist[iList]->GetNeutrino().Lepton().Px()/(mclist[iList]->GetNeutrino().Lepton().P()) );
	      flep_dcosy_truth.push_back(  mclist[iList]->GetNeutrino().Lepton().Py()/(mclist[iList]->GetNeutrino().Lepton().P()) );
	      flep_dcosz_truth.push_back(  mclist[iList]->GetNeutrino().Lepton().Pz()/(mclist[iList]->GetNeutrino().Lepton().P()) );
	    } //angular info for lepton 
	 
	    neutrino_i++;
	  }// mclist is NeutrinoSet()
	}  // loop over MC records 
      } //at least one MC record 

    } // if SaveGenieInfo

  
    // ----------------------- GEANT :: truth list of MC particles -----------------------------------

    art::ValidHandle<std::vector <simb::MCParticle> > mcParticles = e.getValidHandle<std::vector <simb::MCParticle> >(fTruthLabel);
    if (mcParticles.isValid()){
      fnGeantParticles = mcParticles->size();
      
      //Get MC particle list from G4 
      for (unsigned int t = 0; t < mcParticles->size(); ++t){ 
	const simb::MCParticle trueParticle = mcParticles->at(t);

	fTrackId.push_back( trueParticle.TrackId());
	fMother.push_back( trueParticle.Mother());
	fEng.push_back( trueParticle.E() );
	fPdg.push_back( trueParticle.PdgCode());
	fEkin.push_back( trueParticle.E() - trueParticle.Mass() );
	fMass.push_back( trueParticle.Mass() );
	fP.push_back( trueParticle.P()) ;
	fPx.push_back( trueParticle.Px());
	fPy.push_back( trueParticle.Py());
	fPz.push_back( trueParticle.Pz());
	fND.push_back( trueParticle.NumberDaughters());
	fstartX.push_back( trueParticle.Vx());
	fstartY.push_back( trueParticle.Vy());
	fstartZ.push_back( trueParticle.Vz());
	fendX.push_back( trueParticle.EndPosition()[0]);
	fendY.push_back( trueParticle.EndPosition()[1]);
	fendZ.push_back( trueParticle.EndPosition()[2]);

      }
    } //if MC particles are valid 

    //----------------------- LarSoft ::  hit list -------------------------- 
    //art::ValidHandle<std::vector <recob::Hit> > recoHits = e.getValidHandle<std::vector <recob::Hit> >(fHitLabel);
    std::vector< art::Ptr<recob::Hit> >recoHits; 
    auto RecoHitHandle = e.getValidHandle<std::vector<recob::Hit>>(fHitLabel); // art handle for the reco hits 
    if(RecoHitHandle){ art::fill_ptr_vector(recoHits, RecoHitHandle); }

    //Run over the reco hits 
    if ( RecoHitHandle ){
      fnHits = recoHits.size();
      fnColHits = 0;


      //  Build mapping between Hits and MCParticles, starting from ART event record 
      
      //output mapping between MC particles to reco hits <-> 
      lar_pandora::MCParticlesToHits particlesToHits;
      lar_pandora::HitsToMCParticles  hitsToParticles;
      
      lar_pandora::LArPandoraHelper::BuildMCParticleHitMaps(e, fTruthLabel, recoHits,
							    particlesToHits, hitsToParticles,  
							    lar_pandora::LArPandoraHelper::kAddDaughters); 
      
      
      for (unsigned int hit = 0; hit < recoHits.size(); ++hit){
	//const recob::Hit thisHit = recoHits->at(hit);
	const art::Ptr<recob::Hit> thisHit = recoHits.at(hit);

	lar_pandora::HitsToMCParticles::const_iterator vIter = hitsToParticles.find( thisHit ); 
	if (hitsToParticles.end() != vIter) { 
	  const art::Ptr<simb::MCParticle> MCparticle = vIter->second; 
	  fhit_MCparticle.push_back(MCparticle->PdgCode());
	}
	else fhit_MCparticle.push_back(-1); 

	fhit_tpc.push_back(thisHit->WireID().TPC);
	fhit_channel.push_back(thisHit->Channel());
	fhit_time.push_back(thisHit->PeakTime());
	fhit_SADC.push_back(thisHit->SummedADC());
	fhit_wire.push_back(thisHit->WireID().Wire);
	fhit_charge.push_back(thisHit->Integral());
	fhit_plane.push_back(thisHit->View());
	fhit_width.push_back( std::abs( thisHit->PeakTimePlusRMS() - thisHit->PeakTimeMinusRMS()) ); 
	fhit_TOT.push_back(thisHit->EndTick() - thisHit->StartTick()); // time over threshold in ticks 
	// 1 tick = 500 ns (2 MHz sampling frequency)
	
	// extracting physical position of the hit in detector geometry  {
	std::vector<const sim::IDE*> ides;
	try{ ides = bt_serv->HitToSimIDEs_Ps(clockData,thisHit); }
	catch(...){}
	if (ides.size() > 0){
	  //getting physical coordinates of hit from ides 
	  std::vector<double> xyz = bt_serv->SimIDEsToXYZ(ides);
	  fhit_trueX.push_back(xyz[0]);
	  fhit_trueY.push_back(xyz[1]);
	  fhit_trueZ.push_back(xyz[2]);
	}


	if (thisHit->View() == 2){  //collection plane hits only
	  fnColHits +=1;   
	  fcolhit_tpc.push_back(thisHit->WireID().TPC);
	  fcolhit_channel.push_back(thisHit->Channel());
	  fcolhit_time.push_back(thisHit->PeakTime());
	  fcolhit_SADC.push_back(thisHit->SummedADC());
	  fcolhit_wire.push_back(thisHit->WireID().Wire);
	  fcolhit_charge.push_back(thisHit->Integral());
	  fcolhit_width.push_back( std::abs( thisHit->PeakTimePlusRMS() - thisHit->PeakTimeMinusRMS()) );
	  if (isMC){
	    std::vector<const sim::IDE*> ColIDES;
	    try{
	      ColIDES = bt_serv->HitToSimIDEs_Ps(clockData,thisHit);
	    }
	    catch(...){}
	    if (ColIDES.size() > 0){
	      std::vector<double> colxyz = bt_serv->SimIDEsToXYZ(ColIDES);
	      fcolhit_trueX.push_back(colxyz[0]);
	      fcolhit_trueY.push_back(colxyz[1]);
	      fcolhit_trueZ.push_back(colxyz[2]);
	    } // hit XYZ : ides > 0 
	  } // hit XYZ: is MC
	} // run over coll hits 
      } // run over all hits
    } //if reco hits valid 

    
  }//if is MC
  
  //------------------------Pandora :: Collect PFParticles ---------------------------
  
  art::Handle< std::vector<recob::PFParticle> > pfparticleHandle; //art handle for the PFParticles
  std::vector< art::Ptr<recob::PFParticle> > PFParticleVect; // initializiing my PFP vector 

  // fill vector with PFParticles 
  if (e.getByLabel(fParticleLabel, pfparticleHandle)) art::fill_ptr_vector(PFParticleVect, pfparticleHandle);

  // #  PFParticles for current event 
  fnPFParticles_pandora = PFParticleVect.size();
    
  //run over the PFPs
  if (fnPFParticles_pandora != 0){

    for (const art::Ptr<recob::PFParticle> &pfp : PFParticleVect){

      fPFPpdg_pandora.push_back( pfp->PdgCode() );
      fPFPID_pandora.push_back( pfp->Self() );
      fPFPparent_pandora.push_back( pfp->Parent() );
      
      if (!( pfp->IsPrimary()) ) continue; 
      
      fPrimaryID_pandora.push_back( pfp->Self() );
      fNDPrimary_pandora.push_back( pfp->NumDaughters() ); 
      fnPrimaries_pandora ++;
    } // running over the PFPs
  } // if number of PFParticles is non-zero
  
  // ------------------------------- Pandora:: Vertices for this event --------------------------------------------

  //vector holding reco vertices for this event
  lar_pandora::VertexVector vertexVector; 
  //PFParticles associated with each vector 
  lar_pandora::PFParticlesToVertices particlesToVertices;

  //Fill the vectors
  lar_pandora::LArPandoraHelper::CollectVertices(e, fPandoraVertexLabel, vertexVector, particlesToVertices); 

  //get the number of vertices/ primary vertices reconstructed 
  //run over the vertices 
  fnVertices_pandora = vertexVector.size(); 
  for (unsigned int n = 0; n < PFParticleVect.size(); ++n){
    const art::Ptr<recob::PFParticle> particle = PFParticleVect.at(n);
    if ( particle->IsPrimary() ){
      fnPrimaryVertices_pandora++; 
      lar_pandora::PFParticlesToVertices::const_iterator vIter = particlesToVertices.find(particle);
      if (particlesToVertices.end() != vIter) {
	const lar_pandora::VertexVector &vertexVector = vIter->second; 
	if (!vertexVector.empty()){
	  if (vertexVector.size() == 1) {
	    const art::Ptr<recob::Vertex> vertex = *(vertexVector.begin()); 
	    double xyz[3] = {0.0, 0.0, 0.0};
	    vertex->XYZ(xyz);
	    fPrimaryVrtxX_pandora.push_back(xyz[0]);
	    fPrimaryVrtxY_pandora.push_back(xyz[1]);
	    fPrimaryVrtxZ_pandora.push_back(xyz[2]);
	    fPrimaryVrtxPDG_pandora.push_back( particle->PdgCode() ); 
	  }
	}
      }
    } // if particle is primary 
  } //run over PFParticle list 
  if (!vertexVector.empty()){
    for (unsigned int v = 0; v < vertexVector.size(); ++v){ 
      const art::Ptr<recob::Vertex> currentVertex = vertexVector.at(v); 
      double xyz[3] = {0.0, 0.0, 0.0};
      currentVertex->XYZ(xyz); 
      fVrtxX_pandora.push_back(xyz[0]);
      fVrtxY_pandora.push_back(xyz[1]);
      fVrtxZ_pandora.push_back(xyz[2]);
    }
  } //run over the vertices

  // ------------------------------ Pandora :: Clusters + associations  --------------------------------------------------
 
  // 20/02/2023
  lar_pandora::ClusterVector clusterVector; 
  lar_pandora::ClustersToHits clustersToHits; //map between clusters and pandora reco hits 
  lar_pandora::LArPandoraHelper::CollectClusters(e, fSpacePointLabel, clusterVector, clustersToHits); 
  fnClusters_pandora = clusterVector.size(); 

  for (unsigned int c = 0; c < fnClusters_pandora; c++){
    const art::Ptr<recob::Cluster>  thisCluster = clusterVector.at(c);

    fClusterID_pandora.push_back( thisCluster->ID() );
    fClusterView_pandora.push_back( thisCluster->View() );
    fClusterNHits_pandora.push_back( thisCluster->NHits() );
    fClusterIntegral_pandora.push_back( thisCluster->Integral() ); //total charge of cluster from hit shape
    fClusterIntegralAverage_pandora.push_back( thisCluster->IntegralAverage() ); //average charge of cluster hits 
    fClusterSummedADC_pandora.push_back( thisCluster->SummedADC() );
    fClusterSummedADCaverage_pandora.push_back( thisCluster->SummedADCaverage() );
    fClusterStartAngle_pandora.push_back( thisCluster->StartAngle() );
    fClusterEndAngle_pandora.push_back( thisCluster->EndAngle() );
    fClusterStartOpeningAngle_pandora.push_back( thisCluster->StartOpeningAngle() );
    fClusterEndOpeningAngle_pandora.push_back( thisCluster->EndOpeningAngle() );


  }

  // ------------------------------ Pandora :: Reco Track + Hits  --------------------------------------------------

  art::ValidHandle<std::vector<recob::PFParticle>> recoParticles = e.getValidHandle<std::vector<recob::PFParticle>>(fParticleLabel);
  
  //Get track list so they can be associted to hits and calorimetric info
  art::ValidHandle< std::vector<recob::Track> > recoTracks = e.getValidHandle< std::vector<recob::Track> >(fTrackLabel);

  
  if (recoParticles.isValid()){

    //Find vector of tracks for each associated PFParticle 
    const art::FindManyP<recob::Track> findTracks(recoParticles, e, fTrackLabel);

    //Find reco hits for this track
    const art::FindManyP<recob::Hit> findHits(recoTracks, e, fTrackLabel);


    //Get calorimetric information for this track
    const art::FindManyP<anab::Calorimetry> calorimetryAssoc(recoTracks, e, fCaloLabel); //02.12.21

    //Looping over individual particles 
    for (unsigned int p = 0 ; p < recoParticles->size(); ++p){    
  
      //current particle
      const recob::PFParticle particle = recoParticles->at(p); 

      //Get the association between current particle and track(s)
      const std::vector<art::Ptr<recob::Track>> pfpTracks = findTracks.at(p);

      //if 'track-like' particle -> loop over tracks 
      if (pfpTracks.size() != 0){
	for ( unsigned int k = 0; k < pfpTracks.size(); ++k){

	  //current track
	  const art::Ptr<recob::Track> thisTrack = pfpTracks.at(k);
     
	  int ntraj = thisTrack->NumberTrajectoryPoints(); //# traj points in track 
	  if ( ntraj > 0){
	    
	    fTrackStartx_pandora.push_back( thisTrack->Vertex().X() );
	    fTrackStarty_pandora.push_back( thisTrack->Vertex().Y() );
	    fTrackStartz_pandora.push_back( thisTrack->Vertex().Z() );
	    fTrackEndx_pandora.push_back( thisTrack->End().X() );
	    fTrackEndy_pandora.push_back( thisTrack->End().Y() );
	    fTrackEndz_pandora.push_back( thisTrack->End().Z() );
	  }  
	   
	  fTrackLength_pandora.push_back( thisTrack->Length() );
	  fTrackID_pandora.push_back( thisTrack->ID() );
	  fnRecoTracks_pandora +=1;


	  int track_ID = thisTrack->ID();


	  //Get calorimetric information for this track 02.12 
	  
	  const std::vector< art::Ptr<anab::Calorimetry> > trackCalo = calorimetryAssoc.at(track_ID );
	  if (trackCalo.size() !=0){
	    for (unsigned int ll = 0; ll < trackCalo.size(); ++ll){
	      const art::Ptr<anab::Calorimetry> cal = trackCalo.at(ll);

	      //looking only at the collection plane: 
	      int planenum = cal->PlaneID().Plane;
	      if (planenum !=2) continue; 

	      fCalo_TrackID.push_back( track_ID );
	      fCalo_TrackResidualRange.push_back( cal->ResidualRange() );
	      fCalo_TrackdEdx.push_back( cal->dEdx() );

	    } //loop over trackCalo
	  }// Calorimetry
		  
	
	  //Get hits associated to this track
	  const std::vector< art::Ptr<recob::Hit> > trackHits = findHits.at(track_ID);
	  fnTrkHits_pandora.push_back( trackHits.size());


	  //build recoTrkHit --> MCParticle map                                                                   
	  lar_pandora::MCParticlesToHits particlesToTrkHits;
	  lar_pandora::HitsToMCParticles TrkhitsToParticles;

	  lar_pandora::LArPandoraHelper::BuildMCParticleHitMaps(e, fTruthLabel, trackHits,
								particlesToTrkHits, TrkhitsToParticles,
								lar_pandora::LArPandoraHelper::kIgnoreDaughters); //kIgnoreDaughters


	  //build recoTrkHit -->PFParticle map
	  bool useClusters = false;
	  lar_pandora::PFParticlesToHits pfpsToTrkHits;
	  lar_pandora::HitsToPFParticles TrkhitsToPfps;
	  lar_pandora::LArPandoraHelper::BuildPFParticleHitMaps(e, fParticleLabel, 
								pfpsToTrkHits, TrkhitsToPfps, 
								lar_pandora::LArPandoraHelper::kIgnoreDaughters, useClusters); //kIgnoreDaughters

	  //running over track reco hits 
	  if (trackHits.size() != 0){
	    for (unsigned int jj = 0; jj < trackHits.size(); ++jj){
	      const art::Ptr<recob::Hit> thisHit = trackHits.at(jj);

	      //find cluster associated with this hit
	      bool cluster_assoc = false;
	      for (auto const& x : clustersToHits) 
		{

		  lar_pandora::HitVector clusterhits = x.second; 
		  if (std::find(clusterhits.begin(), clusterhits.end(), thisHit) != clusterhits.end()){
		    const art::Ptr<recob::Cluster> cluster = x.first; 
		    fTrkHit_clusterID.push_back( cluster->ID() );
		    cluster_assoc = true; 
		  }
		}
	      if (!cluster_assoc) fTrkHit_clusterID.push_back(-1);


	      //find MC particle associated with this hit
	      lar_pandora::HitsToMCParticles::const_iterator vIter = TrkhitsToParticles.find( thisHit );
              if (TrkhitsToParticles.end() != vIter) {
                const art::Ptr<simb::MCParticle> MCparticle = vIter->second;
                fTrkHit_MCparticle.push_back(MCparticle->PdgCode());
              }
              else fTrkHit_MCparticle.push_back(-1);

	      //find track associated with this hit
	      lar_pandora::HitsToPFParticles::const_iterator pIter = TrkhitsToPfps.find( thisHit );
              if (TrkhitsToPfps.end() != pIter) {
                const art::Ptr<recob::PFParticle> PFparticle = pIter->second;
                fTrkHit_PFParticle.push_back(PFparticle->PdgCode());
              }
              else fTrkHit_PFParticle.push_back(-1);

	      fTrkHit_trkID.push_back( track_ID );
	      fTrkHit_tpc.push_back( thisHit->WireID().TPC );
	      fTrkHit_channel.push_back( thisHit->Channel() );
	      fTrkHit_time.push_back( thisHit->PeakTime() );
	      fTrkHit_SADC.push_back( thisHit->SummedADC() );
	      fTrkHit_wire.push_back( thisHit->WireID().Wire );
	      fTrkHit_charge.push_back( thisHit->Integral() );
	      fTrkHit_plane.push_back( thisHit->View() );
	      if (isMC){
		std::vector<const sim::IDE*> HitIDES;
		try{
		  HitIDES = bt_serv->HitToSimIDEs_Ps(clockData,thisHit);
		}
		catch(...){}
		if (HitIDES.size() > 0){
		  std::vector<double> reco_track_xyz = bt_serv->SimIDEsToXYZ(HitIDES);
		  fTrkHit_trueX.push_back(reco_track_xyz[0]);
		  fTrkHit_trueY.push_back(reco_track_xyz[1]);
		  fTrkHit_trueZ.push_back(reco_track_xyz[2]);
		  
		}// hit XYZ : ides > 0 
	      } // hit XYZ: is MC
	      
	    } // trk hit for-loop
	  } // if hit vect > 0
	} // run over pfp tracks 
      } // if pfp.size !=0
    } //loop over reco particles
  } // if reco particles is Valid


  // -------------------------------- Pandora :: Reco Shower + Hits -----------------------------------------------

  //Shower list
  art::ValidHandle< std::vector<recob::Shower> > recoShowers = e.getValidHandle<std::vector<recob::Shower>>(fShowerLabel);

  if (recoParticles.isValid()){
    //Find vector of showers for each associated PFParticle
    const art::FindManyP<recob::Shower> findShowers( recoParticles, e, fShowerLabel);

    //Find reco hits for this shower 
    const art::FindManyP<recob::Hit> findShowerHits( recoShowers, e, fShowerLabel); 
    
    //Looping over individual PFParticles
    for (unsigned int p = 0; p < recoParticles->size(); ++p){

      //curent particle
      const recob::PFParticle particle = recoParticles->at(p);

      //Get association between current particle and shower(s)
      const std::vector<art::Ptr<recob::Shower>> pfpShowers = findShowers.at(p);

      //if 'shower-like' particle --> loop over reconstructed showers :
      if (pfpShowers.size() !=0){
	for (unsigned int k = 0; k < pfpShowers.size(); ++k){

	  const art::Ptr<recob::Shower> thisShower = pfpShowers.at(k);
	  fnRecoShowers_pandora +=1;
	  fShowerID_pandora.push_back( thisShower->ID() );
	  fShowerLength_pandora.push_back( thisShower->Length() );
	  fShowerStartx_pandora.push_back( thisShower->ShowerStart().X() );
	  fShowerStarty_pandora.push_back( thisShower->ShowerStart().Y() );
	  fShowerStartz_pandora.push_back( thisShower->ShowerStart().Z() );

	  //Get hits associated to this shower
	  const std::vector< art::Ptr<recob::Hit> > showerHits = findShowerHits.at(k);	  
	  if (showerHits.size() !=0){

	    int shower_ID = thisShower->ID();
	    fnShowerHits_pandora.push_back( showerHits.size() ); 

	    //build recoShwrtHit --> MCParticle map 
	    lar_pandora::MCParticlesToHits particlesToShwrHits;
	    lar_pandora::HitsToMCParticles ShwrhitsToParticles;

	    lar_pandora::LArPandoraHelper::BuildMCParticleHitMaps(e, fTruthLabel, showerHits,
								  particlesToShwrHits, ShwrhitsToParticles,
								  lar_pandora::LArPandoraHelper::kIgnoreDaughters); // kIgnoreDaughters

	    //build recoShwrHit -->PFParticle map                               
	    bool useClusters = false;
	    lar_pandora::PFParticlesToHits pfpsToShwrHits;
	    lar_pandora::HitsToPFParticles ShwrhitsToPfps;
	    lar_pandora::LArPandoraHelper::BuildPFParticleHitMaps(e, fParticleLabel,
								  pfpsToShwrHits, ShwrhitsToPfps,
								  lar_pandora::LArPandoraHelper::kIgnoreDaughters, useClusters); //kIgnoreDaughters


	    //Run over the shower hits 
	    for (unsigned int kk = 0; kk < showerHits.size(); ++kk){
	      const art::Ptr<recob::Hit> thisShwHit = showerHits.at(kk);

	      //find cluster associated with this hit
	      bool cluster_assoc = false;
	      for (auto const& x : clustersToHits) 
		{

		  lar_pandora::HitVector clusterhits = x.second; 
		  if (std::find(clusterhits.begin(), clusterhits.end(), thisShwHit) != clusterhits.end()){
		    const art::Ptr<recob::Cluster> cluster = x.first; 
		    fShowerHit_clusterID.push_back( cluster->ID() );
		    cluster_assoc = true; 
		  }
		}
	      if (!cluster_assoc) fShowerHit_clusterID.push_back(-1);


	      //push back MC pdg
	      lar_pandora::HitsToMCParticles::const_iterator vIter = ShwrhitsToParticles.find( thisShwHit );
	      if (ShwrhitsToParticles.end() != vIter) {
		const art::Ptr<simb::MCParticle> MCparticle = vIter->second;
		fShowerHit_MCparticle.push_back(MCparticle->PdgCode());
	      }
	      else fShowerHit_MCparticle.push_back(-1);

	      //push back PFP pdg
	      lar_pandora::HitsToPFParticles::const_iterator pIter = ShwrhitsToPfps.find( thisShwHit );
              if (ShwrhitsToPfps.end() != pIter) {
                const art::Ptr<recob::PFParticle> PFparticle = pIter->second;
                fShowerHit_PFParticle.push_back(PFparticle->PdgCode());
              }
              else fShowerHit_PFParticle.push_back(-1);


	      fShowerHit_ShowerID.push_back(shower_ID);
	      fShowerHit_tpc.push_back( thisShwHit->WireID().TPC );
	      fShowerHit_channel.push_back( thisShwHit->Channel() );
	      fShowerHit_time.push_back( thisShwHit->PeakTime() );
	      fShowerHit_SADC.push_back( thisShwHit->SummedADC() );
	      fShowerHit_wire.push_back( thisShwHit->WireID().Wire );
	      fShowerHit_charge.push_back( thisShwHit->Integral() );
	      fShowerHit_plane.push_back( thisShwHit->View() );
	  
	      if (isMC){
		std::vector<const sim::IDE*> HitIDES;
		try{
		  HitIDES = bt_serv->HitToSimIDEs_Ps(clockData,thisShwHit);
		}
		catch(...){}
		if (HitIDES.size() > 0){
		  std::vector<double> reco_shower_xyz = bt_serv->SimIDEsToXYZ(HitIDES);
		  fShowerHit_trueX.push_back(reco_shower_xyz[0]);
		  fShowerHit_trueY.push_back(reco_shower_xyz[1]);
		  fShowerHit_trueZ.push_back(reco_shower_xyz[2]);
		  
		} // hitIDES > 0
	      } // is MC 
	    }//loop over shower hits 
	  } // if size of shower hits != 0 	
	} //loop over showers for current pfparticle 
      }
    }//loop over PFParticles
  } // if reco particles is valid 

  fTree->Fill();
  
} // analyze function end 

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///
///                                       BEGIN JOB
///
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void dune::EventAnaTree::beginJob()
{

  art::ServiceHandle< art::TFileService> tfs;
  fTree = tfs->make<TTree>("tree", "Analyser Output Tree"); 

  fTree->Branch("event",&fEvent,"event/i");
  fTree->Branch("run",&fRun,"run/i");
  fTree->Branch("subrun",&fSubRun,"subrun/i");

  std::cout << std::string(80, '=') << std::endl;
  std::cout << "Truth information" << std::endl;
  std::cout << std::string(80, '=') << std::endl;


  if (fSaveGenieInfo){

    fTree->Branch("genie_no_primaries",&fgenie_no_primaries,"genie_no_primaries/i");
    fTree->Branch("genie_primaries_pdg",&fgenie_primaries_pdg);
    fTree->Branch("genie_Eng",&fgenie_Eng);
    fTree->Branch("genie_Px",&fgenie_Px);
    fTree->Branch("genie_Py",&fgenie_Py);
    fTree->Branch("genie_Pz",&fgenie_Pz);
    fTree->Branch("genie_P",&fgenie_P);
    fTree->Branch("genie_status_code",&fgenie_status_code); 

    fTree->Branch("nuPDG_truth",&fnuPDG_truth);
    fTree->Branch("ccnc_truth",&fccnc_truth);
    fTree->Branch("mode_truth",&fmode_truth);
    fTree->Branch("hitnuc_truth",&fhitnuc_truth);
    fTree->Branch("nu_vx_truth",&fnu_vx_truth);
    fTree->Branch("nu_vy_truth",&fnu_vy_truth);
    fTree->Branch("nu_vz_truth",&fnu_vz_truth);
    fTree->Branch("nu_dcosx_truth",&fnu_dcosx_truth);
    fTree->Branch("nu_dcosy_truth",&fnu_dcosy_truth);
    fTree->Branch("nu_dcosz_truth",&fnu_dcosz_truth);
    fTree->Branch("lep_mom_truth",&flep_mom_truth);
    fTree->Branch("lep_dcosx_truth",&flep_dcosx_truth);
    fTree->Branch("lep_dcosy_truth",&flep_dcosy_truth);
    fTree->Branch("lep_dcosz_truth",&flep_dcosz_truth);


  }

  fTree->Branch("nGeantParticles",&fnGeantParticles,"nGeantParticles/i"); 
  fTree->Branch("nPrimaries",&fnPrimaries,"nPrimaries/i");
  fTree->Branch("TrackId",&fTrackId);
  fTree->Branch("Mother",&fMother);
  fTree->Branch("Pdg",&fPdg);
  fTree->Branch("Eng",&fEng);
  fTree->Branch("Ekin",&fEkin);
  fTree->Branch("Mass",&fMass);
  fTree->Branch("P",&fP);
  fTree->Branch("Px",&fPx);
  fTree->Branch("Py",&fPy);
  fTree->Branch("Pz",&fPz);
  fTree->Branch("ND",&fND);
  fTree->Branch("startX",&fstartX);
  fTree->Branch("startY",&fstartY);
  fTree->Branch("startZ",&fstartZ);
  fTree->Branch("endX",&fendX);
  fTree->Branch("endY",&fendY);
  fTree->Branch("endZ",&fendZ);

  fTree->Branch("nHits",&fnHits,"nHits/i");
  fTree->Branch("nColHits",&fnColHits,"nColHits/i");
  fTree->Branch("hit_MCparticle",&fhit_MCparticle);
  fTree->Branch("hit_tpc",&fhit_tpc);
  fTree->Branch("hit_channel",&fhit_channel);
  fTree->Branch("hit_time",&fhit_time);
  fTree->Branch("hit_SADC",&fhit_SADC);
  fTree->Branch("hit_wire",&fhit_wire);
  fTree->Branch("hit_charge",&fhit_charge);
  fTree->Branch("hit_plane",&fhit_plane);
  fTree->Branch("hit_width",&fhit_width);
  fTree->Branch("hit_TOT",&fhit_TOT);
  fTree->Branch("hit_trueX",&fhit_trueX);
  fTree->Branch("hit_trueY",&fhit_trueY);
  fTree->Branch("hit_trueZ",&fhit_trueZ);
  fTree->Branch("colhit_tpc",&fcolhit_tpc);
  fTree->Branch("colhit_channel",&fcolhit_channel);
  fTree->Branch("colhit_time",&fcolhit_time);
  fTree->Branch("colhit_SADC",&fcolhit_SADC);
  fTree->Branch("colhit_wire",&fcolhit_wire);
  fTree->Branch("colhit_charge",&fcolhit_charge);
  fTree->Branch("colhit_width",&fcolhit_width);
  fTree->Branch("colhit_trueX",&fcolhit_trueX);
  fTree->Branch("colhit_trueY",&fcolhit_trueY);
  fTree->Branch("colhit_trueZ",&fcolhit_trueZ);
  
  
  fTree->Branch("nPFParticles_pandora",&fnPFParticles_pandora,"nPFParticles/i");
  fTree->Branch("nPrimaries_pandora",&fnPrimaries_pandora,"nPrimaries_pandora/i");
  fTree->Branch("NDPrimary_pandora",&fNDPrimary_pandora);
  fTree->Branch("PrimaryID_pandora",&fPrimaryID_pandora);
  fTree->Branch("PFPpdg_pandora",&fPFPpdg_pandora);
  fTree->Branch("PFPparent_pandora",&fPFPparent_pandora);
  fTree->Branch("PFPID_pandora",&fPFPID_pandora);

  fTree->Branch("nVertices_pandora",&fnVertices_pandora,"nVertices_pandora/i");
  fTree->Branch("nPrimaryVertices_pandora",&fnPrimaryVertices_pandora,"nPrimaryVertices_pandora/i");
  fTree->Branch("VrtxX_pandora",&fVrtxX_pandora); 
  fTree->Branch("VrtxY_pandora",&fVrtxY_pandora); 
  fTree->Branch("VrtxZ_pandora",&fVrtxZ_pandora); 
  fTree->Branch("PrimaryVrtxX_pandora",&fPrimaryVrtxX_pandora); 
  fTree->Branch("PrimaryVrtxY_pandora",&fPrimaryVrtxY_pandora); 
  fTree->Branch("PrimaryVrtxZ_pandora",&fPrimaryVrtxZ_pandora); 
  fTree->Branch("PrimaryVrtxPDG_pandora",&fPrimaryVrtxPDG_pandora); 


  fTree->Branch("nClusters_pandora",&fnClusters_pandora,"nClusters_pandora/i");
  fTree->Branch("ClusterID_pandora",&fClusterID_pandora);
  fTree->Branch("ClusterView_pandora",&fClusterView_pandora);
  fTree->Branch("ClusterNHits_pandora",&fClusterNHits_pandora);
  fTree->Branch("ClusterIntegral_pandora",&fClusterIntegral_pandora);
  fTree->Branch("ClusterIntegralAverage_pandora",&fClusterIntegralAverage_pandora);
  fTree->Branch("ClusterSummedADC_pandora",&fClusterSummedADC_pandora);
  fTree->Branch("ClusterSummedADCaverage_pandora",&fClusterSummedADCaverage_pandora);
  fTree->Branch("ClusterStartAngle_pandora",&fClusterStartAngle_pandora);
  fTree->Branch("ClusterEndAngle_pandora",&fClusterEndAngle_pandora);
  fTree->Branch("ClusterStartOpeningAngle_pandora",&fClusterStartOpeningAngle_pandora);
  fTree->Branch("ClusterEndOpeningAngle_pandora",&fClusterEndOpeningAngle_pandora);



  
  fTree->Branch("nRecoTracks_pandora",&fnRecoTracks_pandora,"nRecoTracks_pandora/i");
  fTree->Branch("TrackID_pandora",&fTrackID_pandora);
  fTree->Branch("TrackLength_pandora",&fTrackLength_pandora); 
  fTree->Branch("TrackStartx_pandora",&fTrackStartx_pandora);
  fTree->Branch("TrackStarty_pandora",&fTrackStarty_pandora);
  fTree->Branch("TrackStartz_pandora",&fTrackStartz_pandora);
  fTree->Branch("TrackEndx_pandora",&fTrackEndx_pandora);
  fTree->Branch("TrackEndy_pandora",&fTrackEndy_pandora);
  fTree->Branch("TrackEndz_pandora",&fTrackEndz_pandora);

  fTree->Branch("nTrkHits_pandora",&fnTrkHits_pandora);
  fTree->Branch("TrkHit_trkID",&fTrkHit_trkID);
  fTree->Branch("TrkHit_MCparticle",&fTrkHit_MCparticle);
  fTree->Branch("TrkHit_PFParticle",&fTrkHit_PFParticle);
  fTree->Branch("TrkHit_tpc",&fTrkHit_tpc);
  fTree->Branch("TrkHit_channel",&fTrkHit_channel);
  fTree->Branch("TrkHit_time",&fTrkHit_time);
  fTree->Branch("TrkHit_SADC",&fTrkHit_SADC);
  fTree->Branch("TrkHit_wire",&fTrkHit_wire);
  fTree->Branch("TrkHit_charge",&fTrkHit_charge);
  fTree->Branch("TrkHit_plane",&fTrkHit_plane);
  fTree->Branch("TrkHit_trueX",&fTrkHit_trueX);
  fTree->Branch("TrkHit_trueY",&fTrkHit_trueY);
  fTree->Branch("TrkHit_trueZ",&fTrkHit_trueZ);
  fTree->Branch("TrkHit_clusterID",&fTrkHit_clusterID);
  
  
  fTree->Branch("nRecoShowers_pandora",&fnRecoShowers_pandora,"nRecoShowers_pandora/i");
  fTree->Branch("ShowerID_pandora",&fShowerID_pandora);
  fTree->Branch("ShowerLength_pandora",&fShowerLength_pandora);
  fTree->Branch("ShowerStartx_pandora",&fShowerStartx_pandora);
  fTree->Branch("ShowerStarty_pandora",&fShowerStarty_pandora);
  fTree->Branch("ShowerStartz_pandora",&fShowerStartz_pandora);


  fTree->Branch("nShowerHits_pandora",&fnShowerHits_pandora);
  fTree->Branch("ShowerHit_ShowerID",&fShowerHit_ShowerID);
  fTree->Branch("ShowerHit_MCparticle",&fShowerHit_MCparticle);
  fTree->Branch("ShowerHit_PFParticle",&fShowerHit_PFParticle);
  fTree->Branch("ShowerHit_tpc",&fShowerHit_tpc);
  fTree->Branch("ShowerHit_channel",&fShowerHit_channel);
  fTree->Branch("ShowerHit_time",&fShowerHit_time);
  fTree->Branch("ShowerHit_SADC",&fShowerHit_SADC);
  fTree->Branch("ShowerHit_wire",&fShowerHit_wire);
  fTree->Branch("ShowerHit_charge",&fShowerHit_charge);
  fTree->Branch("ShowerHit_plane",&fShowerHit_plane);
  fTree->Branch("ShowerHit_trueX",&fShowerHit_trueX);
  fTree->Branch("ShowerHit_trueY",&fShowerHit_trueY);
  fTree->Branch("ShowerHit_trueZ",&fShowerHit_trueZ);
  fTree->Branch("ShowerHit_clusterID",&fShowerHit_clusterID);
  

  fTree->Branch("Calo_TrackID",&fCalo_TrackID);
  fTree->Branch("Calo_TrackResidualRange",&fCalo_TrackResidualRange);
  fTree->Branch("Calo_TrackdEdx",&fCalo_TrackdEdx);
  
  
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///
///                                       END JOB
///
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void dune::EventAnaTree::endJob()
{
    
}



DEFINE_ART_MODULE(dune::EventAnaTree)

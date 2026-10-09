/* Copyright(c) 1998-2026, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice */

#ifndef AliAnalysisTaskSEFMClusterTreeESD_H
#define AliAnalysisTaskSEFMClusterTreeESD_H

//////////////////////////////////////////////////////////////////////////////
//
// AliAnalysisTaskSEFMClusterTreeESD
//
// ESD + AliESDfriends version of the FM input dumper (TPC clusters only).
//
// Per event it stores
//   * TPC clusters (AliTPCclusterMI) attached to TPC seeds stored in the
//     ESD friends: global x,y,z, sector/ROC, pad, time bin, Q, Qmax,
//     cluster widths, applied distortion corrections and the 3 MC labels
//     of each cluster  -> per-cluster truth for track finding / PID
//   * reconstructed ESD tracks
//   * MC particles (AliMCEvent) that are referenced by clusters/tracks,
//     their full ancestry, and physical primaries in acceptance
//   * optional MC track references in the TPC (true crossing points)
//
// Friend / cluster access follows TPC/TPCcalib/AliTPCAnalysisTaskcalib:
//   AliESDEvent::FindFriend() -> AliESDtrack::GetFriendTrack()
//   -> AliESDfriendTrack::GetTPCseed() -> AliTPCseed::GetClusterPointer(row)
// and skips rows flagged with 0x8000 in the cluster index as done in
// AliTPCseed / AliTPCcalibAlign.
//
// A cluster attached to more than one seed is written once; hit_trk is the
// first track, hit_nTrk counts the attaching tracks and hit_trk2 stores the
// second one.
//
// Requirements in the train / run macro:
//   AliESDInputHandler::SetReadFriends(kTRUE)  (+ SetFriendFileName)
//   AliMCEventHandler  (MC);  SetReadTR(kTRUE) only if track refs wanted
//
//////////////////////////////////////////////////////////////////////////////

#include "AliAnalysisTaskSE.h"

#include <map>
#include <vector>

class TTree;
class TList;
class TH1F;
class AliESDEvent;
class AliESDfriend;
class AliESDtrack;
class AliESDfriendTrack;
class AliMCEvent;
class AliPIDResponse;
class AliTPCseed;

class AliAnalysisTaskSEFMClusterTreeESD : public AliAnalysisTaskSE
{
 public:
  // bits packed into hit_flag
  enum EHitFlag {
    kHitSharedSeed  = 1 << 0,  // cluster attached to >1 stored seed
    kHitUsed        = 1 << 1,  // AliTPCclusterMI::IsUsed(10)
    kHitSectorChgd  = 1 << 2,  // AliTPCclusterMI::IsSectorChanged() (distortion moved it)
    kHitEdge        = 1 << 3   // cluster type < 0 : edge cluster
  };

  enum EMCFlag {
    kMCPhysPrim   = 1 << 0,
    kMCSecWeak    = 1 << 1,
    kMCSecMat     = 1 << 2,
    kMCHasCluster = 1 << 3,  // referenced by >=1 stored cluster label
    kMCHasRecoTrk = 1 << 4,  // label of >=1 stored ESD track
    kMCAncestor   = 1 << 5   // kept only as an ancestor
  };

  AliAnalysisTaskSEFMClusterTreeESD();
  AliAnalysisTaskSEFMClusterTreeESD(const char *name);
  virtual ~AliAnalysisTaskSEFMClusterTreeESD();

  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *);
  virtual void Terminate(Option_t *);

  // ---- configuration -------------------------------------------------------
  void SetMC(Bool_t ismc)                   { fIsMC = ismc; }
  void SetUsePhysicsSelection(Bool_t b)     { fUsePhysSel = b; }
  void SetTriggerMask(UInt_t m)             { fTriggerMask = m; }
  void SetMaxVertexZ(Float_t z)             { fMaxVtxZ = z; }
  void SetMinTPCClusters(Int_t n)           { fMinTPCNcls = n; }
  void SetMinPt(Float_t pt)                 { fMinPt = pt; }
  void SetMaxAbsEta(Float_t eta)            { fMaxAbsEta = eta; }
  void SetStorePID(Bool_t b)                { fStorePID = b; }
  void SetStoreTrackRefs(Bool_t b)          { fStoreTrackRefs = b; }  // needs AliMCEventHandler::SetReadTR(kTRUE)
  void SetWriteEventsWithoutFriends(Bool_t b) { fWriteNoFriend = b; }  // no effect on the tree: events without TPC clusters are never written
  void SetMCPrimMaxAbsEta(Float_t eta)      { fMCPrimMaxAbsEta = eta; }
  void SetMaxEvents(Int_t n)                { fMaxEvents = n; }  // max WRITTEN events (= events with TPC clusters)

  static const Int_t kNTPCRows = 159;

 private:
  AliAnalysisTaskSEFMClusterTreeESD(const AliAnalysisTaskSEFMClusterTreeESD &);            // not implemented
  AliAnalysisTaskSEFMClusterTreeESD &operator=(const AliAnalysisTaskSEFMClusterTreeESD &); // not implemented

  void   ClearEvent();
  Bool_t AcceptTrack(const AliESDtrack *trk) const;
  void   FillTrack(AliESDtrack *trk, Int_t trkIdx);
  void   FillTPCClusters(const AliTPCseed *seed, Int_t trkIdx);
  void   FillMCParticles();
  const AliESDfriendTrack *GetFriendTrk(const AliESDtrack *trk) const;

  // ---- settings ------------------------------------------------------------
  Bool_t  fIsMC;
  Bool_t  fUsePhysSel;
  UInt_t  fTriggerMask;
  Float_t fMaxVtxZ;
  Int_t   fMinTPCNcls;
  Float_t fMinPt;
  Float_t fMaxAbsEta;
  Bool_t  fStorePID;
  Bool_t  fStoreTrackRefs;
  Bool_t  fWriteNoFriend;
  Float_t fMCPrimMaxAbsEta;
  Int_t   fMaxEvents;

  // ---- transient -----------------------------------------------------------
  AliESDEvent    *fESD;          //!
  AliESDfriend   *fESDfriend;    //!
  AliMCEvent     *fMCEvt;        //!
  AliPIDResponse *fPIDResponse;  //!
  TList          *fOutput;       //!
  TH1F           *fHEvents;      //!
  TH1F           *fHTracks;      //!
  TTree          *fTree;         //!
  Int_t           fNEventsWritten; //!
  std::map<Int_t, Int_t> fClIndexToHit; //! TPC cluster index -> hit position (dedup)

  // ---- event branches ------------------------------------------------------
  Int_t    fRunNumber;     //!
  UInt_t   fPeriod;        //!
  UInt_t   fOrbit;         //!
  UShort_t fBC;            //!
  Int_t    fEvInFile;      //!
  Float_t  fBz;            //!
  Float_t  fVtx[3];        //!
  Int_t    fVtxNContrib;   //!
  Float_t  fMCVtx[3];      //!
  Float_t  fImpactPar;     //!
  Int_t    fHasFriend;     //!

  // ---- hits ----------------------------------------------------------------
  std::vector<float> fHitX, fHitY, fHitZ;      //! global (cm)
  std::vector<int>   fHitLayer;                //! TPC global padrow 0-158
  std::vector<int>   fHitROC;                  //! TPC ROC 0-71 (IROC 0-35, OROC 36-71)
  std::vector<int>   fHitSector;               //! TPC sector 0-35 (A: 0-17, C: 18-35)
  std::vector<float> fHitPad, fHitTime;        //! TPC pad, time bin
  std::vector<float> fHitQ, fHitQmax;          //! TPC Qtot, Qmax (ADC)
  std::vector<float> fHitSY2, fHitSZ2;         //! cluster sigma^2 (cm^2)
  std::vector<float> fHitDX, fHitDY, fHitDZ;   //! applied distortion corrections (cm)
  std::vector<int>   fHitFlag;                 //! EHitFlag
  std::vector<int>   fHitTrk;                  //! first track (index into trk_*)
  std::vector<int>   fHitTrk2;                 //! second track or -1
  std::vector<int>   fHitNTrk;                 //! number of stored seeds using the cluster
  std::vector<int>   fHitMC0, fHitMC1, fHitMC2;//! cluster MC labels (-1 = none)

  // ---- tracks --------------------------------------------------------------
  std::vector<float> fTrkPt, fTrkEta, fTrkPhi; //!
  std::vector<int>   fTrkCharge;               //!
  std::vector<float> fTrkDCAxy, fTrkDCAz;      //!
  std::vector<int>   fTrkTPCNcls, fTrkTPCNclsF, fTrkTPCNCR; //!
  std::vector<float> fTrkTPCchi2, fTrkTPCdEdx, fTrkTPCmom;  //!
  std::vector<unsigned long long> fTrkStatus;  //!
  std::vector<int>   fTrkLabel, fTrkTPCLabel;  //!
  std::vector<int>   fTrkHasSeed;              //! 1 if a TPC seed was found in the friends
  std::vector<int>   fTrkNClStored;            //! clusters of this seed written (incl. shared)
  std::vector<float> fTrkNSigTPC[4];           //! e, pi, K, p

  // ---- MC particles --------------------------------------------------------
  std::vector<int>   fMCLabel, fMCPdg, fMCMother, fMCDauFirst, fMCDauLast, fMCCharge; //!
  std::vector<float> fMCPx, fMCPy, fMCPz, fMCVx, fMCVy, fMCVz, fMCVt;               //!
  std::vector<int>   fMCFlag, fMCGenIdx;       //!

  // ---- MC track references in the TPC (optional) --------------------------
  std::vector<float> fTRX, fTRY, fTRZ, fTRPx, fTRPy, fTRPz; //!
  std::vector<int>   fTRLabel;                 //!

  ClassDef(AliAnalysisTaskSEFMClusterTreeESD, 1);
};

#endif

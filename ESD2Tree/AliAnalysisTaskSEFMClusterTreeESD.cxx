/**************************************************************************
 * Copyright(c) 1998-2026, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

//////////////////////////////////////////////////////////////////////////////
//
// AliAnalysisTaskSEFMClusterTreeESD
//
// TPC clusters from ESD friends + ESD tracks + MC truth
// for Foundation-Model training. See header for details.
//
//////////////////////////////////////////////////////////////////////////////

#include "AliAnalysisTaskSEFMClusterTreeESD.h"

// ROOT / AliRoot headers (TF1.h, TFormula.h, AliComplexCluster.h, AliITSPIDResponse.h)
// trigger -Weffc++ / -Wnon-virtual-dtor / -Wignored-qualifiers that we cannot fix here.
// Suppress them for these includes only; warnings in this file stay enabled.
#if defined(__GNUC__) || defined(__clang__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Weffc++"
#pragma GCC diagnostic ignored "-Wnon-virtual-dtor"
#pragma GCC diagnostic ignored "-Wignored-qualifiers"
#endif
#include "AliAnalysisManager.h"
#include "AliCollisionGeometry.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliESDVertex.h"
#include "AliESDfriend.h"
#include "AliESDfriendTrack.h"
#include "AliESDtrack.h"
#include "AliExternalTrackParam.h"
#include "AliGenEventHeader.h"
#include "AliInputEventHandler.h"
#include "AliLog.h"
#include "AliMCEvent.h"
#include "AliMCParticle.h"
#include "AliPID.h"
#include "AliPIDResponse.h"
#include "AliTPCclusterMI.h"
#include "AliTPCseed.h"
#include "AliTrackReference.h"

#include "TH1F.h"
#include "TList.h"
#include "TMath.h"
#include "TParticle.h"
#include "TTree.h"
#if defined(__GNUC__) || defined(__clang__)
#pragma GCC diagnostic pop
#endif

ClassImp(AliAnalysisTaskSEFMClusterTreeESD)

namespace {
  const AliPID::EParticleType kPIDSpecies[4] = {AliPID::kElectron, AliPID::kPion, AliPID::kKaon, AliPID::kProton};
  const Int_t kClIndexNotUsed = 0x8000; // flag in AliTPCseed::GetClusterIndex2(row): cluster not used / not assigned
}

// Full member-initialiser list in declaration order (shared by both constructors).
// Keeps -Weffc++ (used by AliRoot/AliPhysics builds) and -Wreorder quiet.
#define FMTREE_MEMBER_INIT \
    fIsMC(kTRUE), fUsePhysSel(kFALSE), fTriggerMask(AliVEvent::kAny), fMaxVtxZ(10.), \
    fMinTPCNcls(1), fMinPt(0.), fMaxAbsEta(1.5), fStorePID(kFALSE), fStoreTrackRefs(kFALSE), \
    fWriteNoFriend(kFALSE), fMCPrimMaxAbsEta(1.5), fMaxEvents(-1), fESD(0), fESDfriend(0), \
    fMCEvt(0), fPIDResponse(0), fOutput(0), fHEvents(0), fHTracks(0), fTree(0), fNEventsWritten(0), \
    fClIndexToHit(), fRunNumber(0), fPeriod(0), fOrbit(0), fBC(0), fEvInFile(0), fBz(0), fVtx(), \
    fVtxNContrib(0), fMCVtx(), fImpactPar(-1), fHasFriend(0), fHitX(), fHitY(), fHitZ(), \
    fHitLayer(), fHitROC(), fHitSector(), fHitPad(), fHitTime(), fHitQ(), fHitQmax(), fHitSY2(), \
    fHitSZ2(), fHitDX(), fHitDY(), fHitDZ(), fHitFlag(), fHitTrk(), fHitTrk2(), fHitNTrk(), \
    fHitMC0(), fHitMC1(), fHitMC2(), fTrkPt(), fTrkEta(), fTrkPhi(), fTrkCharge(), fTrkDCAxy(), \
    fTrkDCAz(), fTrkTPCNcls(), fTrkTPCNclsF(), fTrkTPCNCR(), fTrkTPCchi2(), fTrkTPCdEdx(), \
    fTrkTPCmom(), fTrkStatus(), fTrkLabel(), fTrkTPCLabel(), fTrkHasSeed(), fTrkNClStored(), \
    fTrkNSigTPC(), fMCLabel(), fMCPdg(), fMCMother(), fMCDauFirst(), fMCDauLast(), fMCCharge(), \
    fMCPx(), fMCPy(), fMCPz(), fMCVx(), fMCVy(), fMCVz(), fMCVt(), fMCFlag(), fMCGenIdx(), fTRX(), \
    fTRY(), fTRZ(), fTRPx(), fTRPy(), fTRPz(), fTRLabel()

//____________________________________________________________________________
AliAnalysisTaskSEFMClusterTreeESD::AliAnalysisTaskSEFMClusterTreeESD()
  : AliAnalysisTaskSE(),
FMTREE_MEMBER_INIT
{
}

//____________________________________________________________________________
AliAnalysisTaskSEFMClusterTreeESD::AliAnalysisTaskSEFMClusterTreeESD(const char *name)
  : AliAnalysisTaskSE(name),
FMTREE_MEMBER_INIT
{
  DefineOutput(1, TList::Class());
  DefineOutput(2, TTree::Class());
}

#undef FMTREE_MEMBER_INIT

//____________________________________________________________________________
AliAnalysisTaskSEFMClusterTreeESD::~AliAnalysisTaskSEFMClusterTreeESD()
{
  if (AliAnalysisManager::GetAnalysisManager() &&
      AliAnalysisManager::GetAnalysisManager()->IsProofMode()) return;
  delete fOutput;
  delete fTree;
}

//____________________________________________________________________________
void AliAnalysisTaskSEFMClusterTreeESD::UserCreateOutputObjects()
{
  fOutput = new TList();
  fOutput->SetOwner(kTRUE);

  const char *evLab[10] = {"All", "ESD ok", "Trigger", "Vertex", "|Zvtx|", "Friend ok", "MC ok",
                           "Has TPC clusters", "Written", "Max reached"};
  fHEvents = new TH1F("hEvents", "Event counter", 10, 0.5, 10.5);
  for (Int_t i = 0; i < 10; ++i) fHEvents->GetXaxis()->SetBinLabel(i + 1, evLab[i]);
  fOutput->Add(fHEvents);

  const char *trLab[5] = {"Accepted", "Friend track", "TPC seed", "Clusters written", "Clusters shared"};
  fHTracks = new TH1F("hTracks", "Track / cluster counter", 5, 0.5, 5.5);
  for (Int_t i = 0; i < 5; ++i) fHTracks->GetXaxis()->SetBinLabel(i + 1, trLab[i]);
  fOutput->Add(fHTracks);

  OpenFile(2);
  fTree = new TTree("fmTree", "FM input (ESD+friends): TPC clusters, tracks, MC truth");

  // event
  fTree->Branch("run",         &fRunNumber,   "run/I");
  fTree->Branch("period",      &fPeriod,      "period/i");
  fTree->Branch("orbit",       &fOrbit,       "orbit/i");
  fTree->Branch("bc",          &fBC,          "bc/s");
  fTree->Branch("evInFile",    &fEvInFile,    "evInFile/I");
  fTree->Branch("bz",          &fBz,          "bz/F");
  fTree->Branch("vtx",         fVtx,          "vtx[3]/F");
  fTree->Branch("vtxNContrib", &fVtxNContrib, "vtxNContrib/I");
  fTree->Branch("mcVtx",       fMCVtx,        "mcVtx[3]/F");
  fTree->Branch("impactPar",   &fImpactPar,   "impactPar/F");
  fTree->Branch("hasFriend",   &fHasFriend,   "hasFriend/I");

  // hits
  fTree->Branch("hit_x",      &fHitX);
  fTree->Branch("hit_y",      &fHitY);
  fTree->Branch("hit_z",      &fHitZ);
  fTree->Branch("hit_layer",  &fHitLayer);
  fTree->Branch("hit_roc",    &fHitROC);
  fTree->Branch("hit_sector", &fHitSector);
  fTree->Branch("hit_pad",    &fHitPad);
  fTree->Branch("hit_time",   &fHitTime);
  fTree->Branch("hit_q",      &fHitQ);
  fTree->Branch("hit_qmax",   &fHitQmax);
  fTree->Branch("hit_sy2",    &fHitSY2);
  fTree->Branch("hit_sz2",    &fHitSZ2);
  fTree->Branch("hit_dx",     &fHitDX);
  fTree->Branch("hit_dy",     &fHitDY);
  fTree->Branch("hit_dz",     &fHitDZ);
  fTree->Branch("hit_flag",   &fHitFlag);
  fTree->Branch("hit_trk",    &fHitTrk);
  fTree->Branch("hit_trk2",   &fHitTrk2);
  fTree->Branch("hit_nTrk",   &fHitNTrk);
  fTree->Branch("hit_mc0",    &fHitMC0);
  fTree->Branch("hit_mc1",    &fHitMC1);
  fTree->Branch("hit_mc2",    &fHitMC2);

  // tracks
  fTree->Branch("trk_pt",        &fTrkPt);
  fTree->Branch("trk_eta",       &fTrkEta);
  fTree->Branch("trk_phi",       &fTrkPhi);
  fTree->Branch("trk_charge",    &fTrkCharge);
  fTree->Branch("trk_dcaXY",     &fTrkDCAxy);
  fTree->Branch("trk_dcaZ",      &fTrkDCAz);
  fTree->Branch("trk_tpcNcls",   &fTrkTPCNcls);
  fTree->Branch("trk_tpcNclsF",  &fTrkTPCNclsF);
  fTree->Branch("trk_tpcNCR",    &fTrkTPCNCR);
  fTree->Branch("trk_tpcChi2",   &fTrkTPCchi2);
  fTree->Branch("trk_tpcdEdx",   &fTrkTPCdEdx);
  fTree->Branch("trk_tpcMom",    &fTrkTPCmom);
  fTree->Branch("trk_status",    &fTrkStatus);
  fTree->Branch("trk_label",     &fTrkLabel);
  fTree->Branch("trk_tpcLabel",  &fTrkTPCLabel);
  fTree->Branch("trk_hasSeed",   &fTrkHasSeed);
  fTree->Branch("trk_nClStored", &fTrkNClStored);
  if (fStorePID) {
    fTree->Branch("trk_nSigTPC_e",  &fTrkNSigTPC[0]);
    fTree->Branch("trk_nSigTPC_pi", &fTrkNSigTPC[1]);
    fTree->Branch("trk_nSigTPC_K",  &fTrkNSigTPC[2]);
    fTree->Branch("trk_nSigTPC_p",  &fTrkNSigTPC[3]);
  }

  // MC
  if (fIsMC) {
    fTree->Branch("mc_label",    &fMCLabel);
    fTree->Branch("mc_pdg",      &fMCPdg);
    fTree->Branch("mc_mother",   &fMCMother);
    fTree->Branch("mc_dauFirst", &fMCDauFirst);
    fTree->Branch("mc_dauLast",  &fMCDauLast);
    fTree->Branch("mc_charge",   &fMCCharge);
    fTree->Branch("mc_px",       &fMCPx);
    fTree->Branch("mc_py",       &fMCPy);
    fTree->Branch("mc_pz",       &fMCPz);
    fTree->Branch("mc_vx",       &fMCVx);
    fTree->Branch("mc_vy",       &fMCVy);
    fTree->Branch("mc_vz",       &fMCVz);
    fTree->Branch("mc_vt",       &fMCVt);
    fTree->Branch("mc_flag",     &fMCFlag);
    fTree->Branch("mc_genIdx",   &fMCGenIdx);
    if (fStoreTrackRefs) {
      fTree->Branch("tr_x",     &fTRX);
      fTree->Branch("tr_y",     &fTRY);
      fTree->Branch("tr_z",     &fTRZ);
      fTree->Branch("tr_px",    &fTRPx);
      fTree->Branch("tr_py",    &fTRPy);
      fTree->Branch("tr_pz",    &fTRPz);
      fTree->Branch("tr_label", &fTRLabel);
    }
  }

  PostData(1, fOutput);
  PostData(2, fTree);
}

//____________________________________________________________________________
void AliAnalysisTaskSEFMClusterTreeESD::ClearEvent()
{
  fClIndexToHit.clear();

  fHitX.clear(); fHitY.clear(); fHitZ.clear(); fHitLayer.clear(); fHitROC.clear();
  fHitSector.clear(); fHitPad.clear(); fHitTime.clear(); fHitQ.clear(); fHitQmax.clear();
  fHitSY2.clear(); fHitSZ2.clear(); fHitDX.clear(); fHitDY.clear(); fHitDZ.clear();
  fHitFlag.clear(); fHitTrk.clear(); fHitTrk2.clear(); fHitNTrk.clear();
  fHitMC0.clear(); fHitMC1.clear(); fHitMC2.clear();

  fTrkPt.clear(); fTrkEta.clear(); fTrkPhi.clear(); fTrkCharge.clear(); fTrkDCAxy.clear(); fTrkDCAz.clear();
  fTrkTPCNcls.clear(); fTrkTPCNclsF.clear(); fTrkTPCNCR.clear(); fTrkTPCchi2.clear(); fTrkTPCdEdx.clear();
  fTrkTPCmom.clear(); fTrkStatus.clear(); fTrkLabel.clear(); fTrkTPCLabel.clear();
  fTrkHasSeed.clear(); fTrkNClStored.clear();
  for (Int_t i = 0; i < 4; ++i) fTrkNSigTPC[i].clear();

  fMCLabel.clear(); fMCPdg.clear(); fMCMother.clear(); fMCDauFirst.clear(); fMCDauLast.clear(); fMCCharge.clear();
  fMCPx.clear(); fMCPy.clear(); fMCPz.clear(); fMCVx.clear(); fMCVy.clear(); fMCVz.clear(); fMCVt.clear();
  fMCFlag.clear(); fMCGenIdx.clear();

  fTRX.clear(); fTRY.clear(); fTRZ.clear(); fTRPx.clear(); fTRPy.clear(); fTRPz.clear(); fTRLabel.clear();

  for (Int_t i = 0; i < 3; ++i) { fVtx[i] = -999.; fMCVtx[i] = -999.; }
  fVtxNContrib = 0; fImpactPar = -1.; fHasFriend = 0;
}

//____________________________________________________________________________
Bool_t AliAnalysisTaskSEFMClusterTreeESD::AcceptTrack(const AliESDtrack *trk) const
{
  if (!trk) return kFALSE;
  if (!(trk->GetStatus() & AliVTrack::kTPCin)) return kFALSE;
  if (trk->GetTPCNcls() < fMinTPCNcls) return kFALSE;
  if (trk->Pt() < fMinPt) return kFALSE;
  if (fMaxAbsEta > 0 && TMath::Abs(trk->Eta()) > fMaxAbsEta) return kFALSE;
  return kTRUE;
}

//____________________________________________________________________________
void AliAnalysisTaskSEFMClusterTreeESD::FillTrack(AliESDtrack *trk, Int_t /*trkIdx*/)
{
  fTrkPt.push_back(trk->Pt());
  fTrkEta.push_back(trk->Eta());
  fTrkPhi.push_back(trk->Phi());
  fTrkCharge.push_back(trk->Charge());

  Float_t b[2] = {-999., -999.}, bCov[3] = {0., 0., 0.};
  trk->GetImpactParameters(b, bCov);
  fTrkDCAxy.push_back(b[0]);
  fTrkDCAz.push_back(b[1]);

  fTrkTPCNcls.push_back(trk->GetTPCNcls());
  fTrkTPCNclsF.push_back(trk->GetTPCNclsF());
  fTrkTPCNCR.push_back((Int_t)trk->GetTPCCrossedRows());
  fTrkTPCchi2.push_back(trk->GetTPCchi2());
  fTrkTPCdEdx.push_back(trk->GetTPCsignal());
  const AliExternalTrackParam *tpcIn = trk->GetTPCInnerParam();
  fTrkTPCmom.push_back(tpcIn ? tpcIn->GetP() : trk->P());
  fTrkStatus.push_back(trk->GetStatus());
  fTrkLabel.push_back(trk->GetLabel());
  fTrkTPCLabel.push_back(trk->GetTPCLabel());
  fTrkHasSeed.push_back(0);
  fTrkNClStored.push_back(0);

  if (fStorePID) {
    for (Int_t s = 0; s < 4; ++s)
      fTrkNSigTPC[s].push_back(fPIDResponse ? fPIDResponse->NumberOfSigmasTPC(trk, kPIDSpecies[s]) : -999.);
  }
}

//____________________________________________________________________________
void AliAnalysisTaskSEFMClusterTreeESD::FillTPCClusters(const AliTPCseed *seed, Int_t trkIdx)
{
  Int_t nStored = 0;
  for (Int_t irow = 0; irow < kNTPCRows; ++irow) {
    const Int_t index = seed->GetClusterIndex2(irow);
    if (index < 0 || (index & kClIndexNotUsed)) continue;
    const AliTPCclusterMI *cl = seed->GetClusterPointer(irow);
    if (!cl) continue;
    ++nStored;

    // the same physical cluster can be attached to several seeds -> write once
    std::map<Int_t, Int_t>::iterator it = fClIndexToHit.find(index);
    if (it != fClIndexToHit.end()) {
      const Int_t h = it->second;
      if (fHitNTrk[h] == 1) fHitTrk2[h] = trkIdx;
      fHitNTrk[h] += 1;
      fHitFlag[h] |= kHitSharedSeed;
      fHTracks->Fill(5);
      continue;
    }

    // cluster x,y are in the sector tracking frame (already transformed:
    // alignment + distortion corrections); rotate by the sector angle.
    // (AliCluster::GetGlobalXYZ() would need gGeoManager.)
    const Int_t roc    = cl->GetDetector();      // 0-35 IROC, 36-71 OROC
    const Int_t sector = roc % 36;               // 0-17 A side, 18-35 C side
    const Double_t alpha = ((sector % 18) * 20. + 10.) * TMath::DegToRad();
    const Double_t ca = TMath::Cos(alpha), sa = TMath::Sin(alpha);
    const Double_t lx = cl->GetX(), ly = cl->GetY();

    Int_t flag = 0;
    if (cl->IsUsed(10))        flag |= kHitUsed;
    if (cl->IsSectorChanged()) flag |= kHitSectorChgd;
    if (cl->GetType() < 0)     flag |= kHitEdge;

    fClIndexToHit[index] = fHitX.size();
    fHitX.push_back(lx * ca - ly * sa);
    fHitY.push_back(lx * sa + ly * ca);
    fHitZ.push_back(cl->GetZ());
    fHitLayer.push_back(irow);                   // global pad row 0-158
    fHitROC.push_back(roc);
    fHitSector.push_back(sector);
    fHitPad.push_back(cl->GetPad());
    fHitTime.push_back(cl->GetTimeBin());
    fHitQ.push_back(cl->GetQ());
    fHitQmax.push_back(cl->GetMax());
    fHitSY2.push_back(cl->GetSigmaY2());
    fHitSZ2.push_back(cl->GetSigmaZ2());
    fHitDX.push_back(cl->GetDistortionX());
    fHitDY.push_back(cl->GetDistortionY());
    fHitDZ.push_back(cl->GetDistortionZ());
    fHitFlag.push_back(flag);
    fHitTrk.push_back(trkIdx);
    fHitTrk2.push_back(-1);
    fHitNTrk.push_back(1);
    fHitMC0.push_back(cl->GetLabel(0));
    fHitMC1.push_back(cl->GetLabel(1));
    fHitMC2.push_back(cl->GetLabel(2));
    fHTracks->Fill(4);
  }
  fTrkNClStored[trkIdx] = nStored;
}

//____________________________________________________________________________
const AliESDfriendTrack *AliAnalysisTaskSEFMClusterTreeESD::GetFriendTrk(const AliESDtrack *trk) const
{
  // as in AliTPCAnalysisTaskcalib::Exec(): AliESDtrack::GetFriendTrack()
  const AliESDfriendTrack *ft = trk->GetFriendTrack();
  if (!ft && fESDfriend) {
    const Int_t fid = trk->GetFriendTrackID();
    if (fid >= 0 && fid < fESDfriend->GetNumberOfTracks()) ft = fESDfriend->GetTrack(fid);
  }
  return ft;
}

//____________________________________________________________________________
void AliAnalysisTaskSEFMClusterTreeESD::FillMCParticles()
{
  const Int_t nMC = fMCEvt->GetNumberOfTracks();
  std::vector<int> keep(nMC, 0);

  // particles seen by stored clusters / tracks
  for (size_t i = 0; i < fHitMC0.size(); ++i) {
    const Int_t lab[3] = {fHitMC0[i], fHitMC1[i], fHitMC2[i]};
    for (Int_t k = 0; k < 3; ++k) if (lab[k] >= 0 && lab[k] < nMC) keep[lab[k]] |= kMCHasCluster;
  }
  for (size_t i = 0; i < fTrkLabel.size(); ++i) {
    const Int_t lab = TMath::Abs(fTrkLabel[i]);
    if (lab < nMC) keep[lab] |= kMCHasRecoTrk;
  }

  // charged physical primaries in acceptance (efficiency denominators)
  for (Int_t i = 0; i < nMC; ++i) {
    if (!fMCEvt->IsPhysicalPrimary(i)) continue;
    AliMCParticle *p = (AliMCParticle *)fMCEvt->GetTrack(i);
    if (!p || p->Charge() == 0) continue;
    if (fMCPrimMaxAbsEta > 0 && TMath::Abs(p->Eta()) > fMCPrimMaxAbsEta) continue;
    keep[i] |= kMCPhysPrim;
  }

  // full ancestry of everything kept (needed later for decay-topology labels, e.g. b/c hadrons)
  for (Int_t i = 0; i < nMC; ++i) {
    if (!keep[i] || (keep[i] == kMCAncestor)) continue;
    Int_t m = fMCEvt->GetLabelOfParticleMother(i);
    while (m >= 0 && m < nMC && !keep[m]) {
      keep[m] = kMCAncestor;
      m = fMCEvt->GetLabelOfParticleMother(m);
    }
  }

  for (Int_t i = 0; i < nMC; ++i) {
    if (!keep[i]) continue;
    AliMCParticle *p = (AliMCParticle *)fMCEvt->GetTrack(i);
    if (!p) continue;
    Int_t flag = keep[i] & (kMCHasCluster | kMCHasRecoTrk | kMCAncestor);
    if (fMCEvt->IsPhysicalPrimary(i))        flag |= kMCPhysPrim;
    if (fMCEvt->IsSecondaryFromWeakDecay(i)) flag |= kMCSecWeak;
    if (fMCEvt->IsSecondaryFromMaterial(i))  flag |= kMCSecMat;

    fMCLabel.push_back(i);
    fMCPdg.push_back(p->PdgCode());
    fMCMother.push_back(p->GetMother());
    fMCDauFirst.push_back(p->GetDaughterFirst());
    fMCDauLast.push_back(p->GetDaughterLast());
    fMCCharge.push_back(p->Charge() / 3);      // TParticlePDG charge is in |e|/3
    fMCPx.push_back(p->Px());
    fMCPy.push_back(p->Py());
    fMCPz.push_back(p->Pz());
    fMCVx.push_back(p->Xv());
    fMCVy.push_back(p->Yv());
    fMCVz.push_back(p->Zv());
    fMCVt.push_back(p->Particle() ? p->Particle()->T() : 0.);
    fMCFlag.push_back(flag);
    fMCGenIdx.push_back(p->GetGeneratorIndex());

    // true TPC crossing points (requires AliMCEventHandler::SetReadTR(kTRUE))
    if (fStoreTrackRefs && (keep[i] & (kMCHasCluster | kMCPhysPrim))) {
      const Int_t nRef = p->GetNumberOfTrackReferences();
      for (Int_t r = 0; r < nRef; ++r) {
        AliTrackReference *ref = p->GetTrackReference(r);
        if (!ref || ref->DetectorId() != AliTrackReference::kTPC) continue;
        fTRX.push_back(ref->X());   fTRY.push_back(ref->Y());   fTRZ.push_back(ref->Z());
        fTRPx.push_back(ref->Px()); fTRPy.push_back(ref->Py()); fTRPz.push_back(ref->Pz());
        fTRLabel.push_back(i);
      }
    }
  }
}

//____________________________________________________________________________
void AliAnalysisTaskSEFMClusterTreeESD::UserExec(Option_t *)
{
  fHEvents->Fill(1);
  if (fMaxEvents > 0 && fNEventsWritten >= fMaxEvents) { fHEvents->Fill(10); return; }

  fESD = dynamic_cast<AliESDEvent *>(InputEvent());
  if (!fESD) { AliError("No ESD event"); return; }
  fHEvents->Fill(2);

  AliInputEventHandler *inputHandler =
    (AliInputEventHandler *)AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler();
  if (fUsePhysSel && (!inputHandler || !(inputHandler->IsEventSelected() & fTriggerMask))) return;
  fHEvents->Fill(3);

  ClearEvent();

  // ---- vertex (tracks, fall back to SPD) -----------------------------------
  const AliESDVertex *vtx = fESD->GetPrimaryVertexTracks();
  if (!vtx || vtx->GetNContributors() < 1) vtx = fESD->GetPrimaryVertexSPD();
  if (!vtx || vtx->GetNContributors() < 1) return;
  fHEvents->Fill(4);
  if (fMaxVtxZ > 0 && TMath::Abs(vtx->GetZ()) > fMaxVtxZ) return;
  fHEvents->Fill(5);

  // ---- friends ----------------------------------------------------------------
  fESDfriend = fESD->FindFriend();
  if (!fESDfriend) {
    AliESDInputHandler *esdH = dynamic_cast<AliESDInputHandler *>(inputHandler);
    if (esdH) fESDfriend = esdH->GetESDfriend();
  }
  fHasFriend = (fESDfriend && !fESDfriend->TestSkipBit()) ? 1 : 0;
  // no friend -> no TPC seed -> no clusters: such an event is never written (see below),
  // so it is rejected here already (fWriteNoFriend only keeps it counted in "Friend ok")
  if (!fHasFriend && !fWriteNoFriend) return;
  fHEvents->Fill(6);

  // ---- MC -------------------------------------------------------------------
  fMCEvt = 0;
  if (fIsMC) {
    fMCEvt = MCEvent();
    if (!fMCEvt) { AliError("No MC event"); return; }
    const AliVVertex *mcVtx = fMCEvt->GetPrimaryVertex();
    if (mcVtx) { fMCVtx[0] = mcVtx->GetX(); fMCVtx[1] = mcVtx->GetY(); fMCVtx[2] = mcVtx->GetZ(); }
    AliCollisionGeometry *cg = dynamic_cast<AliCollisionGeometry *>(fMCEvt->GenEventHeader());
    if (cg) fImpactPar = cg->ImpactParameter();
  }
  fHEvents->Fill(7);

  fRunNumber   = fESD->GetRunNumber();
  fPeriod      = fESD->GetPeriodNumber();
  fOrbit       = fESD->GetOrbitNumber();
  fBC          = fESD->GetBunchCrossNumber();
  fEvInFile    = fESD->GetEventNumberInFile();
  fBz          = fESD->GetMagneticField();
  fVtx[0] = vtx->GetX(); fVtx[1] = vtx->GetY(); fVtx[2] = vtx->GetZ();
  fVtxNContrib = vtx->GetNContributors();
  if (fStorePID && inputHandler) fPIDResponse = inputHandler->GetPIDResponse();

  // ---- tracks / clusters ----------------------------------------------------
  const Int_t nTracks = fESD->GetNumberOfTracks();
  Int_t nStored = 0;
  for (Int_t it = 0; it < nTracks; ++it) {
    AliESDtrack *trk = fESD->GetTrack(it);
    if (!AcceptTrack(trk)) continue;
    fHTracks->Fill(1);
    FillTrack(trk, nStored);
    const AliESDfriendTrack *ft = fHasFriend ? GetFriendTrk(trk) : 0;
    if (ft) {
      fHTracks->Fill(2);
      // pointer access to the stored calib object avoids the per-track copy
      // done by AliESDfriendTrack::GetTPCseed(AliTPCseed&)
      const AliTPCseed *seed = dynamic_cast<const AliTPCseed *>(ft->GetTPCseed());
      if (seed) {
        fHTracks->Fill(3);
        fTrkHasSeed[nStored] = 1;
        FillTPCClusters(seed, nStored);
      }
    }
    ++nStored;
  }

  // ---- write only events with at least one stored TPC cluster -----------------
  // (no friend, no TPC seeds, or seeds without clusters -> event skipped;
  //  MC particles are not even collected for such events)
	cout << "Number of TPC clusters: " << fHitX.size() << endl;
  if (fHitX.empty()) return;
  fHEvents->Fill(8);

  if (fIsMC) FillMCParticles();

  fTree->Fill();
  ++fNEventsWritten;
  fHEvents->Fill(9);

  PostData(1, fOutput);
  PostData(2, fTree);
}

//____________________________________________________________________________
void AliAnalysisTaskSEFMClusterTreeESD::Terminate(Option_t *)
{
}

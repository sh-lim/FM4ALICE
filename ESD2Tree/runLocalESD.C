// runLocalESD.C -- local test on MC ESDs with friends
//
//   root -l -b -q 'runLocalESD.C("esdlist.txt", 1, 50)'
//
// esdlist.txt : one AliESDs.root path per line. AliESDfriends.root, galice.root,
//               Kinematics.root (and TrackRefs.root if track refs are wanted)
//               must sit in the same directory as each AliESDs.root.

#include <fstream>

void runLocalESD(TString esd = "AliESDs.root", Bool_t isMC = kTRUE, Long64_t nEvents = 1000000000,
                 Bool_t storeTrackRefs = kFALSE)
{
  gInterpreter->ProcessLine(".include $ROOTSYS/include");
  gInterpreter->ProcessLine(".include $ALICE_ROOT/include");
  gInterpreter->ProcessLine(".include $ALICE_PHYSICS/include");

  AliAnalysisManager *mgr = new AliAnalysisManager("FMClusterTreeESD");

  AliESDInputHandler *esdH = new AliESDInputHandler();
  esdH->SetReadFriends(kTRUE);                   // TPC seeds + clusters live in the friends
  esdH->SetFriendFileName("AliESDfriends.root");
  mgr->SetInputEventHandler(esdH);

  if (isMC) {
    AliMCEventHandler *mcH = new AliMCEventHandler();
    mcH->SetReadTR(storeTrackRefs);
    mgr->SetMCtruthEventHandler(mcH);
  }

  // PID response only if SetStorePID(kTRUE):
  // gInterpreter->ExecuteMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C(kTRUE)");

  gInterpreter->LoadMacro("AliAnalysisTaskSEFMClusterTreeESD.cxx++g");
  gInterpreter->ExecuteMacro(Form("AddTaskFMClusterTreeESD.C(%d, 1, %d, kFALSE, \"\")", (Int_t)isMC, (Int_t)storeTrackRefs));

  if (!mgr->InitAnalysis()) return;
  mgr->PrintStatus();

  TChain *chain = new TChain("esdTree");
	chain->Add(esd.Data());

  mgr->StartAnalysis("local", chain, nEvents);
}

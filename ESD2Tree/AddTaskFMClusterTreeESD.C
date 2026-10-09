// AddTaskFMClusterTreeESD.C
// ESD + friends version. The input handler must read friends:
//   AliESDInputHandler *esdH = new AliESDInputHandler();
//   esdH->SetReadFriends(kTRUE);
//   esdH->SetFriendFileName("AliESDfriends.root");
// For MC: AliMCEventHandler (SetReadTR(kTRUE) only if storeTrackRefs).

AliAnalysisTaskSEFMClusterTreeESD *AddTaskFMClusterTreeESD(Bool_t isMC = kTRUE,
                                                           Int_t minTPCNcls = 1,
                                                           Bool_t storeTrackRefs = kFALSE,
                                                           Bool_t storePID = kFALSE,
                                                           TString suffix = "")
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) { ::Error("AddTaskFMClusterTreeESD", "No analysis manager"); return 0x0; }
  AliESDInputHandler *esdH = dynamic_cast<AliESDInputHandler *>(mgr->GetInputEventHandler());
  if (!esdH) { ::Error("AddTaskFMClusterTreeESD", "ESD input handler required"); return 0x0; }
  if (!esdH->GetReadFriends())
    ::Warning("AddTaskFMClusterTreeESD", "Friends are not read: no TPC clusters will be written. Call SetReadFriends(kTRUE).");

  TString name = Form("FMClusterTreeESD%s", suffix.Data());
  AliAnalysisTaskSEFMClusterTreeESD *task = new AliAnalysisTaskSEFMClusterTreeESD(name.Data());
  task->SetMC(isMC);
  task->SetMinTPCClusters(minTPCNcls);
  task->SetMaxVertexZ(10.);
  task->SetMaxAbsEta(1.5);
  task->SetStoreTrackRefs(storeTrackRefs);
  task->SetStorePID(storePID);
  task->SetUsePhysicsSelection(kFALSE);
  task->SetWriteEventsWithoutFriends(kTRUE);
  mgr->AddTask(task);

  TString qaFile = AliAnalysisManager::GetCommonFileName();
  qaFile += Form(":FMClusterTreeESD%s", suffix.Data());
  AliAnalysisDataContainer *cQA   = mgr->CreateContainer(Form("FMQAESD%s", suffix.Data()), TList::Class(),
                                                         AliAnalysisManager::kOutputContainer, qaFile.Data());
  AliAnalysisDataContainer *cTree = mgr->CreateContainer(Form("fmTreeESD%s", suffix.Data()), TTree::Class(),
                                                         AliAnalysisManager::kOutputContainer,
                                                         Form("FMClusterTreeESD%s.root", suffix.Data()));
  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 1, cQA);
  mgr->ConnectOutput(task, 2, cTree);
  return task;
}

// plotHitsGen3D.C -- 3D TPC cluster maps of ONE event, split by mc_genIdx
//
// For the chosen event one TH3F (x, y, z in cm) is filled per generator index:
//   h3Gen<G>_ev<N>  : hits whose track's MC particle has mc_genIdx == G,  G = genMin ... genMin+nGen-1
// All nGen histograms are drawn on one canvas divided into nCol x nRow pads (default 10 x 5 = 50).
//
// Track-based matching only (hit_mc0/1/2 are NOT used), as in plotHits3D.C:
//   hit --hit_trk--> track index --trk_label--> MC particle with the same mc_label --> mc_genIdx
//
// Usage:
//   root -l 'plotHitsGen3D.C("FMClusterTreeESD.root", 3)'                 // event 3, genIdx 0..49
//   root -l 'plotHitsGen3D.C("FMClusterTreeESD.root", -1)'                // first event that has hits
//   root -l -b -q 'plotHitsGen3D.C("in.root", 3, 0, 50, 1)'               // trk_tpcLabel
//
// Arguments
//   inFile
//   event        : entry number of the event to plot (-1 = first event with nHit > 0)
//   genMin, nGen : genIdx range [genMin, genMin + nGen); one pad per value
//   labelType    : 0 = trk_label (global), 1 = trk_tpcLabel (TPC clusters only)
//   acceptFake   : kTRUE -> fake tracks (label < 0) use |label| (ALICE convention), kFALSE -> skipped
//   nBins        : bins per axis (x, y, z all in [-250, 250] cm)
//                  memory ~ nGen * (nBins+2)^3 * 4 B : 50 hists x 50^3 ~ 28 MB, x 100^3 ~ 210 MB
//   nCol, nRow   : canvas division (nCol * nRow must be >= nGen)
//   outFile      : histograms and the canvas are written here; the PNG goes next to it

#include <map>
#include <unordered_map>
#include <vector>

#include "TCanvas.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TH3F.h"
#include "TKey.h"
#include "TROOT.h"
#include "TString.h"
#include "TStyle.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"

namespace {
// find the tree by name at the top level or one directory down
TTree *FindTree(TFile *f, const char *name)
{
  TTree *t = f->Get<TTree>(name);
  if (t) return t;
  TIter next(f->GetListOfKeys());
  while (TKey *k = (TKey *)next()) {
    if (!TClass::GetClass(k->GetClassName())->InheritsFrom(TDirectory::Class())) continue;
    TDirectory *d = (TDirectory *)k->ReadObj();
    if ((t = d->Get<TTree>(name))) return t;
  }
  return 0;
}
}

void plotHitsGen3D(const char *inFile = "treeout_fm4alice_tpc_00008.root", Long64_t event = 1,
                   Int_t genMin = 0, Int_t nGen = 50,
                   Int_t labelType = 0, Bool_t acceptFake = kTRUE, Int_t nBins = 50,
                   Int_t nCol = 10, Int_t nRow = 5,
                   const char *outFile = "hitsGen3D.root", const char *treeName = "fmTree")
{
  if (nCol * nRow < nGen) {
    ::Error("plotHitsGen3D", "canvas %d x %d has fewer pads than nGen = %d", nCol, nRow, nGen);
    return;
  }

  TFile *fin = TFile::Open(inFile);
  if (!fin || fin->IsZombie()) { ::Error("plotHitsGen3D", "cannot open %s", inFile); return; }
  TTree *tree = FindTree(fin, treeName);
  if (!tree) { ::Error("plotHitsGen3D", "tree '%s' not found in %s", treeName, inFile); return; }

  const char *trkLabBr = (labelType == 1) ? "trk_tpcLabel" : "trk_label";
  const char *need[] = {"hit_x", "hit_y", "hit_z", "hit_trk", trkLabBr, "mc_label", "mc_genIdx"};
  for (const char *b : need)
    if (!tree->GetBranch(b)) { ::Error("plotHitsGen3D", "branch %s missing (MC tree required)", b); return; }

  TTreeReader reader(tree);
  TTreeReaderValue<std::vector<float>> hx(reader, "hit_x"), hy(reader, "hit_y"), hz(reader, "hit_z");
  TTreeReaderValue<std::vector<int>>   hTrk(reader, "hit_trk");
  TTreeReaderValue<std::vector<int>>   trkLab(reader, trkLabBr);
  TTreeReaderValue<std::vector<int>>   mcLab(reader, "mc_label"), mcGen(reader, "mc_genIdx");

  const Long64_t nTot = reader.GetEntries(true);

  // ---- load the requested event -------------------------------------------------
  Long64_t iev = event;
  if (event < 0) {                                   // first event with hits
    iev = -1;
    while (reader.Next())
      if (!hx->empty()) { iev = reader.GetCurrentEntry(); break; }
    if (iev < 0) { ::Error("plotHitsGen3D", "no event with hits in %s", inFile); return; }
  } else {
    if (event >= nTot) { ::Error("plotHitsGen3D", "event %lld >= entries %lld", event, nTot); return; }
    if (reader.SetEntry(event) != TTreeReader::kEntryValid) {
      ::Error("plotHitsGen3D", "cannot read entry %lld", event); return;
    }
  }
  const size_t nHit = hx->size(), nTrk = trkLab->size();
  ::Info("plotHitsGen3D", "%s: event %lld, %zu hits, %zu tracks, genIdx %d..%d via %s",
         inFile, iev, nHit, nTrk, genMin, genMin + nGen - 1, trkLabBr);
  if (nHit == 0) ::Warning("plotHitsGen3D", "event %lld has no hits; pads will be empty", iev);

  // ---- MC label -> genIdx, and number of MC particles per genIdx -----------------
  std::unordered_map<int, int> lab2gen;
  std::map<int, int> nPartPerGen;
  lab2gen.reserve(mcLab->size());
  for (size_t i = 0; i < mcLab->size(); ++i) {
    lab2gen[(*mcLab)[i]] = (*mcGen)[i];
    ++nPartPerGen[(*mcGen)[i]];
  }

  // ---- one TH3F per genIdx ------------------------------------------------------
  TFile *fout = TFile::Open(outFile, "RECREATE");
  const Double_t R = 250.;
  std::vector<TH3F *> h(nGen);
  std::vector<TH2F *> h2(nGen);
  for (Int_t k = 0; k < nGen; ++k) {
    const Int_t g = genMin + k;
    h[k] = new TH3F(Form("h3Gen%d_ev%lld", g, iev),
                    Form("genIdx = %d;x (cm);y (cm);z (cm)", g),
                    nBins, -R, R, nBins, -R, R, nBins, -R, R);
    h[k]->SetDirectory(nullptr);       // keep them alive after fout->Close() for the canvas

    h2[k] = new TH2F(Form("h2Gen%d_ev%lld", g, iev),
                    Form("genIdx = %d;z (cm);r (cm)", g),
										nBins, -R, R, nBins, 0, R);
    h2[k]->SetDirectory(nullptr);       // keep them alive after fout->Close() for the canvas
  }

  // ---- fill ---------------------------------------------------------------------
  std::vector<Long64_t> nSel(nGen, 0);
  Long64_t nNoTrk = 0, nFake = 0, nNoMC = 0, nOutRange = 0;
  for (size_t ih = 0; ih < nHit; ++ih) {
    // 1) track of this hit
    const int it = (*hTrk)[ih];
    if (it < 0 || it >= (int)nTrk) { ++nNoTrk; continue; }

    // 2) MC label of that track; fake tracks have label < 0
    int lab = (*trkLab)[it];
    if (lab < 0) { ++nFake; if (!acceptFake) continue; lab = -lab; }

    // 3) MC particle with mc_label == trk_label -> its mc_genIdx
    auto m = lab2gen.find(lab);
    if (m == lab2gen.end()) { ++nNoMC; continue; }

    const Int_t k = m->second - genMin;
    if (k < 0 || k >= nGen) { ++nOutRange; continue; }

		float xx = (*hx)[ih];
		float yy = (*hy)[ih];
		float zz = (*hz)[ih];

    h[k]->Fill(xx, yy, zz);
		h2[k]->Fill(zz, sqrt(xx*xx + yy*yy));
    ++nSel[k];
  }


  // ---- summary ------------------------------------------------------------------
  Long64_t nInRange = 0;
  for (Int_t k = 0; k < nGen; ++k) {
    nInRange += nSel[k];
    const Int_t g = genMin + k;
    const int nPart = nPartPerGen.count(g) ? nPartPerGen[g] : 0;
    if (nSel[k] > 0 || nPart > 0)
      ::Info("plotHitsGen3D", "  genIdx %3d : %6lld hits  (%d MC particles)", g, nSel[k], nPart);
  }
  ::Info("plotHitsGen3D", "hits %zu | in genIdx range %lld | other genIdx %lld | no-track %lld  "
         "fake-track hits %lld%s  label-not-in-MC %lld",
         nHit, nInRange, nOutRange, nNoTrk, nFake, acceptFake ? " (used |label|)" : " (skipped)", nNoMC);
  if (nPartPerGen.size() == 1 && nPartPerGen.begin()->first == -1)
    ::Warning("plotHitsGen3D", "all mc_genIdx are -1 (single-generator MC, no cocktail): use genMin = -1, nGen = 1");

	//return;

  // ---- draw: nCol x nRow pads ------------------------------------------------------
  gStyle->SetOptStat(0);
  gStyle->SetTitleFontSize(0.09);
  TCanvas *c = new TCanvas(Form("cGen_ev%lld", iev), Form("event %lld: hits per genIdx", iev),
                           150 * nCol, 150 * nRow);
  c->Divide(nCol, nRow, 0.001, 0.001);

  TCanvas *c2 = new TCanvas(Form("cGen2_ev%lld", iev), Form("event %lld: hits per genIdx", iev),
                           150 * nCol, 150 * nRow);
  c2->Divide(nCol, nRow, 0.001, 0.001);

  for (Int_t k = 0; k < nGen; ++k) {
    c->cd(k + 1);
    gPad->SetMargin(0.02, 0.02, 0.05, 0.10);
    h[k]->SetTitle(Form("genIdx = %d  (%lld hits)", genMin + k, nSel[k]));
    for (TAxis *ax : {h[k]->GetXaxis(), h[k]->GetYaxis(), h[k]->GetZaxis()}) {
      ax->SetLabelSize(0.04);
      ax->SetTitleSize(0.05);
      ax->SetNdivisions(503);
    }
    h[k]->SetMarkerStyle(1);
    h[k]->SetMarkerColor(kRed + 1);
    h[k]->Draw("BOX2");
  }
  c->Update();

  for (Int_t k = 0; k < nGen; ++k) {
    c2->cd(k + 1);
    gPad->SetMargin(0.02, 0.02, 0.05, 0.10);
    h2[k]->SetTitle(Form("genIdx = %d  (%lld hits)", genMin + k, nSel[k]));
    for (TAxis *ax : {h2[k]->GetXaxis(), h2[k]->GetYaxis(), h2[k]->GetZaxis()}) {
      ax->SetLabelSize(0.04);
      ax->SetTitleSize(0.05);
      ax->SetNdivisions(503);
    }
    h2[k]->SetMarkerStyle(1);
    h2[k]->SetMarkerColor(kRed + 1);
    h2[k]->Draw("BOX");
  }
  c2->Update();

	return;

  TString png(outFile);
  png.ReplaceAll(".root", "");
  png += Form("_ev%lld.png", iev);
  c->SaveAs(png);

  fout->cd();
  for (TH3F *hk : h) hk->Write();
  c->Write();
  fout->Close();
  ::Info("plotHitsGen3D", "histograms + canvas written to %s, picture %s", outFile, png.Data());
}

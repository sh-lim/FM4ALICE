// plotHits3D.C -- per-event 3D cluster maps from FMClusterTreeESD.root
//
// For every selected event two TH3F (x, y, z in cm) are filled:
//   h3All_ev<N>      : all hits of the event
//   h3Gen<G>_ev<N>   : hits whose track's MC particle has mc_genIdx == G (default G = 0)
//
// Track-based matching only (hit_mc0/1/2 are NOT used):
//   hit --hit_trk--> track index --trk_label--> MC particle with the same mc_label --> mc_genIdx
//   hit is selected if that mc_genIdx == genIdxSel
//
// Usage:
//   root -l 'plotHits3D.C("FMClusterTreeESD.root")'
//   root -l -b -q 'plotHits3D.C("in.root", "hits3D.root", 0, 20, 0, 1)'   // 20 events, trk_tpcLabel
//
// Arguments
//   inFile, outFile
//   firstEvent, nEvents  : event range (nEvents < 0 : all)
//   genIdxSel            : mc_genIdx value for the second histogram
//   labelType            : 0 = trk_label (global), 1 = trk_tpcLabel (TPC clusters only)
//   acceptFake           : kTRUE -> fake tracks (label < 0) use |label| (ALICE convention), kFALSE -> skipped
//   nBins                : bins per axis (x, y, z all in [-250, 250] cm)
//   nDraw                : number of events WITH HITS (nHit > 0) to draw, one PNG each (0 = none)
//                          counting starts at the first event with hits; empty events are skipped
//                          all processed events are always written to outFile

#include <algorithm>
#include <unordered_map>
#include <vector>

#include "TCanvas.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TH3F.h"
#include "TKey.h"
#include "TMath.h"
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

void plotHits3D(const char *inFile = "FMClusterTreeESD.root", const char *outFile = "hits3D.root",
                Long64_t firstEvent = 0, Long64_t nEvents = 10, Int_t genIdxSel = 1,
                Int_t labelType = 0, Bool_t acceptFake = kTRUE, Int_t nBins = 100,
                Int_t nDraw = 10, const char *treeName = "fmTree")
{
  TFile *fin = TFile::Open(inFile);
  if (!fin || fin->IsZombie()) { ::Error("plotHits3D", "cannot open %s", inFile); return; }
  TTree *tree = FindTree(fin, treeName);
  if (!tree) { ::Error("plotHits3D", "tree '%s' not found in %s", treeName, inFile); return; }

  const char *trkLabBr = (labelType == 1) ? "trk_tpcLabel" : "trk_label";
  const char *need[] = {"hit_x", "hit_y", "hit_z", "hit_trk", trkLabBr, "mc_label", "mc_genIdx"};
  for (const char *b : need)
    if (!tree->GetBranch(b)) { ::Error("plotHits3D", "branch %s missing (MC tree required)", b); return; }

  TTreeReader reader(tree);
  TTreeReaderValue<std::vector<float>> hx(reader, "hit_x"), hy(reader, "hit_y"), hz(reader, "hit_z");
  TTreeReaderValue<std::vector<int>>   hTrk(reader, "hit_trk");
  TTreeReaderValue<std::vector<int>>   trkLab(reader, trkLabBr);
  TTreeReaderValue<std::vector<int>>   mcLab(reader, "mc_label"), mcGen(reader, "mc_genIdx");

  const Long64_t nTot = reader.GetEntries(true);
  const Long64_t last = (nEvents < 0) ? nTot : TMath::Min(nTot, firstEvent + nEvents);
  if (firstEvent >= nTot) { ::Error("plotHits3D", "firstEvent %lld >= entries %lld", firstEvent, nTot); return; }
  ::Info("plotHits3D", "%s: %lld events, processing [%lld, %lld), selection mc_genIdx == %d via %s",
         inFile, nTot, firstEvent, last, genIdxSel, trkLabBr);

  TFile *fout = TFile::Open(outFile, "RECREATE");
  const Double_t R = 250.;

  // totals over all processed events
  Long64_t tHits = 0, tSel = 0, tNoTrk = 0, tFake = 0, tNoMC = 0;
  Int_t nDrawn = 0;           // events drawn so far (only events with nHit > 0 count)
  std::vector<int> genSeen;   // distinct genIdx values encountered (for the warning below)

  reader.SetEntriesRange(firstEvent, last);
  while (reader.Next()) {
    const Long64_t iev = reader.GetCurrentEntry();
    const size_t nHit = hx->size(), nTrk = trkLab->size();

		int nGenIdxSel = 0;

    // MC label -> genIdx for this event
    std::unordered_map<int, int> lab2gen;
    lab2gen.reserve(mcLab->size());
    for (size_t i = 0; i < mcLab->size(); ++i) {
      lab2gen[(*mcLab)[i]] = (*mcGen)[i];
			if ( genIdxSel==(*mcGen)[i] ) nGenIdxSel++;
      if (std::find(genSeen.begin(), genSeen.end(), (*mcGen)[i]) == genSeen.end()) genSeen.push_back((*mcGen)[i]);
    }

		cout << "CHECK nGenIdxSel: " << nGenIdxSel << ", genSeen size: " << genSeen.size() << endl;

    TH3F *hAll = new TH3F(Form("h3All_ev%lld", iev), Form("all hits, event %lld;x (cm);y (cm);z (cm)", iev),
                          nBins, -R, R, nBins, -R, R, nBins, -R, R);
    TH3F *hSel = new TH3F(Form("h3Gen%d_ev%lld", genIdxSel, iev),
                          Form("hits with track MC genIdx = %d, event %lld;x (cm);y (cm);z (cm)", genIdxSel, iev),
                          nBins, -R, R, nBins, -R, R, nBins, -R, R);
    hAll->SetDirectory(fout);
    hSel->SetDirectory(fout);

    Long64_t nSel = 0, nNoTrk = 0, nFake = 0, nNoMC = 0;
    for (size_t h = 0; h < nHit; ++h) {
      const float x = (*hx)[h], y = (*hy)[h], z = (*hz)[h];
      hAll->Fill(x, y, z);

      // 1) track of this hit: index hit_trk into the trk_* vectors
      const int it = (*hTrk)[h];
      if (it < 0 || it >= (int)nTrk) { ++nNoTrk; continue; }

      // 2) MC label of that track (trk_label); fake tracks have label < 0
      int lab = (*trkLab)[it];
      if (lab < 0) { ++nFake; if (!acceptFake) continue; lab = -lab; }

      // 3) MC particle with mc_label == trk_label -> its mc_genIdx
      auto m = lab2gen.find(lab);
      if (m == lab2gen.end()) { ++nNoMC; continue; }
			//if (m->second<10 ) cout << "m->second " << m->second << endl;
      if (m->second == genIdxSel) { hSel->Fill(x, y, z); ++nSel; }
    }

    if (nHit == 0) ::Warning("plotHits3D", "ev %lld: no hits (hasFriend=0 or no TPC seeds in friends)", iev);
    ::Info("plotHits3D", "ev %lld: hits %zu | selected %lld | no-track %lld  fake-track hits %lld  label-not-in-MC %lld",
           iev, nHit, nSel, nNoTrk, nFake, nNoMC);

    tHits += nHit; tSel += nSel; tNoTrk += nNoTrk; tFake += nFake; tNoMC += nNoMC;

    const Bool_t doDraw = (nHit > 0 && nDrawn < nDraw);
    if (doDraw) {
      gStyle->SetOptStat(0);
      TCanvas *c = new TCanvas(Form("c_ev%lld", iev), Form("event %lld", iev), 1400, 650);
      c->Divide(2, 1);
      c->cd(1); hAll->SetMarkerStyle(1); hAll->SetMarkerColor(kBlue + 1); hAll->Draw("BOX2");
      c->cd(2); hSel->SetMarkerStyle(1); hSel->SetMarkerColor(kRed + 1);  hSel->Draw("BOX2");
      c->SaveAs(Form("hits3D_ev%lld.png", iev));
      if (nDrawn == 0) ::Info("plotHits3D", "first event with hits: %lld", iev);
      ++nDrawn;
    }

    fout->cd();
    hAll->Write();
    hSel->Write();
    if (!doDraw) { delete hAll; delete hSel; }   // drawn ones stay alive for the canvases
  }

  ::Info("plotHits3D", "TOTAL hits %lld | selected %lld (%.3f) | no-track %lld  fake-track hits %lld%s  label-not-in-MC %lld",
         tHits, tSel, tHits ? double(tSel) / tHits : 0., tNoTrk, tFake,
         acceptFake ? " (used |label|)" : " (skipped)", tNoMC);
  TString gens;
  for (int g : genSeen) gens += Form(" %d", g);
  ::Info("plotHits3D", "mc_genIdx values present:%s", gens.Data());
  ::Info("plotHits3D", "PNG files written: %d (events with hits only)", nDrawn);
  if (genSeen.size() == 1 && genSeen[0] == -1)
    ::Warning("plotHits3D", "all mc_genIdx are -1 (single-generator MC, no cocktail): use genIdxSel = -1");

  fout->Close();
  ::Info("plotHits3D", "histograms written to %s", outFile);
}

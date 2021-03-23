#include <iostream>
#include <random>
#include <algorithm>
#include <iterator>
#include <iostream>
#include <vector>

#include <Rtypes.h>
#include <TChain.h>
#include <TFile.h>
#include <TClonesArray.h>
#include <TString.h>
#include <TMath.h>
#include <TStopwatch.h>
#include <TRandom3.h>

// Define whether new or old mpdroot version is used
#define _NEW_MCSTACK_
//#define _OLD_MCSTACK_

#include <FairMCEventHeader.h>
#ifdef _OLD_MCSTACK_
#include <FairMCTrack.h>
#endif
#ifdef _NEW_MCSTACK_
#include <MpdMCTrack.h>
#endif
#include <MpdEvent.h>
#include <MpdZdcDigi.h>
#include <MpdPid.h>
#include <MpdTrack.h>
#include <MpdKalmanTrack.h>
#include <MpdVertex.h>

void MpdDstReader(TChain *inChain, TString outFileName, Bool_t primaryOnly = false, Bool_t eventShuffle = false)
{
  TStopwatch timer;
  timer.Start();

  if (!inChain) return;

  const Double_t cut_pt = 0.15; // default: 0.15 GeV/c
  const Double_t cut_eta = 0.5; // default: 0.5
  const Int_t cut_nhits = 16;   // default: 16
  const Double_t dca_cut = 0.2; // dca cut

  TH1F *hRefMult = new TH1F("hRefMultSTAR","hRefMultSTAR",2500,0,2500);
  TH2F *hBvsRefMult = new TH2F("hBvsRefMult","hBvsRefMult",2500,0,2500,200,0.,20.);

  FairMCEventHeader *MCHeader;
  TClonesArray *MCTracks;
  MpdEvent *MPDEvent;
  //TClonesArray *FHCalHits;
  TClonesArray *MpdGlobalTracks;
  //MpdZdcDigi *FHCalHit;
  //TClonesArray *mpdKalmanTracks;
  //TClonesArray *vertexes;

  MCHeader = nullptr;
  MCTracks = nullptr;
  MPDEvent = nullptr;
  //FHCalHits = nullptr;
  //mpdKalmanTracks = nullptr;
  //vertexes = nullptr;

  inChain->SetBranchAddress("MCEventHeader.", &MCHeader);
  inChain->SetBranchAddress("MCTrack", &MCTracks);
  inChain->SetBranchAddress("MPDEvent.", &MPDEvent);
  //inChain->SetBranchAddress("ZdcDigi", &FHCalHits);
  //inChain->SetBranchAddress("TpcKalmanTrack", &mpdKalmanTracks);
  //inChain->SetBranchAddress("Vertex", &vertexes);
  
  std::vector<Long64_t> vEvents;
  std::random_device rd;
  std::mt19937 g(rd());

  Long64_t Nentries = 5e5; // (Long64_t) inChain->GetEntries();
  // Shuffle events - for PHSD only!
  /*if (eventShuffle)
  {
    std::cout << "Applying event shuffle..." << std::endl;
    for (Long64_t jentry=0; jentry<Nentries;jentry++)
    {
      vEvents.push_back(jentry);
    }
    std::shuffle(vEvents.begin(), vEvents.end(), g);
  }*/

  // Starting event loop
  Long64_t Nevents;
  if (!eventShuffle) Nevents = (Nentries < 5e5) ? Nentries : 5e5;
  if (eventShuffle) Nevents = Nentries;
  Int_t refMult, nhits;
  Double_t pt, eta;
  Long64_t iEvent;
  TRandom3 *rnd = new TRandom3();
  Double_t prob_skip;
  for (Long64_t jentry=0; jentry<Nevents;jentry++)
  {
    if (!eventShuffle) iEvent = jentry;
    //if (eventShuffle) iEvent  = vEvents.at(jentry);
    if (eventShuffle)
    {
      prob_skip = rnd->Rndm();
      if ( prob_skip > (1. - 5e5/Nentries) ) continue;
      iEvent = jentry;
    }

    inChain->GetEntry(iEvent);
    if (jentry%1000 == 0) std::cout << "Event ["
      << jentry << "/" << Nevents << "]: reading event " << iEvent << std::endl;
    refMult = 0;
    MpdGlobalTracks = (TClonesArray*) MPDEvent->GetGlobalTracks();
    Int_t ntracks = MpdGlobalTracks->GetEntriesFast();
    for (int iTr=0; iTr<ntracks; iTr++)
    {
      auto mpdtrack = (MpdTrack*) MpdGlobalTracks->UncheckedAt(iTr);
#ifdef _OLD_MCSTACK_
      auto mctrack = (FairMCTrack*) MCTracks->UncheckedAt(mpdtrack->GetID());
#endif
#ifdef _NEW_MCSTACK_
      auto mctrack = (MpdMCTrack*) MCTracks->UncheckedAt(mpdtrack->GetID());
#endif
    
      pt = mpdtrack->GetPt();
      eta = mpdtrack->GetEta();
      nhits = mpdtrack->GetNofHits();

      if (pt < cut_pt) continue;
      if (TMath::Abs(eta) > cut_eta) continue;
      if (nhits < cut_nhits) continue;

      // Skip secondary particles if the flag primaryOnly is set
      //if (primaryOnly && mctrack->GetMotherId() != -1) continue;
      if (primaryOnly && TMath::Sqrt(TMath::Power(mpdtrack->GetDCAX(),2) + TMath::Power(mpdtrack->GetDCAY(),2) + TMath::Power(mpdtrack->GetDCAZ(),2)) > dca_cut) continue;

      refMult++;
    }
    hRefMult->Fill(refMult);
    hBvsRefMult->Fill(refMult,MCHeader->GetB());
  }
  TFile *fo = new TFile(outFileName.Data(),"recreate");
  fo->cd();
  hRefMult->Write();
  hBvsRefMult->Write();
  fo->Close();

  timer.Stop();
  timer.Print();
}

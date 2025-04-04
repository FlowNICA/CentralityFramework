/** @file   FitterHelper.h
    @author Viktor Klochkov (klochkov44@gmail.com)
    @author Ilya Selyuzhenkov (ilya.selyuzhenkov@gmail.com)
    @brief  Methods for fit QA
*/

#ifndef GlauberFitterHelper_H
#define GlauberFitterHelper_H 1

#include "TCanvas.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TLegend.h"
#include "TPad.h"

#include "Fitter.h"

namespace Glauber {
inline void DrawHistos(const Fitter &fit, Bool_t isSim, Bool_t isData,
                       Bool_t isGlauber, Bool_t isNBD) {
  std::unique_ptr<TCanvas> c1{new TCanvas("c1", "canvas", 1500, 900)};

  c1->Divide(2, 2);

  std::unique_ptr<TPad> c1_1{
      (TPad *)c1->GetListOfPrimitives()->FindObject("c1_1")};
  std::unique_ptr<TPad> c1_2{
      (TPad *)c1->GetListOfPrimitives()->FindObject("c1_2")};
  std::unique_ptr<TPad> c1_4{
      (TPad *)c1->GetListOfPrimitives()->FindObject("c1_4")};

  c1_1->SetLogy(1);
  c1_2->SetLogy(1);
  c1_4->SetLogy(1);

  /*const*/ TH1F hGlaub = fit.GetGlauberFitHisto();
  /*const*/ TH1F hGlPlp = fit.GetGlauberPlpHisto();
  /*const*/ TH1F hGlSng = fit.GetGlauberSngHisto();
  /*const*/ TH1F hData = fit.GetDataHisto();
  /*const*/ TH1F hNBD = fit.GetNBDHisto();
  /*const*/ TH1F hNcoll = fit.GetNcollHisto();
  /*const*/ TH1F hNpart = fit.GetNpartHisto();
  /*const*/ TH1F hBestFit = fit.GetBestFitHisto();
  /*const*/ TH1F hBestPlp = fit.GetBestPlpHisto();
  /*const*/ TH1F hBestSng = fit.GetBestSngHisto();

  TH2F hBestPlpEv1Ev2 = fit.GetBestPlpEv1Ev2();
  TH2F hBestB_VS_Multiplicity = fit.GetBestB_VS_Multiplicity();
  TH2F hBestNpart_VS_Multiplicity = fit.GetBestNpart_VS_Multiplicity();
  TH2F hBestNcoll_VS_Multiplicity = fit.GetBestNcoll_VS_Multiplicity();
  TH2F hBestEcc1_VS_Multiplicity = fit.GetBestEcc1_VS_Multiplicity();
  TH2F hBestPsi1_VS_Multiplicity = fit.GetBestPsi1_VS_Multiplicity();
  TH2F hBestEcc2_VS_Multiplicity = fit.GetBestEcc2_VS_Multiplicity();
  TH2F hBestPsi2_VS_Multiplicity = fit.GetBestPsi2_VS_Multiplicity();
  TH2F hBestEcc3_VS_Multiplicity = fit.GetBestEcc3_VS_Multiplicity();
  TH2F hBestPsi3_VS_Multiplicity = fit.GetBestPsi3_VS_Multiplicity();
  TH2F hBestEcc4_VS_Multiplicity = fit.GetBestEcc4_VS_Multiplicity();
  TH2F hBestPsi4_VS_Multiplicity = fit.GetBestPsi4_VS_Multiplicity();
  TH2F hBestEcc5_VS_Multiplicity = fit.GetBestEcc5_VS_Multiplicity();
  TH2F hBestPsi5_VS_Multiplicity = fit.GetBestPsi5_VS_Multiplicity();

  std::unique_ptr<TFile> fOut{TFile::Open("glauber_qa.root", "recreate")};

  if (isSim) {
    c1->cd(1);
    hNcoll.SetLineColor(2);

    hNcoll.Draw();
    hNpart.Draw("same");

    std::unique_ptr<TLegend> legSim{new TLegend(0.6, 0.75, 0.75, 0.83)};
    legSim->AddEntry(&hNpart, "Npart", "l");
    legSim->AddEntry(&hNcoll, "hNcoll", "l");
    legSim->Draw("same");

    hNcoll.Write();
    hNpart.Write();
  }

  if (isData) {
    c1->cd(2);
    hData.Draw();
    hData.Write();
    if (isGlauber) {
      hBestFit.SetLineColor(kRed);
      hBestFit.Draw("same");
      hBestPlp.SetLineColor(kBlack);
      hBestPlp.Draw("same");
      hBestSng.SetLineColor(kGreen + 1);
      hBestSng.Draw("same");

      std::unique_ptr<TLegend> legData{new TLegend(0.6, 0.75, 0.75, 0.83)};
      legData->AddEntry(&hBestFit, "Fit", "l");
      legData->AddEntry(&hBestPlp, "Pile-up", "l");
      legData->AddEntry(&hBestSng, "Single", "l");
      legData->AddEntry(&hData, "Data", "l");
      legData->Draw("same");
      hBestFit.Write();
      hBestPlp.Write();
      hBestSng.Write();
      hBestPlpEv1Ev2.Write();
      hBestB_VS_Multiplicity.Write();
      hBestNpart_VS_Multiplicity.Write();
      hBestNcoll_VS_Multiplicity.Write();
      hBestEcc1_VS_Multiplicity.Write();
      hBestEcc1_VS_Multiplicity.Write();
      hBestEcc2_VS_Multiplicity.Write();
      hBestEcc2_VS_Multiplicity.Write();
      hBestEcc3_VS_Multiplicity.Write();
      hBestEcc3_VS_Multiplicity.Write();
      hBestEcc4_VS_Multiplicity.Write();
      hBestEcc4_VS_Multiplicity.Write();
      hBestEcc5_VS_Multiplicity.Write();
      hBestEcc5_VS_Multiplicity.Write();
    }
  }

  if (isNBD) {
    c1->cd(3);
    hNBD.Draw();
    hNBD.Write();
  }

  if (isGlauber) {
    c1->cd(4);
    hBestFit.Draw();
    hBestPlp.Draw("same");
    hBestSng.Draw("same");
  }

  TTree *BestResult = new TTree("BestResult", "BestResult");
  Float_t mu, f, k, p, chi2_error, chi2;
  BestResult->Branch("mu", &mu);
  BestResult->Branch("f", &f);
  BestResult->Branch("k", &k);
  BestResult->Branch("p", &p);
  BestResult->Branch("chi2", &chi2);
  BestResult->Branch("chi2_error", &chi2_error);

  mu = fit.GetOptimalMu();
  f = fit.GetOptimalF();
  k = fit.GetOptimalK();
  chi2_error = fit.GetOptimalChi2Error();
  p = fit.GetOptimalP();
  chi2 = fit.GetOptimalChi2();

  BestResult->Fill();
  BestResult->Write();

  c1->Write();
  c1->SaveAs("glauber.pdf");
  fOut->Close();
}

} // namespace Glauber

#endif

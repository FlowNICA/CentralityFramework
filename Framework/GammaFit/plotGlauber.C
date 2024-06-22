#include "PlotFunc.C"

void plotGlauber() //Результаты фитирования множественности с ратио и классами центральности для Глаубера
{
  // gROOT->SetStyle("Pub");
  gROOT->ForceStyle();
  // gStyle->SetErrorX(0);

  const int MStyleData[3] = {25, 28, 27};
  const int MStyleFit[3] = {20, 20, 20};
  const int MStyleGl[3] = {22, 22, 22};

  const double MSizeData[3] = {1.4, 1.7, 2.2};
  const double MSizeFit[3] = {1.4, 1.4, 1.4};
  const double MSizeGl[3] = {1.4, 1.4, 1.4};

  Color_t MColor[3] = {1, kBlue + 1, kRed + 2};

  const char *model0 = "GSM";
  const char *modelT = "DCM-QGSM-SMM, GEANT4";
  const char *energ = "4.5";
  const char *sys = "AgAg";
  const char *systemT = Form("Ag+Ag, #sqrt{s_{NN}}=%s GeV", energ);
  const char *energ0 = "4";

  // glaub->Rebin(2);Ideal->Rebin(2);glaub->Rebin(2);

  TFile *file3 = new TFile("file:///home/dim/Work/FIT/dataNEW/FITGlaub/DCA_AgAg_4GeV/glauber_qa_dca_AgAg4.0GeV.root");
  TFile *fileFin = new TFile("file:///home/dim/Work/FIT/dataNEW/FITGlaub/DCA_AgAg_4GeV/FINAL_dca_AgAg4.0GeV.root");
  TFile *fileIn = new TFile("file:///home/dim/Work/FIT/dataNEW/GSM_AgAg4.root");
  // TFile *file3 = new TFile("/home/dim/Work/FIT/FIToutGamma/OUT/GSM610_5_fitGamma.root");

  TH1F *urqmd = (TH1F *)file3->Get("hRefMultSTAR");
  int Nn = 1.05 * urqmd->FindLastBinAbove();
  TH1F *hNEWfc2 = (TH1F *)file3->Get("glaub_fit_histo");
  TH2F *hBvsRefMultFit = (TH2F *)file3->Get("B_VS_Multiplicity_Histo");//B_VS_Multiplicity");
  TH2F *hBvsRefMult = (TH2F *)fileIn->Get("hBvsRefMult");

  TH1F *Data0_10 = (TH1F *)hNEWfc2->Clone();
  TH1F *Data10_40 = (TH1F *)hNEWfc2->Clone();
  TH1F *Data40_80 = (TH1F *)hNEWfc2->Clone();
  Data0_10->Scale(2);
  Data10_40->Scale(2);
  Data40_80->Scale(2);

  urqmd->Rebin(2);
  hNEWfc2->Rebin(2);
  TGraphErrors *Gurqmd, *NEWfc2;
  Gurqmd = HistGraph(urqmd);
  NEWfc2 = HistGraph(hNEWfc2);

  TTree *Result = (TTree *)fileFin->Get("Result");
  if (!Result)
    return;
  Int_t MinBorder;
  Int_t MaxBorder;
  Int_t fMinBorder;
  Int_t fMaxBorder;
  Float_t MinPercent;
  Float_t MaxPercent;
  Result->SetBranchAddress("MinBorder", &MinBorder);
  Result->SetBranchAddress("MaxBorder", &MaxBorder);
  Result->SetBranchAddress("MinPercent", &MinPercent);
  Result->SetBranchAddress("MaxPercent", &MaxPercent);
  Int_t const Nclasses = Result->GetEntries();
  Int_t binNch[Nclasses];
  Float_t GlaubV[Nclasses];

  TH1D *HFit[Nclasses];
  TH1D *HData[Nclasses];
  // Color_t color[11] = {kRed + 2, kBlue + 1, 14, kGreen + 3, kMagenta + 3, kGreen + 1, kYellow + 2, 46, kBlue - 9, kViolet + 8, kRed};
  Color_t color[10] = {kRed + 2, kBlue + 1, 14, kGreen + 3, kMagenta + 3, kGreen + 1, kYellow + 2, 46, kBlue - 9, kViolet + 8};

  int IntMinC, IntMaxC;
  double BmeanF[Nclasses], BsigmaF[Nclasses], BmeanErF[Nclasses], BsigmaErF[Nclasses];
  double BmeanD[Nclasses], BsigmaD[Nclasses], BmeanErD[Nclasses], BsigmaErD[Nclasses], centbin[Nclasses], centEr[Nclasses];
  for (int cent = 0; cent < Nclasses; cent++)
  {
    Result->GetEntry(cent);
    IntMaxC = (float)MaxPercent;
    IntMinC = (float)MinPercent;
    fMinBorder = MinBorder;
    fMaxBorder = MaxBorder;

    binNch[cent] = MinBorder;
    GlaubV[cent] = Gurqmd->Eval(binNch[cent]);

    // cout<<endl<<Form("%s%f_%f", "B_VS_CentralityClass_", MinPercent, MaxPercent);
    // HFit[cent] = (TH1F *)fileFin->Get(Form("B_VS_CentralityClass %0.1f%s-%0.1f%s", IntMinC,"%", IntMaxC,"%"));
    // cout<<endl<<Form("B_VS_CentralityClass %0.1f%s-%0.1f%s", IntMinC,"%", IntMaxC,"%");

    if (cent == Nclasses - 1)
    {
      Result->GetEntry(Nclasses - 2);
      fMaxBorder = MinBorder;
      fMinBorder = 1;
    }

    HFit[cent] = hBvsRefMultFit->ProjectionY(Form("%d_hist_Impact", cent), fMinBorder, fMaxBorder);
    HFit[cent]->Scale(1 / HFit[cent]->Integral(1, HFit[cent]->GetNbinsX(), "width"));
    HFit[cent]->SetLineColorAlpha(color[cent], 1);
    BmeanF[cent] = HFit[cent]->GetMean();
    BsigmaF[cent] = HFit[cent]->GetStdDev();
    BmeanErF[cent] = HFit[cent]->GetMeanError();
    BsigmaErF[cent] = HFit[cent]->GetStdDevError();

    HData[cent] = hBvsRefMult->ProjectionY(Form("%d_hist_Impact", cent), fMinBorder, fMaxBorder);
    cout << endl
         << "min=" << fMinBorder << " ,max=" << fMaxBorder;
    // HData[cent] = (TH1F *)file3->Get(Form("%s%d_%d_Data", "B_VS_CentralityClass_", IntMinC, IntMaxC));
    HData[cent]->Scale(1 / HData[cent]->Integral(1, HData[cent]->GetNbinsX(), "width"));
    HData[cent]->SetLineColorAlpha(color[cent], 1);
    BmeanD[cent] = HData[cent]->GetMean();
    BsigmaD[cent] = HData[cent]->GetStdDev();
    BmeanErD[cent] = HData[cent]->GetMeanError();
    BsigmaErD[cent] = HData[cent]->GetStdDevError();
    centbin[cent] = 0.5 * (IntMinC + IntMaxC);
    centEr[cent] = 0;
  }

  TGraphErrors *GrData = new TGraphErrors(Nclasses, centbin, BmeanD, centEr, BmeanErD);
  GrData->SetMarkerStyle(25);
  GrData->SetMarkerSize(1.5);
  // GrData->SetMarkerColor(4);
  // GrData->SetLineColor(4);
  GrData->GetXaxis()->SetTitle("Centrality, %");
  GrData->GetYaxis()->SetTitle("<b>, fm");
  GrData->SetTitle("");
  GrData->GetYaxis()->SetTitleOffset(1.55);
  GrData->GetYaxis()->SetRangeUser(-0.4, 18);
  GrData->GetXaxis()->SetLimits(0, 100);
  GrData->GetXaxis()->SetNdivisions(210);
  GrData->GetYaxis()->SetTitleSize(0.05);
  GrData->GetYaxis()->SetLabelSize(0.05);
  GrData->GetYaxis()->SetTitleOffset(1.3);
  GrData->GetYaxis()->SetLabelOffset(0.018);

  TGraphErrors *GrFit = new TGraphErrors(Nclasses, centbin, BmeanF, centEr, BmeanErF);
  GrFit->SetMarkerStyle(MStyleGl[0]);
  GrFit->SetMarkerSize(1.4);
  GrFit->SetMarkerColorAlpha(MColor[1], 1);
  GrFit->SetLineColorAlpha(MColor[1], 1);
  GrFit->SetLineWidth(2.);
  GrFit->SetTitle("");
  GrFit->GetXaxis()->SetTitle("Centrality, %");
  GrFit->GetYaxis()->SetTitle("<b>, fm");
  TGraphErrors *GrRatio = RatioGr(GrFit, GrData, 0, 0.85, 100, 1.15);
  GrRatio->SetTitle("");
  GrRatio->GetYaxis()->SetTitle("MC-Gl/Model");
  GrRatio->GetXaxis()->SetTitle("Centrality, %");
  GrRatio->GetXaxis()->SetTitleSize(0.1);
  GrRatio->GetXaxis()->SetLabelSize(0.1);
  GrRatio->GetYaxis()->SetTitleSize(0.1);
  GrRatio->GetYaxis()->SetLabelSize(0.1);
  GrRatio->GetYaxis()->SetTitleOffset(0.65);
  GrRatio->GetXaxis()->SetTitleOffset(1.5);
  GrRatio->GetYaxis()->SetLabelOffset(0.013);
  GrRatio->GetXaxis()->SetLabelOffset(0.04);
  GrRatio->GetXaxis()->SetNdivisions(210);
  GrRatio->GetYaxis()->SetNdivisions(505);
  GrRatio->GetXaxis()->SetTickLength(0.1);
  GrRatio->GetYaxis()->CenterTitle(true);
  TCanvas *CRatio2 = new TCanvas("CRatio_of_BMean2", "CRatio2", 800, 800);
  TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
  pad1->SetLeftMargin(0.13);
  pad1->SetBottomMargin(0);
  pad1->SetTopMargin(0.05);
  pad1->Draw();
  pad1->cd();
  GrData->Draw("AP");
  GrFit->Draw("SAME P");
  // TLegend *legdif = new TLegend(0.2, 0.7, 0.4, 0.9);
  TLegend *legdif = new TLegend(0.65, 0.15, 0.85, 0.3);
  legdif->SetBorderSize(0);
  legdif->AddEntry(GrData, "Model", "p");
  legdif->AddEntry(GrFit, "MC-Gl", "p");
  legdif->SetTextSize(0.04);
  legdif->Draw();
  TPaveText *ptB = new TPaveText(0.35, 0.65, 0.7, 0.9, "NDC NB");
  ptB->AddText(systemT);
  ptB->AddText(modelT);
 ptB->AddText("N_{hits} > 16, DCA > 0.5 cm");
  //ptB->AddText("Primary only");
  ptB->SetBorderSize(0);
  ptB->SetFillColor(0);
  ptB->SetTextSize(0.04);
  ptB->Draw();
  CRatio2->cd(); // Go back to the main canvas before defining pad2
  TPad *pad2 = new TPad("pad2", "pad2", 0, 0, 1, 0.3);
  pad2->SetLeftMargin(0.13);
  pad2->SetBottomMargin(0.3);
  pad2->SetTopMargin(0.0);
  pad2->Draw();
  pad2->cd();

  GrRatio->Draw("AP");
  TLine lineB;
  lineB.SetLineWidth(2);
  lineB.SetLineStyle(2);
  lineB.DrawLine(0, 1, 100, 1);
  CRatio2->SaveAs(Form("/home/dim/Work/FIT/dataNEW/FITplot/Glauber/GlauberFitmean_%s_%s_%s_DCA.png", model0, energ0, sys));
  //CRatio2->SaveAs(Form("/home/dim/Work/FIT/dataNEW/FITplot/GlauberFitmean_%s_%s_%sTEST.png", model0, energ0, sys));

  NEWfc2->SetMarkerStyle(MStyleGl[0]);
  NEWfc2->SetMarkerColor(MColor[1]);
  NEWfc2->SetMarkerSize(MSizeGl[0]);
  NEWfc2->SetLineColorAlpha(MColor[1], 1);
  NEWfc2->SetLineWidth(2);
  // NEWfc2->SetNpx(n1bins);
  Gurqmd->SetMarkerStyle(MStyleData[0]);
  Gurqmd->SetMarkerSize(MSizeData[0]);
  Gurqmd->SetLineWidth(2);
  Gurqmd->SetLineColorAlpha(1, 1);
  Gurqmd->GetXaxis()->SetLimits(0, Nn);
  Gurqmd->GetXaxis()->SetNdivisions(509);
  Gurqmd->GetYaxis()->SetRangeUser(0.5, 5000000);
  Gurqmd->SetTitle("");
  Gurqmd->GetXaxis()->SetTitle("N_{ch}");
  Gurqmd->GetYaxis()->SetTitle("dN/dN_{ch}");

  Gurqmd->GetYaxis()->SetTitleSize(25);
  Gurqmd->GetYaxis()->SetTitleFont(43);
  Gurqmd->GetYaxis()->SetTitleOffset(1.9);
  // Gurqmd->GetYaxis()->CenterTitle(true);
  Gurqmd->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
  Gurqmd->GetYaxis()->SetLabelSize(25);
  Gurqmd->GetYaxis()->SetLabelOffset(0.018);

  // Gurqmd->GetXaxis()->CenterTitle(true);
  Gurqmd->GetXaxis()->SetTitleSize(25);
  Gurqmd->GetXaxis()->SetTitleFont(43);
  Gurqmd->GetXaxis()->SetTitleOffset(2.65);
  Gurqmd->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
  Gurqmd->GetXaxis()->SetLabelSize(25);
  Gurqmd->GetXaxis()->SetLabelOffset(0.035);
  // Gurqmd->GetXaxis()->SetTickLength(0.08);

  TGraphErrors *RatioUrFit;
  RatioUrFit = RatioGr(NEWfc2, Gurqmd, 0, 0.65, Nn, 1.35);
  RatioUrFit->GetXaxis()->SetTitle("N_{ch}");
  RatioUrFit->GetYaxis()->SetTitle("MC-Gl/Model");
  RatioUrFit->SetTitle("");
  RatioUrFit->GetYaxis()->SetTitleSize(25);
  RatioUrFit->GetYaxis()->SetTitleFont(43);
  RatioUrFit->GetYaxis()->SetTitleOffset(2);
  RatioUrFit->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
  RatioUrFit->GetYaxis()->SetLabelSize(25);
  RatioUrFit->GetYaxis()->SetLabelOffset(0.03);
  RatioUrFit->GetYaxis()->CenterTitle(true);
  //  RatioUrFit->GetXaxis()->CenterTitle(true);
  RatioUrFit->GetYaxis()->SetNdivisions(504);
  RatioUrFit->GetXaxis()->SetNdivisions(509);
  RatioUrFit->GetXaxis()->SetTitleSize(25);
  RatioUrFit->GetXaxis()->SetTitleFont(43);
  RatioUrFit->GetXaxis()->SetTitleOffset(3.9);
  RatioUrFit->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
  RatioUrFit->GetXaxis()->SetLabelSize(25);
  RatioUrFit->GetXaxis()->SetLabelOffset(0.035);
  RatioUrFit->GetXaxis()->SetTickLength(0.08);
  /*RatioUrFit->SetMarkerStyle(MStyleFit[0]);
  RatioUrFit->SetMarkerSize(MSizeFit[0]);
  RatioUrFit->SetLineWidth(2);
  RatioUrFit->SetLineColorAlpha(MColor[2], 1);
  RatioUrFit->SetMarkerColor(MColor[2]);*/
  RatioUrFit->SetLineColorAlpha(MColor[1], 1);
  RatioUrFit->SetLineWidth(2);

  TCanvas *c = new TCanvas("Mean", "Mean", 900, 800);
  float ratY = 0.35, ratx = 1;
  TPad *pad[4];

  pad[0] = new TPad("pad22", "pad22", 0.0, 0, ratx, ratY);
  pad[1] = new TPad("pad22", "pad22", 0.0, ratY, ratx, 1);
  pad[2] = new TPad("pad22", "pad22", 0, 0, 1, ratY);
  pad[3] = new TPad("pad22", "pad22", 0, ratY, 1, 1);

  c->cd();
  pad[0]->SetBottomMargin(0.3);
  pad[0]->SetTopMargin(0.0);
  pad[0]->SetLeftMargin(0.12);
  pad[0]->SetRightMargin(0.03);
  pad[0]->Draw();
  pad[0]->cd();
  RatioUrFit->Draw("AP");
  TLine line;
  line.SetLineWidth(2);
  line.SetLineStyle(2);
  line.DrawLine(0, 1, Nn, 1);

  pad[1]->SetBottomMargin(0.0);
  pad[1]->SetTopMargin(0.05);
  pad[1]->SetLeftMargin(0.12);
  pad[1]->SetRightMargin(0.03);
  c->cd();
  pad[1]->Draw();
  pad[1]->cd();
  pad[1]->cd()->SetLogy();

  CutHisto(Data0_10, binNch[1], binNch[0]);
  Data0_10->SetFillColor(kWhite);
  // Data0_10->SetFillStyle(3144);
  Data0_10->SetLineColor(kGray + 1);

  CutHisto(Data10_40, binNch[7], binNch[1]);
  Data10_40->SetFillColor(kGray + 1);
  Data10_40->SetFillStyle(3003);
  Data10_40->SetLineColor(kGray + 1);

  CutHisto(Data40_80, binNch[11], binNch[7]);
  Data40_80->SetFillColor(kGray + 5);
  Data40_80->SetFillStyle(3002);
  Data40_80->SetLineColor(kGray + 1);

  Gurqmd->Draw("AP");

  Data0_10->Draw("SAME H");
  Data10_40->Draw("SAME H");
  Data40_80->Draw("SAME H");
  Gurqmd->Draw("SAME P");

  NEWfc2->Draw("SAME P");
  // TLegend *legdif1 = new TLegend(0.25, 0.7, 0.45, 0.91);
   TLegend *legdif1 = new TLegend(0.8, 0.6, 0.95, 0.91);
  legdif1->SetBorderSize(0);
  legdif1->AddEntry(Gurqmd, "Model", "p");
  legdif1->AddEntry(NEWfc2, "#Gamma-fit", "p");
  legdif1->AddEntry(Data0_10, "0-10%","f");
  legdif1->AddEntry(Data10_40, "10-40%", "f");
  legdif1->AddEntry(Data40_80, "40-80%", "f");
  legdif1->SetTextSize(0.04);
  legdif1->Draw();

  TPaveText *pt = new TPaveText(0.45, 0.65, 0.65, 0.91, "NDC NB"); // 0.56,0.72,0.89,0.89 right corner
  pt->SetBorderSize(0);
  pt->SetFillColor(0);
  pt->AddText(systemT);
  pt->AddText(modelT);
  // pt->AddText("N_{hits} > 16, DCA > 0.5 cm");
  pt->AddText("Primary only");
  pt->SetTextSize(0.04);
  pt->Draw();

line.SetLineWidth(2);
  line.SetLineStyle(1);

  // Read TTree
  TLatex latex;
  latex.SetTextAlign(12);
  latex.SetTextSize(0.04);
  latex.SetTextFont(42);
  latex.SetTextAngle(90);

  for (int i = 0; i < Nclasses; i++)
  {
    Result->GetEntry(i);
    line.DrawLine(binNch[i], 0, binNch[i], GlaubV[i]);
    IntMaxC = (int)MaxPercent;
    IntMinC = (int)MinPercent;

    if (i == 0)
    {
      latex.DrawLatex(MinBorder + 0.2 * (MaxBorder - MinBorder), 0.5* sqrt(GlaubV[0]), Form("%d-%d%s", IntMinC, IntMaxC, "%"));
    }

    else if (i > 0 && i < 5)
    {
      latex.DrawLatex(0.5 * (MinBorder + MaxBorder), 0.5 * sqrt(GlaubV[0]), Form("%d-%d%s", IntMinC, IntMaxC, "%"));
    }

    else if (i > 7 && i < 9)
    {
      latex.DrawLatex(0.5 * (MinBorder + MaxBorder), 0.5 * sqrt(GlaubV[0]), Form("%d-%d%s", IntMinC, IntMaxC, "%"));
    }
    if (i > 6)
    {
      line.SetLineStyle(2);
    }
  }



 c->SaveAs(Form("/home/dim/Work/FIT/dataNEW/FITplot/Glauber/GlauberFit_%s%s_%s_DCA.png", model0, energ0, sys));
  //c->SaveAs(Form("/home/dim/Work/FIT/dataNEW/FITplot/GlauberFit_%s%s_%s_TEST.png", model0, energ0, sys));
}

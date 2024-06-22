#include "PlotFunc.C"

void plotGamma2() //Результаты фитирования множественности с ратио и классами центральности Гамма фита
{
  // gROOT->SetStyle("Pub");
  gROOT->ForceStyle();
  // gStyle->SetErrorX(0);
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);

  const int MStyleData[3] = {25, 28, 27};
  const int MStyleFit[3] = {20, 20, 20};
  const int MStyleGl[3] = {22, 22, 22};

  const double MSizeData[3] = {1.3, 1.7, 2.2};
  const double MSizeFit[3] = {1.3, 1.4, 1.4};
  const double MSizeGl[3] = {1.3, 1.4, 1.4};

  Color_t MColor[3] = {1, kBlue + 1, kRed + 2};

   const char *model0 = "UrQMD_recoMPD_FXT_10bin";
  const char *modelT = "UrQMD+GEAN4";
  const char *energ = "2.5";
  const char *sys = "BiBi";
  const char *systemT = Form("Bi+Bi, #sqrt{s_{NN}}=%s GeV", energ);
  //const char *sys = "Xe+CsI";
  //const char *systemT = Form("Xe+CsI(protons), %s AGeV", energ);
  const char *energ0 = "4";
  const char *outDir = "~/Documents/Work/Dataset/QnToolsMPD/pics-cent/";

  TPaveText *ptB = new TPaveText(0.45, 0.7, 0.65, 0.92, "NDC NB");
  ptB->AddText(systemT);
  ptB->AddText(modelT);
  //ptB->AddText("N_{hits} > 4, |DCA_{x,y}| < 3 cm");
  //  ptB->AddText("Reconstructed");
  // ptB->AddText("Primary only");
  ptB->SetBorderSize(0);
  ptB->SetFillColor(0);
  ptB->SetTextSize(0.035);

  TFile *file3 = new TFile("~/Documents/Work/Dataset/QnToolsMPD/files-new2/cent_urqmd_bibi_2.5gev_fxt.root");

  TH1F *urqmd = (TH1F *)file3->Get("h1_refMult");//h_refmult_Event_Mult");//

  TH1F *hNEWfc2 = (TH1F *)file3->Get("fit_hist");
  int Nn = (int)(1.05 * hNEWfc2->FindLastBinAbove() / 20) * 20 + 10;

  TTree *ResultGamma = (TTree *)file3->Get("Result");
  Int_t MinBorderGamma;
  Int_t MaxBorderGamma;
  Float_t MinPercent;
  Float_t MaxPercent;
  ResultGamma->SetBranchAddress("MinBorder", &MinBorderGamma);
  ResultGamma->SetBranchAddress("MaxBorder", &MaxBorderGamma);
  ResultGamma->SetBranchAddress("MinPercent", &MinPercent);
  ResultGamma->SetBranchAddress("MaxPercent", &MaxPercent);
  Int_t const NclassesGamma = ResultGamma->GetEntries();

  TTree *FitResult = (TTree *)file3->Get("FitResult");
  if (!FitResult)
    return;

  Float_t a1, a2, a3, teta, n_knee, chi2, NDF;
  Int_t minNch;
  FitResult->SetBranchAddress("teta", &teta);
  FitResult->SetBranchAddress("n_knee", &n_knee);
  FitResult->SetBranchAddress("a1", &a1);
  FitResult->SetBranchAddress("a2", &a2);
  FitResult->SetBranchAddress("a3", &a3);
  FitResult->SetBranchAddress("chi2", &chi2);
  FitResult->SetBranchAddress("NDF", &NDF);
  FitResult->SetBranchAddress("minNch", &minNch);
  FitResult->GetEntry(0);

  Int_t binNchGamma[NclassesGamma];
  Int_t binNchMax[NclassesGamma];
  Float_t GammaV[NclassesGamma];

  TH1F *HFit[NclassesGamma];
  TH1F *HData[NclassesGamma];
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
  // Color_t color[11] = {kRed + 2, kBlue + 1, 14, kGreen + 3, kMagenta + 3, kGreen + 1, kYellow + 2, 46, kBlue - 9, kViolet + 8, kRed};

  int IntMinC, IntMaxC;
  double BmeanF[NclassesGamma], BsigmaF[NclassesGamma], BmeanErF[NclassesGamma], BsigmaErF[NclassesGamma];
  double BmeanD[NclassesGamma], BsigmaD[NclassesGamma], BmeanErD[NclassesGamma], BsigmaErD[NclassesGamma], centbin[NclassesGamma], centEr[NclassesGamma];
  for (int cent = 0; cent < NclassesGamma; cent++)
  {
    ResultGamma->GetEntry(cent);
    IntMaxC = (int)MaxPercent;
    IntMinC = (int)MinPercent;
    binNchGamma[cent] = MinBorderGamma;
    binNchMax[cent] = MaxBorderGamma;
    GammaV[cent] = Gurqmd->Eval(binNchGamma[cent]);
    // cout<<endl<<Form("%s%f_%f", "B_VS_CentralityClass_", MinPercent, MaxPercent);
    HFit[cent] = (TH1F *)file3->Get(Form("%s%d_%d", "B_VS_CentralityClass_", IntMinC, IntMaxC));
    HFit[cent]->Scale(1 / HFit[cent]->Integral(1, HFit[cent]->GetNbinsX(), "width"));
    HFit[cent]->SetLineColorAlpha(kRed + 2, 1);
    // HFit[cent]->SetFillColor(kRed + 2 + 1);
    // HFit[cent]->SetFillStyle(3003);
    HFit[cent]->SetLineWidth(2.);

    BmeanF[cent] = HFit[cent]->GetMean();
    BsigmaF[cent] = HFit[cent]->GetStdDev();
    BmeanErF[cent] = HFit[cent]->GetMeanError();
    BsigmaErF[cent] = HFit[cent]->GetStdDevError();

    HData[cent] = (TH1F *)file3->Get(Form("%s%d_%d_Data", "B_VS_CentralityClass_", IntMinC, IntMaxC));
    HData[cent]->Scale(1 / HData[cent]->Integral(1, HData[cent]->GetNbinsX(), "width"));
    HData[cent]->SetLineColorAlpha(1, 1);
    HData[cent]->SetMarkerStyle(25);
    HData[cent]->SetFillColor(kGray + 1);
    HData[cent]->SetFillStyle(3003);
    HData[cent]->SetLineWidth(2.);
    BmeanD[cent] = HData[cent]->GetMean();
    BsigmaD[cent] = HData[cent]->GetStdDev();
    BmeanErD[cent] = HData[cent]->GetMeanError();
    BsigmaErD[cent] = HData[cent]->GetStdDevError();
    centbin[cent] = 0.5 * (IntMinC + IntMaxC);
    centEr[cent] = 0;
  }

  /*TCanvas *cimp = new TCanvas("cimp", "cimp", 900, 800);
  HData[0]->GetYaxis()->SetRangeUser(0., 1.2);
  HData[0]->GetXaxis()->SetRangeUser(0., 13.);
  HData[0]->GetYaxis()->SetTitle("P(b)");
  HData[0]->GetXaxis()->SetTitle("b, fm");
  // HData[0]->Draw("P");
  HData[0]->Draw("Hist");

  HFit[0]->Draw("HIST SAME");

  HData[3]->Draw("HistSame");
  // HData[3]->Draw("PSame");
  HFit[3]->Draw("HIST SAME");

  // HData[7]->Draw("P SAME");
  HData[7]->Draw("Hist SAME");
  HFit[7]->Draw("HIST SAME");

  // HData[10]->Draw("P SAME");
  HData[10]->Draw("Hist SAME");
  HFit[10]->Draw("HIST SAME");

  TLegend *leg = new TLegend(0.15, 0.75, 0.35, 0.85);
  // leg->SetHeader(systemT,"C");
  leg->SetBorderSize(0);
  leg->SetTextSize(0.025);
  leg->AddEntry(HData[0], "Model", "f");
  leg->AddEntry(HFit[10], "#Gamma-fit", "l");
  leg->Draw();

  TLatex text;
  text.SetTextAlign(12);
  text.SetTextSize(0.028);
  text.SetTextFont(62);
  text.DrawLatex(1.8, 0.5, "0-5%");
  text.DrawLatex(4, 0.8, "15-20%");
  text.DrawLatex(6.5, 0.9, "35-40%");
  text.DrawLatex(9., 0.9, "60-70%");

  text.DrawLatex(4.5, 1.13, systemT);
  // text.DrawLatex(4.1,1.05,Form("%s, GEANT4",modelT));
  text.DrawLatex(4.5, 1.05, Form("%s", modelT));
  cimp->SaveAs(Form("/home/dim/Work/FIT/dataNEW/FITGlaub/gen_AgAg9/GammaImpactDis_%s_%s_%s_NewBin.png", model0, energ0, sys));*/

  TGraphErrors *GrData = new TGraphErrors(NclassesGamma, centbin, BmeanD, centEr, BsigmaD);
  GrData->SetMarkerStyle(25);
  GrData->SetMarkerSize(1.3);
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

  TGraphErrors *GrFit = new TGraphErrors(NclassesGamma, centbin, BmeanF, centEr, BsigmaF);
  GrFit->SetMarkerStyle(20);
  GrFit->SetMarkerSize(1.3);
  GrFit->SetMarkerColor(MColor[2]);
  GrFit->SetMarkerColorAlpha(MColor[2], 1);
  GrFit->SetLineWidth(2.);
  GrFit->SetTitle("");
  GrFit->GetXaxis()->SetTitle("Centrality, %");
  GrFit->GetYaxis()->SetTitle("<b>, fm");
  TGraphErrors *GrRatio = RatioGr(GrFit, GrData, 0, 0.85, 100, 1.15);
  GrRatio->SetTitle("");
  GrRatio->GetYaxis()->SetTitle("#Gamma-fit/Model");
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
  GrRatio->SetLineColor(1);
  TCanvas *CRatio2 = new TCanvas("CRatio_of_BMean2", "CRatio2", 800, 800);
  TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
  pad1->SetLeftMargin(0.13);
  pad1->SetBottomMargin(0);
  pad1->SetTopMargin(0.05);
  pad1->SetRightMargin(0.03);
  pad1->Draw();
  pad1->cd();
  GrData->Draw("AP");
  GrFit->Draw("SAME P");
  // TLegend *legdif = new TLegend(0.2, 0.7, 0.4, 0.9);
  TLegend *legdif = new TLegend(0.75, 0.15, 0.9, 0.3);
  legdif->SetBorderSize(0);
  legdif->AddEntry(GrData, "Model", "p");
  legdif->AddEntry(GrFit, "#Gamma-Fit", "p");
  legdif->SetTextSize(0.035);
  legdif->Draw();

  ptB->Draw();
  CRatio2->cd(); // Go back to the main canvas before defining pad2
  TPad *pad2 = new TPad("pad2", "pad2", 0, 0, 1, 0.3);
  pad2->SetLeftMargin(0.13);
  pad2->SetBottomMargin(0.3);
  pad2->SetTopMargin(0.0);
  pad2->SetRightMargin(0.03);
  pad2->Draw();
  pad2->cd();

  GrRatio->Draw("AP");
  TLine lineB;
  lineB.SetLineWidth(2);
  lineB.SetLineStyle(2);
  lineB.DrawLine(0, 1, 100, 1);
  // CRatio2->SaveAs(Form("/home/dim/Work/FIT/dataNEW/FITplot/GammaFitmean_%s_%s_%s_MotherID.png", model0, energ0, sys));
  CRatio2->SaveAs(Form("%sGammaFitmean_%s_%s_%s_NewBin.png", outDir, model0, energ0, sys));

  NEWfc2->SetMarkerStyle(MStyleFit[0]);
  NEWfc2->SetMarkerColor(MColor[2]);
  NEWfc2->SetMarkerSize(MSizeFit[0]);
  NEWfc2->SetLineColorAlpha(MColor[2], 1);
  // NEWfc2->SetNpx(n1bins);
  Gurqmd->SetMarkerStyle(MStyleData[0]);
  Gurqmd->SetMarkerSize(MSizeData[0]);
  Gurqmd->SetLineWidth(2);
  Gurqmd->SetLineColorAlpha(1, 1);
  Gurqmd->GetXaxis()->SetLimits(0, Nn);
  Gurqmd->GetXaxis()->SetNdivisions(509);
  Gurqmd->GetYaxis()->SetRangeUser(0.5, 2e8);
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
  RatioUrFit = RatioGr(NEWfc2, Gurqmd, 0, 0.65, Nn, 1.35, true);
  RatioUrFit->GetXaxis()->SetTitle("N_{ch}");
  RatioUrFit->GetYaxis()->SetTitle("#Gamma-fit/Model");
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
  RatioUrFit->GetXaxis()->SetTitleOffset(1.4);
  RatioUrFit->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
  RatioUrFit->GetXaxis()->SetLabelSize(25);
  RatioUrFit->GetXaxis()->SetLabelOffset(0.035);
  RatioUrFit->GetXaxis()->SetTickLength(0.08);
  RatioUrFit->SetMarkerStyle(MStyleFit[0]);
  RatioUrFit->SetMarkerSize(MSizeFit[0]);
  RatioUrFit->SetLineWidth(2);
  RatioUrFit->SetLineColorAlpha(MColor[2], 1);
  RatioUrFit->SetMarkerColor(MColor[2]);

  TCanvas *c = new TCanvas("Mean", "Mean", 800, 800);
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

  CutHisto(Data0_10, binNchGamma[0], binNchMax[1]);
  Data0_10->SetFillColor(kWhite);
  // Data0_10->SetFillStyle(3144);
  Data0_10->SetLineColor(kGray + 1);

  CutHisto(Data10_40, binNchGamma[7], binNchGamma[1]);
  Data10_40->SetFillColor(kGray + 1);
  Data10_40->SetFillStyle(3003);
  Data10_40->SetLineColor(kGray + 1);

  CutHisto(Data40_80, binNchGamma[11], binNchGamma[7]);
  Data40_80->SetFillColor(kGray + 5);
  Data40_80->SetFillStyle(3002);
  Data40_80->SetLineColor(kGray + 1);

  Gurqmd->Draw("AP");
  Data0_10->Draw("SAME H");
  Data10_40->Draw("SAME H");
  Data40_80->Draw("SAME H");
  Gurqmd->Draw("SAME P");
  NEWfc2->SetLineColorAlpha(kRed, 1);
  NEWfc2->Draw("SAME P");

  TLatex FitResultText0;
  FitResultText0.SetTextAlign(12);
  FitResultText0.SetTextSize(0.035);
  FitResultText0.SetTextFont(62);
  FitResultText0.DrawLatex(0.05 * Nn, 5e7, Form("#chi^{2}/ndf=%0.2f", chi2 / NDF));
  TLatex FitResultText;
  FitResultText.SetTextAlign(12);
  FitResultText.SetTextSize(0.035);
  FitResultText.SetTextFont(42);
  // FitResultText.DrawLatex(0.05 * Nn, 5e7, Form("#chi^{2}/ndf=%0.2f", chi2 / NDF));
  FitResultText.DrawLatex(0.05 * Nn, 2e7, Form("#theta=%1.1f, N_{knee}=%0.0f", teta, n_knee));
  FitResultText.DrawLatex(0.05 * Nn, 8e6, Form("a_{1}=%0.1f, a_{2}=%0.1f, a_{3}=%0.1f", a1, a2, a3));
  FitResultText.DrawLatex(0.05 * Nn, 2e6, "#frac{N_{knee}}{#theta}#times#sum_{i=1}^{3}a_{i}c_{b}^{i}");

  // TLegend *legdif1 = new TLegend(0.25, 0.7, 0.45, 0.91);
  TLegend *legdif1 = new TLegend(0.8, 0.6, 0.95, 0.91);
  legdif1->SetBorderSize(0);
  legdif1->AddEntry(Gurqmd, "Model", "p");
  legdif1->AddEntry(NEWfc2, "#Gamma-fit", "p");
  legdif1->AddEntry(Data0_10, "0-10%", "f");
  legdif1->AddEntry(Data10_40, "10-40%", "f");
  legdif1->AddEntry(Data40_80, "40-80%", "f");
  legdif1->SetTextSize(0.03);
  legdif1->Draw();

  /*TPaveText *pt = new TPaveText(0.45, 0.75, 0.65, 0.91, "NDC NB"); // 0.56,0.72,0.89,0.89 right corner
  pt->SetBorderSize(0);
  pt->SetFillColor(0);
  pt->AddText(systemT);
  pt->AddText(modelT);
  pt->AddText("N_{hits} > 16, |DCA| < 0.5 cm");
  //  pt->AddText("Primary only");
  pt->SetTextSize(0.035);*/
  ptB->Draw();

  line.SetLineWidth(2);
  line.SetLineStyle(1);

  // Read TTree
  TLatex latex;
  latex.SetTextAlign(12);
  latex.SetTextSize(0.04);
  latex.SetTextFont(42);
  latex.SetTextAngle(90);

  for (int i = 0; i < NclassesGamma; i++)
  {
    ResultGamma->GetEntry(i);
    line.DrawLine(binNchGamma[i], 0, binNchGamma[i], GammaV[i]);
    IntMaxC = (int)MaxPercent;
    IntMinC = (int)MinPercent;

    if (i == 0)
    {
      latex.DrawLatex(MinBorderGamma + 0.2 * (MaxBorderGamma - MinBorderGamma), 0.5 * sqrt(GammaV[0]), Form("%d-%d%s", IntMinC, IntMaxC, "%"));
    }

    else if (i > 0 && i < 5)
    {
      latex.DrawLatex(0.5 * (MinBorderGamma + MaxBorderGamma), 0.5 * sqrt(GammaV[0]), Form("%d-%d%s", IntMinC, IntMaxC, "%"));
    }

    else if (i > 4 && i < 6)
    {
      latex.DrawLatex(0.5 * (MinBorderGamma + MaxBorderGamma), 0.5 * sqrt(GammaV[0]), Form("%d-%d%s", IntMinC, IntMaxC, "%"));
    }
    if (i > 6)
    {
      line.SetLineStyle(2);
    }
  }

  /*TPaveText *TextCent = new TPaveText(0.6, 0.15, 0.7, 0.25, "NDC NB"); // 0.56,0.72,0.89,0.89 right corner
  // TPaveText *TextCent = new TPaveText(0.5, 0.15, 0.65, 0.25, "NDC NB");
  TextCent->SetBorderSize(0);
  TextCent->SetFillColor(0);
  TextCent->AddText("0-10%");
  TextCent->SetTextAngle(90);
  TextCent->SetTextSize(0.04);
  TextCent->Draw();*/

  // c->SaveAs(Form("/home/dim/Work/FIT/dataNEW/FITplot/GammaFit_%s%s_%s_MotherId.png", model0, energ0, sys));
  c->SaveAs(Form("%sGammaFit_%s%s_%s_NewBin.png", outDir, model0, energ0, sys));
}

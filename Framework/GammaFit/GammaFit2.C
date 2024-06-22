#include <TH2.h>
#include <TTree.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TROOT.h>
#include <TF1.h>
#include <TLine.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TGraph.h>
#include <Math/PdfFuncMathCore.h>
#include <Math/IntegratorOptions.h>
#include <TLegend.h>
#include <TPaveText.h>

bool GetSigma = true;
Float_t sigma = 607.084;
Float_t pi = TMath::Pi(), bmax = 18;
Color_t color[10] = {kRed + 2, kBlue + 1, 14, kGreen + 3, kMagenta + 3, kGreen + 1, kYellow + 2, 46, kBlue - 9, kViolet + 8};
//  Float_t teta = 1.367447, n_knee = 140.194, a1 = -2.54682, a2 = 0.993259, a3 = -2.94017, chi2, NDF, chi2_NDF;
 //Float_t teta = 1, n_knee = 70.194, a1 = -4, a2 = 2.4, a3 = -3.94017, chi2, NDF, chi2_NDF;
// Au+Au
// Float_t teta = 0.6, n_knee = 137, a1 = -2.4, a2 = 0.58, a3 = -2.2, chi2, NDF, chi2_NDF; // 2.4, 2.7 GeV
Float_t teta = 1.1, n_knee = 220, a1 = -5., a2 = 2., a3 = -5., chi2, NDF, chi2_NDF;
// Float_t teta = 1.1, n_knee = 900, a1 = -4., a2 = -4., a3 = 5., chi2, NDF, chi2_NDF;
// Xe+Cs
//Float_t teta = 0.6, n_knee = 100, a1 = -2.4, a2 = 0.58, a3 = -2.2, chi2, NDF, chi2_NDF;
const int CentClass = 14;
int bin_cent[CentClass];
int range_cent[CentClass] = {0, 5, 10, 15, 20, 25, 30, 35, 40, 50, 60, 70, 80, 100};

// const int CentClass = 11;
// int bin_cent[CentClass];
// int range_cent[CentClass] = {0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100};

TFile *d_outfile;

double Scale(TH1D *HistTrue, TH1D *Hist)
{
    TF1 *f1 = new TF1("f1", "gaus", 0.75 * HistTrue->FindLastBinAbove(), 0.98 * HistTrue->FindLastBinAbove());
    HistTrue->Fit(f1, "RMN");
    double MaxX1 = f1->GetParameter(1) + 3 * f1->GetParameter(2);

    TF1 *f2 = new TF1("f2", "gaus", 0.75 * Hist->FindLastBinAbove(), 0.98 * Hist->FindLastBinAbove());
    Hist->Fit(f2, "RMN");
    double MaxX2 = f2->GetParameter(1) + 3 * f2->GetParameter(2);
    cout << "MaxX2 / MaxX1 " << MaxX2 / MaxX1 << endl;

    return MaxX2 / MaxX1;
}

TGraphErrors *RatioGr(TGraphErrors *const &gr1, TGraphErrors *const &gr2, double Xmin, double Ymin, double Xmax, double Ymax)
{
    // Read points
    Double_t *vx_gr1 = gr1->GetX();
    Double_t *vy_gr1 = gr1->GetY();
    Double_t *vx_gr2 = gr2->GetX();
    Double_t *vy_gr2 = gr2->GetY();
    // Read errors
    Double_t *ex_gr1 = gr1->GetEX();
    Double_t *ey_gr1 = gr1->GetEY();
    Double_t *ex_gr2 = gr2->GetEX();
    Double_t *ey_gr2 = gr2->GetEY();

    int n1bins = gr1->GetN();
    int n2bins = gr2->GetN();

    Color_t GrColor = gr1->GetMarkerColor();
    Style_t GrStyle = gr1->GetMarkerStyle();
    Size_t GrSize = gr1->GetMarkerSize();

    if (n2bins < n1bins)
    {
        n1bins = n2bins;
    }
    Double_t vx_gr3[n1bins], vy_gr3[n1bins], ex_gr3[n1bins], ey_gr3[n1bins];
    for (int i = 0; i < n1bins; i++)
    {
        vx_gr3[i] = vx_gr1[i];
        ex_gr3[i] = ex_gr1[i];
        vy_gr3[i] = vy_gr1[i] / vy_gr2[i];
        ey_gr3[i] = sqrt(pow(ey_gr1[i] / vy_gr2[i], 2) + pow(vy_gr1[i] * ey_gr2[i] / (vy_gr2[i] * vy_gr2[i]), 2));
    }

    TGraphErrors *grRatio = new TGraphErrors(n1bins, vx_gr3, vy_gr3, ex_gr3, ey_gr3);
    grRatio->GetXaxis()->SetLimits(Xmin, Xmax);
    grRatio->GetYaxis()->SetRangeUser(Ymin, Ymax);
    grRatio->SetMarkerColor(GrColor);
    grRatio->SetMarkerStyle(GrStyle);
    grRatio->SetMarkerSize(GrSize);

    return grRatio;
}

//Уссловная вероятность;
double PnbGamma(double *x, double *par)
{
    double Cb = x[0], n = par[0], teta = par[1], n_knee = par[2], a1 = par[3], a2 = par[4], a3 = par[5];
    double fn = n_knee * exp(a1 * Cb + a2 * pow(Cb, 2) + a3 * pow(Cb, 3)) / teta;
    return ROOT::Math::gamma_pdf(n, fn, teta);
};

double PnbGamma2(double *x, double *par)
{
    double n = x[0], Cb = par[0], teta = par[1], n_knee = par[2], a1 = par[3], a2 = par[4], a3 = par[5];
    double fn = n_knee * exp(a1 * Cb + a2 * pow(Cb, 2) + a3 * pow(Cb, 3)) / teta;
    return ROOT::Math::gamma_pdf(n, fn, teta);
};

double ftPn(double *x, double *par)
{
    // Parameters
    double teta = par[0];
    double n_knee = par[1];
    double a1 = par[2];
    double a2 = par[3];
    double a3 = par[4];

    // Variables
    double n = x[0];
    // Function
    TF1 *f = new TF1("f", PnbGamma, 0, 10, 6);
    f->SetParameters(n, teta, n_knee, a1, a2, a3);
    double func = f->Integral(0, 1);
    return func;
}
double PbFit(double *x, double *par)
{
    double b = x[0], sigma = par[0], b0 = par[1], db = par[2], alpha = par[3];
    double PB = 2 * 3.1415926 * b / sigma;
    double var = (b - b0) / db;
    double Pinel = 1 / pow(1 + exp(var), alpha);
    return PB * Pinel;
};

double GetSigmaFit(const char *fileadres = "", const char *current_B_vs_mult = "hBvsRefMult")
{
    TFile *file = new TFile(fileadres);
    TH2D *BvsN;
    TH1D *Hist;
    BvsN = (TH2D *)file->Get(current_B_vs_mult);
    int NbinY = BvsN->GetNbinsY();
    int NbinX = BvsN->GetNbinsX();
    Hist = BvsN->ProjectionY(Form("hist_mult_%d", 1), 1, NbinX);
    Hist->Rebin(2);
    Hist->Scale(1 / Hist->Integral(1, Hist->FindLastBinAbove(), "width"));
    // Hist->Scale(1 /1.02853);

    // TF1 *fPbFit = new TF1("fPbFit", PbFit, 0.2, 15.5, 4);
    TF1 *fPbFit = new TF1("fPbFit", PbFit, 0.2, 15, 4);
    fPbFit->SetParameters(370, 11, 0.3, 1);
    Hist->Fit(fPbFit, "RMQB");
    Hist->Fit(fPbFit, "RM");
    float sig = fPbFit->GetParameter(0);
    // file->Close();
    return sig;
}

void Start(const char *fileadres, const char *current_mult, const char *outadres, int minNch, const char *current_B_vs_mult , bool efficiencyFit, const char *fileadres2, const char *current_mult2)
{
    if (GetSigma)
    {
        sigma = GetSigmaFit(fileadres,current_B_vs_mult);
    }
    int Nch0 = minNch;
    TFile *file = new TFile(fileadres);

    TH1D *hRefmult = (TH1D *)file->Get(current_mult);
    // file->Close();
    //TH1D *Gev = (TH1D *)hRefmult->Clone();
    bin_cent[0] = hRefmult->FindLastBinAbove();
    TH1D *Gev= new TH1D("RefMult","RefMult",bin_cent[0],0,bin_cent[0]);
    for(int i=0;i<bin_cent[0];i++)
    {
        float v=hRefmult->GetBinContent(i);
        Gev->SetBinContent(i,v);
    }
    Gev->Scale(1 / Gev->Integral(1, Gev->GetNbinsX(), "width"));
   

    if (efficiencyFit == true)
    {
        int mediumNn = 0.2 * Gev->FindLastBinAbove();
        TFile *fileGl = new TFile(fileadres2);
        TH1D *Ideal = (TH1D *)fileGl->Get(current_mult2);
        Ideal->Scale(1 / Ideal->Integral(1, Ideal->GetNbinsX(), "width"));
        float Integr = Ideal->Integral(mediumNn, Ideal->GetNbinsX(), "width");
        int EfmediumNn = mediumNn * Scale(Ideal, Gev);
        Gev->Scale(Integr / Gev->Integral(EfmediumNn, Gev->GetNbinsX(), "width"));
        // cout << "start3" << endl;
    }

    Gev->SetTitle("");
    Gev->GetYaxis()->SetTitle("1/N dN_{ch}/dN");
    Gev->GetXaxis()->SetTitle("N_{ch}");
    Gev->GetXaxis()->SetRangeUser(0, bin_cent[0]);
    Gev->GetYaxis()->SetRangeUser(0.5 * 1e-6, 0.5);

    // Y axis plot settings
    Gev->GetYaxis()->SetTitleSize(14);
    Gev->GetYaxis()->SetLabelSize(14);
    Gev->GetYaxis()->SetTitleFont(43);
    Gev->GetYaxis()->SetLabelFont(43);
    Gev->GetYaxis()->SetTitleOffset(1.5);
    Gev->GetYaxis()->SetLabelOffset(0.006);
    // X axis plot settings
    Gev->GetXaxis()->SetTitleSize(14);
    Gev->GetXaxis()->SetLabelSize(14);
    Gev->GetXaxis()->SetTitleFont(43);
    Gev->GetXaxis()->SetLabelFont(43);
    Gev->GetXaxis()->SetTitleOffset(3.6);
    Gev->GetXaxis()->SetLabelOffset(0.025);

    TH1F *GevC = (TH1F *)Gev->Clone();
    GevC->GetYaxis()->SetTitle("Data/Fit");
    GevC->GetYaxis()->CenterTitle(true);
    GevC->GetXaxis()->SetTitle("N_{ch}");
    GevC->GetXaxis()->SetTickLength(0.08);
    GevC->GetYaxis()->SetLabelOffset(0.015);

    //Задаем функцию для фитирования Gaus и строим верхнию половину графика
    TF1 *fc22 = new TF1("fit_func_full", ftPn, 0.1, bin_cent[0], 5);
    TF1 *fc2 = new TF1("fit_func", ftPn, Nch0, bin_cent[0], 5);
     fc2->SetParameters(teta, n_knee, a1, a2, a3);
   // fc2->SetParameters(teta, 0.9 * bin_cent[0], a1, a2, a3);
    /*fc2->SetParLimits(0,0.1,3);
    fc2->SetParLimits(1,10,bin_cent[0]);
    fc2->SetParLimits(2,-10,10);
    fc2->SetParLimits(3,-10,10);
    fc2->SetParLimits(4,-10,10);*/


    TCanvas *c2 = new TCanvas("Canvas_0f_fit_result", "ratio data/fit", 650, 500);
    TPad *pad12 = new TPad("pad12", "pad12", 0, 0.3, 1, 1.0);
    pad12->SetBottomMargin(0.0);
    pad12->Draw();
    pad12->cd()->SetLogy();
    auto fit_result1 = (TFitResultPtr) Gev->Fit(fc2, "RLQS");
    if (!fit_result1) std::cerr << "First approach fit results an error!" << std::endl;
    auto fit_result2 = (TFitResultPtr) Gev->Fit(fc2, "RMQS");
    if (!fit_result2) std::cerr << "Second approach fit results an error!" << std::endl;
    auto fit_result3 = (TFitResultPtr) Gev->Fit(fc2, "RMS");
    if (!fit_result3) std::cerr << "Second approach fit results an error!" << std::endl;
    Gev->Draw();

    teta = fc2->GetParameter(0);
    n_knee = fc2->GetParameter(1);
    a1 = fc2->GetParameter(2);
    a2 = fc2->GetParameter(3);
    a3 = fc2->GetParameter(4);
    chi2 = fc2->GetChisquare();
    NDF = fc2->GetNDF();
    chi2_NDF = chi2 / NDF;

    fc22->SetParameters(teta, n_knee, a1, a2, a3);
    fc22->SetNpx(2000);
    fc22->Draw("SAME");
    cout << endl
         << "Fit results" << endl;
    cout << "teta=" << teta << ", n_knee=" << n_knee << ", a1=" << a1 << " ,a2=" << a2 << ", a3=" << a3 << ", Sigma=" << sigma << endl
         << ", chi^2/NDF=" << chi2_NDF << endl
         << endl;
    TLine *line2 = new TLine(n_knee, 0.5 * 1e-6, n_knee, 0.5);
    line2->Draw("SAME");

    TLegend *legdif = new TLegend(0.18, 0.15, 0.33, 0.35);
    legdif->SetBorderSize(0);
    legdif->AddEntry(Gev, "Data", "pl");
    legdif->AddEntry(fc2, "Fit", "pl");
    legdif->Draw();
    c2->cd(); // Go back to the main canvas before defining pad2
    TPad *pad22 = new TPad("pad22", "pad22", 0, 0, 1, 0.3);
    pad22->SetBottomMargin(0.3);
    pad22->SetTopMargin(0.0);
    pad22->Draw();
    pad22->cd();

    // TF1 *NEWfc2 = new TF1;
    // NEWfc2 = Gev->GetFunction("fit_func");
    GevC->Divide(fc22);
    GevC->SetMinimum(0.45); // Define Y ..
    GevC->SetMaximum(1.55); // .. range
    GevC->SetStats(0);      // No statistics on lower plot
    GevC->Draw();
    TLine *line = new TLine(0, 1, bin_cent[0], 1);
    line->Draw("SAME");
    d_outfile = new TFile(outadres, "recreate");
    c2->Write();
    fc2->Write();
    fc22->Write();

    TTree *treeFit = new TTree("FitResult", "FitResult");
    treeFit->Branch("teta", &teta, "teta/F");
    treeFit->Branch("n_knee", &n_knee, "n_knee/F");
    treeFit->Branch("a1", &a1, "a1/F");
    treeFit->Branch("a2", &a2, "a2/F");
    treeFit->Branch("a3", &a3, "a3/F");
    treeFit->Branch("chi2", &chi2, "chi2/F");
    treeFit->Branch("NDF", &NDF, "NDF/F");
    treeFit->Branch("minNch", &minNch, "minNch/I");
    treeFit->Branch("sigma", &sigma, "sigma/F");
    treeFit->Fill();

    treeFit->Write();

    TH1F *hFitResult = new TH1F("fit_hist", "fit_hist", 2000, 0, 2000);
    hFitResult->FillRandom("fit_func_full", 10 * hRefmult->GetEntries());

    float IntSc = hRefmult->Integral(minNch, hRefmult->GetNbinsX(), "width");
    hFitResult->Scale(IntSc / hFitResult->Integral(minNch, hFitResult->GetNbinsX(), "width"));
    hFitResult->SetLineColor(2);
    hFitResult->Write();
    hRefmult->Write();
}
double Pb(double *b, double *Nch)
{
    // Parameters
    double n0 = Nch[0];
    double nn = Nch[1];
    // Variables
    double cb = pi * b[0] * b[0] / sigma;
    // Function
    TF1 *Pnb = new TF1("Pnb", PnbGamma2, 0, bin_cent[0], 6);
    Pnb->SetParameters(cb, teta, n_knee, a1, a2, a3);
    double IPnb = Pnb->Integral(n0, nn);
    TF1 *Pn = new TF1("Pn", ftPn, 0, bin_cent[0], 5);
    Pn->SetParameters(teta, n_knee, a1, a2, a3);
    double IPn = Pn->Integral(n0, nn);

    return 2 * pi * b[0] * IPnb / (sigma * IPn);
}
double bmean(double n0, double nn)
{
    TF1 *fPb = new TF1("fPb", Pb, 0, bmax, 2);
    fPb->SetParameters(n0, nn);
    return fPb->Mean(0, bmax);
}

double bsigma(double n0, double nn)
{
    double mean_b = 0, mean_b2 = 0;
    TF1 *fPb = new TF1("fPb", Pb, 0, bmax, 2);
    fPb->SetParameters(n0, nn);
    mean_b = fPb->Mean(0, bmax);
    mean_b2 = fPb->Moment(2, 0, bmax);
    return sqrt(mean_b2 - (mean_b * mean_b));
}

double integ(double n0, double nn)
{
    TF1 *Pn = new TF1("Pn", ftPn, 0, bin_cent[0], 5);
    Pn->SetParameters(teta, n_knee, a1, a2, a3);
    double IPn = Pn->Integral(n0, nn);
    return IPn;
}

void Rebin(double N0)
{
    int n0 = bin_cent[0] - 1;
    double integr = 0, norm;
    norm = integ(N0, bin_cent[0]);
    bin_cent[CentClass - 1] = N0;
    for (int i = 1; i < (CentClass - 1); i++)
    {
        while (integr < range_cent[i] * 0.01)
        {
            n0 = n0 - 1;
            integr = integ(n0, bin_cent[0]) / norm;
            // cout<<"cent="<<range_cent[i]*0.01<<" ,n="<<n0<<" ,integr="<<integr<<endl;
        }
        bin_cent[i] = n0+1;
    }
    cout << "Centrality classes in multiplicity" << endl;
    for (int i = 0; i < (CentClass - 1); i++)
    {
        cout << Form("%d%s%d%s", range_cent[i], "% - ", range_cent[i + 1], "%, ") << bin_cent[i] << " - " << bin_cent[i + 1] << " ,centrality " << integ(bin_cent[i], bin_cent[0]) / norm << " - " << integ(bin_cent[i + 1], bin_cent[0]) / norm << endl;
    }
}

void PlotMeanb()
{
    d_outfile->cd();
    Rebin(1);

    /*const int CentClass = 12;
    int bin_cent[CentClass];
    int range_cent[CentClass] = {0, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100};*/

    TF1 *fPb[CentClass - 1];
    TH1F *Himpact[CentClass - 1];
    double Bmean[CentClass - 1], Bsigma[CentClass - 1], centbin[CentClass - 1], centEr[CentClass - 1];

    for (int cent = 0; cent < CentClass - 1; cent++)
    {
        fPb[cent] = new TF1(Form("fPb%d", cent), Pb, 0, bmax, 2);
        fPb[cent]->SetParameters(bin_cent[cent + 1], bin_cent[cent]);
        Himpact[cent] = new TH1F(Form("%s%d_%d", "B_VS_CentralityClass_", range_cent[cent], range_cent[cent + 1]), "ImpactParametDist", 200, 0, bmax);
        Himpact[cent]->FillRandom(Form("fPb%d", cent), 20000);
        Himpact[cent]->Scale(1 / Himpact[cent]->Integral(1, Himpact[cent]->GetNbinsX(), "width"));
        Himpact[cent]->Write();

        Bmean[cent] = Himpact[cent]->GetMean();
        Bsigma[cent] = Himpact[cent]->GetRMS();
        centbin[cent] = 0.5 * (range_cent[cent] + range_cent[cent + 1]);
        centEr[cent] = 0;
    }

    TGraphErrors *GrFit = new TGraphErrors(CentClass - 1, centbin, Bmean, centEr, Bsigma);
    GrFit->SetMarkerStyle(20);
    GrFit->SetMarkerSize(1.2);
    GrFit->SetMarkerColor(2);
    GrFit->SetLineColor(2);
    GrFit->SetLineWidth(2.);
    GrFit->SetName("Fit_B_Mean");
    GrFit->SetTitle("");
    GrFit->GetXaxis()->SetTitle("Centrality, %");
    GrFit->GetYaxis()->SetTitle("<b>, fm");
    GrFit->SetTitle("B vs Centraliry");
    GrFit->SetName("Fit_B_Mean");
    GrFit->Write();
    TCanvas *CRatio = new TCanvas("CRatio_of_BMean", "CRatio", 500, 500);
    GrFit->Draw("AP");
    TTree *tree = new TTree("Result", "Result");
    Float_t MinPercent;
    Float_t MaxPercent;
    Int_t MinBorder;
    Int_t MaxBorder;
    tree->Branch("MinPercent", &MinPercent, "MinPercent/F");
    tree->Branch("MaxPercent", &MaxPercent, "MaxPercent/F");
    tree->Branch("MinBorder", &MinBorder, "MinBorder/I");
    tree->Branch("MaxBorder", &MaxBorder, "MaxBorder/I");

    for (int i = 0; i < CentClass - 1; i++)
    {
        MinPercent = range_cent[i];
        MaxPercent = range_cent[i + 1];
        MinBorder = bin_cent[i + 1];
        MaxBorder = bin_cent[i];
        tree->Fill();
    }

    tree->Write();
    d_outfile->Close();

    cout << endl
         << "Centrality classes in impact parameter" << endl;
    cout << "Centr. class, bmin-bmax , bmean" << endl;
    cout << Form("%d%s%d", range_cent[0], "% - ", range_cent[1]) << "%, " << 0. << " - " << GrFit->Eval(range_cent[1]) << " fm, " << GrFit->Eval(0.5 * (range_cent[0] + range_cent[1])) << " fm" << endl;
    for (int i = 1; i < CentClass - 1; i++)
    {
        cout << Form("%d%s%d", range_cent[i], "% - ", range_cent[i + 1]) << "%, " << GrFit->Eval(range_cent[i]) << " - " << GrFit->Eval(range_cent[i + 1]) << " fm, " << GrFit->Eval(0.5 * (range_cent[i] + range_cent[i + 1])) << " fm" << endl;
    }
}

void PlotMeanb2(const char *fileadres = "/home/dim/FIT/data/NEWurqmd.root", const char *current_B_vs_mult = "hNpartImpact")
{

    d_outfile->cd();
    Rebin(1);
    TF1 *fPb[CentClass - 1];
    TH1F *Himpact[CentClass - 1];
    double Bmean[CentClass - 1], Bsigma[CentClass - 1], centbin[CentClass - 1], centEr[CentClass - 1];

    for (int cent = 0; cent < CentClass - 1; cent++)
    {
        fPb[cent] = new TF1(Form("fPb%d", cent), Pb, 0, bmax, 2);
        fPb[cent]->SetParameters(bin_cent[cent + 1], bin_cent[cent]);
        Himpact[cent] = new TH1F(Form("%s%d_%d", "B_VS_CentralityClass_", range_cent[cent], range_cent[cent + 1]), "ImpactParametDist", 200, 0, bmax);
        Himpact[cent]->FillRandom(Form("fPb%d", cent), 50000);
        Himpact[cent]->Scale(1 / Himpact[cent]->Integral(1, Himpact[cent]->GetNbinsX(), "width"));
        Himpact[cent]->Write();

        Bmean[cent] = Himpact[cent]->GetMean();
        Bsigma[cent] = Himpact[cent]->GetRMS();
        centbin[cent] = 0.5 * (range_cent[cent] + range_cent[cent + 1]);
        centEr[cent] = 0;
    }

    TGraphErrors *GrFit = new TGraphErrors(CentClass - 1, centbin, Bmean, centEr, Bsigma);
    GrFit->SetMarkerStyle(20);
    GrFit->SetMarkerSize(1.2);
    GrFit->SetMarkerColor(2);
    GrFit->SetLineColor(2);
    GrFit->SetLineWidth(2.);
    GrFit->SetName("Fit_B_Mean");
    GrFit->SetTitle("");
    GrFit->GetXaxis()->SetTitle("Centrality, %");
    GrFit->GetYaxis()->SetTitle("<b>, fm");
    GrFit->SetTitle("B vs Centraliry");
    GrFit->SetName("Fit_B_Mean");
    GrFit->Write();
    TTree *tree = new TTree("Result", "Result");
    Float_t MinPercent;
    Float_t MaxPercent;
    Int_t MinBorder;
    Int_t MaxBorder;
    tree->Branch("MinPercent", &MinPercent, "MinPercent/F");
    tree->Branch("MaxPercent", &MaxPercent, "MaxPercent/F");
    tree->Branch("MinBorder", &MinBorder, "MinBorder/I");
    tree->Branch("MaxBorder", &MaxBorder, "MaxBorder/I");

    for (int i = 0; i < CentClass - 1; i++)
    {
        MinPercent = range_cent[i];
        MaxPercent = range_cent[i + 1];
        MinBorder = bin_cent[i + 1];
        MaxBorder = bin_cent[i];
        tree->Fill();
    }

    tree->Write();
    // d_outfile->Close();

    TF1 *fitB = new TF1("fitB", "pol5", 0., 100.);
   
    GrFit->Fit(fitB, "R0Q");

    cout << endl
         << "Centrality classes in impact parameter" << endl;
    cout << "Centr. class, bmin-bmax , bmean" << endl;
    cout << Form("%d%s%d", range_cent[0], "% - ", range_cent[1]) << "%, " << 0. << " - " << fitB->Eval(range_cent[1]) << " fm, " << fitB->Eval(0.5 * (range_cent[0] + range_cent[1])) << " fm" << endl;
    for (int i = 1; i < CentClass - 1; i++)
    {
        cout << Form("%d%s%d", range_cent[i], "% - ", range_cent[i + 1]) << "%, " << fitB->Eval(range_cent[i]) << " - " << fitB->Eval(range_cent[i + 1]) << " fm, " << fitB->Eval(0.5 * (range_cent[i] + range_cent[i + 1])) << " fm" << endl;
    }

    TFile *file1 = new TFile(fileadres);
    TH2F *MeanB;
    MeanB = (TH2F *)file1->Get(current_B_vs_mult);
    TH1D *Hist[CentClass - 1];
    d_outfile->cd();
    for (int i = 0; i < (CentClass - 1); i++)
    {
        Hist[i] = MeanB->ProjectionY(Form("%s%d_%d_Data", "B_VS_CentralityClass_", range_cent[i], range_cent[i + 1]), bin_cent[i + 1] + 1, bin_cent[i] + 1);
        Hist[i]->Scale(1 / Hist[i]->Integral(1, Hist[i]->GetNbinsX(), "width"));
        Hist[i]->Write();
    }

    double y[CentClass - 1], ey[CentClass - 1];

    for (int i = 0; i < (CentClass - 1); i++)
    {
        y[i] = Hist[i]->GetMean();
        ey[i] = Hist[i]->GetStdDev();
    }

    TGraphErrors *GrData = new TGraphErrors(CentClass - 1, centbin, y, centEr, ey);
    GrData->SetMarkerStyle(21);
    GrData->SetMarkerSize(1.5);
    GrData->SetMarkerColor(4);
    GrData->SetLineColor(4);
    GrData->GetYaxis()->SetTitle("<b>, fm");
    GrData->SetTitle("");
    GrData->SetName("Data_B_Mean");
    d_outfile->cd();
    GrData->Write();

    TCanvas *CRatio = new TCanvas("CRatio_of_BMean", "CRatio", 1000, 500);
    CRatio->Divide(2, 1);
    CRatio->cd(1);
    GrData->Draw("AP");
    GrFit->Draw("SAME P");
    fitB->SetLineColor(1);
    fitB->SetLineStyle(1);
    fitB->Draw("SAME");
    CRatio->cd(2);
    TGraphErrors *GrRatio = RatioGr(GrFit, GrData, 0, 0.89, 100, 1.11);
    GrRatio->SetTitle("");
    GrRatio->GetYaxis()->SetTitle("Fit/Data");
    GrRatio->GetXaxis()->SetTitle("Cent., %");
    GrRatio->Draw();
    TLine line;
    line.SetLineWidth(2);
    line.SetLineStyle(2);
    line.DrawLine(0, 1, 100, 1);
    CRatio->Write();
    d_outfile->Close();
}

// const char *fileadres = "/home/dim/Work/FIT/data/refMult_DCMQGSMSMM_11gev_500k.root", const char *current_mult = "hRefMultSTAR",
// const char *fileadres = "/home/dim/Work/FIT/data/refMult_ampt15_11gev_500k.root", const char *current_mult = "hRefMultSTAR",
// const char *fileadres = "file:///home/dim/Work/FIT/data/refMult_UrQMD_11gev_500k.root", const char *current_mult = "hRefMultSTAR",
// const char *fileadres = "file:///home/dim/Work/FIT/dataNEW/15ampt200_1M.root", const char *current_mult = "hRefMultSTAR",

void GammaFit2(const char *fileadres = "file:///home/dim/Work/FIT/dataNEW/FITGlaub/gen_AuAu9/gen_GSM_AuAu9.root", const char *current_mult = "hRefMultSTAR",
               const char *outadres = "/home/dim/Work/FIT/dataNEW/FITGlaub/gen_AuAu9/gamma.root", int minNch = 10, const char *current_B_vs_mult = "hBvsRefMult",bool efficiencyFit = false,
               const char *fileadres2 = "", const char *current_mult2 = "")
{
    // Enable implicit parallelization
    ROOT::EnableImplicitMT();
    Start(fileadres, current_mult, outadres, minNch,current_B_vs_mult, efficiencyFit, fileadres2, current_mult2);
    // PlotMeanb();

    PlotMeanb2(fileadres, current_B_vs_mult);
    // PlotMeanb();
}

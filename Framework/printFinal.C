#include <iostream>
#include <vector>

#include <Rtypes.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TString.h>

void printFinal(TString inFileName="", TString outFileName="")
{
    if (inFileName == "") return;

    Bool_t isTypeSimple = false;
    Bool_t isTypeLatex = false;
    Bool_t isTypeSvc = false;

    // Define output file format
    if (outFileName == "") isTypeSimple = true;
    else if (outFileName.Contains(".tex")) isTypeLatex = true;
    else if (outFileName.Contains(".cvs")) isTypeSvc = true;

    TFile *fi = new TFile(inFileName.Data(),"read");
    if (!fi) return;

    TH1D *hBavg = (TH1D*) fi->Get("B_average_VS_Centrality");
    TH1D *hNpartavg = (TH1D*) fi->Get("Npart_average_VS_Centrality");
    TH1D *hNcollavg = (TH1D*) fi->Get("Ncoll_average_VS_Centrality");

    TTree *Result = (TTree*) fi->Get("Result");
    if (!Result) return;

    Int_t Ncc;
    Float_t MinPercent;
    Float_t MaxPercent;
    Int_t MinBorder;
    Int_t MaxBorder;

    Result->SetBranchAddress("Ncc", &Ncc);
    Result->SetBranchAddress("MinPercent", &MinPercent);
    Result->SetBranchAddress("MaxPercent", &MaxPercent);
    Result->SetBranchAddress("MinBorder", &MinBorder);
    Result->SetBranchAddress("MaxBorder", &MaxBorder);

    std::vector<std::pair<Float_t, Float_t>> vCent;
    std::vector<std::pair<Int_t, Int_t>> vBorders;
    std::vector<std::pair<Float_t, Float_t>> vBimp;
    std::vector<std::pair<Float_t, Float_t>> vNpart;
    std::vector<std::pair<Float_t, Float_t>> vNcoll;
    std::vector<Float_t> vBavg;
    std::vector<Float_t> vNpartavg;
    std::vector<Float_t> vNcollavg;

    Int_t Nclasses = Result->GetEntries();

    // Read TTree
    for (int i=0; i<Nclasses; i++)
    {
        Result->GetEntry(i);
        vCent.push_back({MinPercent, MaxPercent});
        vBorders.push_back({MinBorder, MaxBorder});
    }

    // Read averaged values
    for (int i=0; i<hBavg->GetNbinsX(); i++)
    {
        if (hBavg->GetBinContent(i+1) != 0)
            vBavg.push_back(hBavg->GetBinContent(i+1));
    }
    for (int i=0; i<hNpartavg->GetNbinsX(); i++)
    {
        if (hNpartavg->GetBinContent(i+1) != 0)
            vNpartavg.push_back(hNpartavg->GetBinContent(i+1));
    }
    for (int i=0; i<hNcollavg->GetNbinsX(); i++)
    {
        if (hNcollavg->GetBinContent(i+1) != 0)
            vNcollavg.push_back(hNcollavg->GetBinContent(i+1));
    }

    // Fitting mean values with pol5 function
    TF1 *fitB = new TF1("fitB","pol5",0.,100.);
    TF1 *fitNpart = new TF1("fitNpart","pol5",0.,100.);
    TF1 *fitNcoll = new TF1("fitNcoll","pol5",0.,100.);

    hBavg->Fit(fitB,"R0Q");
    hNpartavg->Fit(fitNpart,"R0Q");
    hNcollavg->Fit(fitNcoll,"R0Q");

    // Extract min/max values for B, Npart, Ncoll
    Float_t fmin, fmax;
    for (int i=0; i<vBavg.size(); i++)
    {
        fmin = fitB->Eval(vCent.at(i).first);
        fmax = fitB->Eval(vCent.at(i).second);
        vBimp.push_back({fmin, fmax});
    }
    for (int i=0; i<vNpartavg.size(); i++)
    {
        fmin = fitNpart->Eval(vCent.at(i).first);
        fmax = fitNpart->Eval(vCent.at(i).second);
        vNpart.push_back({fmin, fmax});
    }
    for (int i=0; i<vNcollavg.size(); i++)
    {
        fmin = fitNcoll->Eval(vCent.at(i).first);
        fmax = fitNcoll->Eval(vCent.at(i).second);
        vNcoll.push_back({fmin, fmax});
    }

    Int_t NreasonableClasses = vBavg.size();

    // Output for isTypeSimple case
    if (isTypeSimple)
    {
        std::cout << "File: " << inFileName.Data() << "." << std::endl;
        std::cout << "Cent, %   | Mult_min | Mult_max | <b>, fm | bmin, fm | bmax, fm | <Npart> | Npart_min | Npart_max | <Ncoll> | Ncoll_min | Ncoll_max |" << std::endl;
        std::cout << "----------|----------|----------|---------|----------|----------|---------|-----------|-----------|---------|-----------|-----------|" << std::endl;
        for (int i=0; i<NreasonableClasses; i++)
        {
            std::cout << Form("%3.0f - %3.0f | %8i | %8i | %7.2f | %8.2f | %8.2f | %7.2f | %9.2f | %9.2f | %7.2f | %9.2f | %9.2f |", 
                vCent.at(i).first, vCent.at(i).second,
                vBorders.at(i).first, vBorders.at(i).second,
                vBavg.at(i), vBimp.at(i).first, vBimp.at(i).second,
                vNpartavg.at(i), vNpart.at(i).second, vNpart.at(i).first,
                vNcollavg.at(i), vNcoll.at(i).second, vNcoll.at(i).first) << std::endl;
        }
        std::cout << "-------------------------------------------------------------------------------------------------------------------------------------" << std::endl;
    }
    
    // Output for isTypeLatex case
    if (isTypeLatex)
    {
        std::cout << "File: " << inFileName.Data() << "." << std::endl;

    }
}
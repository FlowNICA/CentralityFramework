#include <iostream>
#include <fstream>
#include <vector>
#include <string>

#include <Rtypes.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>

void printFinal(TString inFileName="", TString outFileName="")
{
    if (inFileName == "") return;

    Bool_t isTypeSimple = false;
    Bool_t isTypeLatex = false;
    Bool_t isTypeScv = false;
    Bool_t isTypeCpp = false;

    // Define output file format
    if (outFileName == "") isTypeSimple = true;
    else if (outFileName.Contains(".tex")) isTypeLatex = true;
    else if (outFileName.Contains(".csv")) isTypeScv = true;
    else if (outFileName.Contains(".C")) isTypeCpp = true;
    else if (outFileName.Contains(".cpp")) isTypeCpp = true;
    else if (outFileName.Contains(".h")) isTypeCpp = true;
    else if (outFileName.Contains(".hpp")) isTypeCpp = true;

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
    Double_t Bavg;
    Double_t Bwidth;
    Double_t NpartAvg;
    Double_t NpartWidth;
    Double_t NcollAvg;
    Double_t NcollWidth;

    Result->SetBranchAddress("Ncc", &Ncc);
    Result->SetBranchAddress("MinPercent", &MinPercent);
    Result->SetBranchAddress("MaxPercent", &MaxPercent);
    Result->SetBranchAddress("MinBorder", &MinBorder);
    Result->SetBranchAddress("MaxBorder", &MaxBorder);
    if (!hBavg){
      Result->SetBranchAddress("BAverage", &Bavg);
      Result->SetBranchAddress("BWidth", &Bwidth);
    }
    if (!hNpartavg){
      Result->SetBranchAddress("NpartAverage", &NpartAvg);
      Result->SetBranchAddress("NpartWidth", &NpartWidth);
    }
    if (!hNcollavg){
      Result->SetBranchAddress("NcollAverage", &NcollAvg);
      Result->SetBranchAddress("NcollWidth", &NcollWidth);
    }

    std::vector<std::pair<Float_t, Float_t>> vCent;
    std::vector<std::pair<Int_t, Int_t>> vBorders;
    std::vector<std::pair<Float_t, Float_t>> vBimp;
    std::vector<std::pair<Float_t, Float_t>> vNpart;
    std::vector<std::pair<Float_t, Float_t>> vNcoll;
    std::vector<Float_t> vBavg;
    std::vector<Float_t> vBavgRMS;
    std::vector<Float_t> vNpartavg;
    std::vector<Float_t> vNpartavgRMS;
    std::vector<Float_t> vNcollavg;
    std::vector<Float_t> vNcollavgRMS;

    Int_t Nclasses = Result->GetEntriesFast();

    // Read TTree
    for (int i=0; i<Nclasses; i++)
    {
        Result->GetEntry(i);
        vCent.push_back({MinPercent, MaxPercent});
        vBorders.push_back({MinBorder, MaxBorder});
        if (!hBavg){
          vBavg.push_back(Bavg);
          vBavgRMS.push_back(Bwidth);
        }
        if (!hNpartavg){
          vNpartavg.push_back(NpartAvg);
          vNpartavgRMS.push_back(NpartWidth);
        }
        if (!hNcollavg){
          vNcollavg.push_back(NcollAvg);
          vNcollavgRMS.push_back(NcollWidth);
        }
    }
    
    // Read averaged values
    if (hBavg){
        for (int i=0; i<hBavg->GetNbinsX(); i++)
        {
            if (hBavg->GetBinContent(i+1) != 0){
                vBavg.push_back(hBavg->GetBinContent(i+1));
                vBavgRMS.push_back(hBavg->GetBinError(i+1));
            }
        }
    }
    if (hNpartavg){
        for (int i=0; i<hNpartavg->GetNbinsX(); i++)
        {
            if (hNpartavg->GetBinContent(i+1) != 0)
            {
                vNpartavg.push_back(hNpartavg->GetBinContent(i+1));
                vNpartavgRMS.push_back(hNpartavg->GetBinError(i+1));
            }
        }
    }
    if (hNcollavg){
        for (int i=0; i<hNcollavg->GetNbinsX(); i++)
        {
            if (hNcollavg->GetBinContent(i+1) != 0)
            {
                vNcollavg.push_back(hNcollavg->GetBinContent(i+1));
                vNcollavgRMS.push_back(hNcollavg->GetBinError(i+1));
            }
        }
    }

    if (!hBavg){
        hBavg = new TH1D("B_average_VS_Centrality", "B_average_VS_Centrality", Nclasses, vCent.at(0).first, vCent.at(vCent.size()-1).second);
        for (int i=0; i<Nclasses; i++){
            hBavg->SetBinContent(i+1, vBavg.at(i));
            hBavg->SetBinError(i+1, vBavgRMS.at(i));
        }
    }

    if (!hNpartavg){
        hNpartavg = new TH1D("Npart_average_VS_Centrality", "Npart_average_VS_Centrality", Nclasses, vCent.at(0).first, vCent.at(vCent.size()-1).second);
        for (int i=0; i<Nclasses; i++){
            hNpartavg->SetBinContent(i+1, vNpartavg.at(i));
            hNpartavg->SetBinError(i+1, vNpartavgRMS.at(i));
        }
    }

    if (!hNcollavg){
        hNcollavg = new TH1D("Ncoll_average_VS_Centrality", "Ncoll_average_VS_Centrality", Nclasses, vCent.at(0).first, vCent.at(vCent.size()-1).second);
        for (int i=0; i<Nclasses; i++){
            hNcollavg->SetBinContent(i+1, vNcollavg.at(i));
            hNcollavg->SetBinError(i+1, vNcollavgRMS.at(i));
        }
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

    if (!isTypeSimple && !isTypeLatex && !isTypeScv && !isTypeCpp)
    {
        std::cerr << "Output format is not known!" << std::endl;
        return;
    }

    // Output for isTypeSimple case
    if (isTypeSimple)
    {
        std::cout << "File: " << inFileName.Data() << "." << std::endl;
        std::cout << "Cent, %   | Mult_min | Mult_max | <b>, fm |   RMS   | bmin, fm | bmax, fm | <Npart> |   RMS   | Npart_min | Npart_max | <Ncoll> |   RMS   | Ncoll_min | Ncoll_max |" << std::endl;
        std::cout << "----------|----------|----------|---------|---------|----------|----------|---------|---------|-----------|-----------|---------|---------|-----------|-----------|" << std::endl;
        for (int i=0; i<NreasonableClasses; i++)
        {
            std::cout << Form("%3.0f - %3.0f | %8i | %8i | %7.2f | %7.2f | %8.2f | %8.2f | %7.2f | %7.2f | %9.2f | %9.2f | %7.2f | %7.2f | %9.2f | %9.2f |", 
                vCent.at(i).first, vCent.at(i).second,
                vBorders.at(i).first, vBorders.at(i).second,
                vBavg.at(i), vBavgRMS.at(i), vBimp.at(i).first, vBimp.at(i).second,
                vNpartavg.at(i), vNpartavgRMS.at(i), vNpart.at(i).second, vNpart.at(i).first,
                vNcollavg.at(i), vNcollavgRMS.at(i), vNcoll.at(i).second, vNcoll.at(i).first) << std::endl;
        }
        std::cout << "-------------------------------------------------------------------------------------------------------------------------------------------------------------------" << std::endl;
    }
    
    // Output for isTypeLatex case
    ofstream myfile;
    if (isTypeLatex)
    {
        std::cout << "File: " << inFileName.Data() << "." << std::endl;
        myfile.open(outFileName.Data());
        myfile << "\\documentclass[11pt]{article}\n";
        myfile << "\\usepackage[utf8]{inputenc}\n";
        myfile << "\\usepackage{geometry}";
        // myfile << "\\geometry{a4paper, total={170mm,257mm}, left=20mm, top=20mm}\n" << std::endl;
        myfile << "\\geometry{legalpaper, landscape, margin=2in}\n" << std::endl;
        myfile << "\\begin{document}\n" << std::endl;
        myfile << "Generated from a file:\n";
        myfile << "\\begin{verbatim*}\n";
        myfile << inFileName.Data() << std::endl;
        myfile << "\\end{verbatim*}\n";
        myfile << "\\begin{center}\n";
        myfile << "\\begin{tabular}{ |c|c|c|c|c|c|c|c|c|c|c|c|c|c|c| }\n";
        myfile << "\t\\hline\n";
        // myfile << "\t Centrality, \\% & $N_{ch}^{min}$ & $N_{ch}^{max}$ & $\\langle b \\rangle$, fm & RMS & $b_{min}$, fm & $b_{max}$, fm & $\\langle N_{part} \\rangle$ & RMS & $N_{part}^{min}$ & $N_{part}^{max}$ & $\\langle N_{coll} \\rangle$ & RMS & $N_{coll}^{min}$ & $N_{coll}^{max}$ \\\\\n";
        myfile << "\t Centrality, \\% & $N_{ch}^{min}$ & $N_{ch}^{max}$ & $\\langle b \\rangle$, fm & RMS & $\\langle N_{part} \\rangle$ & RMS & $\\langle N_{coll} \\rangle$ & RMS \\\\\n";
        for (int i=0; i<NreasonableClasses; i++)
        {
            myfile << "\t\\hline\n";
            // myfile << Form("\t%.0f - %.0f & %i & %i & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f \\\\\n",
            //     vCent.at(i).first, vCent.at(i).second,
            //     vBorders.at(i).first, vBorders.at(i).second,
            //     vBavg.at(i), vBavgRMS.at(i), vBimp.at(i).first, vBimp.at(i).second,
            //     vNpartavg.at(i), vNpartavgRMS.at(i), vNpart.at(i).second, vNpart.at(i).first,
            //     vNcollavg.at(i), vNcollavgRMS.at(i), vNcoll.at(i).second, vNcoll.at(i).first);
            myfile << Form("\t%.0f - %.0f & %i & %i & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f \\\\\n",
                vCent.at(i).first, vCent.at(i).second,
                vBorders.at(i).first, vBorders.at(i).second,
                vBavg.at(i), vBavgRMS.at(i), //vBimp.at(i).first, vBimp.at(i).second,
                vNpartavg.at(i), vNpartavgRMS.at(i), //vNpart.at(i).second, vNpart.at(i).first,
                vNcollavg.at(i), vNcollavgRMS.at(i)); //vNcoll.at(i).second, vNcoll.at(i).first);
        }
        myfile << "\t\\hline\n";
        myfile << "\\end{tabular}\n";
        myfile << "\\end{center}\n";

        myfile << "\\end{document}";
        std::cout << "Output file " << outFileName.Data() << " is created." << std::endl;
        myfile.close();
    }

    // Outout for isTypeScv case
    if (isTypeScv)
    {
        std::cout << "File: " << inFileName.Data() << "." << std::endl;
        myfile.open(outFileName.Data());

        myfile << "Generated from a file: " << inFileName.Data() << "\n";
        myfile << "Centrality class,Mult_min,Mult_max,Mean b,RMS,b_min,b_max,Mean Npart,RMS,Npart_min,Npart_max,Mean Ncoll,RMS,Ncoll_min,Ncoll_max,\n";
        for (int i=0; i<NreasonableClasses; i++)
        {
            myfile << Form("%.0f - %.0f,%i,%i,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,\n",
                vCent.at(i).first, vCent.at(i).second,
                vBorders.at(i).first, vBorders.at(i).second,
                vBavg.at(i), vBavgRMS.at(i), vBimp.at(i).first, vBimp.at(i).second,
                vNpartavg.at(i), vNpartavgRMS.at(i), vNpart.at(i).second, vNpart.at(i).first,
                vNcollavg.at(i), vNcollavgRMS.at(i), vNcoll.at(i).second, vNcoll.at(i).first);
        }
        std::cout << "Output file " << outFileName.Data() << " is created." << std::endl;
        myfile.close();
    }

    // Output for isTypeCpp case
    if (isTypeCpp)
    {
        std::cout << "File: " << inFileName.Data() << "." << std::endl;
        myfile.open(outFileName.Data());

        myfile << "// Generated from a file: " << inFileName.Data() << "\n";
        myfile << "Float_t minCentPercent [" << NreasonableClasses << "] = { ";
        for (int i=0; i<NreasonableClasses; i++)
        {
            if (i!=NreasonableClasses-1) myfile << vCent.at(i).first << ", ";
            if (i==NreasonableClasses-1) myfile << vCent.at(i).first << "};\n";
        }
        myfile << "Float_t maxCentPercent [" << NreasonableClasses << "] = { ";
        for (int i=0; i<NreasonableClasses; i++)
        {
            if (i!=NreasonableClasses-1) myfile << vCent.at(i).second << ", ";
            if (i==NreasonableClasses-1) myfile << vCent.at(i).second << "};\n";
        }
        myfile << "Int_t minMult [" << NreasonableClasses << "] = { ";
        for (int i=0; i<NreasonableClasses; i++)
        {
            if (i!=NreasonableClasses-1) myfile << vBorders.at(i).first << ", ";
            if (i==NreasonableClasses-1) myfile << vBorders.at(i).first << "};\n";
        }
        myfile << "Int_t maxMult [" << NreasonableClasses << "] = { ";
        for (int i=0; i<NreasonableClasses; i++)
        {
            if (i!=NreasonableClasses-1) myfile << vBorders.at(i).second << ", ";
            if (i==NreasonableClasses-1) myfile << vBorders.at(i).second << "};\n";
        }
        myfile << "Float_t meanB [" << NreasonableClasses << "] = { ";
        for (int i=0; i<NreasonableClasses; i++)
        {
            if (i!=NreasonableClasses-1) myfile << vBavg.at(i) << ", ";
            if (i==NreasonableClasses-1) myfile << vBavg.at(i) << "};\n";
        }
        myfile << "Float_t rmsB [" << NreasonableClasses << "] = { ";
        for (int i=0; i<NreasonableClasses; i++)
        {
            if (i!=NreasonableClasses-1) myfile << vBavgRMS.at(i) << ", ";
            if (i==NreasonableClasses-1) myfile << vBavgRMS.at(i) << "};\n";
        }
        myfile << "Float_t minB [" << NreasonableClasses << "] = { ";
        for (int i=0; i<NreasonableClasses; i++)
        {
            if (i!=NreasonableClasses-1) myfile << vBimp.at(i).first << ", ";
            if (i==NreasonableClasses-1) myfile << vBimp.at(i).first << "};\n";
        }
        myfile << "Float_t maxB [" << NreasonableClasses << "] = { ";
        for (int i=0; i<NreasonableClasses; i++)
        {
            if (i!=NreasonableClasses-1) myfile << vBimp.at(i).second << ", ";
            if (i==NreasonableClasses-1) myfile << vBimp.at(i).second << "};\n";
        }
        myfile << "Float_t meanNpart [" << NreasonableClasses << "] = { ";
        for (int i=0; i<NreasonableClasses; i++)
        {
            if (i!=NreasonableClasses-1) myfile << vNpartavg.at(i) << ", ";
            if (i==NreasonableClasses-1) myfile << vNpartavg.at(i) << "};\n";
        }
        myfile << "Float_t rmsNpart [" << NreasonableClasses << "] = { ";
        for (int i=0; i<NreasonableClasses; i++)
        {
            if (i!=NreasonableClasses-1) myfile << vNpartavgRMS.at(i) << ", ";
            if (i==NreasonableClasses-1) myfile << vNpartavgRMS.at(i) << "};\n";
        }
        myfile << "Float_t minNpart [" << NreasonableClasses << "] = { ";
        for (int i=0; i<NreasonableClasses; i++)
        {
            if (i!=NreasonableClasses-1) myfile << vNpart.at(i).first << ", ";
            if (i==NreasonableClasses-1) myfile << vNpart.at(i).first << "};\n";
        }
        myfile << "Float_t maxNpart [" << NreasonableClasses << "] = { ";
        for (int i=0; i<NreasonableClasses; i++)
        {
            if (i!=NreasonableClasses-1) myfile << vNpart.at(i).second << ", ";
            if (i==NreasonableClasses-1) myfile << vNpart.at(i).second << "};\n";
        }
        myfile << "Float_t meanNcoll [" << NreasonableClasses << "] = { ";
        for (int i=0; i<NreasonableClasses; i++)
        {
            if (i!=NreasonableClasses-1) myfile << vNcollavg.at(i) << ", ";
            if (i==NreasonableClasses-1) myfile << vNcollavg.at(i) << "};\n";
        }
        myfile << "Float_t rmsNcoll [" << NreasonableClasses << "] = { ";
        for (int i=0; i<NreasonableClasses; i++)
        {
            if (i!=NreasonableClasses-1) myfile << vNcollavgRMS.at(i) << ", ";
            if (i==NreasonableClasses-1) myfile << vNcollavgRMS.at(i) << "};\n";
        }
        myfile << "Float_t minNcoll [" << NreasonableClasses << "] = { ";
        for (int i=0; i<NreasonableClasses; i++)
        {
            if (i!=NreasonableClasses-1) myfile << vNcoll.at(i).first << ", ";
            if (i==NreasonableClasses-1) myfile << vNcoll.at(i).first << "};\n";
        }
        myfile << "Float_t maxNcoll [" << NreasonableClasses << "] = { ";
        for (int i=0; i<NreasonableClasses; i++)
        {
            if (i!=NreasonableClasses-1) myfile << vNcoll.at(i).second << ", ";
            if (i==NreasonableClasses-1) myfile << vNcoll.at(i).second << "};\n";
        }

        myfile << std::endl;
        myfile << "Float_t GetCentMult(Int_t mult)\n";
        myfile << "//Function returns centrality (in percent) for a given multiplicity\n";
        myfile << "{\n";
        myfile << "\tInt_t centBin = -1;\n";
        myfile << "\tfor (Int_t i = 0; i < " << NreasonableClasses << "; i++)\n";
        myfile << "\t{\n";
        myfile << "\t\tif (mult >= minMult[i] && mult < maxMult[i])\n";
        myfile << "\t\t\tcentBin = i;\n";
        myfile << "\t}\n";
        myfile << "\tif (centBin == -1) return -1.;\n";
        myfile << "\treturn (Double_t)((maxCentPercent[centBin] - minCentPercent[centBin]) / 2. + minCentPercent[centBin]);\n";
        myfile << "}";
        
        std::cout << "Output file " << outFileName.Data() << " is created." << std::endl;
        myfile.close();
    }
    
}

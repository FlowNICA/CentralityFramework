#include "Fitter.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TMath.h"
#include "TRandom.h"

ClassImp(Glauber::Fitter)

    // -----   Default constructor   -------------------------------------------
    Glauber::Fitter::Fitter(std::unique_ptr<TTree> tree)
{
    fSimTree = std::move(tree);
    std::cout << fSimTree->GetEntries() << std::endl;

    if (!fSimTree)
    {
        std::cout << "SetSimHistos: *** Error - " << std::endl;
        exit(EXIT_FAILURE);
    }

    fSimTree->SetBranchAddress("B", &fB);
    fSimTree->SetBranchAddress("Npart", &fNpart);
    fSimTree->SetBranchAddress("Ncoll", &fNcoll);
    fSimTree->SetBranchAddress("Ecc1", &fEcc1);
    fSimTree->SetBranchAddress("Psi1", &fPsi1);
    fSimTree->SetBranchAddress("Ecc2", &fEcc2);
    fSimTree->SetBranchAddress("Psi2", &fPsi2);
    fSimTree->SetBranchAddress("Ecc3", &fEcc3);
    fSimTree->SetBranchAddress("Psi3", &fPsi3);
    fSimTree->SetBranchAddress("Ecc4", &fEcc4);
    fSimTree->SetBranchAddress("Psi4", &fPsi4);
    fSimTree->SetBranchAddress("Ecc5", &fEcc5);
    fSimTree->SetBranchAddress("Psi5", &fPsi5);
}

#ifdef __THREADS_ON__
Glauber::Fitter::Fitter(std::unique_ptr<TTree> tree, unsigned int Nthreads)
{
    fNthreads = Nthreads;

    fSimTree = std::move(tree);
    std::cout << fSimTree->GetEntries() << std::endl;

    if (!fSimTree)
    {
        std::cout << "SetSimHistos: *** Error - " << std::endl;
        exit(EXIT_FAILURE);
    }

    fSimTree->SetBranchAddress("B", &fB);
    fSimTree->SetBranchAddress("Npart", &fNpart);
    fSimTree->SetBranchAddress("Ncoll", &fNcoll);
    fSimTree->SetBranchAddress("Ecc1", &fEcc1);
    fSimTree->SetBranchAddress("Psi1", &fPsi1);
    fSimTree->SetBranchAddress("Ecc2", &fEcc2);
    fSimTree->SetBranchAddress("Psi2", &fPsi2);
    fSimTree->SetBranchAddress("Ecc3", &fEcc3);
    fSimTree->SetBranchAddress("Psi3", &fPsi3);
    fSimTree->SetBranchAddress("Ecc4", &fEcc4);
    fSimTree->SetBranchAddress("Psi4", &fPsi4);
    fSimTree->SetBranchAddress("Ecc5", &fEcc5);
    fSimTree->SetBranchAddress("Psi5", &fPsi5);
}
#endif

void Glauber::Fitter::Init(int nEntries, TString fmode)
{

    if (nEntries < 0 || nEntries > fSimTree->GetEntries())
    {
        std::cout << "Init: *** ERROR - number of entries < 0 or less that number of entries in input tree" << std::endl;
        std::cout << "Init: *** number of entries in input tree = " << fSimTree->GetEntries() << std::endl;
        exit(EXIT_FAILURE);
    }

    fvB.clear();
    fvNpart.clear();
    fvNcoll.clear();
    fvEcc1.clear();
    fvPsi1.clear();
    fvEcc2.clear();
    fvPsi2.clear();
    fvEcc3.clear();
    fvPsi3.clear();
    fvEcc4.clear();
    fvPsi4.clear();
    fvEcc5.clear();
    fvPsi5.clear();

    const int BMax = int(fSimTree->GetMaximum("B"));
    const int NpartMax = int(fSimTree->GetMaximum("Npart"));
    const int NcollMax = int(fSimTree->GetMaximum("Ncoll"));

    const float Ecc1Min = fSimTree->GetMinimum("Ecc1");
    const float Psi1Min = fSimTree->GetMinimum("Psi1");
    const float Ecc2Min = fSimTree->GetMinimum("Ecc2");
    const float Psi2Min = fSimTree->GetMinimum("Psi2");
    const float Ecc3Min = fSimTree->GetMinimum("Ecc3");
    const float Psi3Min = fSimTree->GetMinimum("Psi3");
    const float Ecc4Min = fSimTree->GetMinimum("Ecc4");
    const float Psi4Min = fSimTree->GetMinimum("Psi4");
    const float Ecc5Min = fSimTree->GetMinimum("Ecc5");
    const float Psi5Min = fSimTree->GetMinimum("Psi5");
    const float Ecc1Max = fSimTree->GetMaximum("Ecc1");
    const float Psi1Max = fSimTree->GetMaximum("Psi1");
    const float Ecc2Max = fSimTree->GetMaximum("Ecc2");
    const float Psi2Max = fSimTree->GetMaximum("Psi2");
    const float Ecc3Max = fSimTree->GetMaximum("Ecc3");
    const float Psi3Max = fSimTree->GetMaximum("Psi3");
    const float Ecc4Max = fSimTree->GetMaximum("Ecc4");
    const float Psi4Max = fSimTree->GetMaximum("Psi4");
    const float Ecc5Max = fSimTree->GetMaximum("Ecc5");
    const float Psi5Max = fSimTree->GetMaximum("Psi5");

    fBHisto = TH1F("fBHisto", "B", BMax / fBinSize, 0, BMax);
    fNpartHisto = TH1F("fNpartHisto", "Npart", NpartMax / fBinSize, 0, NpartMax);
    fNcollHisto = TH1F("fNcollHisto", "Ncoll", NcollMax / fBinSize, 0, NcollMax);
    fEcc1Histo = TH1F("fEcc1Histo", "#epsilon1", (Ecc1Max - Ecc1Min) / 0.01, Ecc1Min, Ecc1Max);
    fPsi1Histo = TH1F("fPsi1Histo", "#psi1", (Psi1Max - Psi1Min) / 0.01, Psi1Min, Psi1Max);
    fEcc2Histo = TH1F("fEcc2Histo", "#epsilon2", (Ecc2Max - Ecc2Min) / 0.01, Ecc2Min, Ecc2Max);
    fPsi2Histo = TH1F("fPsi2Histo", "#psi2", (Psi2Max - Psi2Min) / 0.01, Psi2Min, Psi2Max);
    fEcc3Histo = TH1F("fEcc3Histo", "#epsilon3", (Ecc3Max - Ecc3Min) / 0.01, Ecc3Min, Ecc3Max);
    fPsi3Histo = TH1F("fPsi3Histo", "#psi3", (Psi3Max - Psi3Min) / 0.01, Psi3Min, Psi3Max);
    fEcc4Histo = TH1F("fEcc4Histo", "#epsilon4", (Ecc4Max - Ecc4Min) / 0.01, Ecc4Min, Ecc4Max);
    fPsi4Histo = TH1F("fPsi4Histo", "#psi4", (Psi4Max - Psi4Min) / 0.01, Psi4Min, Psi4Max);
    fEcc5Histo = TH1F("fEcc5Histo", "#epsilon5", (Ecc5Max - Ecc5Min) / 0.01, Ecc5Min, Ecc5Max);
    fPsi5Histo = TH1F("fPsi5Histo", "#psi5", (Psi5Max - Psi5Min) / 0.01, Psi5Min, Psi5Max);

    for (int i = 0; i < nEntries; i++)
    {
        fSimTree->GetEntry(i);
        fBHisto.Fill(fB);
        fNcollHisto.Fill(fNcoll);
        fNpartHisto.Fill(fNpart);
        fEcc1Histo.Fill(fEcc1);
        fPsi1Histo.Fill(fPsi1);
        fEcc2Histo.Fill(fEcc2);
        fPsi2Histo.Fill(fPsi2);
        fEcc3Histo.Fill(fEcc3);
        fPsi3Histo.Fill(fPsi3);
        fEcc4Histo.Fill(fEcc4);
        fPsi4Histo.Fill(fPsi4);
        fEcc5Histo.Fill(fEcc5);
        fPsi5Histo.Fill(fPsi5);

        fvB.push_back(fB);
        fvNpart.push_back(fNpart);
        fvNcoll.push_back(fNcoll);
        fvEcc1.push_back(fEcc1);
        fvPsi1.push_back(fPsi1);
        fvEcc2.push_back(fEcc2);
        fvPsi2.push_back(fPsi2);
        fvEcc3.push_back(fEcc3);
        fvPsi3.push_back(fPsi3);
        fvEcc4.push_back(fEcc4);
        fvPsi4.push_back(fPsi4);
        fvEcc5.push_back(fEcc5);
        fvPsi5.push_back(fPsi5);
    }
    std::cout << fSimTree->GetEntries() << std::endl;

    fNbins = fDataHisto.GetNbinsX();

    while (fDataHisto.GetBinContent(fNbins - 1) == 0)
        fNbins--;

    fNbins++;

    const float min = fDataHisto.GetXaxis()->GetXmin();
    const float max = fDataHisto.GetXaxis()->GetXmax();

    fMaxValue = min + (max - min) * fNbins / fDataHisto.GetNbinsX();

    std::cout << "fNbins = " << fNbins << std::endl;
    std::cout << "fMaxValue = " << fMaxValue << std::endl;

#ifdef __THREADS_ON__
    std::cout << std::endl;
    std::cout << "Number of threads found: " << fNthreads << std::endl;
    std::cout << std::endl;
#endif
}

float Glauber::Fitter::Nancestors(float f) const
{
    if (fMode == "Default")
        return f * fNpart + (1 - f) * fNcoll;
    else if (fMode == "PSD")
        return f - fNpart;
    else if (fMode == "Npart")
        return pow(fNpart, f);
    else if (fMode == "Ncoll")
        return pow(fNcoll, f);
    else if (fMode == "NpartFast")
        return pow(fNpart, f) / pow(10, f);
    else if (fMode == "NcollFast")
        return pow(fNcoll, f) / pow(100, f);
    else if (fMode == "STAR")
        return (1 - f) * fNpart / 2. + f * fNcoll;
    else if (fMode == "HADES")
        return (1. - f * fNpart * fNpart) * fNpart;

    return -1.;
}

float Glauber::Fitter::Nancestors(float f, float npart, float ncoll) const
{
    if (fMode == "Default")
        return f * npart + (1 - f) * ncoll;
    else if (fMode == "PSD")
        return f - npart;
    else if (fMode == "Npart")
        return pow(npart, f);
    else if (fMode == "Ncoll")
        return pow(ncoll, f);
    else if (fMode == "NpartFast")
        return pow(npart, f) / pow(10, f);
    else if (fMode == "NcollFast")
        return pow(ncoll, f) / pow(100, f);
    else if (fMode == "STAR")
        return (1 - f) * npart / 2. + f * ncoll;
    else if (fMode == "HADES")
        return (1. - f * npart * npart) * npart;

    return -1.;
}

float Glauber::Fitter::NancestorsMax(float f) const
{
    const int NpartMax = fNpartHisto.GetXaxis()->GetXmax(); // some magic
    const int NcollMax = fNcollHisto.GetXaxis()->GetXmax(); // TODO

    if (fMode == "Default")
        return f * NpartMax + (1 - f) * NcollMax;
    else if (fMode == "PSD")
        return f;
    else if (fMode == "Npart")
        return pow(NpartMax, f);
    else if (fMode == "Ncoll")
        return pow(NcollMax, f);
    else if (fMode == "NpartFast")
        return pow(NpartMax, f) / pow(10, f);
    else if (fMode == "NcollFast")
        return pow(NcollMax, f) / pow(100, f);
    else if (fMode == "STAR")
        return (1 - f) * NpartMax / 2. + f * NcollMax;
    else if (fMode == "HADES")
        return (1. - f * NpartMax * NpartMax) * NpartMax;

    return -1.;
}

/*
 * take Glauber MC data from fSimTree
 * Populate fGlauberFitHisto with NBD x Na
 */

void Glauber::Fitter::SetGlauberFitHisto(float f, float mu, float k, float p, int n, Bool_t Norm2Data)
{
    fGlauberFitHisto = TH1F("glaub", "", fNbins * 1.3, 0, 1.3 * fMaxValue);
    fGlauberPlpHisto = TH1F("glplp", "", fNbins * 1.3, 0, 1.3 * fMaxValue);
    fGlauberSngHisto = TH1F("glsng", "", fNbins * 1.3, 0, 1.3 * fMaxValue);
    fGlauberPlpEv1Ev2 = TH2F("", "Multiplicity ev1 vs Multiplicity ev2;nHits 1;nHits 2", fNbins * 1.3, 0, 1.3 * fMaxValue, fNbins * 1.3, 0, 1.3 * fMaxValue);
    fB_VS_Multiplicity = TH2F("", "B VS Multiplicity;nHits;B, fm", fNbins * 1.3, 0, 1.3 * fMaxValue, 200, 0, 20);
    fNpart_VS_Multiplicity = TH2F("", "N_{part} VS Multiplicity;nHits;N_{part}", fNbins * 1.3, 0, 1.3 * fMaxValue, 10000, 0, 10000);
    fNcoll_VS_Multiplicity = TH2F("", "N_{coll} VS Multiplicity;nHits;N_{coll}", fNbins * 1.3, 0, 1.3 * fMaxValue, 10000, 0, 10000);
    fEcc1_VS_Multiplicity = TH2F("", "#epsilon1 VS Multiplicity;nHits;#epsilon1", fNbins * 1.3, 0, 1.3 * fMaxValue, 100, 0, 1);
    fPsi1_VS_Multiplicity = TH2F("", "#psi1 VS Multiplicity;nHits;#psi1", fNbins * 1.3, 0, 1.3 * fMaxValue, 2 * 3.14 / 0.01, 0, 2 * 3.14);
    fEcc2_VS_Multiplicity = TH2F("", "#epsilon2 VS Multiplicity;nHits;#epsilon2", fNbins * 1.3, 0, 1.3 * fMaxValue, 100, 0, 1);
    fPsi2_VS_Multiplicity = TH2F("", "#psi2 VS Multiplicity;nHits;#psi2", fNbins * 1.3, 0, 1.3 * fMaxValue, 2 * 3.14 / 0.01, 0, 2 * 3.14);
    fEcc3_VS_Multiplicity = TH2F("", "#epsilon3 VS Multiplicity;nHits;#epsilon3", fNbins * 1.3, 0, 1.3 * fMaxValue, 100, 0, 1);
    fPsi3_VS_Multiplicity = TH2F("", "#psi3 VS Multiplicity;nHits;#psi3", fNbins * 1.3, 0, 1.3 * fMaxValue, 2 * 3.14 / 0.01, 0, 2 * 3.14);
    fEcc4_VS_Multiplicity = TH2F("", "#epsilon4 VS Multiplicity;nHits;#epsilon4", fNbins * 1.3, 0, 1.3 * fMaxValue, 100, 0, 1);
    fPsi4_VS_Multiplicity = TH2F("", "#psi4 VS Multiplicity;nHits;#psi4", fNbins * 1.3, 0, 1.3 * fMaxValue, 2 * 3.14 / 0.01, 0, 2 * 3.14);
    fEcc5_VS_Multiplicity = TH2F("", "#epsilon5 VS Multiplicity;nHits;#epsilon5", fNbins * 1.3, 0, 1.3 * fMaxValue, 100, 0, 1);
    fPsi5_VS_Multiplicity = TH2F("", "#psi5 VS Multiplicity;nHits;#psi5", fNbins * 1.3, 0, 1.3 * fMaxValue, 2 * 3.14 / 0.01, 0, 2 * 3.14);

    fGlauberFitHisto.SetName("glaub_fit_histo");
    fGlauberPlpHisto.SetName("glaub_plp_histo");
    fGlauberSngHisto.SetName("glaub_sng_histo");
    fGlauberPlpEv1Ev2.SetName("glaub_plp_ev1ev2");
    fB_VS_Multiplicity.SetName("B_VS_Multiplicity");
    fNpart_VS_Multiplicity.SetName("Npart_VS_Multiplicity");
    fNcoll_VS_Multiplicity.SetName("Ncoll_VS_Multiplicity");
    fEcc1_VS_Multiplicity.SetName("Ecc1_VS_Multiplicity");
    fPsi1_VS_Multiplicity.SetName("Psi1_VS_Multiplicity");
    fEcc2_VS_Multiplicity.SetName("Ecc2_VS_Multiplicity");
    fPsi2_VS_Multiplicity.SetName("Psi2_VS_Multiplicity");
    fEcc3_VS_Multiplicity.SetName("Ecc3_VS_Multiplicity");
    fPsi3_VS_Multiplicity.SetName("Psi3_VS_Multiplicity");
    fEcc4_VS_Multiplicity.SetName("Ecc4_VS_Multiplicity");
    fPsi4_VS_Multiplicity.SetName("Psi4_VS_Multiplicity");
    fEcc5_VS_Multiplicity.SetName("Ecc5_VS_Multiplicity");
    fPsi5_VS_Multiplicity.SetName("Psi5_VS_Multiplicity");

    int nentries = (int)(n * (1. - p));
    int plp_counter = nentries;
#ifndef __THREADS_ON__
    BuildMultiplicity(f, mu, k, p, 0, nentries, plp_counter, n, nentries);
#endif
#ifdef __THREADS_ON__
    std::vector<std::thread> v_thr;
    std::deque<std::atomic<int>> v_progress;

    for (unsigned int i = 0; i < fNthreads; ++i)
        v_progress.emplace_back(0);

    for (unsigned int i = 0; i < fNthreads; ++i)
    {
        int n_part = (int)(n / fNthreads);
        int i_start = i * n_part;
        int i_stop = (int)((i + 1) * n_part * (1. - p));
        int p_start = i_stop;
        int p_stop = (i + 1) * n_part;
        v_thr.emplace_back([&]
                           { Glauber::Fitter::BuildMultiplicity(f, mu, k, p, i_start, i_stop, p_start, p_stop, std::ref(v_progress[i])); });
    }

    bool isOver = false;
    int tot_progress, tot_denum;
    while (not isOver)
    {
        isOver = true;
        tot_progress = 0;
        tot_denum = 0;
        if (fFirstIteration)
            std::cout << "\tGlauber::Fitter::SetGlauberFitHisto: Initialization, progress: ";
        else
            std::cout << "\tGlauber::Fitter::SetGlauberFitHisto: Constructing multiplicity, progress: ";
        for (int j = 0; j < fNthreads; ++j)
        {
            if (v_progress[j].load() == 0)
                continue;
            tot_progress += v_progress[j].load();
            tot_denum++;
        }
        if (tot_denum > 0)
            tot_progress /= tot_denum;
        std::cout << tot_progress << "% \r" << std::flush;
        if (tot_progress < 100 && tot_progress > 0)
            isOver = false;
        std::chrono::milliseconds dura(200);
        std::this_thread::sleep_for(dura);
    }

    for (auto &thread : v_thr)
        thread.join();
#endif
    if (Norm2Data)
        NormalizeGlauberFit();

    std::cout << "\t                                                                                                \r" << std::flush;
}

#ifndef __THREADS_ON__
bool Glauber::Fitter::BuildMultiplicity(float f, float mu, float k, float p, int i_start, int i_stop, int plp_start, int plp_stop, int n)
#endif
#ifdef __THREADS_ON__
    bool Glauber::Fitter::BuildMultiplicity(float f, float mu, float k, float p, int i_start, int i_stop, int plp_start, int plp_stop, std::atomic<int> &_progress)
#endif
{
    std::random_device rd;
    std::mt19937 rngnum(rd());
    std::uniform_real_distribution<float> unirnd(0., 1.);
    std::negative_binomial_distribution<> nbddist(k, (float)(k / (k + mu)));
    std::gamma_distribution<> gammadist((float)((mu * k) / (mu + k)), (float)((k + mu) / k));
    int plp_counter = plp_start;
#ifdef __THREADS_ON__
    std::lock_guard<std::mutex> guard(fMtx);
#endif
    for (int i = i_start; i < i_stop; i++)
    {
#ifndef __THREADS_ON__
        std::cout << "\tGlauber::Fitter::SetGlauberFitHisto: Constructing multiplicity, event [" << i
                  << "/" << n << "]\r" << std::flush;
#endif
#ifdef __THREADS_ON__
        _progress.store((i - i_start) * 100 / (i_stop - i_start) + 1);
#endif
        const int Na = int(Nancestors(f, fvNpart.at(i), fvNcoll.at(i)));

        float nHits{0.}, nPlp{0.};
        if (fUseNbd)
            for (int j = 0; j < Na; j++)
                nHits += nbddist(rngnum);
        else
            for (int j = 0; j < Na; j++)
                nHits += gammadist(rngnum);
        if (p > 1e-10 && unirnd(rngnum) <= p)
        {
            const int Na1 = int(Nancestors(f, fvNpart.at(plp_counter), fvNcoll.at(plp_counter)));
            if (fUseNbd)
                for (int j = 0; j < Na1; j++)
                    nPlp += nbddist(rngnum);
            else
                for (int j = 0; j < Na1; j++)
                    nPlp += gammadist(rngnum);
            plp_counter++;

            fGlauberPlpHisto.Fill(nHits + nPlp);
            fGlauberPlpEv1Ev2.Fill(nHits, nPlp);
            nHits += nPlp;
        }
        else
        {
            fGlauberSngHisto.Fill(nHits);
        }
        fGlauberFitHisto.Fill(nHits);
        fB_VS_Multiplicity.Fill(nHits, fvB.at(i));
        fNpart_VS_Multiplicity.Fill(nHits, fvNpart.at(i));
        fNcoll_VS_Multiplicity.Fill(nHits, fvNcoll.at(i));
        fEcc1_VS_Multiplicity.Fill(nHits, fvEcc1.at(i));
        fPsi1_VS_Multiplicity.Fill(nHits, fvPsi1.at(i));
        fEcc2_VS_Multiplicity.Fill(nHits, fvEcc2.at(i));
        fPsi2_VS_Multiplicity.Fill(nHits, fvPsi2.at(i));
        fEcc3_VS_Multiplicity.Fill(nHits, fvEcc3.at(i));
        fPsi3_VS_Multiplicity.Fill(nHits, fvPsi3.at(i));
        fEcc4_VS_Multiplicity.Fill(nHits, fvEcc4.at(i));
        fPsi4_VS_Multiplicity.Fill(nHits, fvPsi4.at(i));
        fEcc5_VS_Multiplicity.Fill(nHits, fvEcc5.at(i));
        fPsi5_VS_Multiplicity.Fill(nHits, fvPsi5.at(i));
    }
    return true;
}

void Glauber::Fitter::NormalizeGlauberFit()
{

    int fGlauberFitHistoInt{0};
    int fDataHistoInt{0};

    const int lowchibin = fFitMinBin;
    const int highchibin = fFitMaxBin < fNbins ? fFitMaxBin : fNbins;

    for (int i = lowchibin; i < highchibin; i++)
    {
        fGlauberFitHistoInt += fGlauberFitHisto.GetBinContent(i + 1);
        fDataHistoInt += fDataHisto.GetBinContent(i + 1);
    }

    const float ScaleFactor = (float)fDataHistoInt / fGlauberFitHistoInt;

    //     std::cout << "Scale = " << Scale << std::endl;
    fGlauberFitHisto.Scale(ScaleFactor);
    fGlauberPlpHisto.Scale(ScaleFactor);
    fGlauberSngHisto.Scale(ScaleFactor);
    fGlauberPlpEv1Ev2.Scale(ScaleFactor);
}

/**
 *
 * @param mu mean value of negative binominal distribution (we are looking for it)
 * @param chi2 return value (indicates good match)
 * @param mu_min lower search edge for mean value NBD
 * @param mu_max upper search edge for mean value NBD
 * @param f parameter of Na
 * @param k NBD parameter
 * @param nEvents
 * @param nIter
 */
void Glauber::Fitter::FindMuGoldenSection(float *mu, float *chi2, float *chi2_error, float mu_min, float mu_max, float f, float k, float p, int nEvents, int nIter, int n)
{
    const float phi = (float)((1 + TMath::Sqrt(5)) / 2);

    /* left */
    float mu_1 = mu_max - (mu_max - mu_min) / phi;

    /* right */
    float mu_2 = mu_min + (mu_max - mu_min) / phi;
    fFirstIteration = true;

    SetGlauberFitHisto(f, mu_1, k, p, nEvents);
    float chi2_mu1 = GetChi2();
    float chi2_mu1_error = GetChi2Error();

    SetGlauberFitHisto(f, mu_2, k, p, nEvents);
    float chi2_mu2 = GetChi2();
    float chi2_mu2_error = GetChi2Error();

    fFirstIteration = false;
    for (int j = 0; j < nIter; j++)
    {
        if (chi2_mu1 >= chi2_mu2)
        {
            mu_min = mu_1;
            mu_1 = mu_2;
            mu_2 = mu_min + (mu_max - mu_min) / phi;
            chi2_mu1 = chi2_mu2;
            SetGlauberFitHisto(f, mu_2, k, p, nEvents);
            chi2_mu2 = GetChi2();
            chi2_mu2_error = GetChi2Error();
        }
        else
        {
            mu_max = mu_2;
            mu_2 = mu_1;
            mu_1 = mu_max - (mu_max - mu_min) / phi;
            chi2_mu2 = chi2_mu1;
            SetGlauberFitHisto(f, mu_1, k, p, nEvents);
            chi2_mu1 = GetChi2();
            chi2_mu1_error = GetChi2Error();
        }

        std::cout << "n = " << n + j << " f = " << f << " k = " << k << " p = " << p << " mu1 = " << mu_1 << " mu2 = " << mu_2 << " chi2_mu1 = " << chi2_mu1 << " chi2_mu2 = " << chi2_mu2 << std::endl;
    }

    /* take min(mu) */
    *mu = (chi2_mu1 < chi2_mu2) ? mu_1 : mu_2;
    /* take min(chi2) */
    *chi2 = (chi2_mu1 < chi2_mu2) ? chi2_mu1 : chi2_mu2;
    /* take min(chi2_error) */
    *chi2_error = (chi2_mu1 < chi2_mu2) ? chi2_mu1_error : chi2_mu2_error;
}

/**
 * Find the best match
 *
 * @param return value of best fit parameters
 * @param f0 lower search edge for parameter of Na, for which chi2 will be calculated
 * @param f1 upper search edge for parameter of Na, for which chi2 will be calculated
 * @param k0 lower search edge for NBD parameter
 * @param k1 upper search edge for NBD parameter
 * @param nEvents
 */
float Glauber::Fitter::FitGlauber(Float_t f0, Float_t f1, Float_t k0, Float_t k1, Float_t p0, Float_t p1, Int_t nEvents)
{
    float f_fit{-1};
    float mu_fit{-1};
    float k_fit{-1};
    float p_fit{-1};
    float Chi2Min{1e10};
    float Chi2Min_error{0};

    const TString filename = Form("%s/fit_%4.2f_%4.2f_%4.2f_%4.2f_%d.root", fOutDirName.Data(), f0, k0, k1, p0, fFitMinBin);

    //     std::unique_ptr<TFile> file {TFile::Open(filename, "recreate")};
    //     std::unique_ptr<TTree> tree {new TTree("test_tree", "tree" )};

    TFile *file{TFile::Open(filename, "recreate")};
    TTree *tree{new TTree("test_tree", "tree")};

    float f, mu, k, p, chi2, chi2_error, sigma;

    tree->Branch("f", &f, "f/F");
    tree->Branch("mu", &mu, "mu/F");
    tree->Branch("k", &k, "k/F");
    tree->Branch("p", &p, "k/F");
    tree->Branch("chi2", &chi2, "chi2/F");
    tree->Branch("chi2_error", &chi2_error, "chi2_error/F");
    tree->Branch("sigma", &sigma, "sigma/F");

    int n = 1;
    for (float i = f0; i <= f1; i = i + fFstep)
    {
        f = i;
        for (float j = k0; j <= k1; j = j + fKstep)
        {
            k = j;
            for (float h = p0; h <= p1; h = h + fPstep)
            {
                p = h;
                mu = fMaxValue / NancestorsMax(f);
                const float mu_min = 0.0 * mu;
                const float mu_max = 1.0 * mu;

                if (fNiter == 0)
                    fNiter = 2;

                FindMuGoldenSection(&mu, &chi2, &chi2_error, mu_min, mu_max, f, k, p, nEvents, fNiter, n);
                n = n + fNiter;
                sigma = (mu / k + 1) * mu;

                tree->Fill();

                if (chi2 < Chi2Min)
                {
                    f_fit = f;
                    mu_fit = mu;
                    k_fit = k;
                    p_fit = p;
                    Chi2Min = chi2;
                    Chi2Min_error = chi2_error;
                    fBestFitHisto = fGlauberFitHisto;
                    fBestPlpHisto = fGlauberPlpHisto;
                    fBestSngHisto = fGlauberSngHisto;
                    fBestPlpEv1Ev2 = fGlauberPlpEv1Ev2;
                    fBestB_VS_Multiplicity = fB_VS_Multiplicity;
                    fBestNpart_VS_Multiplicity = fNpart_VS_Multiplicity;
                    fBestNcoll_VS_Multiplicity = fNcoll_VS_Multiplicity;
                    fBestEcc1_VS_Multiplicity = fEcc1_VS_Multiplicity;
                    fBestPsi1_VS_Multiplicity = fPsi1_VS_Multiplicity;
                    fBestEcc2_VS_Multiplicity = fEcc2_VS_Multiplicity;
                    fBestPsi2_VS_Multiplicity = fPsi2_VS_Multiplicity;
                    fBestEcc3_VS_Multiplicity = fEcc3_VS_Multiplicity;
                    fBestPsi3_VS_Multiplicity = fPsi3_VS_Multiplicity;
                    fBestEcc4_VS_Multiplicity = fEcc4_VS_Multiplicity;
                    fBestPsi4_VS_Multiplicity = fPsi4_VS_Multiplicity;
                    fBestEcc5_VS_Multiplicity = fEcc5_VS_Multiplicity;
                    fBestPsi5_VS_Multiplicity = fPsi5_VS_Multiplicity;
                }
            }
        }
    }

    SetNBDhist(mu_fit, k_fit);

    tree->Write();
    file->Write();
    file->Close();

    fOptimalF = f_fit;
    fOptimalK = k_fit;
    fOptimalMu = mu_fit;
    fOptimalP = p_fit;
    fOptimalChi2Ndf = Chi2Min;
    fOptimalChi2NdfError = Chi2Min_error;

    return Chi2Min;
}

/**
 * Compare fGlauberFitHisto with fDataHisto
 * @return chi2 value
 */
float Glauber::Fitter::GetChi2() const
{
    float chi2{0.0};

    const int lowchibin = fFitMinBin;
    const int highchibin = fFitMaxBin < fNbins ? fFitMaxBin : fNbins;

    for (int i = lowchibin; i <= highchibin; ++i)
    {
        if (fDataHisto.GetBinContent(i) < 1.0)
            continue;
        const float error2 = pow(fDataHisto.GetBinError(i), 2) + pow(fGlauberFitHisto.GetBinError(i), 2);
        //	std::cout<<"i="<<i<<"  DataError="<<fDataHisto.GetBinError(i)<<std::endl;
        //	std::cout<<"i="<<i<<"  GlauberError="<<fGlauberFitHisto.GetBinError(i)<<std::endl;
        const float diff2 = pow((fGlauberFitHisto.GetBinContent(i) - fDataHisto.GetBinContent(i)), 2) / error2;
        //	std::cout<<"i="<<i<<"  DataContent="<<fDataHisto.GetBinContent(i)<<std::endl;
        //	std::cout<<"i="<<i<<"  GlauberContent="<<fGlauberFitHisto.GetBinContent(i)<<std::endl;

        chi2 += diff2;
    }

    chi2 = chi2 / (highchibin - lowchibin + 1);
    return chi2;
}

/**
 * Compare fGlauberFitHisto with fDataHisto
 * @return chi2 error value
 */
float Glauber::Fitter::GetChi2Error() const
{
    float chi2_error{0.0};

    const int lowchibin = fFitMinBin;
    const int highchibin = fFitMaxBin < fNbins ? fFitMaxBin : fNbins;

    for (int i = lowchibin; i <= highchibin; ++i)
    {
        if (fDataHisto.GetBinContent(i) < 1.0)
            continue;
        const float error2 = pow(fDataHisto.GetBinError(i), 2) + pow(fGlauberFitHisto.GetBinError(i), 2);
        //	std::cout<<"i="<<i<<"  DataError="<<fDataHisto.GetBinError(i)<<std::endl;
        //	std::cout<<"i="<<i<<"  GlauberError="<<fGlauberFitHisto.GetBinError(i)<<std::endl;
        const float diff = (fGlauberFitHisto.GetBinContent(i) - fDataHisto.GetBinContent(i));
        //	std::cout<<"i="<<i<<"  DataContent="<<fDataHisto.GetBinContent(i)<<std::endl;
        //	std::cout<<"i="<<i<<"  GlauberContent="<<fGlauberFitHisto.GetBinContent(i)<<std::endl;
        const float error_diff = (fGlauberFitHisto.GetBinError(i) - fDataHisto.GetBinError(i));

        chi2_error += pow(diff * error_diff / error2, 2);
    }

    chi2_error = 2 * pow(chi2_error, 0.5) / (highchibin - lowchibin + 1);
    return chi2_error;
}

/**
 * Populates histogram nbd_<mean>_<k> with values of NBD
 * @param mu
 * @param k
 */
void Glauber::Fitter::SetNBDhist(float mu, float k)
{
    // Interface for TH1F.
    const int nBins = (mu + 1.) * 3 < 10 ? 10 : (mu + 1.) * 3;

    fNbdHisto = TH1F("fNbdHisto", "", nBins, 0, nBins);
    if (fUseNbd)
        fNbdHisto.SetName("nbd");
    else
        fNbdHisto.SetName("gamma");

    std::random_device rd;
    std::mt19937 rngnum(rd());
    std::uniform_real_distribution<float> unirnd(0., 1.);
    std::negative_binomial_distribution<> nbddist(k, (float)(k / (k + mu)));
    std::gamma_distribution<> gammadist((float)((mu * k) / (mu + k)), (float)((k + mu) / k));

    for (int i = 0; i < 1e5; ++i)
    {
        if (fUseNbd)
            fNbdHisto.Fill(nbddist(rngnum));
        else
            fNbdHisto.Fill(gammadist(rngnum));
    }
}

/**
 * Negative Binomial Distribution (by definition)
 * @param n argument
 * @param mu mean
 * @param k argument
 * @return NBD for a given parameters
 */
float Glauber::Fitter::NBD(float n, float mu, float k) const
{
    // Compute NBD.
    float F;
    float f;

    if (n + k > 100.0)
    {
        // log method for handling large numbers
        F = TMath::LnGamma(n + k) - TMath::LnGamma(n + 1.) - TMath::LnGamma(k);
        f = n * TMath::Log(mu / k) - (n + k) * TMath::Log(1.0 + mu / k);
        F += f;
        F = TMath::Exp(F);
    }
    else
    {
        F = TMath::Gamma(n + k) / (TMath::Gamma(n + 1.) * TMath::Gamma(k));
        f = n * TMath::Log(mu / k) - (n + k) * TMath::Log(1.0 + mu / k);
        f = TMath::Exp(f);
        F *= f;
    }

    return F;
}
/**
 * Creates histo with a given model parameter distribution
 * @param range observable range
 * @param name name of the MC-Glauber model parameter
 * @param par array with fit parameters
 * @param Nevents
 * @return pointer to the histogram
 */
std::unique_ptr<TH1F> Glauber::Fitter::GetModelHisto(const float range[2], TString name, int nEvents)
{
    fvModelInput.clear();

    const float p = fOptimalP;

    float modelpar{-999.};
    fSimTree->SetBranchAddress(name, &modelpar);
    for (int i = 0; i < nEvents; ++i)
    {
        fSimTree->GetEntry(i);
        fvModelInput.push_back(modelpar);
    }

    std::unique_ptr<TH1F> hModel(new TH1F("hModel", "name", 100, fSimTree->GetMinimum(name), fSimTree->GetMaximum(name)));

    int nentries = (int)(nEvents * (1. - p));
    int plp_counter = nentries;
#ifndef __THREADS_ON__
    BuildModel(range, 0, nentries, plp_counter, nEvents, nentries);
#endif
#ifdef __THREADS_ON__
    std::vector<std::thread> v_thr;
    std::deque<std::atomic<int>> v_progress;

    for (unsigned int i = 0; i < fNthreads; ++i)
        v_progress.emplace_back(0);

    for (unsigned int i = 0; i < fNthreads; ++i)
    {
        int n_part = (int)(nEvents / fNthreads);
        int i_start = i * n_part;
        int i_stop = (int)((i + 1) * n_part * (1. - p));
        int p_start = i_stop;
        int p_stop = (i + 1) * n_part;
        v_thr.emplace_back([&]
                           { Glauber::Fitter::BuildModel(range, i_start, i_stop, p_start, p_stop, std::ref(v_progress[i])); });
    }

    bool isOver = false;
    int tot_progress, tot_denum;
    while (not isOver)
    {
        isOver = true;
        tot_progress = 0;
        tot_denum = 0;
        std::cout << "\tGlauber::Fitter::GetModelHisto: Constructing multiplicity, progress: ";
        for (int j = 0; j < fNthreads; ++j)
        {
            if (v_progress[j].load() == 0)
                continue;
            tot_progress += v_progress[j].load();
            tot_denum++;
        }
        if (tot_denum > 0)
            tot_progress /= tot_denum;
        std::cout << tot_progress << "% \r" << std::flush;
        if (tot_progress < 100)
            isOver = false;
        std::chrono::milliseconds dura(200);
        std::this_thread::sleep_for(dura);
    }

    for (auto &thread : v_thr)
        thread.join();
#endif

    for (auto &mult : fvModel)
    {
        hModel->Fill(mult);
    }
    std::cout << "\t                                                                                                \r" << std::flush;

    fvModel.clear();
    fvModelInput.clear();

    // return std::move(hModel);
    return hModel;
}

#ifndef __THREADS_ON__
bool Glauber::Fitter::BuildModel(const float range[2], int i_start, int i_stop, int plp_start, int plp_stop, int n)
#endif
#ifdef __THREADS_ON__
    bool Glauber::Fitter::BuildModel(const float range[2], int i_start, int i_stop, int plp_start, int plp_stop, std::atomic<int> &_progress)
#endif
{
    const float f = fOptimalF;
    const float mu = fOptimalMu;
    const float k = fOptimalK;
    const float p = fOptimalP;

    fvModel.clear();

    std::random_device rd;
    std::mt19937 rngnum(rd());
    std::uniform_real_distribution<float> unirnd(0., 1.);
    std::negative_binomial_distribution<> nbddist(k, (float)(k / (k + mu)));
    std::gamma_distribution<> gammadist((float)((mu * k) / (mu + k)), (float)((k + mu) / k));
    int plp_counter = plp_start;
#ifdef __THREADS_ON__
    std::lock_guard<std::mutex> guard(fMtx);
#endif
    for (int i = i_start; i < i_stop; i++)
    {
#ifndef __THREADS_ON__
        std::cout << "\tGlauber::Fitter::BuildModel: Constructing multiplicity, event [" << i
                  << "/" << n << "]\r" << std::flush;
#endif
#ifdef __THREADS_ON__
        _progress.store((i - i_start) * 100 / (i_stop - i_start) + 1);
#endif
        const int Na = int(Nancestors(f, fvNpart.at(i), fvNcoll.at(i)));

        float nHits{0.};
        if (fUseNbd)
            for (int j = 0; j < Na; j++)
                nHits += nbddist(rngnum);
        else
            for (int j = 0; j < Na; j++)
                nHits += gammadist(rngnum);
        if (p > 1e-10 && unirnd(rngnum) <= p)
        {
            const int Na1 = int(Nancestors(f, fvNpart.at(plp_counter), fvNcoll.at(plp_counter)));
            if (fUseNbd)
                for (int j = 0; j < Na1; j++)
                    nHits += nbddist(rngnum);
            else
                for (int j = 0; j < Na1; j++)
                    nHits += gammadist(rngnum);
            plp_counter++;
        }
        if (nHits > range[0] && nHits < range[1])
        {
            fvModel.push_back(fvModelInput.at(i));
        }
    }
    return true;
}
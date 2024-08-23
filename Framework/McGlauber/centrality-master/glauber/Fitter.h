/** @file   Fitter.h
    @class  Glauber::Fitter
    @author Viktor Klochkov (klochkov44@gmail.com)
    @author Ilya Selyuzhenkov (ilya.selyuzhenkov@gmail.com)
    @brief  Class to fit histo with Glauber based function
*/

#ifndef GlauberFitter_H
#define GlauberFitter_H 1

#include <vector>
#include <random>
#include "TString.h"
#include "TNamed.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TTree.h"
// #include "TMinuit.h"

#ifdef __THREADS_ON__
#include <thread>
#include <mutex>
#include <atomic>
#include <deque>
#endif

namespace Glauber
{
    class Fitter
    {

    public:
        /**   Default constructor   **/
        Fitter() {};
        Fitter(std::unique_ptr<TTree> tree);
#ifdef __THREADS_ON__
        Fitter(std::unique_ptr<TTree> tree, unsigned int Nthreads);
#endif
        /**   Destructor   **/
        virtual ~Fitter() {};

        void Init(int nEntries, TString fmode);
        void SetGlauberFitHisto(Float_t f, Float_t mu, Float_t k, Float_t p, Int_t n = 10000, Bool_t Norm2Data = true);
        void NormalizeGlauberFit();
        void DrawHistos(Bool_t isSim = true, Bool_t isData = true, Bool_t isGlauber = false, Bool_t isNBD = false);

        float FitGlauber(Float_t f0, Float_t f1, Float_t k0, Float_t k1, Float_t p0, Float_t p1, Int_t nEvents);
        void FindMuGoldenSection(Float_t *mu, Float_t *chi2, float *chi2_error, Float_t mu_min, Float_t mu_max, Float_t f, Float_t k, Float_t p, Int_t nEvents = 10000, Int_t nIter = 5, int n = 0);

        Float_t GetChi2(void) const;
        Float_t GetChi2Error(void) const;

        Float_t NBD(Float_t n, Float_t mu, Float_t k) const;
        void SetNBDhist(Float_t mu, Float_t k);

        float Nancestors(float f) const;
        float Nancestors(float f, float npart, float ncoll) const;
        float NancestorsMax(float f) const;
#ifndef __THREADS_ON__
        bool BuildMultiplicity(float f, float mu, float k, float p, int i_start, int i_stop, int plp_start, int plp_stop, int n);
#endif
#ifdef __THREADS_ON__
        bool BuildMultiplicity(float f, float mu, float k, float p, int i_start, int i_stop, int plp_start, int plp_stop, std::atomic<int> &_progress);
#endif
#ifndef __THREADS_ON__
        bool BuildModel(const float range[2], int i_start, int i_stop, int plp_start, int plp_stop, int n);
#endif
#ifdef __THREADS_ON__
        bool BuildModel(const float range[2], int i_start, int i_stop, int plp_start, int plp_stop, std::atomic<int> &_progress);
#endif

#ifdef __THREADS_ON__
        void SetNthreads(unsigned int n) { fNthreads = n; }
        unsigned int GetNthreads() const { return fNthreads; }
#endif

        void UseGamma() { fUseNbd = false; }
        void UseNbd() { fUseNbd = true; }

        std::unique_ptr<TH1F> GetModelHisto(const Float_t range[2], TString name, Int_t nEvents);

        //
        //         Setters
        //
        void SetInputHisto(const TH1F &h) { fDataHisto = h; }
        void SetFitMinBin(Int_t min) { fFitMinBin = min; }
        void SetFitMaxBin(Int_t min) { fFitMaxBin = min; }
        void SetNormMinBin(Int_t min) { fNormMinBin = min; }
        void SetBinSize(Float_t size) { fBinSize = size; }
        void SetOutDirName(TString name) { fOutDirName = name; }
        void SetMode(const TString mode) { fMode = mode; }
        void SetMassNumber(Float_t A) { fA = A; }
        void SetNiter(Int_t a) { fNiter = a; }
        void SetFstepSize(Float_t a) { fFstep = a; }
        void SetKstepSize(Float_t a) { fKstep = a; }
        void SetPstepSize(Float_t a) { fPstep = a; }

        //
        //         Getters
        //
        TH1F GetGlauberFitHisto() const { return fGlauberFitHisto; }
        TH1F GetGlauberPlpHisto() const { return fGlauberPlpHisto; }
        TH1F GetGlauberSngHisto() const { return fGlauberSngHisto; }
        TH1F GetDataHisto() const { return fDataHisto; }
        TH1F GetNBDHisto() const { return fNbdHisto; }
        TH1F GetBHisto() const { return fBHisto; }
        TH1F GetNpartHisto() const { return fNpartHisto; }
        TH1F GetNcollHisto() const { return fNcollHisto; }
        TH1F GetBestFitHisto() const { return fBestFitHisto; }
        TH1F GetBestPlpHisto() const { return fBestPlpHisto; }
        TH1F GetBestSngHisto() const { return fBestSngHisto; }
        TH1F GetEcc1Histo() const { return fEcc1Histo; }
        TH1F GetPsi1Histo() const { return fPsi1Histo; }
        TH1F GetEcc2Histo() const { return fEcc2Histo; }
        TH1F GetPsi2Histo() const { return fPsi2Histo; }
        TH1F GetEcc3Histo() const { return fEcc3Histo; }
        TH1F GetPsi3Histo() const { return fPsi3Histo; }
        TH1F GetEcc4Histo() const { return fEcc4Histo; }
        TH1F GetPsi4Histo() const { return fPsi4Histo; }
        TH1F GetEcc5Histo() const { return fEcc5Histo; }
        TH1F GetPsi5Histo() const { return fPsi5Histo; }

        TH2F GetGlauberPlpEv1Ev2() const { return fGlauberPlpEv1Ev2; }
        TH2F GetB_VS_Multiplicity() const { return fB_VS_Multiplicity; }
        TH2F GetNpart_VS_Multiplicity() const { return fNpart_VS_Multiplicity; }
        TH2F GetNcoll_VS_Multiplicity() const { return fNcoll_VS_Multiplicity; }
        TH2F GetEcc1_VS_Multiplicity() const { return fEcc1_VS_Multiplicity; }
        TH2F GetPsi1_VS_Multiplicity() const { return fPsi1_VS_Multiplicity; }
        TH2F GetEcc2_VS_Multiplicity() const { return fEcc2_VS_Multiplicity; }
        TH2F GetPsi2_VS_Multiplicity() const { return fPsi2_VS_Multiplicity; }
        TH2F GetEcc3_VS_Multiplicity() const { return fEcc3_VS_Multiplicity; }
        TH2F GetPsi3_VS_Multiplicity() const { return fPsi3_VS_Multiplicity; }
        TH2F GetEcc4_VS_Multiplicity() const { return fEcc4_VS_Multiplicity; }
        TH2F GetPsi4_VS_Multiplicity() const { return fPsi4_VS_Multiplicity; }
        TH2F GetEcc5_VS_Multiplicity() const { return fEcc5_VS_Multiplicity; }
        TH2F GetPsi5_VS_Multiplicity() const { return fPsi5_VS_Multiplicity; }

        TH2F GetBestPlpEv1Ev2() const { return fBestPlpEv1Ev2; }
        TH2F GetBestB_VS_Multiplicity() const { return fBestB_VS_Multiplicity; }
        TH2F GetBestNpart_VS_Multiplicity() const { return fBestNpart_VS_Multiplicity; }
        TH2F GetBestNcoll_VS_Multiplicity() const { return fBestNcoll_VS_Multiplicity; }
        TH2F GetBestEcc1_VS_Multiplicity() const { return fBestEcc1_VS_Multiplicity; }
        TH2F GetBestPsi1_VS_Multiplicity() const { return fBestPsi1_VS_Multiplicity; }
        TH2F GetBestEcc2_VS_Multiplicity() const { return fBestEcc2_VS_Multiplicity; }
        TH2F GetBestPsi2_VS_Multiplicity() const { return fBestPsi2_VS_Multiplicity; }
        TH2F GetBestEcc3_VS_Multiplicity() const { return fBestEcc3_VS_Multiplicity; }
        TH2F GetBestPsi3_VS_Multiplicity() const { return fBestPsi3_VS_Multiplicity; }
        TH2F GetBestEcc4_VS_Multiplicity() const { return fBestEcc4_VS_Multiplicity; }
        TH2F GetBestPsi4_VS_Multiplicity() const { return fBestPsi4_VS_Multiplicity; }
        TH2F GetBestEcc5_VS_Multiplicity() const { return fBestEcc5_VS_Multiplicity; }
        TH2F GetBestPsi5_VS_Multiplicity() const { return fBestPsi5_VS_Multiplicity; }

        float GetOptimalF() const { return fOptimalF; }
        float GetOptimalK() const { return fOptimalK; }
        float GetOptimalMu() const { return fOptimalMu; }
        float GetOptimalP() const { return fOptimalP; }
        float GetOptimalChi2() const { return fOptimalChi2Ndf; }
        float GetOptimalChi2Error() const { return fOptimalChi2NdfError; }

    private:
        /**   Data members  **/
        TH1F fBHisto;
        TH1F fNpartHisto;
        TH1F fNcollHisto;
        TH1F fEcc1Histo;
        TH1F fPsi1Histo;
        TH1F fEcc2Histo;
        TH1F fPsi2Histo;
        TH1F fEcc3Histo;
        TH1F fPsi3Histo;
        TH1F fEcc4Histo;
        TH1F fPsi4Histo;
        TH1F fEcc5Histo;
        TH1F fPsi5Histo;

        TH1F fDataHisto;
        TH1F fNbdHisto;
        TH1F fGlauberFitHisto;
        TH1F fGlauberPlpHisto;
        TH1F fGlauberSngHisto;
        TH1F fBestFitHisto;
        TH1F fBestPlpHisto;
        TH1F fBestSngHisto;
        Int_t fNiter;
        Float_t fFstep;
        Float_t fKstep;
        Float_t fPstep;

        TH2F fGlauberPlpEv1Ev2;
        TH2F fB_VS_Multiplicity;
        TH2F fNpart_VS_Multiplicity;
        TH2F fNcoll_VS_Multiplicity;
        TH2F fEcc1_VS_Multiplicity;
        TH2F fPsi1_VS_Multiplicity;
        TH2F fEcc2_VS_Multiplicity;
        TH2F fPsi2_VS_Multiplicity;
        TH2F fEcc3_VS_Multiplicity;
        TH2F fPsi3_VS_Multiplicity;
        TH2F fEcc4_VS_Multiplicity;
        TH2F fPsi4_VS_Multiplicity;
        TH2F fEcc5_VS_Multiplicity;
        TH2F fPsi5_VS_Multiplicity;

        TH2F fBestPlpEv1Ev2;
        TH2F fBestB_VS_Multiplicity;
        TH2F fBestNpart_VS_Multiplicity;
        TH2F fBestNcoll_VS_Multiplicity;
        TH2F fBestEcc1_VS_Multiplicity;
        TH2F fBestPsi1_VS_Multiplicity;
        TH2F fBestEcc2_VS_Multiplicity;
        TH2F fBestPsi2_VS_Multiplicity;
        TH2F fBestEcc3_VS_Multiplicity;
        TH2F fBestPsi3_VS_Multiplicity;
        TH2F fBestEcc4_VS_Multiplicity;
        TH2F fBestPsi4_VS_Multiplicity;
        TH2F fBestEcc5_VS_Multiplicity;
        TH2F fBestPsi5_VS_Multiplicity;

        /* MC data */
        std::unique_ptr<TTree> fSimTree{nullptr};

        Float_t fA{-1.}; // mass number
        Float_t fB{-1.};
        Float_t fNpart{-1.};
        Float_t fNcoll{-1.};
        Float_t fEcc1{-1.};
        Float_t fPsi1{-1.};
        Float_t fEcc2{-1.};
        Float_t fPsi2{-1.};
        Float_t fEcc3{-1.};
        Float_t fPsi3{-1.};
        Float_t fEcc4{-1.};
        Float_t fPsi4{-1.};
        Float_t fEcc5{-1.};
        Float_t fPsi5{-1.};

        std::vector<float> fvB{};
        std::vector<float> fvNpart{};
        std::vector<float> fvNcoll{};
        std::vector<float> fvEcc1{};
        std::vector<float> fvPsi1{};
        std::vector<float> fvEcc2{};
        std::vector<float> fvPsi2{};
        std::vector<float> fvEcc3{};
        std::vector<float> fvPsi3{};
        std::vector<float> fvEcc4{};
        std::vector<float> fvPsi4{};
        std::vector<float> fvEcc5{};
        std::vector<float> fvPsi5{};

        std::vector<float> fvModel{};
        std::vector<float> fvModelInput{};

        Float_t fMaxValue{-1.};

        Int_t fNbins{-1};
        Float_t fBinSize{1};

        Int_t fFitMinBin{-1};
        Int_t fFitMaxBin{-1};

        Int_t fNormMinBin{-1};

        TString fMode{"Default"};

        TString fOutDirName{""};

        float fOptimalF{0.};
        float fOptimalK{0.};
        float fOptimalMu{0.};
        float fOptimalP{0.};
        float fOptimalChi2Ndf{0.};
        float fOptimalChi2NdfError{0.};

        bool fUseNbd{false};
        bool fFirstIteration{false};

#ifdef __THREADS_ON__
        std::mutex fMtx;
        unsigned int fNthreads{std::thread::hardware_concurrency()};
#endif

        ClassDef(Fitter, 2);
    };
}

#endif

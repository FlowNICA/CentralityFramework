#include "Fitter.h"
#include "TAxis.h"
#include "TFile.h"
#include "TMath.h"
#include "TTree.h"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <random>
#include <string>
#ifdef __THREADS_ON__
#include <chrono>
#include <thread>
#endif

namespace {

/* Functional forms of the number of ancestors, see Fitter::SetMode */
enum class NaMode {
  kDefault,
  kPSD,
  kNpart,
  kNcoll,
  kNpartFast,
  kNcollFast,
  kSTAR,
  kHADES,
  kUnknown
};

NaMode ParseMode(const TString &mode) {
  if (mode == "Default")
    return NaMode::kDefault;
  if (mode == "PSD")
    return NaMode::kPSD;
  if (mode == "Npart")
    return NaMode::kNpart;
  if (mode == "Ncoll")
    return NaMode::kNcoll;
  if (mode == "NpartFast")
    return NaMode::kNpartFast;
  if (mode == "NcollFast")
    return NaMode::kNcollFast;
  if (mode == "STAR")
    return NaMode::kSTAR;
  if (mode == "HADES")
    return NaMode::kHADES;
  return NaMode::kUnknown;
}

inline float Nancestors(NaMode mode, float f, float npart, float ncoll) {
  switch (mode) {
  case NaMode::kDefault:
    return f * npart + (1 - f) * ncoll;
  case NaMode::kPSD:
    return f - npart;
  case NaMode::kNpart:
    return pow(npart, f);
  case NaMode::kNcoll:
    return pow(ncoll, f);
  case NaMode::kNpartFast:
    return pow(npart, f) / pow(10, f);
  case NaMode::kNcollFast:
    return pow(ncoll, f) / pow(100, f);
  case NaMode::kSTAR:
    return (1 - f) * npart / 2. + f * ncoll;
  case NaMode::kHADES:
    return (1. - f * npart * npart) * npart;
  default:
    return -1.;
  }
}

/* Sum of n i.i.d. Gamma(alpha, theta) draws is one Gamma(n*alpha, theta) draw */
float SumOfGammas(int n, std::gamma_distribution<> &gammadist,
                  std::mt19937 &rngnum) {
  if (n <= 0)
    return 0.;
  const auto &par = gammadist.param();
  return gammadist(rngnum, std::gamma_distribution<>::param_type(
                               n * par.alpha(), par.beta()));
}

/* Gamma distribution with mean mu and NBD-like parameter k */
std::gamma_distribution<> MakeGamma(float mu, float k) {
  return std::gamma_distribution<>((float)((mu * k) / (mu + k)),
                                   (float)((k + mu) / k));
}

/*
 * Scan points v0, v0 + step, ... <= v1 (computed by index, not by
 * accumulation). A non-positive step or v1 <= v0 gives the single point v0,
 * too many points give an empty scan.
 */
std::vector<float> ScanPoints(float v0, float v1, float step,
                              const char *name) {
  const long kMaxPoints{100000};
  if (v1 <= v0)
    return {v0};
  if (!(step > 0.)) {
    std::cout << "FitGlauber: *** Warning - " << name << " step is " << step
              << ", using only " << name << " = " << v0 << std::endl;
    return {v0};
  }
  const double n_points = std::floor((v1 - v0) / step + 1e-4) + 1;
  if (!(n_points <= kMaxPoints)) {
    std::cout << "FitGlauber: *** Error - " << n_points << " points for "
              << name << " (step " << step << "), maximum is " << kMaxPoints
              << std::endl;
    return {};
  }
  const long n = (long)n_points;
  std::vector<float> points;
  points.reserve(n);
  for (long i = 0; i < n; i++)
    points.push_back(v0 + i * step);
  return points;
}

/*
 * Histogram over [min, max] with bins of the given width; a degenerate range
 * (e.g. a variable which is always 0) still gives one bin of that width
 */
TH1F MakeRangeHisto(const char *name, const char *title, float min, float max,
                    float width) {
  if (!(max > min))
    max = min + width;
  const int nbins = std::max(1, (int)((max - min) / width));
  return TH1F(name, title, nbins, min, max);
}

/* Text progress bar, redrawn in place only when the percentage changes */
class ProgressBar {
public:
  ProgressBar(std::string label, long total)
      : fLabel(std::move(label)), fTotal(total > 0 ? total : 1) {}

  void Print(long done) {
    if (done > fTotal)
      done = fTotal;
    const int percent = (int)(100 * done / fTotal);
    if (percent == fLastPercent)
      return;
    fLastPercent = percent;
    const int filled = (int)(kWidth * done / fTotal);
    std::cout << "\t" << fLabel << " [" << std::string(filled, '#')
              << std::string(kWidth - filled, '.') << "] " << std::setw(3)
              << percent << "%\r" << std::flush;
  }

  void Finish() {
    Print(fTotal);
    std::cout << std::endl;
  }

private:
  static constexpr int kWidth{50};
  std::string fLabel;
  long fTotal;
  int fLastPercent{-1};
};

/*
 * Progress shared by workers: printed by the worker itself in single-thread
 * mode and by the main thread (RunWorkers) in multi-thread mode
 */
class Progress {
public:
  static constexpr int kStep{1024}; // events between progress updates

  Progress(std::string label, long total) : fBar(std::move(label), total) {}

  void Add(long n) {
    fDone += n;
#ifndef __THREADS_ON__
    fBar.Print(fDone);
#endif
  }
  long Done() const { return fDone; }
  ProgressBar &Bar() { return fBar; }

private:
  ProgressBar fBar;
  std::atomic<long> fDone{0};
};

/* Run work(i_worker, n_workers) on n_workers threads and show progress */
template <class Work>
void RunWorkers(unsigned int n_workers, Progress &progress, Work work) {
#ifndef __THREADS_ON__
  (void)n_workers;
  work(0u, 1u);
#endif
#ifdef __THREADS_ON__
  std::atomic<unsigned int> n_finished{0};
  std::vector<std::thread> threads;
  for (unsigned int i = 0; i < n_workers; i++)
    threads.emplace_back([&work, &n_finished, i, n_workers] {
      work(i, n_workers);
      n_finished++;
    });
  while (n_finished < n_workers) {
    progress.Bar().Print(progress.Done());
    std::this_thread::sleep_for(std::chrono::milliseconds(100));
  }
  for (auto &thread : threads)
    thread.join();
#endif
  progress.Bar().Finish();
}

/* Pile-up partner index from [plp_start, plp_stop), cycling; -1 if empty */
inline int NextPileUp(int &plp_counter, int plp_start, int plp_stop) {
  if (plp_start >= plp_stop)
    return -1;
  if (plp_counter < plp_start || plp_counter >= plp_stop)
    plp_counter = plp_start;
  return plp_counter++;
}

/* Multiplicity of a single simulated event */
struct SimEvent {
  int i;       // index of the Glauber event
  float nHits; // total multiplicity (with pile-up)
  float nPlp;  // multiplicity of the pile-up partner
  bool isPlp;  // pile-up happened
};

} // namespace

ClassImp(Glauber::Fitter)

    // -----   Default constructor   -------------------------------------------
    Glauber::Fitter::Fitter(std::unique_ptr<TTree> tree) {
  fSimTree = std::move(tree);

  if (!fSimTree) {
    std::cout << "SetSimHistos: *** Error - " << std::endl;
    exit(EXIT_FAILURE);
  }
  std::cout << fSimTree->GetEntries() << std::endl;

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
    : Fitter(std::move(tree)) {
  fNthreads = Nthreads;
}
#endif

void Glauber::Fitter::Init(int nEntries, TString fmode) {

  if (nEntries < 0 || nEntries > fSimTree->GetEntries()) {
    std::cout << "Init: *** ERROR - number of entries < 0 or less that number "
                 "of entries in input tree"
              << std::endl;
    std::cout << "Init: *** number of entries in input tree = "
              << fSimTree->GetEntries() << std::endl;
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

  fBHisto = MakeRangeHisto("fBHisto", "B", 0, BMax, fBinSize);
  fNpartHisto = MakeRangeHisto("fNpartHisto", "Npart", 0, NpartMax, fBinSize);
  fNcollHisto = MakeRangeHisto("fNcollHisto", "Ncoll", 0, NcollMax, fBinSize);
  fEcc1Histo =
      MakeRangeHisto("fEcc1Histo", "#epsilon1", Ecc1Min, Ecc1Max, 0.01);
  fPsi1Histo = MakeRangeHisto("fPsi1Histo", "#psi1", Psi1Min, Psi1Max, 0.01);
  fEcc2Histo =
      MakeRangeHisto("fEcc2Histo", "#epsilon2", Ecc2Min, Ecc2Max, 0.01);
  fPsi2Histo = MakeRangeHisto("fPsi2Histo", "#psi2", Psi2Min, Psi2Max, 0.01);
  fEcc3Histo =
      MakeRangeHisto("fEcc3Histo", "#epsilon3", Ecc3Min, Ecc3Max, 0.01);
  fPsi3Histo = MakeRangeHisto("fPsi3Histo", "#psi3", Psi3Min, Psi3Max, 0.01);
  fEcc4Histo =
      MakeRangeHisto("fEcc4Histo", "#epsilon4", Ecc4Min, Ecc4Max, 0.01);
  fPsi4Histo = MakeRangeHisto("fPsi4Histo", "#psi4", Psi4Min, Psi4Max, 0.01);
  fEcc5Histo =
      MakeRangeHisto("fEcc5Histo", "#epsilon5", Ecc5Min, Ecc5Max, 0.01);
  fPsi5Histo = MakeRangeHisto("fPsi5Histo", "#psi5", Psi5Min, Psi5Max, 0.01);

  for (int i = 0; i < nEntries; i++) {
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

  while (fNbins > 1 && fDataHisto.GetBinContent(fNbins - 1) == 0)
    fNbins--;

  if (fNbins <= 1) {
    std::cout << "Init: *** ERROR - input data histogram is empty" << std::endl;
    exit(EXIT_FAILURE);
  }

  fNbins++;

  const float min = fDataHisto.GetXaxis()->GetXmin();
  const float max = fDataHisto.GetXaxis()->GetXmax();

  fMaxValue = min + (max - min) * fNbins / fDataHisto.GetNbinsX();

  std::cout << "fNbins = " << fNbins << std::endl;
  std::cout << "fMaxValue = " << fMaxValue << std::endl;

#ifdef __THREADS_ON__
  std::cout << std::endl;
  std::cout << "Number of threads found: " << NumThreads() << std::endl;
  std::cout << std::endl;
#endif
}

float Glauber::Fitter::Nancestors(float f, float npart, float ncoll) const {
  return ::Nancestors(ParseMode(fMode), f, npart, ncoll);
}

float Glauber::Fitter::NancestorsMax(float f) const {
  const int NpartMax = fNpartHisto.GetXaxis()->GetXmax(); // some magic
  const int NcollMax = fNcollHisto.GetXaxis()->GetXmax(); // TODO

  const NaMode mode = ParseMode(fMode);
  if (mode == NaMode::kPSD)
    return f;
  return ::Nancestors(mode, f, NpartMax, NcollMax);
}

/*
 * take Glauber MC data from fSimTree
 * Populate fGlauberFitHisto with NBD x Na
 */

void Glauber::Fitter::SetGlauberFitHisto(float f, float mu, float k, float p,
                                         int n, Bool_t Norm2Data) {
  fGlauberFitHisto = TH1F("glaub", "", fNbins * 1.3, 0, 1.3 * fMaxValue);
  fGlauberPlpHisto = TH1F("glplp", "", fNbins * 1.3, 0, 1.3 * fMaxValue);
  fGlauberSngHisto = TH1F("glsng", "", fNbins * 1.3, 0, 1.3 * fMaxValue);
  fGlauberPlpEv1Ev2 =
      TH2F("", "Multiplicity ev1 vs Multiplicity ev2;nHits 1;nHits 2",
           fNbins * 1.3, 0, 1.3 * fMaxValue, fNbins * 1.3, 0, 1.3 * fMaxValue);
  fB_VS_Multiplicity = TH2F("", "B VS Multiplicity;nHits;B, fm", fNbins * 1.3,
                            0, 1.3 * fMaxValue, 200, 0, 20);
  fNpart_VS_Multiplicity =
      TH2F("", "N_{part} VS Multiplicity;nHits;N_{part}", fNbins * 1.3, 0,
           1.3 * fMaxValue, 10000, 0, 10000);
  fNcoll_VS_Multiplicity =
      TH2F("", "N_{coll} VS Multiplicity;nHits;N_{coll}", fNbins * 1.3, 0,
           1.3 * fMaxValue, 10000, 0, 10000);
  fEcc1_VS_Multiplicity = TH2F("", "#epsilon1 VS Multiplicity;nHits;#epsilon1",
                               fNbins * 1.3, 0, 1.3 * fMaxValue, 100, 0, 1);
  fPsi1_VS_Multiplicity =
      TH2F("", "#psi1 VS Multiplicity;nHits;#psi1", fNbins * 1.3, 0,
           1.3 * fMaxValue, 2 * 3.14 / 0.01, 0, 2 * 3.14);
  fEcc2_VS_Multiplicity = TH2F("", "#epsilon2 VS Multiplicity;nHits;#epsilon2",
                               fNbins * 1.3, 0, 1.3 * fMaxValue, 100, 0, 1);
  fPsi2_VS_Multiplicity =
      TH2F("", "#psi2 VS Multiplicity;nHits;#psi2", fNbins * 1.3, 0,
           1.3 * fMaxValue, 2 * 3.14 / 0.01, 0, 2 * 3.14);
  fEcc3_VS_Multiplicity = TH2F("", "#epsilon3 VS Multiplicity;nHits;#epsilon3",
                               fNbins * 1.3, 0, 1.3 * fMaxValue, 100, 0, 1);
  fPsi3_VS_Multiplicity =
      TH2F("", "#psi3 VS Multiplicity;nHits;#psi3", fNbins * 1.3, 0,
           1.3 * fMaxValue, 2 * 3.14 / 0.01, 0, 2 * 3.14);
  fEcc4_VS_Multiplicity = TH2F("", "#epsilon4 VS Multiplicity;nHits;#epsilon4",
                               fNbins * 1.3, 0, 1.3 * fMaxValue, 100, 0, 1);
  fPsi4_VS_Multiplicity =
      TH2F("", "#psi4 VS Multiplicity;nHits;#psi4", fNbins * 1.3, 0,
           1.3 * fMaxValue, 2 * 3.14 / 0.01, 0, 2 * 3.14);
  fEcc5_VS_Multiplicity = TH2F("", "#epsilon5 VS Multiplicity;nHits;#epsilon5",
                               fNbins * 1.3, 0, 1.3 * fMaxValue, 100, 0, 1);
  fPsi5_VS_Multiplicity =
      TH2F("", "#psi5 VS Multiplicity;nHits;#psi5", fNbins * 1.3, 0,
           1.3 * fMaxValue, 2 * 3.14 / 0.01, 0, 2 * 3.14);

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

  n = std::min(n, (int)fvNpart.size());
  const NaMode mode = ParseMode(fMode);
  const unsigned int n_workers = NumThreads();

  /* Multiplicities are simulated in parallel, histograms are filled after */
  std::vector<std::vector<SimEvent>> sim(n_workers);
  Progress progress("Glauber::Fitter::SetGlauberFitHisto",
                    (long)(n * (1. - p)));

  RunWorkers(n_workers, progress,
             [&](unsigned int i_worker, unsigned int n_workers_) {
               /* events [i_start, plp_stop) of the worker: main events
                * first, then the pile-up pool */
               const int i_start = (int)((long long)i_worker * n / n_workers_);
               const int plp_stop =
                   (int)((long long)(i_worker + 1) * n / n_workers_);
               const int i_stop =
                   i_start + (int)((plp_stop - i_start) * (1. - p));
               int plp_counter = i_stop;

               std::random_device rd;
               std::mt19937 rngnum(rd());
               std::uniform_real_distribution<float> unirnd(0., 1.);
               auto gammadist = MakeGamma(mu, k);

               auto &out = sim[i_worker];
               out.reserve(i_stop - i_start);
               for (int i = i_start; i < i_stop; i++) {
                 if ((i - i_start + 1) % Progress::kStep == 0)
                   progress.Add(Progress::kStep);
                 const int Na =
                     int(::Nancestors(mode, f, fvNpart[i], fvNcoll[i]));

                 SimEvent ev{i, SumOfGammas(Na, gammadist, rngnum), 0., false};
                 if (p > 1e-10 && unirnd(rngnum) <= p) {
                   const int j = NextPileUp(plp_counter, i_stop, plp_stop);
                   if (j >= 0) {
                     const int Na1 =
                         int(::Nancestors(mode, f, fvNpart[j], fvNcoll[j]));
                     ev.nPlp = SumOfGammas(Na1, gammadist, rngnum);
                     ev.isPlp = true;
                   }
                 }
                 ev.nHits += ev.nPlp;
                 out.push_back(ev);
               }
               progress.Add((i_stop - i_start) % Progress::kStep);
             });

  for (const auto &out : sim)
    for (const auto &ev : out) {
      const int i = ev.i;
      const float nHits = ev.nHits;
      if (ev.isPlp) {
        fGlauberPlpHisto.Fill(nHits);
        fGlauberPlpEv1Ev2.Fill(nHits - ev.nPlp, ev.nPlp);
      } else {
        fGlauberSngHisto.Fill(nHits);
      }
      fGlauberFitHisto.Fill(nHits);
      fB_VS_Multiplicity.Fill(nHits, fvB[i]);
      fNpart_VS_Multiplicity.Fill(nHits, fvNpart[i]);
      fNcoll_VS_Multiplicity.Fill(nHits, fvNcoll[i]);
      fEcc1_VS_Multiplicity.Fill(nHits, fvEcc1[i]);
      fPsi1_VS_Multiplicity.Fill(nHits, fvPsi1[i]);
      fEcc2_VS_Multiplicity.Fill(nHits, fvEcc2[i]);
      fPsi2_VS_Multiplicity.Fill(nHits, fvPsi2[i]);
      fEcc3_VS_Multiplicity.Fill(nHits, fvEcc3[i]);
      fPsi3_VS_Multiplicity.Fill(nHits, fvPsi3[i]);
      fEcc4_VS_Multiplicity.Fill(nHits, fvEcc4[i]);
      fPsi4_VS_Multiplicity.Fill(nHits, fvPsi4[i]);
      fEcc5_VS_Multiplicity.Fill(nHits, fvEcc5[i]);
      fPsi5_VS_Multiplicity.Fill(nHits, fvPsi5[i]);
    }

  if (Norm2Data)
    NormalizeGlauberFit();
}

void Glauber::Fitter::NormalizeGlauberFit() {

  int fGlauberFitHistoInt{0};
  int fDataHistoInt{0};

  const int lowchibin = fFitMinBin;
  const int highchibin = fFitMaxBin < fNbins ? fFitMaxBin : fNbins;

  for (int i = lowchibin; i < highchibin; i++) {
    fGlauberFitHistoInt += fGlauberFitHisto.GetBinContent(i + 1);
    fDataHistoInt += fDataHisto.GetBinContent(i + 1);
  }

  if (fGlauberFitHistoInt == 0) {
    std::cout << "NormalizeGlauberFit: *** Warning - empty model histogram in "
                 "the fit range"
              << std::endl;
    return;
  }

  const float ScaleFactor = (float)fDataHistoInt / fGlauberFitHistoInt;

  fGlauberFitHisto.Scale(ScaleFactor);
  fGlauberPlpHisto.Scale(ScaleFactor);
  fGlauberSngHisto.Scale(ScaleFactor);
  fGlauberPlpEv1Ev2.Scale(ScaleFactor);
}

/**
 * Find the best match
 *
 * All (f, k, p) grid points are fitted simultaneously: in each golden section
 * iteration the multiplicity distributions for all grid points are built in
 * one pass over the Glauber events.
 *
 * @param f0 lower search edge for parameter of Na, for which chi2 will be
 * calculated
 * @param f1 upper search edge for parameter of Na, for which chi2 will be
 * calculated
 * @param k0 lower search edge for NBD parameter
 * @param k1 upper search edge for NBD parameter
 * @param p0 lower search edge for pile-up probability
 * @param p1 upper search edge for pile-up probability
 * @param nEvents number of Glauber events used to build multiplicity
 * @return best chi2/ndf
 */
float Glauber::Fitter::FitGlauber(Float_t f0, Float_t f1, Float_t k0,
                                  Float_t k1, Float_t p0, Float_t p1,
                                  Int_t nEvents) {
  float f_fit{-1};
  float mu_fit{-1};
  float k_fit{-1};
  float p_fit{-1};
  float Chi2Min{1e10};
  float Chi2Min_error{0};
  bool isFitted{false};

  const NaMode mode = ParseMode(fMode);
  if (mode == NaMode::kUnknown) {
    std::cout << "FitGlauber: *** Error - unknown mode " << fMode << std::endl;
    return Chi2Min;
  }
  if (nEvents > (int)fvNpart.size()) {
    std::cout << "FitGlauber: *** Warning - only " << fvNpart.size()
              << " Glauber events are loaded by Init, using them" << std::endl;
    nEvents = fvNpart.size();
  }
  if (fNiter <= 0)
    fNiter = 2;

  const TString filename = Form("%s/fit_%4.2f_%4.2f_%4.2f_%4.2f_%d.root",
                                fOutDirName.Data(), f0, k0, k1, p0, fFitMinBin);

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

  /* Golden section state of a single (f, k, p) grid point */
  struct GridPoint {
    float f, k, p;
    float mu_min, mu_max, mu_1, mu_2;
    float chi2_mu1, chi2_mu2, chi2_mu1_error, chi2_mu2_error;
    bool valid;
  };
  /* Multiplicity to be built with a given mu for a given grid point */
  struct Evaluation {
    int g;
    float mu;
    bool isMu2;
  };

  std::vector<GridPoint> grid;
  for (float fi : ScanPoints(f0, f1, fFstep, "f"))
    for (float ki : ScanPoints(k0, k1, fKstep, "k"))
      for (float pi : ScanPoints(p0, p1, fPstep, "p")) {
        GridPoint gp{};
        gp.f = fi;
        gp.k = ki;
        gp.p = pi;
        gp.mu_min = 0.;
        gp.mu_max = fMaxValue / NancestorsMax(gp.f);
        gp.chi2_mu1 = gp.chi2_mu2 = 1e10;
        gp.valid = std::isfinite(gp.mu_max) && gp.mu_max > 0. && gp.k > 0. &&
                   gp.p >= 0. && gp.p < 1.;
        if (!gp.valid)
          std::cout << "FitGlauber: *** Warning - skipping f = " << gp.f
                    << " k = " << gp.k << " p = " << gp.p
                    << " (mu_max = " << gp.mu_max << ")" << std::endl;
        grid.push_back(gp);
      }

  const float phi = (float)((1 + TMath::Sqrt(5)) / 2);

  /* Model histogram binning, same as in SetGlauberFitHisto */
  const TAxis axis(fNbins * 1.3, 0, 1.3 * fMaxValue);
  const int nBinsModel = axis.GetNbins() + 2;
  const double xmin = axis.GetXmin();
  const double xmax = axis.GetXmax();
  const double binsPerUnit = axis.GetNbins() / (xmax - xmin);
  /* same as TAxis::FindFixBin for fixed bins, without the function call */
  auto FindBin = [&](double x) {
    if (!(x >= xmin))
      return 0;
    if (!(x < xmax))
      return nBinsModel - 1;
    return std::min(1 + (int)((x - xmin) * binsPerUnit), nBinsModel - 2);
  };

  const int lowchibin = fFitMinBin;
  const int highchibin = fFitMaxBin < fNbins ? fFitMaxBin : fNbins;

  /* Same as NormalizeGlauberFit + chi2/ndf and its error, on raw counts */
  auto Chi2FromCounts = [&](const std::vector<float> &counts, float &chi2_out,
                            float &chi2_error_out) {
    int modelInt{0};
    int dataInt{0};
    for (int i = lowchibin; i < highchibin; i++) {
      modelInt += counts.at(i + 1);
      dataInt += fDataHisto.GetBinContent(i + 1);
    }
    if (modelInt == 0) {
      chi2_out = 1e10;
      chi2_error_out = 0.;
      return;
    }
    const float scale = (float)dataInt / modelInt;

    float sum_chi2{0.};
    float sum_error{0.};
    for (int i = lowchibin; i <= highchibin; ++i) {
      const float data = fDataHisto.GetBinContent(i);
      if (data < 1.0)
        continue;
      const float data_error = fDataHisto.GetBinError(i);
      const float model = counts.at(i) * scale;
      const float model_error = sqrt(counts.at(i)) * scale;
      const float error2 = pow(data_error, 2) + pow(model_error, 2);
      const float diff = model - data;
      sum_chi2 += pow(diff, 2) / error2;
      sum_error += pow(diff * (model_error - data_error) / error2, 2);
    }
    chi2_out = sum_chi2 / (highchibin - lowchibin + 1);
    chi2_error_out = 2 * pow(sum_error, 0.5) / (highchibin - lowchibin + 1);
  };

  const unsigned int n_workers = NumThreads();

  /*
   * Build multiplicity distributions for all evaluations in one pass over
   * Glauber events: for each event i, loop over all requested (f, mu, k, p).
   * Work units are (evaluation, chunk of events); events are split into
   * chunks only if there are fewer evaluations than threads.
   */
  auto BuildAll = [&](const std::vector<Evaluation> &evals,
                      std::vector<std::vector<float>> &counts,
                      const std::string &label) {
    const int n_evals = evals.size();
    const int n_chunks = std::max(1, ((int)n_workers + n_evals - 1) / n_evals);
    const int n_units = n_evals * n_chunks;

    struct Unit {
      int e;
      int i_start, i_stop;     // main events
      int plp_start, plp_stop; // pile-up pool
    };
    std::vector<Unit> units;
    for (int c = 0; c < n_chunks; c++)
      for (int e = 0; e < n_evals; e++) {
        const int nentries = (int)(nEvents * (1. - grid[evals[e].g].p));
        const int nplp = nEvents - nentries;
        units.push_back(
            Unit{e, (int)((long long)c * nentries / n_chunks),
                 (int)((long long)(c + 1) * nentries / n_chunks),
                 nentries + (int)((long long)c * nplp / n_chunks),
                 nentries + (int)((long long)(c + 1) * nplp / n_chunks)});
      }

    /* Worker w takes units w, w + n_workers, ... (units are ordered by
     * chunk); units of the same chunk share one loop over the events */
    struct Group {
      int lo, hi;             // events looped over
      std::vector<int> units; // units in this loop
    };
    std::vector<std::vector<Group>> plan(n_workers);
    long total = 0;
    for (int u = 0; u < n_units; u++) {
      auto &groups = plan[u % n_workers];
      if (groups.empty() || groups.back().units.back() / n_evals != u / n_evals)
        groups.push_back(Group{nEvents, 0, {}});
      auto &group = groups.back();
      group.lo = std::min(group.lo, units[u].i_start);
      group.hi = std::max(group.hi, units[u].i_stop);
      group.units.push_back(u);
    }
    for (const auto &groups : plan)
      for (const auto &group : groups)
        total += std::max(0, group.hi - group.lo);

    std::vector<std::vector<float>> unit_counts(
        n_units, std::vector<float>(nBinsModel, 0.));
    Progress progress(label, total);

    RunWorkers(n_workers, progress, [&](unsigned int i_worker, unsigned int) {
      std::random_device rd;
      std::mt19937 rngnum(rd());
      std::uniform_real_distribution<float> unirnd(0., 1.);

      struct Job {
        const Unit *unit;
        float f, p;
        std::gamma_distribution<> gammadist;
        int plp_counter;
        float *counts;
      };

      for (const auto &group : plan[i_worker]) {
        std::vector<Job> jobs;
        for (int u : group.units) {
          const Evaluation &e = evals[units[u].e];
          const GridPoint &gp = grid[e.g];
          jobs.push_back(Job{&units[u], gp.f, gp.p, MakeGamma(e.mu, gp.k),
                             units[u].plp_start, unit_counts[u].data()});
        }

        for (int i = group.lo; i < group.hi; i++) {
          if ((i - group.lo + 1) % Progress::kStep == 0)
            progress.Add(Progress::kStep);
          const float npart = fvNpart[i];
          const float ncoll = fvNcoll[i];
          for (auto &job : jobs) {
            if (i < job.unit->i_start || i >= job.unit->i_stop)
              continue;
            const int Na = int(::Nancestors(mode, job.f, npart, ncoll));
            float nHits = SumOfGammas(Na, job.gammadist, rngnum);
            if (job.p > 1e-10 && unirnd(rngnum) <= job.p) {
              const int j = NextPileUp(job.plp_counter, job.unit->plp_start,
                                       job.unit->plp_stop);
              if (j >= 0) {
                const int Na1 =
                    int(::Nancestors(mode, job.f, fvNpart[j], fvNcoll[j]));
                nHits += SumOfGammas(Na1, job.gammadist, rngnum);
              }
            }
            job.counts[FindBin(nHits)] += 1.;
          }
        }
        progress.Add(std::max(0, group.hi - group.lo) % Progress::kStep);
      }
    });

    counts.assign(n_evals, std::vector<float>(nBinsModel, 0.));
    for (int u = 0; u < n_units; u++)
      for (int b = 0; b < nBinsModel; b++)
        counts[units[u].e][b] += unit_counts[u][b];
  };

  std::vector<Evaluation> evals;
  std::vector<std::vector<float>> counts;

  auto UpdateChi2 = [&]() {
    for (size_t e = 0; e < evals.size(); e++) {
      auto &gp = grid[evals[e].g];
      if (evals[e].isMu2)
        Chi2FromCounts(counts[e], gp.chi2_mu2, gp.chi2_mu2_error);
      else
        Chi2FromCounts(counts[e], gp.chi2_mu1, gp.chi2_mu1_error);
    }
  };

  /* Initial golden section points (mu_1 and mu_2) for all grid points */
  for (int g = 0; g < (int)grid.size(); g++) {
    auto &gp = grid[g];
    if (!gp.valid)
      continue;
    gp.mu_1 = gp.mu_max - (gp.mu_max - gp.mu_min) / phi;
    gp.mu_2 = gp.mu_min + (gp.mu_max - gp.mu_min) / phi;
    evals.push_back(Evaluation{g, gp.mu_1, false});
    evals.push_back(Evaluation{g, gp.mu_2, true});
  }
  std::cout << "FitGlauber: " << grid.size() << " (f, k, p) points"
            << std::endl;

  if (!evals.empty()) {
    BuildAll(evals, counts, "FitGlauber: initialization");
    UpdateChi2();

    /* Golden section iterations, all grid points at once */
    for (int j = 0; j < fNiter; j++) {
      evals.clear();
      for (int g = 0; g < (int)grid.size(); g++) {
        auto &gp = grid[g];
        if (!gp.valid)
          continue;
        if (gp.chi2_mu1 >= gp.chi2_mu2) {
          gp.mu_min = gp.mu_1;
          gp.mu_1 = gp.mu_2;
          gp.mu_2 = gp.mu_min + (gp.mu_max - gp.mu_min) / phi;
          gp.chi2_mu1 = gp.chi2_mu2;
          gp.chi2_mu1_error = gp.chi2_mu2_error;
          evals.push_back(Evaluation{g, gp.mu_2, true});
        } else {
          gp.mu_max = gp.mu_2;
          gp.mu_2 = gp.mu_1;
          gp.mu_1 = gp.mu_max - (gp.mu_max - gp.mu_min) / phi;
          gp.chi2_mu2 = gp.chi2_mu1;
          gp.chi2_mu2_error = gp.chi2_mu1_error;
          evals.push_back(Evaluation{g, gp.mu_1, false});
        }
      }
      BuildAll(evals, counts,
               Form("FitGlauber: iteration [%d/%d]", j + 1, fNiter));
      UpdateChi2();

      for (int g = 0; g < (int)grid.size(); g++) {
        const auto &gp = grid[g];
        if (!gp.valid)
          continue;
        std::cout << "n = " << 1 + g * fNiter + j << " f = " << gp.f
                  << " k = " << gp.k << " p = " << gp.p << " mu1 = " << gp.mu_1
                  << " mu2 = " << gp.mu_2 << " chi2_mu1 = " << gp.chi2_mu1
                  << " chi2_mu2 = " << gp.chi2_mu2 << std::endl;
      }
    }
  }

  for (const auto &gp : grid) {
    f = gp.f;
    k = gp.k;
    p = gp.p;
    /* take min(mu), min(chi2), min(chi2_error) */
    const bool isFirst = gp.chi2_mu1 < gp.chi2_mu2;
    mu = isFirst ? gp.mu_1 : gp.mu_2;
    chi2 = isFirst ? gp.chi2_mu1 : gp.chi2_mu2;
    chi2_error = isFirst ? gp.chi2_mu1_error : gp.chi2_mu2_error;
    sigma = (mu / k + 1) * mu;

    tree->Fill();

    if (gp.valid && (!isFitted || chi2 < Chi2Min)) {
      isFitted = true;
      f_fit = f;
      mu_fit = mu;
      k_fit = k;
      p_fit = p;
      Chi2Min = chi2;
      Chi2Min_error = chi2_error;
    }
  }

  tree->Write();
  file->Write();
  file->Close();

  if (!isFitted) {
    std::cout << "FitGlauber: *** Error - no valid (f, k, p) points to fit"
              << std::endl;
    return Chi2Min;
  }

  /* Build full set of histograms for the best fit parameters */
  SetGlauberFitHisto(f_fit, mu_fit, k_fit, p_fit, nEvents);
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

  SetNBDhist(mu_fit, k_fit);

  fOptimalF = f_fit;
  fOptimalK = k_fit;
  fOptimalMu = mu_fit;
  fOptimalP = p_fit;
  fOptimalChi2Ndf = Chi2Min;
  fOptimalChi2NdfError = Chi2Min_error;

  return Chi2Min;
}

/**
 * Populates histogram nbd_<mean>_<k> with values of NBD
 * @param mu
 * @param k
 */
void Glauber::Fitter::SetNBDhist(float mu, float k) {
  // Interface for TH1F.
  const int nBins = (mu + 1.) * 3 < 10 ? 10 : (mu + 1.) * 3;

  fNbdHisto = TH1F("fNbdHisto", "", nBins, 0, nBins);
  if (fUseNbd)
    fNbdHisto.SetName("nbd");
  else
    fNbdHisto.SetName("gamma");

  std::random_device rd;
  std::mt19937 rngnum(rd());
  auto gammadist = MakeGamma(mu, k);

  for (int i = 0; i < 1e5; ++i) {
    fNbdHisto.Fill(gammadist(rngnum));
  }
}

/**
 * Creates histo with a given model parameter distribution
 * @param range observable range
 * @param name name of the MC-Glauber model parameter
 * @param nEvents
 * @return pointer to the histogram
 */
std::unique_ptr<TH1F> Glauber::Fitter::GetModelHisto(const float range[2],
                                                     TString name,
                                                     int nEvents) {
  nEvents = std::min(nEvents, (int)fvNpart.size());

  const float f = fOptimalF;
  const float mu = fOptimalMu;
  const float k = fOptimalK;
  const float p = fOptimalP;
  const NaMode mode = ParseMode(fMode);

  std::vector<float> modelInput;
  modelInput.reserve(nEvents);
  float modelpar{-999.};
  fSimTree->SetBranchAddress(name, &modelpar);
  for (int i = 0; i < nEvents; ++i) {
    fSimTree->GetEntry(i);
    modelInput.push_back(modelpar);
  }

  std::unique_ptr<TH1F> hModel(new TH1F("hModel", "name", 100,
                                        fSimTree->GetMinimum(name),
                                        fSimTree->GetMaximum(name)));

  const unsigned int n_workers = NumThreads();
  std::vector<std::vector<float>> model(n_workers);
  Progress progress("Glauber::Fitter::GetModelHisto",
                    (long)(nEvents * (1. - p)));

  RunWorkers(n_workers, progress,
             [&](unsigned int i_worker, unsigned int n_workers_) {
               /* events [i_start, plp_stop) of the worker: main events
                * first, then the pile-up pool */
               const int i_start =
                   (int)((long long)i_worker * nEvents / n_workers_);
               const int plp_stop =
                   (int)((long long)(i_worker + 1) * nEvents / n_workers_);
               const int i_stop =
                   i_start + (int)((plp_stop - i_start) * (1. - p));
               int plp_counter = i_stop;

               std::random_device rd;
               std::mt19937 rngnum(rd());
               std::uniform_real_distribution<float> unirnd(0., 1.);
               auto gammadist = MakeGamma(mu, k);

               for (int i = i_start; i < i_stop; i++) {
                 if ((i - i_start + 1) % Progress::kStep == 0)
                   progress.Add(Progress::kStep);
                 const int Na =
                     int(::Nancestors(mode, f, fvNpart[i], fvNcoll[i]));

                 float nHits = SumOfGammas(Na, gammadist, rngnum);
                 if (p > 1e-10 && unirnd(rngnum) <= p) {
                   const int j = NextPileUp(plp_counter, i_stop, plp_stop);
                   if (j >= 0) {
                     const int Na1 =
                         int(::Nancestors(mode, f, fvNpart[j], fvNcoll[j]));
                     nHits += SumOfGammas(Na1, gammadist, rngnum);
                   }
                 }
                 if (nHits > range[0] && nHits < range[1])
                   model[i_worker].push_back(modelInput[i]);
               }
               progress.Add((i_stop - i_start) % Progress::kStep);
             });

  for (const auto &out : model)
    for (float value : out)
      hModel->Fill(value);

  return hModel;
}

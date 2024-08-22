void config()
{
  // Set up input file from MC-Glauber
  std::string inFileGlauberName = "~/Documents/Work/Dataset/Centrality/Bmn/out_xe124csi_277mb_6M.root";
  std::string inTreeGlauberName = "nt_CsI_Xe124";

  // Set up input file with data multiplicity
  std::string inFileDataName = "~/Documents/Work/Dataset/Centrality/Bmn/run8_mult_negch_RunId_8000_8200_500k.root";
  std::string inHistDataName = "hMultRun8_VtxZ020";

  // Set up output directory
  std::string outDir = ".";

  // Set up number of iteration to find optimal mu parameter
  int Niter = 10;

  // Set up parameters for multiplicity fit
  float f_min  = 0.84;
  float f_max  = 0.84;
  float f_step = 0.01;
  float k_min  = 0.1;
  float k_max  = 0.1;
  float k_step = 0.1;
  float p_min  = 0.041;
  float p_max  = 0.041;
  float p_step = 0.001;

  // Set up fit ranges
  int multMin = 30;
  int multMax = 130;

  // Set up bin size in the data histogram
  int bin_size = 1;

  // Set up mode for number of ancestors calculation
  ///****************************************
  ///  |   mode    |   function for Na      |
  ///****************************************
  ///  |  Default  | f*Npart + (1-f)*Ncoll  |
  ///  |    PSD    |       f-Npart          |
  ///  |   Npart   |       Npart^f          |
  ///  |   Ncoll   |       Ncoll^f          |
  std::string mode = "Default";

  // Set up number of threads
  // sets maximum concurent processes by default
  unsigned int n_thr = std::thread::hardware_concurrency();
  
  std::unique_ptr<TFile> glauber_file{ TFile::Open(inFileGlauberName.c_str(), "read") };
  std::unique_ptr<TTree> glauber_tree{ (TTree*) glauber_file->Get(inTreeGlauberName.c_str()) };

  std::unique_ptr<TFile> data_file{ TFile::Open(inFileDataName.c_str(), "read")};
  std::unique_ptr<TH1F>  data_hist{ (TH1F*) data_file->Get(inHistDataName.c_str())};

  const int nevents = 10*(int(data_hist->Integral(multMin,multMax)));

  Glauber::Fitter fitter ( std::move(glauber_tree));
  fitter.SetMode(mode.c_str());
  fitter.SetMassNumber(f_min/2);
  fitter.SetInputHisto(*data_hist);
  fitter.SetBinSize(bin_size);
  fitter.Init(nevents, mode.c_str());
  
  fitter.SetFitMinBin(multMin);
  fitter.SetFitMaxBin(multMax);
  fitter.SetOutDirName(outDir.c_str());
  fitter.SetNiter(Niter);
  fitter.SetFstepSize(f_step);
  fitter.SetKstepSize(k_step);
  fitter.SetPstepSize(p_step);

  fitter.SetNthreads(n_thr);
  fitter.UseGamma();

  float par[5];
  float chi2=1e10;
  chi2 = fitter.FitGlauber(par, f_min, f_max, k_min, k_max, p_min, p_max, nevents);

  std::cout << "f = " << par[0] << "    mu = " << par[1] << "    k = " << par[2] << "    p = " << par[4] << "    chi2 = " << chi2 << "    chi2_error = " << par[3] << std::endl; 
  
  DrawHistos(fitter, par, chi2, true, true, true, true);

  const float range[2] = {(float)multMin, (float)multMax};
  std::unique_ptr<TH1F> hB(fitter.GetModelHisto (range, "B", par, 100000));
  hB->SaveAs( "b_test.root" );
}
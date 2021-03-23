void plot_b_Cent(TString inFile_hist, TString inFile_FINAL, TString outFileName="test_bCent.root")
{
  TFile *fi_hist = new TFile(inFile_hist.Data(),"read");
  TFile *fi_glaub = new TFile(inFile_FINAL.Data(),"read");

  const int Ncent = 10;
  const TString centClass[Ncent] = {"0.0%-10.0%", "10.0%-20.0%", "20.0%-30.0%", "30.0%-40.0%", "40.0%-50.0%", "50.0%-60.0%", "60.0%-70.0%", "70.0%-80.0%", "80.0%-90.0%", "90.0%-100.0%"};

  TH1D *hBCent[Ncent];
  TTree *Result = (TTree*) fi_glaub->Get("Result");
  int minNch, maxNch;
  Result->SetBranchAddress("MinBorder", &minNch);
  Result->SetBranchAddress("MaxBorder", &maxNch);
  for (int i=0; i<Ncent; i++)
  {
    hBCent[i] = (TH1D*) fi_glaub->Get(Form("B_VS_CentralityClass %s", centClass[i].Data()));
  }
  TH2F *hBvsCent = (TH2F*) fi_hist->Get("hBvsRefMult");
  //TH2F *hBvsCent = (TH2F*) fi_hist->Get("hNpartImpact0");

  TFile *fo = new TFile(outFileName.Data(),"recreate");

  TH1F *pBCent_glaub = new TH1F("pBCent_glaub","<b> vs Centrality from Glauber calculations;Centrality, %;<b>#pm#sigma_{b}",Ncent, 0., 100.);
  TH1F *pBCent_model = new TH1F("pBCent_model","<b> vs Centrality from model;Centrality, %;<b>#pm#sigma_{b}",Ncent, 0., 100.);
  std::vector<std::pair<int, int>> vBorders;
  for (int i=0; i<Result->GetEntriesFast(); i++)
  {
    Result->GetEntry(i); 
    vBorders.push_back({minNch, maxNch});
    std::cout << "Nch_min = " << minNch << " Nch_max = " << maxNch << std::endl;
  }
  
  TH1D *hBCent_model[Ncent];
  float meanB, sigmB, ratioB;
  std::cout << "Centrality | <b>_glauber +- sigma_b | <b>_model +- sigma_b | glauber/model ratio |" << std::endl;
  for (int i=0; i<Ncent; i++)
  {
    std::cout << centClass[i].Data() << " | ";
    // Fill results from Glauber
    meanB = hBCent[i]->GetMean();
    sigmB = hBCent[i]->GetRMS();
    pBCent_glaub->SetBinContent(i+1, meanB);
    pBCent_glaub->SetBinError(i+1, sigmB);
    std::cout << meanB << " +- " << sigmB;
    ratioB = meanB;

    std::cout << " | ";
    // Fill results from the model
    hBCent_model[i] = (TH1D*) hBvsCent->ProjectionY(Form("B_VS_CentralityClass %s",centClass[i].Data()),hBvsCent->GetXaxis()->FindBin(vBorders.at(i).first),hBvsCent->GetXaxis()->FindBin(vBorders.at(i).second));
    meanB = hBCent_model[i]->GetMean();
    sigmB = hBCent_model[i]->GetRMS();
    pBCent_model->SetBinContent(i+1, meanB);
    pBCent_model->SetBinError(i+1, sigmB);
    std::cout << meanB << " +- " << sigmB;

    ratioB /= meanB;
    std::cout << " | ";
    std::cout << ratioB;
    std::cout << std::endl;

    //delete tmp;
  }

  fo->cd();
  pBCent_glaub->Write();
  pBCent_model->Write();
  fo->mkdir("Glauber");
  fo->cd("Glauber");
  for (int i=0; i<Ncent; i++)
  {
    hBCent[i]->Write();
  }
  fo->mkdir("Model");
  fo->cd("Model");
  for (int i=0; i<Ncent; i++)
  {
    hBCent_model[i]->Write();
  }
  fo->Close();
}

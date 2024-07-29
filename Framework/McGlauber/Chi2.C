void Chi2(TString InFileName)
{
  TFile *file = new TFile(InFileName);
  TTree *tree = (TTree *)file->Get("test_tree");
  TTree *delta = new TTree("delta", "delta");

  Float_t f, mu, k, p;
  Float_t chi2, chi2_error, chi2_min_error, CHI2, DELTA = 1e10;
  Float_t chi2_min = 1e10;
  Float_t f_min, mu_min, k_min, p_min;
  Float_t f_max = 0, mu_max = 0, k_max = 0, p_max = 0;
  Float_t f_delta = 0, mu_delta = 0, k_delta = 0, p_delta = 0;
  Float_t sigma;

  tree->SetBranchAddress("f", &f);
  tree->SetBranchAddress("mu", &mu);
  tree->SetBranchAddress("k", &k);
  tree->SetBranchAddress("p", &p);
  tree->SetBranchAddress("chi2", &chi2);
  tree->SetBranchAddress("chi2_error", &chi2_error);
  tree->SetBranchAddress("sigma", &sigma);

  delta->Branch("chi2", &chi2_min);
  delta->Branch("delta_chi2", &chi2_min_error);
  delta->Branch("mu", &mu_min);
  delta->Branch("delta_mu", &mu_delta);
  delta->Branch("k", &k_min);
  delta->Branch("delta_k", &k_delta);
  delta->Branch("f", &f_min);
  delta->Branch("delta_f", &f_delta);
  delta->Branch("p", &p_min);
  delta->Branch("delta_p", &p_delta);

  Int_t n = tree->GetEntries();

  std::cout << " GetEntries = " << n << std::endl;

  std::vector<Float_t> vf;
  std::vector<Float_t> vk;
  std::vector<Float_t> vp;
  std::vector<Float_t> vchi2;

  float X, Y, Z;
  for (Int_t i = 0; i < n; i++)
  {
    tree->GetEntry(i);

    std::cout << " f = " << f << " mu = " << mu << " k = " << k << " p = " << p << " sigma = " << sigma << " chi2 = " << chi2 << " chi2_error = " << chi2_error << std::endl;

    vf.push_back(f);
    vk.push_back(k);
    vp.push_back(p);
    vchi2.push_back(chi2);

    if (chi2 < chi2_min)
    {
      chi2_min = chi2;
      chi2_min_error = chi2_error;
      CHI2 = chi2_min + chi2_min_error;
      f_min = f;
      mu_min = mu;
      k_min = k;
      p_min = p;
      X = f_min;
      Y = k_min;
      Z = p_min;
    }
  }

  // TCanvas *c1 = new TCanvas("c1", "c1", 800, 800);

  // TH3F *h = new TH3F("#chi^{2} vs f, k, p", "#chi^{2} vs f, k, p; f; k; p; #chi^{2}",
  //                    vf.size(), *std::min_element(vf.begin(), vf.end()), *std::max_element(vf.begin(), vf.end()),
  //                    vk.size(), *std::min_element(vk.begin(), vk.end()), *std::max_element(vk.begin(), vk.end()),
  //                    vp.size(), *std::min_element(vp.begin(), vp.end()), *std::max_element(vp.begin(), vp.end()));
  // for (int i = 0; i < vchi2.size(); ++i)
  // {
  //   float val_f = vf.at(i);
  //   float val_k = vk.at(i);
  //   float val_p = vp.at(i);
  //   float val_chi2 = vchi2.at(i);
  //   int bin = h->FindBin(val_f, val_k, val_p);
  //   h->SetBinContent(bin, val_chi2);
  // }
  // h->Draw("colz");

  // float F, K, P;
  // for (double xx = X; (xx <= X + 0.2 && xx <= *std::max_element(vf.begin(), vf.end())); xx += 0.001)
  // {
  //   for (double yy = Y; (yy <= Y + 10.0 && yy <= *std::max_element(vk.begin(), vk.end())); yy += 0.001)
  //   {
  //     for (double zz = Z; (zz <= Z + 0.2 && zz <= *std::max_element(vp.begin(), vp.end())); zz += 0.001)
  //       if (TMath::Abs(h->Interpolate(xx, yy, zz) - CHI2) < DELTA)
  //       {
  //         if (TMath::Abs(xx - f_min) > f_delta)
  //           f_delta = TMath::Abs(xx - f_min);
  //         if (TMath::Abs(yy - k_min) > k_delta)
  //           k_delta = TMath::Abs(yy - k_min);
  //         if (TMath::Abs(zz - p_min) > p_delta)
  //           p_delta = TMath::Abs(zz - p_min);
  //         F = 1e10;
  //         K = 1e10;
  //         P = 1e10;
  //         for (Int_t ii = 0; ii < n; ii++)
  //         {
  //           tree->GetEntry(ii);
  //           if (TMath::Abs(f - xx) <= F && TMath::Abs(k - yy) <= K && TMath::Abs(p - zz) <= P && (TMath::Abs(mu - mu_min)) > mu_delta)
  //             mu_delta = TMath::Abs(mu - mu_min);
  //         }
  //         DELTA = TMath::Abs(h->Interpolate(xx, yy, zz) - CHI2);
  //       }
  //     for (double zz = Z; (zz >= Z - 0.2 && zz >= *std::min_element(vp.begin(), vp.end())); zz -= 0.001)
  //       if (TMath::Abs(h->Interpolate(xx, yy, zz) - CHI2) < DELTA)
  //       {
  //         if (TMath::Abs(xx - f_min) > f_delta)
  //           f_delta = TMath::Abs(xx - f_min);
  //         if (TMath::Abs(yy - k_min) > k_delta)
  //           k_delta = TMath::Abs(yy - k_min);
  //         if (TMath::Abs(zz - p_min) > p_delta)
  //           p_delta = TMath::Abs(zz - p_min);
  //         F = 1e10;
  //         K = 1e10;
  //         P = 1e10;
  //         for (Int_t ii = 0; ii < n; ii++)
  //         {
  //           tree->GetEntry(ii);
  //           if (TMath::Abs(f - xx) <= F && TMath::Abs(k - yy) <= K && TMath::Abs(p - zz) <= P && (TMath::Abs(mu - mu_min)) > mu_delta)
  //             mu_delta = TMath::Abs(mu - mu_min);
  //         }
  //         DELTA = TMath::Abs(h->Interpolate(xx, yy, zz) - CHI2);
  //       }
  //   }
  //   for (double yy = Y; (yy >= Y - 10.0 && yy >= *std::min_element(vk.begin(), vk.end())); yy -= 0.001)
  //   {
  //     for (double zz = Z; (zz <= Z + 0.2 && zz <= *std::max_element(vp.begin(), vp.end())); zz += 0.001)
  //       if (TMath::Abs(h->Interpolate(xx, yy, zz) - CHI2) < DELTA)
  //       {
  //         if (TMath::Abs(xx - f_min) > f_delta)
  //           f_delta = TMath::Abs(xx - f_min);
  //         if (TMath::Abs(yy - k_min) > k_delta)
  //           k_delta = TMath::Abs(yy - k_min);
  //         if (TMath::Abs(zz - p_min) > p_delta)
  //           p_delta = TMath::Abs(zz - p_min);
  //         F = 1e10;
  //         K = 1e10;
  //         P = 1e10;
  //         for (Int_t ii = 0; ii < n; ii++)
  //         {
  //           tree->GetEntry(ii);
  //           if (TMath::Abs(f - xx) <= F && TMath::Abs(k - yy) <= K && TMath::Abs(p - zz) <= P && (TMath::Abs(mu - mu_min)) > mu_delta)
  //             mu_delta = TMath::Abs(mu - mu_min);
  //         }
  //         DELTA = TMath::Abs(h->Interpolate(xx, yy, zz) - CHI2);
  //       }
  //     for (double zz = Z; (zz >= Z - 0.2 && zz >= *std::min_element(vp.begin(), vp.end())); zz -= 0.001)
  //       if (TMath::Abs(h->Interpolate(xx, yy, zz) - CHI2) < DELTA)
  //       {
  //         if (TMath::Abs(xx - f_min) > f_delta)
  //           f_delta = TMath::Abs(xx - f_min);
  //         if (TMath::Abs(yy - k_min) > k_delta)
  //           k_delta = TMath::Abs(yy - k_min);
  //         if (TMath::Abs(zz - p_min) > p_delta)
  //           p_delta = TMath::Abs(zz - p_min);
  //         F = 1e10;
  //         K = 1e10;
  //         P = 1e10;
  //         for (Int_t ii = 0; ii < n; ii++)
  //         {
  //           tree->GetEntry(ii);
  //           if (TMath::Abs(f - xx) <= F && TMath::Abs(k - yy) <= K && TMath::Abs(p - zz) <= P && (TMath::Abs(mu - mu_min)) > mu_delta)
  //             mu_delta = TMath::Abs(mu - mu_min);
  //         }
  //         DELTA = TMath::Abs(h->Interpolate(xx, yy, zz) - CHI2);
  //       }
  //   }
  // }

  // for (double xx = X; (xx >= X - 0.2 && xx >= g->GetXmin()); xx -= 0.001)
  // {
  //   for (double yy = Y; (yy <= Y + 2.0 || yy <= g->GetYmax()); yy += 0.001)
  //     if (TMath::Abs(g->Interpolate(xx, yy) - CHI2) < DELTA)
  //     {
  //       if (TMath::Abs(xx - f_min) > f_delta)
  //         f_delta = TMath::Abs(xx - f_min);
  //       if (TMath::Abs(yy - k_min) > k_delta)
  //         k_delta = TMath::Abs(yy - k_min);
  //       F = 1e10;
  //       K = 1e10;
  //       for (Int_t ii = 0; ii < n; ii++)
  //       {
  //         tree->GetEntry(ii);
  //         if (TMath::Abs(f - xx) <= F && TMath::Abs(k - yy) <= K && (TMath::Abs(mu - mu_min)) > mu_delta)
  //           mu_delta = TMath::Abs(mu - mu_min);
  //       }
  //       DELTA = TMath::Abs(g->Interpolate(xx, yy) - CHI2);
  //     }
  //   for (double yy = Y; (yy >= Y - 2.0 && yy >= g->GetYmin()); yy -= 0.001)
  //     if (TMath::Abs(g->Interpolate(xx, yy) - CHI2) < DELTA)
  //     {
  //       if (TMath::Abs(xx - f_min) > f_delta)
  //         f_delta = TMath::Abs(xx - f_min);
  //       if (TMath::Abs(yy - k_min) > k_delta)
  //         k_delta = TMath::Abs(yy - k_min);
  //       F = 1e10;
  //       K = 1e10;
  //       for (Int_t ii = 0; ii < n; ii++)
  //       {
  //         tree->GetEntry(ii);
  //         if (TMath::Abs(f - xx) <= F && TMath::Abs(k - yy) <= K && (TMath::Abs(mu - mu_min)) > mu_delta)
  //           mu_delta = TMath::Abs(mu - mu_min);
  //       }
  //       DELTA = TMath::Abs(g->Interpolate(xx, yy) - CHI2);
  //     }
  // }

  // g->Draw("colz");

  // std::cout << " f = " << f_min << "+/-" << f_delta << " mu = " << mu_min << "+/-" << mu_delta << " k = " << k_min << "+/-" << k_delta << " chi2 = " << chi2_min << "+/-" << chi2_min_error << std::endl;
  std::cout << " f = " << f_min << " mu = " << mu_min << " k = " << k_min << " p = " << p_min << " sigma = " << sigma << " chi2 = " << chi2 << " chi2_error = " << chi2_error << std::endl;

  TFile *f1 = new TFile("Fit_Errors_RPC.root", "recreate");
  delta->Write();
  // h->Write();
  f1->Close();
}

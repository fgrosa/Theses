void PlotFitMass(Double_t nSigma=3, Double_t pTmin=1, Double_t pTmax=24, TString massfilename="MassFit_pt_1-24", TString masscanvasname="cMassFit") {
  
  TFile inmassfile(Form("%s.root",massfilename.Data()),"UPDATE");
  TCanvas* c=(TCanvas*)inmassfile.Get(masscanvasname.Data());
  TH1F* h = (TH1F*)c->GetPrimitive("fhistoInvMass");
  TF1* fs = (TF1*)c->GetPrimitive("funcmass");
  TF1* fb = (TF1*)c->GetPrimitive("funcbkgFullRange");
  TPaveText* info = new TPaveText(0.6,0.5,0.9,0.7,"NDC");
  TPaveText* info2 = new TPaveText(0.1,0.5,0.5,0.8,"NDC");
  h->SetDirectory(0);
  inmassfile.Close();
  
  Double_t mean = fs->GetParameter(3);
  Double_t sigma = fs->GetParameter(4);
  Double_t errmean = fs->GetParError(3);
  Double_t errsigma = fs->GetParError(4);
  Double_t ints=fs->Integral(mean-nSigma*sigma,mean+nSigma*sigma)/h->GetBinWidth(4);
  Double_t intb=fb->Integral(mean-nSigma*sigma,mean+nSigma*sigma)/h->GetBinWidth(2);
  Double_t signal = ints-intb;
  Double_t signalerr = fs->GetParError(fs->GetNpar()-3)/fs->GetParameter(fs->GetNpar()-3)*signal;
  Double_t bkg = intb;
  Double_t bkgerr = fb->GetParError(0)/fb->GetParameter(0)*bkg;
  Double_t significance = signal/TMath::Sqrt(signal+bkg);
  Double_t significanceerr = significance*TMath::Sqrt((signalerr*signalerr+bkgerr*bkgerr)/(4.*(signal+bkg)*(signal+bkg))+(bkg/(signal+bkg))*(signalerr*signalerr)/signal/signal);
  Double_t signaloverbkg = signal/bkg;

  Double_t binwidth = h->GetBinWidth(10)*1000;//in MeV
  
  info->Clear();
  info->SetTextSize(0.04);
  info->SetBorderSize(0);
  info->SetFillStyle(0);
  info->SetTextColor(kBlue);
  info->SetTextFont(132);
  info->AddText(Form("mean = %.3f #pm %.3f",mean,errmean));
  info->AddText(Form("sigma = %.3f #pm %.3f",sigma,errsigma));
  info2->Clear();    
  info2->SetTextSize(0.04);
  info2->SetBorderSize(0);
  info2->SetTextFont(132);
  info2->SetFillStyle(0);
  info2->AddText(Form("Significance (%0.f#sigma) = %.1f #pm %.1f",nSigma,significance,significanceerr));
  info2->AddText(Form("S (%0.f#sigma) = %.0f #pm %.0f",nSigma,signal,signalerr));
  info2->AddText(Form("B (%0.f#sigma) = %.0f #pm %.0f",nSigma,bkg,bkgerr));
  info2->AddText(Form("S/B (%0.f#sigma) = %.4f",nSigma,signaloverbkg));
  
  TCanvas* cMass = new TCanvas(masscanvasname.Data(),"",1200,900);
  h->Draw("E");
  h->SetStats(0);
  h->SetTitle(Form("%0.f < #it{p}_{T} < %0.f GeV/c", pTmin, pTmax));
  h->SetTitleSize(0.04);
  h->GetXaxis()->SetTitleSize(0.045);
  h->GetYaxis()->SetTitleSize(0.045);
  h->GetYaxis()->SetTitle(Form("Entries/(%0.f Mev/c^{2})",binwidth));
  h->GetYaxis()->SetRangeUser(h->GetMaximum()*0.3,h->GetMaximum()*1.1);
  fs->Draw("same");
  info->Draw("same");
  info2->Draw("same");

  cMass->SaveAs(Form("%s_fit.root",massfilename.Data()));
  cMass->SaveAs(Form("%s_fit.eps",massfilename.Data()));
}

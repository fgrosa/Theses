void results()
{
  TFile *f = new TFile("fgrosa_Photonic.root","open");
  TFile *f2 = new TFile("ULLS.root","open");
  TList *l = f->Get("fgrosa_Photonic");

  //All THnSparse (1st dim-> Pt inc, 2nd dim -> mee, 3rd bin-> Pt ass)
  THnSparseF *spLike = (THnSparseF*)l->FindObject("fHistLikeSign_MC");
  THnSparseF *spUnlike = (THnSparseF*)l->FindObject("fHistUnlikeSign_MC");
  THnSparseF *spLS = (THnSparseF*)l->FindObject("fHistLikeSign_esd");
  THnSparseF *spUS = (THnSparseF*)l->FindObject("fHistUnlikeSign_esd");
  THnSparseF *Pi0mee_MC = (THnSparseF*)l->FindObject("fHistelfromPi0Mass_MC");
  THnSparseF *Etamee_MC = (THnSparseF*)l->FindObject("fHistelfromEtaMass_MC");
  THnSparseF *gammamee_MC = (THnSparseF*)l->FindObject("fHistelfromgammaMass_MC");
  THnSparseF *mee_MC = (THnSparseF*)l->FindObject("fHistelMass_MC");
  THnSparseF *Pi0mee_esd = (THnSparseF*)l->FindObject("fHistelfromPi0Mass_esd");
  THnSparseF *Etamee_esd = (THnSparseF*)l->FindObject("fHistelfromEtaMass_esd");
  THnSparseF *gammamee_esd = (THnSparseF*)l->FindObject("fHistelfromgammaMass_esd");
  THnSparseF *mee_esd = (THnSparseF*)l->FindObject("fHistelMass_esd");
  THnSparseF *spelPip = (THnSparseF*)l->FindObject("fHistElPipPt");
  THnSparseF *spposPim = (THnSparseF*)l->FindObject("fHistPosPimPt");

  //Pt cuts
  /*
  Double_t ptcut_incl = 4;
  Double_t ptcut_ass = 3;

  spLS->GetAxis(0)->SetRange(ptcut_incl,31);
  spUS->GetAxis(0)->SetRange(ptcut_incl,31);
  spLike->GetAxis(0)->SetRange(ptcut_incl,31);
  spUnlike->GetAxis(0)->SetRange(ptcut_incl,31);
  
  Pi0mee_MC->GetAxis(0)->SetRange(ptcut_incl,31);
  Etamee_MC->GetAxis(0)->SetRange(ptcut_incl,31);
  gammamee_MC->GetAxis(0)->SetRange(ptcut_incl,31);
  mee_MC->GetAxis(0)->SetRange(ptcut_incl,31);
  Pi0mee_esd->GetAxis(0)->SetRange(ptcut_incl,31);
  Etamee_esd->GetAxis(0)->SetRange(ptcut_incl,31);
  gammamee_esd->GetAxis(0)->SetRange(ptcut_incl,31);
  mee_esd->GetAxis(0)->SetRange(ptcut_incl,31);

  spLS->GetAxis(2)->SetRange(ptcut_ass,31);
  spUS->GetAxis(2)->SetRange(ptcut_ass,31);
  spLike->GetAxis(2)->SetRange(ptcut_ass,31);
  spUnlike->GetAxis(2)->SetRange(ptcut_ass,31);
  
  Pi0mee_MC->GetAxis(2)->SetRange(ptcut_ass,31);
  Etamee_MC->GetAxis(2)->SetRange(ptcut_ass,31);
  gammamee_MC->GetAxis(2)->SetRange(ptcut_ass,31);
  mee_MC->GetAxis(2)->SetRange(ptcut_ass,31);
  Pi0mee_esd->GetAxis(2)->SetRange(ptcut_ass,31);
  Etamee_esd->GetAxis(2)->SetRange(ptcut_ass,31);
  gammamee_esd->GetAxis(2)->SetRange(ptcut_ass,31);
  mee_esd->GetAxis(2)->SetRange(ptcut_ass,31);  
  */

    //Main functions definitions
  TF1 *dalitz_Pi0 = new TF1("dalitz_Pi0", "(x<[1])*[0]*(2./x)*(((1.+(x/[1])^2)^2-4*(x/[1])^2)^(3./2.))*sqrt(1-(4*(0.000511/[1])^2/(x/[1])^2))*(1+(2*(0.000511/[1])^2/(x/[1])^2))+(x>=[1])*0",0.,2);//*(1/(1-5.5*x^2))^2
  dalitz_Pi0->SetLineWidth(2); 

  TF1 *dalitz_Eta = new TF1("dalitz_Eta", "(x<[1])*[0]*(2./x)*(((1.+(x/[1])^2)^2-4*(x/[1])^2)^(3./2.))*sqrt(1-(4*(0.000511/[1])^2/(x/[1])^2))*(1+(2*(0.000511/[1])^2/(x/[1])^2))+(x>=[1])*0",0.,2);//*(1/(1-1.9*x^2))^2
  dalitz_Eta->SetLineWidth(2); 

  TF1 *fGaussExp = new TF1("fGaussExp","(x>=[1])*([0]*(exp(-0.5*((x-[1])/[2])^2)+exp((x-[1])/[3])*(1.-exp(-0.5*((x-[1])/[2])^2))))+(x<[1])*([0]*exp(-0.5*((x-[1])/[2])^2))",0,1);
  fGaussExp->SetLineWidth(2);
  fGaussExp->SetLineColor(8);

  /*------------------------- PT distribution of sources from MC truth ------------------------------*/
  /*
  TCanvas *c1 = new TCanvas("c1","c1",10,10,1200,900);
  c1->Clear();
  c1->Divide(4,2);
  c1->cd(1);
  l->FindObject("fHistPi0Pt_true")->Draw();
  c1->cd(2);
  l->FindObject("fHistEtaPt_true")->Draw();
  c1->cd(3);
  l->FindObject("fHistgammaPt_true")->Draw();
  c1->cd(4);
  l->FindObject("fHistelPt_true")->Draw();
  c1->cd(5);
  l->FindObject("fHistPi0Pt")->Draw();
  c1->cd(6);
  l->FindObject("fHistEtaPt")->Draw();
  c1->cd(7);
  l->FindObject("fHistgammaPt")->Draw();
  c1->cd(8);
  l->FindObject("fHistelPt")->Draw();

  /*--------------------- PT distribution of electrons from MC truth ------------------------------*/
  /*
  TCanvas *c2 = new TCanvas("c2","c2",10,10,1200,900);
  c2->Clear();
  c2->Divide(2,2);
  c2->cd(1);
  l->FindObject("fHistelfromPi0Pt")->Draw();
  c2->cd(2);
  l->FindObject("fHistelfromEtaPt")->Draw();
  c2->cd(3);
  l->FindObject("fHistelfromgammaPt")->Draw();
  c2->cd(4);
  l->FindObject("fHistelPt")->Draw();

  /*---------------------- Invariant mass from MC truth with MC observables -----------------------------*/

  TCanvas *c3 = new TCanvas("c3","c3",10,10,1200,900);
  
  fGaussExp->SetParameters(150000,0.001,0.005,-0.05);

  Double_t n_pi0_MC;
  Double_t m_pi0_MC;
  Double_t n_eta_MC;
  Double_t m_eta_MC;
  Double_t par0_MC;
  Double_t par1_MC;
  Double_t par2_MC;
  Double_t par3_MC;

  c3->Clear();
  c3->Divide(3,1);
  c3->cd(1)->SetLogy();
  TH1F *hPi0mee_MC = (TH1F*)Pi0mee_MC->Projection(1);
  dalitz_Pi0->SetParameters(1,0.13957);
  //hPi0mee_MC->SetStats(0);
  hPi0mee_MC->SetMarkerSize(0.5);
  hPi0mee_MC->SetMarkerColor(1);
  hPi0mee_MC->SetLineColor(1);
  dalitz_Pi0->SetLineColor(1);
  hPi0mee_MC->Fit("dalitz_Pi0","L");
  hPi0mee_MC->Draw("E0E1");
  n_pi0_MC = dalitz_Pi0->GetParameter(0);
  m_pi0_MC = dalitz_Pi0->GetParameter(1);
  c3->cd(2)->SetLogy();
  TH1F *hEtamee_MC = (TH1F*)Etamee_MC->Projection(1);
  //hEtamee_MC->SetStats(0);
  //hEtamee_MC->GetYaxis()->SetRangeUser(0.1,3000);
  hEtamee_MC->SetMarkerSize(0.5);
  hEtamee_MC->SetMarkerColor(4);
  hEtamee_MC->SetLineColor(4);
  dalitz_Eta->SetParameters(1.,0.547);
  dalitz_Eta->SetLineColor(4);
  hEtamee_MC->Fit("dalitz_Eta","L");
  hEtamee_MC->Draw("E0E1");
  n_eta_MC = dalitz_Eta->GetParameter(0);
  m_eta_MC = dalitz_Eta->GetParameter(1);
  c3->cd(3)->SetLogy();
  TH1F *hgammamee_MC = (TH1F*)gammamee_MC->Projection(1);
  hgammamee_MC->GetXaxis()->SetRangeUser(0,0.2);
  //hgammamee_MC->SetStats(0);
  hgammamee_MC->SetMarkerSize(0.5);
  hgammamee_MC->SetMarkerColor(8);
  hgammamee_MC->SetLineColor(8);
  hgammamee_MC->Fit("fGaussExp");
  hgammamee_MC->Draw("E0E1");
  par0_MC = fGaussExp->GetParameter(0);
  par1_MC = fGaussExp->GetParameter(1);
  par2_MC = fGaussExp->GetParameter(2);
  par3_MC = fGaussExp->GetParameter(3);
  
  TF1 *gamma = new TF1("gamma","fGaussExp",0.,1.);
  gamma->SetParameters(par0_MC,par1_MC,par2_MC,par3_MC);
  gamma->SetLineColor(8);
  gamma->SetLineWidth(2);

  TF1 *Eta = new TF1("Eta","dalitz_Eta",0.,1.);
  Eta->SetParameters(n_eta_MC,m_eta_MC);
  Eta->SetLineColor(4);
  Eta->SetLineWidth(2);

  TF1 *Pi0 = new TF1("Pi0","dalitz_Pi0",0.,1.);
  Pi0->SetParameters(n_pi0_MC,m_pi0_MC);
  Pi0->SetLineColor(1);
  Pi0->SetLineWidth(2);

  TF1 *dist_MC = new TF1("dist_MC","Pi0+Eta+gamma",0.,1.);
  dist_MC->SetLineWidth(2);
  dist_MC->SetLineColor(2);
  dist_MC->SetParameters(n_pi0_MC,m_pi0_MC,n_eta_MC,m_eta_MC,par0_MC,par1_MC,par2_MC,par3_MC);
  dist_MC->FixParameter(1,m_pi0_MC);
  dist_MC->FixParameter(3,m_eta_MC);
  //dist_MC->FixParameter(4,par0_MC);
  dist_MC->FixParameter(5,par1_MC);
  //dist_MC->FixParameter(6,par2_MC);
  //dist_MC->FixParameter(7,par3_MC);

  cout << "\n\n***************** REAL RATIOS ******************\n\n" << endl;

  Double_t n_gamma = hgammamee_MC->GetEntries();
  Double_t n_Pi0 = hPi0mee_MC->GetEntries();
  Double_t n_Eta = hEtamee_MC->GetEntries();
  Double_t err_n_gamma =TMath::Sqrt(n_gamma);
  Double_t err_n_Pi0 =  TMath::Sqrt(n_Pi0);
  Double_t err_n_Eta = TMath::Sqrt(n_Eta);
;

 cout << "n of gamma = " << n_gamma << " +/- " << err_n_gamma << endl; 
 cout << "n of Pi0 = " << n_Pi0 << " +/- " << err_n_Pi0 <<endl;
 cout << "n of gamma = " << n_Eta << " +/- " << err_n_Eta << endl;
 cout << "Real Pi0/Eta Ratio = " << n_Pi0/n_Eta << " +/- " << n_Pi0/n_Eta*TMath::Sqrt(TMath::Power(err_n_Pi0/n_Pi0,2)+TMath::Power(err_n_Eta/n_Eta,2)) <<endl;
 cout << "Real gamma/Pi0 Ratio = "<< n_gamma/n_Pi0 << " +/- " << n_gamma/n_Pi0*TMath::Sqrt(TMath::Power(err_n_Pi0/n_Pi0,2)+TMath::Power(err_n_gamma/n_gamma,2))  <<endl;
  cout << "Real gamma/Eta Ratio = "<< n_gamma/n_Eta << " +/- " << n_gamma/n_Eta*TMath::Sqrt(TMath::Power(err_n_Eta/n_Eta,2)+TMath::Power(err_n_gamma/n_gamma,2)) <<endl;
  cout << "Real gamma/(Pi0+Eta) Ratio = " << n_gamma/(n_Eta+n_Pi0) << " +/- " << n_gamma/(n_Eta+n_Pi0)*TMath::Sqrt(TMath::Power(err_n_gamma/n_gamma,2)+TMath::Power(TMath::Sqrt(TMath::Power(err_n_Pi0,2)+TMath::Power(err_n_Eta,2))/(n_Eta+n_Pi0),2)) << "\n" <<endl;

												      	  
  cout << "\n\n***************** MC DIST ******************\n\n" << endl;

  TCanvas *MC_dist = new TCanvas("MC_dist","MC_dist",10,10,1200,900);
  MC_dist->Clear();
  MC_dist->SetLogy();
  
  TH1F *hmee_MC = (TH1F*)mee_MC->Projection(1);
  hmee_MC->GetXaxis()->SetTitle("m_{ee} (GeV/c^{2})");
  hmee_MC->GetYaxis()->SetTitle("dN/dm_{ee} (c^{2}/GeV)");

  hmee_MC->SetStats(0);
  hmee_MC->SetLineColor(2);
  hmee_MC->SetMarkerColor(2);
  hmee_MC->SetMarkerStyle(20);
  hmee_MC->SetMarkerSize(0.7);
  hmee_MC->SetTitle("  ");
  hmee_MC->GetXaxis()->CenterTitle();
  hmee_MC->GetYaxis()->CenterTitle();
  hmee_MC->Fit("dist_MC","L");
  hmee_MC->Draw("E0E1");

  Double_t n_pi0_MC2 = dist_MC->GetParameter(0);
  Double_t m_pi0_MC2 = dist_MC->GetParameter(1);
  Double_t n_eta_MC2 = dist_MC->GetParameter(2);
  Double_t m_eta_MC2 = dist_MC->GetParameter(3);
  Double_t par0_MC2 = dist_MC->GetParameter(4);
  Double_t par1_MC2 = dist_MC->GetParameter(5);
  Double_t par2_MC2 = dist_MC->GetParameter(6);
  Double_t par3_MC2 = dist_MC->GetParameter(7);

  gamma->SetParameters(par0_MC2,par1_MC2,par2_MC2,par3_MC2);
  Eta->SetParameters(n_eta_MC2,m_eta_MC2);
  Pi0->SetParameters(n_pi0_MC2,m_pi0_MC2);

  Pi0->SetLineStyle(10);
  Eta->SetLineStyle(3);
  gamma->SetLineStyle(9);
  Pi0->Draw("same");
  Eta->Draw("same");
  gamma->Draw("same");

  cout<<"NDF = " << dist_MC->GetNDF()<<endl;   
 
  Double_t gamma_area = gamma->Integral(0,1)*hmee_MC->GetNbinsX()/0.6;
  Double_t Pi0_area = Pi0->Integral(0,1)*hmee_MC->GetNbinsX()/0.6;
  Double_t Eta_area = Eta->Integral(0,1)*hmee_MC->GetNbinsX()/0.6;

  cout << "\n\nIntegrate area of gammas = " << gamma_area <<  endl;
  cout << "Integrate area of pions = " << Pi0_area << endl;
  cout << "Integrate area of etas = " << Eta_area << endl;
  cout << "Pi0/Eta ratio = " << Pi0_area/Eta_area << endl;
  cout << "gamma/Pi0 ratio = " << gamma_area/Pi0_area << endl;
  cout << "gamma/Eta ratio = " << gamma_area/Eta_area << endl;
  cout << "gamma/(Pi0+Eta) ratio = " << gamma_area/(Pi0_area+Eta_area) << "\n" << endl;

  TLatex *latex = new TLatex();
  latex->SetTextSize(0.04);
  latex->SetTextAlign(33);  //align at top
  latex->SetTextFont(42);
  latex->DrawLatex(0.35,300000,"pp, #sqrt{s} = 7 TeV, m_{ee} calculated");
  latex->DrawLatex(0.35,90000,"using the MC truth observables");

  TLegend *legend = new TLegend(0.6,0.5,0.89,0.89);
  legend->SetBorderSize(0);
  legend->SetFillColor(kWhite);
  legend->SetTextSizePixels(1000);
  legend->SetTextFont(42);
  legend->SetTextSize(0.05);
  legend->AddEntry(hmee_MC,"N_{e^{#pm}e^{#mp}}","lep");
  legend->AddEntry(dist_MC,"all sources ","l");
  legend->AddEntry(gamma,"#gamma #rightarrow e^{+} e^{-}","l");
  legend->AddEntry(Pi0,"#pi^{0} #rightarrow e^{+} e^{-} #gamma","l");
  legend->AddEntry(Eta,"#eta #rightarrow e^{+} e^{-} #gamma","l");
  legend->Draw();

  /*-------------------------------- Invariant mass distribution MC truth with esd tracks------------------------------------*/
  
  TCanvas *c4 = new TCanvas("c4","c4",10,10,1200,900);

  fGaussExp->SetParameters(120000,0.015,0.012,-10);

  Double_t n_pi0_esd;
  Double_t m_pi0_esd;
  Double_t n_eta_esd;
  Double_t m_eta_esd;
  Double_t par0_esd;
  Double_t par1_esd;
  Double_t par2_esd;
  Double_t par3_esd;

  c4->Clear();
  c4->Divide(3,1);
  c4->cd(1)->SetLogy();
  TH1F *hPi0mee_esd = (TH1F*)Pi0mee_esd->Projection(1);
  dalitz_Pi0->SetParameters(1,0.13957);
  //hPi0mee_esd->SetStats(0);
  hPi0mee_esd->SetMarkerSize(0.5);
  hPi0mee_esd->SetMarkerColor(1);
  hPi0mee_esd->SetLineColor(1);
  dalitz_Pi0->SetLineColor(1);
  hPi0mee_esd->Fit("dalitz_Pi0","WLI", "  ", 0,0.2);
  hPi0mee_esd->Draw("E0E1");
  n_pi0_esd = dalitz_Pi0->GetParameter(0);
  m_pi0_esd = dalitz_Pi0->GetParameter(1);
  c4->cd(2)->SetLogy();
  TH1F *hEtamee_esd = (TH1F*)Etamee_esd->Projection(1);
  //hEtamee_esd->SetStats(0);
  hEtamee_esd->SetMarkerSize(0.5);
  hEtamee_esd->SetMarkerColor(4);
  hEtamee_esd->SetLineColor(4);
  dalitz_Eta->SetLineColor(4);
  dalitz_Eta->SetParameters(1.,0.54785);
  hEtamee_esd->Fit("dalitz_Eta","WLI","  ",0,0.6);
  hEtamee_esd->Draw("E0E1");
  n_eta_esd = dalitz_Eta->GetParameter(0);
  m_eta_esd = dalitz_Eta->GetParameter(1);
  cout << n_eta_esd << "          " << m_eta_esd << endl;
  c4->cd(3)->SetLogy();
  TH1F *hgammamee_esd = (TH1F*)gammamee_esd->Projection(1);
  //hgammamee_esd->SetStats(0);
  hgammamee_esd->SetMarkerSize(0.5);
  hgammamee_esd->SetMarkerColor(8);
  hgammamee_esd->SetLineColor(8);
  hgammamee_esd->Fit("fGaussExp","WLI","  ",0,0.16);
  hgammamee_esd->Draw("E0E1");
  par0_esd = fGaussExp->GetParameter(0);
  par1_esd = fGaussExp->GetParameter(1);
  par2_esd = fGaussExp->GetParameter(2);
  par3_esd = fGaussExp->GetParameter(3);


  TF1 *Eta2 = new TF1("Eta2","dalitz_Eta",0.,2.);
  Eta2->SetParameters(n_eta_esd,m_eta_esd);
  Eta2->SetLineColor(4);
  Eta2->SetLineWidth(2);

  TF1 *Pi02 = new TF1("Pi02","dalitz_Pi0",0.,1.);
  Pi02->SetParameters(n_pi0_esd,m_pi0_esd);
  Pi02->SetLineColor(1);
  Pi02->SetLineWidth(2);
  
  TF1 *gamma2 = new TF1("gamma2","fGaussExp",0.,1.);
  gamma2->SetParameters(par0_esd,par1_esd,par2_esd,par3_esd);
  gamma2->SetLineColor(8);
  gamma2->SetLineWidth(2);

  TF1 *dist_esd = new TF1("dist_esd","Pi02+Eta2+gamma2",0.,1.);
  dist_esd->SetLineWidth(2);
  dist_esd->SetLineColor(2);
  dist_esd->SetParameters(n_pi0_esd,m_pi0_esd,n_eta_esd,m_eta_esd,par0_esd,par1_esd,par2_esd,par3_esd);
  dist_esd->FixParameter(1,m_pi0_esd);
  dist_esd->FixParameter(3,m_eta_esd);
  //dist_esd->FixParameter(0,n_pi0_esd);
  //dist_esd->FixParameter(2,n_eta_esd);
  //dist_esd->FixParameter(4,par0_esd);
  dist_esd->FixParameter(5,par1_esd);
  //dist_esd->FixParameter(6,par2_esd);
  //dist_esd->FixParameter(7,par3_esd);
    
  TH1F *hmee_esd = (TH1F*)mee_esd->Projection(1);

  hmee_esd->GetXaxis()->SetTitle("m_{ee} (GeV/c^{2})");
  hmee_esd->GetYaxis()->SetTitle("dN/dm_{ee} (c^{2}/GeV)");

 cout << "\n\n***************** ESD DIST ******************\n\n" << endl;
	  
  TCanvas *esd_dist = new TCanvas("esd_dist","esd_dist",10,10,1200,900);
  esd_dist->Clear();
  esd_dist->SetLogy();
 
  hmee_esd->SetMarkerSize(0.7);
  hmee_esd->SetTitle("  ");
  hmee_esd->GetXaxis()->CenterTitle();
  hmee_esd->GetYaxis()->CenterTitle();
  hmee_esd->SetStats(0);
  hmee_esd->SetMarkerColor(2);
  hmee_esd->SetLineColor(2);
  hmee_esd->SetMarkerStyle(20);
  TFitResultPtr fitresults2 = hmee_esd->Fit("dist_esd","SWLI"," ",0.005,0.6);
  hmee_esd->Draw("E0E1");

  Double_t n_pi0_esd2 = dist_esd->GetParameter(0);
  Double_t m_pi0_esd2 = dist_esd->GetParameter(1);
  Double_t n_eta_esd2 = dist_esd->GetParameter(2);
  Double_t m_eta_esd2 = dist_esd->GetParameter(3);
  Double_t par0_esd2 = dist_esd->GetParameter(4);
  Double_t par1_esd2 = dist_esd->GetParameter(5);
  Double_t par2_esd2 = dist_esd->GetParameter(6);
  Double_t par3_esd2 = dist_esd->GetParameter(7);

  gamma2->SetParameters(par0_esd2,par1_esd2,par2_esd2,par3_esd2);
  Eta2->SetParameters(n_eta_esd2,m_eta_esd2);
  Pi02->SetParameters(n_pi0_esd2,m_pi0_esd2);
  Pi02->SetLineStyle(10);
  Eta2->SetLineStyle(3);
  gamma2->SetLineStyle(9);
 
  Pi02->Draw("same");
  Eta2->Draw("same");
  gamma2->Draw("same");

  cout<<"NDF = " << dist_esd->GetNDF()<<endl;   

  //Integrals and errors
  TMatrixDSym cov2(8);
  cov2 = fitresults2->GetCovarianceMatrix();    
  cov2.Print();

  //gamma
  Double_t par_gamma2[4];
  TMatrixDSym cov_gamma2(4);

  for(Int_t i = 0; i<4; i++)
    {
      par_gamma2[i]=gamma2->GetParameter(i);
      for(Int_t j = 0; j< 4; j++)
	{
	  cov_gamma2(i,j) = cov2(i+4,j+4);
	  if(i==j && i==3)
	    cov_gamma2(i,j)=0;
	}
    }

  //Pi0
  Double_t par_Pi02[2];
  TMatrixDSym cov_Pi02(2);

  for(Int_t i = 0; i<2; i++)
    {
      par_Pi02[i]=Pi02->GetParameter(i);
      for(Int_t j = 0; j< 2; j++)
	{
	  cov_Pi02(i,j) = cov2(i,j);
	}
    }

  //Eta
  Double_t par_Eta2[2];
  TMatrixDSym cov_Eta2(2);

  for(Int_t i = 0; i<2; i++)
    {
      par_Eta2[i]=Eta2->GetParameter(i);
      for(Int_t j = 0; j< 2; j++)
	{
	  cov_Eta2(i,j) = cov2(i+2,j+2);
	}
    }

  Double_t gamma_area2 = gamma2->Integral(0.,0.6)*hmee_esd->GetNbinsX()/0.6;
  Double_t Pi0_area2 = Pi02->Integral(0.,0.6)*hmee_esd->GetNbinsX()/0.6;
  Double_t Eta_area2 = Eta2->Integral(0.,0.6)*hmee_esd->GetNbinsX()/0.6;
  Double_t gamma_error_area2 = gamma2->IntegralError(0,0.6,par_gamma2,cov_gamma2.GetMatrixArray())*hmee_esd->GetNbinsX()/0.6;
  Double_t Pi0_error_area2 = Pi02->IntegralError(0,0.6,par_Pi02,cov_Pi02.GetMatrixArray())*hmee_esd->GetNbinsX()/0.6;
  Double_t Eta_error_area2 = Eta2->IntegralError(0,0.6,par_Eta2,cov_Eta2.GetMatrixArray())*hmee_esd->GetNbinsX()/0.6;

cout << "\n\nIntegrate area of gammas = " << gamma_area2 << " +/- " << gamma_error_area2 << endl;
  cout << "Integrate area of pions = " << Pi0_area2 << " +/- " << Pi0_error_area2 <<endl;
  cout << "Integrate area of etas = " << Eta_area2 << " +/- " << Eta_error_area2 <<endl;
  cout << "Pi0/Eta ratio = " << Pi0_area2/Eta_area2 << " +/- " << Pi0_area2/Eta_area2*TMath::Sqrt(TMath::Power(Pi0_error_area2/Pi0_area2,2)+TMath::Power(Eta_error_area2/Eta_area2,2)) << endl;
  cout << "gamma/Pi0 ratio = " << gamma_area2/Pi0_area2 <<" +/- " << gamma_area2/Pi0_area2*TMath::Sqrt(TMath::Power(gamma_error_area2/gamma_area2,2)+TMath::Power(Pi0_error_area2/Pi0_area2,2))<< endl;
  cout << "gamma/Eta ratio = " << gamma_area2/Eta_area2 << " +/- " << gamma_area2/Eta_area2*TMath::Sqrt(TMath::Power(gamma_error_area2/gamma_area2,2)+TMath::Power(Eta_error_area2/Eta_area2,2))<< endl;
  cout << "gamma/(Pi0+Eta) ratio = " << gamma_area2/(Eta_area2+Pi0_area2) << " +/- " << gamma_area2/(Eta_area2+Pi0_area2)*TMath::Sqrt(TMath::Power(gamma_error_area2/gamma_area2,2)+TMath::Power(TMath::Sqrt(TMath::Power(Pi0_error_area2,2)+TMath::Power(Eta_error_area2,2))/(Pi0_area2+Eta_area2),2)) << "\n" << endl;

  TLatex *latex2 = new TLatex();
  latex2->SetTextSize(0.04);
  latex2->SetTextAlign(33);  //align at top
  latex2->SetTextFont(42);
  latex2->DrawLatex(0.35,200000,"pp, #sqrt{s} = 7 TeV, m_{ee} calculated");
  latex2->DrawLatex(0.35,80000,"using the ESD observables");

  TLegend *legend2 = new TLegend(0.6,0.5,0.89,0.89);
  legend2->SetBorderSize(0);
  legend2->SetFillColor(kWhite);
  legend2->SetTextSizePixels(1000);
  legend2->SetTextFont(42);
  legend2->SetTextSize(0.05);
  legend2->AddEntry(hmee_esd,"N_{e^{#pm}e^{#mp}}","lep");
  legend2->AddEntry(dist_esd,"all sources ","l");
  legend2->AddEntry(gamma2,"#gamma #rightarrow e^{+} e^{-}","l");
  legend2->AddEntry(Pi02,"#pi^{0} #rightarrow e^{+} e^{-} #gamma","l");
  legend2->AddEntry(Eta2,"#eta #rightarrow e^{+} e^{-} #gamma","l");
  legend2->Draw();

/*------------------------------- Inclusive and associated pt distributions -------------------------------------*/
  
  TCanvas *c5 = new TCanvas("c5","c5",10,10,1200,900);
  c5->Clear();
  c5->Divide(1,2);
  c5->cd(1)->SetLogy();
  TH1F *hAss=(TH1F*)l->FindObject("fHistAssPt");
  hAss->SetLineColor(1);
  hAss->SetMarkerColor(1);
  hAss->SetMarkerStyle(20);
  hAss->SetMarkerSize(0.8);
  hAss->Draw("E0E1");
  c5->cd(2)->SetLogy();
  TH1F *hIncl=(TH1F*)l->FindObject("fHistInclPt");
  hIncl->SetLineColor(2);
  hIncl->SetMarkerColor(2);
  hIncl->SetMarkerSize(0.8);
  hIncl->SetMarkerStyle(20);
  hIncl->Draw("E0E1");
  TCanvas *IncAss = new TCanvas("IncAss", "IncAss", 10,10,1200,900);
  IncAss->Clear();
  IncAss->SetLogy();
  hAss->SetStats(0);
  hIncl->SetStats(0);
  hAss->SetTitle(" ");
  hIncl->SetTitle(" ");
  hAss->GetXaxis()->CenterTitle();
  hAss->GetYaxis()->CenterTitle();
  hIncl->GetXaxis()->CenterTitle();
  hIncl->GetYaxis()->CenterTitle();
  hAss->Draw("E0E1");
  hIncl->Draw("E0E1 same");

  TLegend *legend6 = new TLegend(0.55,0.55,0.89,0.89);
  legend6->SetBorderSize(0);
  legend6->SetFillColor(kWhite);
  legend6->SetTextSizePixels(1000);
  legend6->SetTextFont(42);
  legend6->SetTextSize(0.045);
  legend6->AddEntry(hIncl,"Inclusive e^{#pm}","lep");
  legend6->AddEntry(hAss,"Associated e^{#pm}","lep");
  legend6->Draw();

  /*------------------------------ LS US distribution with MC truth observables ---------------------------------------*/
  
  TCanvas *c6 = new TCanvas("c6","c6",10,10,1200,900);
  c6->Clear();
  c6->Divide(2,3); 
  c6->cd(1);
  c6->cd(1)->SetLogz();
  TH2D *hLike = (TH2D*)spLike->Projection(0,1);
  hLike->SetStats(0);
  hLike->Draw("lego2z");
  c6->cd(2);
  c6->cd(2)->SetLogz();
  TH2D *hUnlike = (TH2D*)spUnlike->Projection(0,1);
  hUnlike->SetStats(0);
  hUnlike->Draw("lego2z");
  TH1D *h1like =(TH1D*)spLike->Projection(1);
  h1like->Sumw2();
  c6->cd(3);
  c6->cd(3)->SetLogy();
  h1like->SetMarkerStyle(20);
  h1like->SetMarkerColor(1);
  h1like->SetMarkerSize(0.8);
  h1like->SetLineColor(1);
  h1like->Draw();
  TH1D *h1unlike =(TH1D*)spUnlike->Projection(1);
  h1unlike->Sumw2();
  c6->cd(4);
  c6->cd(4)->SetLogy();
  h1unlike->SetMarkerStyle(20);
  h1unlike->SetMarkerColor(2);
  h1unlike->SetMarkerSize(0.8);
  h1unlike->SetLineColor(2);
  h1unlike->Draw();
  c6->cd(5);
  c6->cd(5)->SetLogy();
  h1unlike->Draw("E0E1");
  h1like->Draw("E0E1 same"); 
  
  

  /*
  cout << "\n\n***************** MC COMB ******************\n\n" << endl;
  
  TCanvas *MC_comb = new TCanvas("MC_comb","MC_comb",10,10,1200,900);
  MC_comb->SetLogy();
  TH1F *h3 = (TH1F*)h1unlike->Clone(); 
  h3->GetYaxis()->SetRangeUser(0,5000);
  h3->Add(h1like,-1);
  h3->SetStats(0);  
  h3->Sumw2();
  h3->SetTitle(" ");
  
  for(Int_t i=0; i<31; i++)
    {
      if(h3->GetBinContent(i)<0)
	h3->SetBinContent(i,0);
      if(h3->GetBinContent(i)==0)
	h3->SetBinError(i,1);
    }
  */
  /*
 TF1 *dist_MC2 = new TF1("dist_MC2","[0]+Pi0+Eta+gamma",0,1);
  
  dist_MC2->SetLineColor(2);
  dist_MC2->SetParameters(400,n_pi0_MC2,m_pi0_MC2,n_eta_MC2,m_eta_MC2,par0_MC2,par1_MC2,par2_MC2,par3_MC2);
  //dist_MC2->FixParameter(0,5);
  dist_MC2->FixParameter(2,m_pi0_MC2);
  dist_MC2->FixParameter(4,m_eta_MC2);
  //dist_MC2->FixParameter(5,par0_MC2);
  dist_MC2->FixParameter(6,par1_MC2);
  //dist_MC2->FixParameter(0,n_pi0_MC);
  //dist_MC2->FixParameter(2,n_eta_MC2);
    
  h3->Fit("dist_MC2","L"," ",0,0.3);
  h3->Draw("E0E1");

  TF1 *Pi03 = new TF1("Pi03","dalitz_Pi0",0.,1.);
  TF1 *Eta3 = new TF1("Eta3","dalitz_Eta",0.,1.);
  TF1 *gamma3 = new TF1("gamma3","fGaussExp",0.,1);
  TF1 *cost3 = new TF1("cost3","[0]",0,1);

  Pi03->SetLineColor(1);
  Eta3->SetLineColor(4);
  gamma3->SetLineColor(8);
  Pi03->SetLineWidth(2);
  Eta3->SetLineWidth(2);
  gamma3->SetLineWidth(2);
  cost3->SetLineColor(49);
  Pi03->SetLineStyle(10);
  Eta3->SetLineStyle(3);
  gamma3->SetLineStyle(9);
  cost3->SetLineStyle(7);

  Pi03->SetParameters(dist_MC2->GetParameter(1),dist_MC2->GetParameter(2));
  Eta3->SetParameters(dist_MC2->GetParameter(3),dist_MC2->GetParameter(4));
  gamma3->SetParameters(dist_MC2->GetParameter(5),dist_MC2->GetParameter(6),dist_MC2->GetParameter(7),dist_MC2->GetParameter(8));
  cost3->SetParameter(0,dist_MC2->GetParameter(0));
 
  Pi03->Draw("same");
  Eta3->Draw("same");
  gamma3->Draw("same");
  cost3->Draw("same");

  cout<<"NDF = " << dist_MC2->GetNDF()<<endl; 
  
  Double_t gamma_area3 = gamma3->Integral(0,0.6)*h3->GetNbinsX()/0.6;
  Double_t Pi0_area3 = Pi03->Integral(0,1)*h3->GetNbinsX()/0.6;
  Double_t Eta_area3 = Eta3->Integral(0,1)*h3->GetNbinsX()/0.6;

  cout << "\n\nIntegrate area of gammas = " << gamma_area3 << endl;
  cout << "Integrate area of pions = " << Pi0_area3 << endl;
  cout << "Integrate area of etas = " << Eta_area3 << endl;
  cout << "Pi0/Eta ratio = " << Pi0_area3/Eta_area3 << endl;
  cout << "gamma/Pi0 ratio = " << gamma_area3/Pi0_area3 << endl;
  cout << "gamma/Eta ratio = " << gamma_area3/Eta_area3 << endl;
  cout << "gamma/(Pi0+Eta) ratio = " << gamma_area3/(Pi0_area3+Eta_area3) << "\n" << endl;

  TLatex *latex3 = new TLatex();
  latex3->SetTextSize(0.04);
  latex3->SetTextAlign(33);  //align at top
  latex3->SetTextFont(42);
  latex3->DrawLatex(0.21,300000,"pp, #sqrt{s} = 7 TeV,");
  latex3->DrawLatex(0.21,150000,"N_{e^{#pm}e^{#mp}} - N_{e^{#pm}e^{#pm}}");

  TLegend *legend3 = new TLegend(0.65,0.55,0.89,0.89);
  legend3->SetBorderSize(0);
  legend3->SetFillColor(kWhite);
  legend3->SetTextSizePixels(1000);
  legend3->SetTextFont(42);
  legend3->SetTextSize(0.045);
  legend3->AddEntry(h3,"N_{e^{#pm}e^{#mp}}","lep");
  legend3->AddEntry(dist_MC2,"all sources ","l");
  legend3->AddEntry(gamma3,"#gamma #rightarrow e^{+} e^{-}","l");
  legend3->AddEntry(Pi03,"#pi^{0} #rightarrow e^{+} e^{-} #gamma","l");
  legend3->AddEntry(Eta3,"#eta #rightarrow e^{+} e^{-} #gamma","l");
  legend3->AddEntry(cost3,"const back","l");
  legend3->Draw();



  /*----------------------------- LS US distribution with esd tracks -------------------------------------------*/
  
  TCanvas *c7 = new TCanvas("c7","c7",10,10,1200,900);
  c7->Clear();
  c7->Divide(2,2); 
  c7->cd(1);
  c7->cd(1)->SetLogz();

  TH2D *hLS = (TH2D*)spLS->Projection(0,1);
  hLS->SetStats(0);
  hLS->Draw("lego2z");
  c7->cd(2);
  c7->cd(2)->SetLogz();
  TH2D *hUS = (TH2D*)spUS->Projection(0,1);
  hUS->SetStats(0);
  hUS->Draw("lego2z");
  TH1D *h1LS =(TH1D*)spLS->Projection(1);
  h1LS->Sumw2();
  c7->cd(3);
  c7->cd(3)->SetLogy();
  h1LS->SetMarkerStyle(20);
  h1LS->SetMarkerColor(1);
  h1LS->SetMarkerSize(0.8);
  h1LS->SetLineColor(1);
  h1LS->Draw();
  TH1D *h1US =(TH1D*)spUS->Projection(1);
  h1US->Sumw2();
  c7->cd(4);
  c7->cd(4)->SetLogy();
  h1US->SetMarkerStyle(20);
  h1US->SetMarkerColor(2);
  h1US->SetMarkerSize(0.8);
  h1US->SetLineColor(2);
  h1US->Draw();
  c7->cd(5);
  c7->cd(5)->SetLogy();
  h1US->Draw("E0E1");
  h1LS->Draw("E0E1 same"); 

  TCanvas *LsUs = new TCanvas("LsUs","LsUs",10,10,1200,900);
  LsUs->Clear();
  LsUs->SetLogy();
  h1US->SetStats(0);
  h1LS->SetStats(0);
  h1US->SetTitle(" ");
  h1US->GetYaxis()->CenterTitle();
  h1US->GetXaxis()->CenterTitle();
  h1US->GetYaxis()->SetRangeUser(800,1000000);
  h1US->GetYaxis()->SetTitle("dN/dm_{ee} (c^{2}/GeV)");
  h1US->GetXaxis()->SetTitle("m_{ee} (GeV/c^{2})");
  h1US->Draw("E0E1");
  h1LS->Draw("E0E1 Same");

  TLatex *latex5 = new TLatex();
  latex5->SetTextSize(0.04);
  latex5->SetTextAlign(33);  //align at top
  latex5->SetTextFont(42);
  latex5->DrawLatex(0.21,500000.,"pp, #sqrt{s} = 7 TeV,");
  latex5->DrawLatex(0.35,200000.,"0 (GeV/c) < p_{T} < 8 (GeV/c)");

  TLegend *legend5 = new TLegend(0.64,0.6,0.89,0.85);
  legend5->SetBorderSize(0);
  legend5->SetFillColor(kWhite);
  legend5->SetTextSizePixels(1000);
  legend5->SetTextFont(42);
  legend5->SetTextSize(0.06);
  legend5->AddEntry(h1LS," N_{e^{#pm} e^{#pm}} ","lep");
  legend5->AddEntry(h1US," N_{e^{#pm} e^{#mp}}","lep");
  legend5->Draw();

  cout << "\n\n************* Normalization of LS distribution ***************\n\n" << endl;

  TCanvas *LSUS_ratio = new TCanvas("LSUS_ratio","LSUS_ratio",10,10,1200,900);
  LSUS_ratio->SetLogy();

  Double_t ratio[60];
  Double_t err_ratio[60]; 
  Double_t mee[60];
 
  for(Int_t i=2; i<62; i++)
    {
      if(h1LS->GetBinContent(i)!=0)
	{
	  ratio[i-2] = h1US->GetBinContent(i)/h1LS->GetBinContent(i);
	  err_ratio[i-2]=ratio[i-2]*(h1US->GetBinError(i)/h1US->GetBinContent(i)+h1LS->GetBinError(i)/h1LS->GetBinContent(i));
	  mee[i-2]=0.6/h1LS->GetNbinsX()*i;
	}
    }

  TF1 *fit = new TF1("fit","pol1",0,1);
  fit->SetLineWidth(2.5);
  TGraphErrors *r = new TGraphErrors(59,mee,ratio,0,err_ratio);
  r->SetTitle(" ");
  r->SetMarkerStyle(20);
  r->SetMarkerSize(0.8);
  r->Fit("fit", "  ", "  ",0.3,0.6);
  r->Draw("AP");
  r->GetXaxis()->SetTitle("m_{ee} (GeV/c^{2})");
  r->GetYaxis()->SetTitle("dN/dm_{ee}^{US}/dN/dm_{ee}^{LS}");
  // r->GetXaxis()->SetTitleSize(0.4);
  r->GetXaxis()->CenterTitle();
  r->GetYaxis()->CenterTitle();
  r->GetYaxis()->SetRangeUser(0.8,3);

  TLatex *latex8 = new TLatex();
  latex8->SetTextSize(0.06);
  latex8->SetTextAlign(33);  //align at top
  latex8->SetTextFont(42);
  latex8->DrawLatex(0.5,2.7,"f(m_{ee}) = p_{0}+p_{1}#upoint m_{ee}");
  latex8->SetTextSize(0.04);
  latex8->DrawLatex(0.5,2.2,"#chi^{2}/Ndf = 10.79 / 29");
  latex8->DrawLatex(0.5,2,"p_{0} = 1.11 #pm 0.02");
  latex8->DrawLatex(0.585,1.85,"p_{1} = -0.11 #pm 0.05 c^{2}/GeV");


  TH1F *h1LS_norm = new TH1F("h1LS_norm", "  ", 60,0,0.6);

  for(Int_t i=1; i<62; i++)
    {
  h1LS_norm->SetBinContent(i,h1LS->GetBinContent(i)*(fit->GetParameter(0)+fit->GetParameter(1)*0.6/h1LS_norm->GetNbinsX()*i));
    }

  TCanvas *LsUs_norm = new TCanvas("LsUs_norm","LsUs_norm",10,10,1200,900);
  LsUs_norm->Clear();
  LsUs_norm->SetLogy();
  h1US->SetStats(0);
  h1LS_norm->SetStats(0);
  h1LS_norm->SetLineColor(1);
  h1US->SetTitle(" ");
  h1US->GetYaxis()->CenterTitle();
  h1US->GetXaxis()->CenterTitle();
  h1US->GetYaxis()->SetRangeUser(800,1000000);
  h1US->GetYaxis()->SetTitle("dN/dm_{ee} (c^{2}/GeV)");
  h1US->GetXaxis()->SetTitle("m_{ee} (GeV/c^{2})");
  h1US->Draw("E0E1");
  h1LS_norm->Draw("E0E1 Same");

  TLatex *latex7 = new TLatex();
  latex7->SetTextSize(0.04);
  latex7->SetTextAlign(33);  //align at top
  latex7->SetTextFont(42);
  latex7->DrawLatex(0.18,480000.,"pp, #sqrt{s} = 7 TeV,");
  latex7->DrawLatex(0.30,200000.,"0 (GeV/c) < p_{T} < 8 (GeV/c)");

  TLegend *legend7 = new TLegend(0.58,0.6,0.83,0.85);
  legend7->SetBorderSize(0);
  legend7->SetFillColor(kWhite);
  legend7->SetTextSizePixels(1000);
  legend7->SetTextFont(42);
  legend7->SetTextSize(0.06);
  legend7->AddEntry(h1LS_norm," N_{e^{#pm} e^{#pm} normalized} ","lep");
  legend7->AddEntry(h1US," N_{e^{#pm} e^{#mp}}","lep");
  legend7->Draw();

 cout << "\n\n***************** ESD COMB ******************\n\n" << endl;

  TCanvas *esd_comb = new TCanvas("esd_comb", "esd_comb", 10,10,1200,900); 
  esd_comb->SetLogy();
  TH1F *h4 = (TH1F*)h1US->Clone(); 
  h4->GetYaxis()->SetRange(0,5000);
  h4->Add(h1LS_norm,-1);
  h4->SetStats(0);

   TF1 *dist_esd2 = new TF1("dist_esd2","Pi02+Eta2+gamma2",0,1);
   dist_esd2->SetLineColor(2);
   dist_esd2->SetParameters(n_pi0_esd,m_pi0_esd,n_eta_esd,m_eta_esd,par0_esd,par1_esd, par2_esd, par3_esd); 

   dist_esd2->FixParameter(1,m_pi0_esd2);
   dist_esd2->FixParameter(3,m_eta_esd2);
   //dist_esd2->FixParameter(0,n_pi0_esd2);
   //dist_esd2->FixParameter(2,n_eta_esd2);
   //dist_esd2->FixParameter(4,par0_esd2);
   dist_esd2->FixParameter(5,par1_esd2);
   //dist_esd2->FixParameter(6,par2_esd2);   
   //dist_esd2->FixParameter(7,par3_esd2);

   TFitResultPtr fitresults4 = h4->Fit("dist_esd2","SWLI","  ",0.005,0.5);
   h4->Draw("E1E0");
   
   TF1 *Pi04 = new TF1("Pi04","dalitz_Pi0",0.,1.);
   TF1 *Eta4 = new TF1("Eta4","dalitz_Eta",0.,1.);
   TF1 *gamma4 = new TF1("gamma4","fGaussExp",0.,1);

   Pi04->SetLineColor(1);
   Eta4->SetLineColor(4);
   gamma4->SetLineColor(3);
   Pi04->SetLineWidth(2);
   Eta4->SetLineWidth(2);
   gamma4->SetLineWidth(2);
   Pi04->SetLineStyle(10);
   Eta4->SetLineStyle(3);
   gamma4->SetLineStyle(9);

   Pi04->SetParameters(dist_esd2->GetParameter(0),dist_esd2->GetParameter(1));
   Eta4->SetParameters(dist_esd2->GetParameter(2),dist_esd2->GetParameter(3));  
   gamma4->SetParameters(dist_esd2->GetParameter(4),dist_esd2->GetParameter(5),dist_esd2->GetParameter(6),dist_esd2->GetParameter(7));

   Pi04->Draw("same");
   Eta4->Draw("same");
   gamma4->Draw("same");
   
   cout<<"NDF = " << dist_esd2->GetNDF()<<endl;   

  //Integrals and errors
  TMatrixDSym cov4(8);
  cov4 = fitresults4->GetCovarianceMatrix();    
  cov4.Print();

  //gamma
  Double_t par_gamma4[4];
  TMatrixDSym cov_gamma4(4);

  for(Int_t i = 0; i<4; i++)
    {
      par_gamma4[i]=gamma4->GetParameter(i);
      for(Int_t j = 0; j< 4; j++)
	{
	  cov_gamma4(i,j) = cov4(i+4,j+4);
	  if(i==j && i==3)
	    cov_gamma4(i,j)=0;
	}
    }

  //Pi0
  Double_t par_Pi04[2];
  TMatrixDSym cov_Pi04(2);

  for(Int_t i = 0; i<2; i++)
    {
      par_Pi04[i]=Pi04->GetParameter(i);
      for(Int_t j = 0; j< 2; j++)
	{
	  cov_Pi04(i,j) = cov4(i,j);
	}
    }

  //Eta
  Double_t par_Eta4[2];
  TMatrixDSym cov_Eta4(2);

  for(Int_t i = 0; i<2; i++)
    {
      par_Eta4[i]=Eta4->GetParameter(i);
      for(Int_t j = 0; j< 2; j++)
	{
	  cov_Eta4(i,j) = cov4(i+2,j+2);
	}
    }

  Double_t gamma_area4 = gamma4->Integral(0,0.6)*hmee_esd->GetNbinsX()/0.6;
  Double_t Pi0_area4 = Pi04->Integral(0,0.6)*hmee_esd->GetNbinsX()/0.6;
  Double_t Eta_area4 = Eta4->Integral(0,0.6)*hmee_esd->GetNbinsX()/0.6;
  Double_t gamma_error_area4 = gamma4->IntegralError(0,0.6,par_gamma4,cov_gamma4.GetMatrixArray())*hmee_esd->GetNbinsX()/0.06;
  Double_t Pi0_error_area4 = Pi04->IntegralError(0,0.6,par_Pi04,cov_Pi04.GetMatrixArray())*hmee_esd->GetNbinsX()/0.6;
  Double_t Eta_error_area4 = Eta4->IntegralError(0,0.6,par_Eta4,cov_Eta4.GetMatrixArray())*hmee_esd->GetNbinsX()/0.6;

  cout << "\n\nIntegrate area of gammas = " << gamma_area4 << " +/- " << gamma_error_area4 << endl;
  cout << "Integrate area of pions = " << Pi0_area4 << " +/- " << Pi0_error_area4 <<endl;
  cout << "Integrate area of etas = " << Eta_area4 << " +/- " << Eta_error_area4 <<endl;
  cout << "Pi0/Eta ratio = " << Pi0_area4/Eta_area4 << " +/- " << Pi0_area4/Eta_area4*TMath::Sqrt(TMath::Power(Pi0_error_area4/Pi0_area4,2)+TMath::Power(Eta_error_area4/Eta_area4,2)) << endl;
  cout << "gamma/Pi0 ratio = " << gamma_area4/Pi0_area4 <<" +/- " << gamma_area4/Pi0_area4*TMath::Sqrt(TMath::Power(gamma_error_area4/gamma_area4,2)+TMath::Power(Pi0_error_area4/Pi0_area4,2))<< endl;
  cout << "gamma/Eta ratio = " << gamma_area4/Eta_area4 << " +/- " << gamma_area4/Eta_area4*TMath::Sqrt(TMath::Power(gamma_error_area4/gamma_area4,2)+TMath::Power(Eta_error_area4/Eta_area4,2))<< endl;
  cout << "gamma/(Pi0+Eta) ratio = " << gamma_area4/(Eta_area4+Pi0_area4) << " +/- " << gamma_area4/(Eta_area4+Pi0_area4)*TMath::Sqrt(TMath::Power(gamma_error_area4/gamma_area4,2)+TMath::Power(TMath::Sqrt(TMath::Power(Pi0_error_area4,2)+TMath::Power(Eta_error_area4,2))/(Pi0_area4+Eta_area4),2)) << "\n" << endl; 
					 
						     
  TLatex *latex4 = new TLatex();
  latex4->SetTextSize(0.04);
  latex4->SetTextAlign(33);  //align at top
  latex4->SetTextFont(42);
  latex4->DrawLatex(0.25,300000,"pp, #sqrt{s} = 7 TeV,");
  latex4->DrawLatex(0.25,150000,"N_{e^{#pm}e^{#mp}} - N_{e^{#pm}e^{#pm}}");
  
  TLegend *legend4 = new TLegend(0.65,0.55,0.89,0.89);
  legend4->SetBorderSize(0);
  legend4->SetFillColor(kWhite);
  legend4->SetTextSizePixels(1000);
  legend4->SetTextFont(42);
  legend4->SetTextSize(0.045);
  legend4->AddEntry(h4,"N_{e^{#pm}e^{#mp}}","lep");
  legend4->AddEntry(dist_esd2,"all sources ","l");
  legend4->AddEntry(gamma4,"#gamma #rightarrow e^{+} e^{-}","l");
  legend4->AddEntry(Pi04,"#pi^{0} #rightarrow e^{+} e^{-} #gamma","l");
  legend4->AddEntry(Eta4,"#eta #rightarrow e^{+} e^{-} #gamma","l");
  legend4->Draw();

  /*
  //print values
  TGraph *point_gamma = new TGraph();
  TGraph *point_Pi0 = new TGraph();
  TGraph *point_Eta = new TGraph();
  TGraph *point_tot = new TGraph();
  point_gamma->SetMarkerStyle(20);
  point_Pi0->SetMarkerStyle(20);
  point_Eta->SetMarkerStyle(20);
  point_tot->SetMarkerStyle(20);
  point_gamma->SetMarkerColor(8);
  point_Pi0->SetMarkerColor(1);
  point_Eta->SetMarkerColor(4);
  point_tot->SetMarkerColor(2);

  for(Int_t i=0; i<16; i++)
    {
      cout <<"tot = " << dist_esd2->Eval(0.005+i*0.01)<< " sum cont = " <<gamma4->Eval(0.005+i*0.01)+Pi04->Eval(0.005+i*0.01)+Eta4->Eval(0.005+i*0.01)<< " gamma = " << gamma4->Eval(0.005+i*0.01) << " Pi0 = "<<Pi04->Eval(0.005+i*0.01)<<" Eta = " <<Eta4->Eval(0.005+i*0.01) << endl ;
      point_gamma->SetPoint(i,0.005+i*0.01,gamma4->Eval(0.005+i*0.01));
      point_Pi0->SetPoint(i,0.005+i*0.01,Pi04->Eval(0.005+i*0.01));
      point_Eta->SetPoint(i,0.005+i*0.01,Eta4->Eval(0.005+i*0.01));
      point_tot->SetPoint(i,0.005+i*0.01,dist_esd2->Eval(0.005+i*0.01));
    }

  point_gamma->Draw("same P");
  point_Pi0->Draw("same P");
  point_Eta->Draw("same P");
  point_tot->Draw("same P");

  /*------------------------------ LS US distributions for different pt range -----------------------------------------*/
  
  TCanvas *c8 = new TCanvas("c8","c8",10,10,1200,900);
  c8->Clear();
  c8->Divide(3,3);
  
  TH1D *ULproj1 = (TH1D*)hUS->ProjectionX("ULproj1",1,2);
  TH1D *ULproj2 = (TH1D*)hUS->ProjectionX("ULproj2",2,3);
  TH1D *ULproj3 = (TH1D*)hUS->ProjectionX("ULproj3",3,4);
  TH1D *ULproj4 = (TH1D*)hUS->ProjectionX("ULproj4",4,5);
  TH1D *ULproj5 = (TH1D*)hUS->ProjectionX("ULproj5",5,6);
  TH1D *ULproj6 = (TH1D*)hUS->ProjectionX("ULproj6",6,7);
  TH1D *ULproj7 = (TH1D*)hUS->ProjectionX("ULproj7",7,8);
  TH1D *ULproj8 = (TH1D*)hUS->ProjectionX("ULproj8",8,9);
  TH1D *ULproj9 = (TH1D*)hUS->ProjectionX("ULproj9",9,10);
  
  TH1D *Lproj1 = (TH1D*)hLS->ProjectionX("Lproj1",1,2);
  TH1D *Lproj2 = (TH1D*)hLS->ProjectionX("Lproj2",2,3);
  TH1D *Lproj3 = (TH1D*)hLS->ProjectionX("Lproj3",3,4);
  TH1D *Lproj4 = (TH1D*)hLS->ProjectionX("Lproj4",4,5);
  TH1D *Lproj5 = (TH1D*)hLS->ProjectionX("Lproj5",5,6);
  TH1D *Lproj6 = (TH1D*)hLS->ProjectionX("Lproj6",6,7);
  TH1D *Lproj7 = (TH1D*)hLS->ProjectionX("Lproj7",7,8);
  TH1D *Lproj8 = (TH1D*)hLS->ProjectionX("Lproj8",8,9);
  TH1D *Lproj9 = (TH1D*)hLS->ProjectionX("Lproj9",9,10);

  c8->cd(1)->SetLogy();
  ULproj1->SetMarkerColor(2);
  ULproj1->SetLineColor(2);
  Lproj1->SetMarkerColor(1);
  Lproj1->SetLineColor(1);
  ULproj1->Draw("E0E1");
  Lproj1->Draw("E0E1 same");
  c8->cd(2)->SetLogy();
  ULproj2->SetMarkerColor(2);
  ULproj2->SetLineColor(2);
  Lproj2->SetMarkerColor(1);
  Lproj2->SetLineColor(1);
  ULproj2->Draw("E0E1");
  Lproj2->Draw("E0E1 same");
  c8->cd(3)->SetLogy();
  ULproj3->SetMarkerColor(2);
  ULproj3->SetLineColor(2);
  Lproj3->SetMarkerColor(1);
  Lproj3->SetLineColor(1);
  ULproj3->Draw("E0E1");
  Lproj3->Draw("E0E1 same");
  c8->cd(4)->SetLogy();
  ULproj4->SetMarkerColor(2);
  ULproj4->SetLineColor(2);
  Lproj4->SetMarkerColor(1);
  Lproj4->SetLineColor(1);
  ULproj4->Draw("E0E1");
  Lproj4->Draw("E0E1 same");
  c8->cd(5)->SetLogy();
  ULproj5->SetMarkerColor(2);
  ULproj5->SetLineColor(2);
  Lproj5->SetMarkerColor(1);
  Lproj5->SetLineColor(1);
  ULproj5->Draw("E0E1");
  Lproj5->Draw("E0E1 same");
  c8->cd(6)->SetLogy();
  ULproj6->SetMarkerColor(2);
  ULproj6->SetLineColor(2);
  Lproj6->SetMarkerColor(1);
  Lproj6->SetLineColor(1);
  ULproj6->Draw("E0E1");
  Lproj6->Draw("E0E1 same");
  c8->cd(7)->SetLogy();
  ULproj7->SetMarkerColor(2);
  ULproj7->SetLineColor(2);
  Lproj7->SetMarkerColor(1);
  Lproj7->SetLineColor(1);
  ULproj7->Draw("E0E1");
  Lproj7->Draw("E0E1 same");
  c8->cd(8)->SetLogy();
  ULproj8->SetMarkerColor(2);
  ULproj8->SetLineColor(2);
  Lproj8->SetMarkerColor(1);
  Lproj8->SetLineColor(1);
  ULproj8->Draw("E0E1");
  Lproj8->Draw("E0E1 same");
  c8->cd(9)->SetLogy();
  ULproj9->SetMarkerColor(2);
  ULproj9->SetLineColor(2);
  Lproj9->SetMarkerColor(1);
  Lproj9->SetLineColor(1);
  ULproj9->Draw("E0E1");
  Lproj9->Draw("E0E1 same");
 
  cout << "\n\n************* Normalization of LS distribution for different pt range ***************\n\n" << endl;

  TCanvas *c12 = new TCanvas("c12","c12",10,10,1200,900);
  c12->Clear();
  c12->Divide(3,3);
  c12->cd(1)->SetLogy();

  Double_t ratio1[60];
  Double_t err_ratio1[60]; 
  Double_t mee1[60];
 
  for(Int_t i=2; i<61; i++)
    {
      if(Lproj1->GetBinContent(i)!=0)
	{
	  ratio1[i-2] = ULproj1->GetBinContent(i)/Lproj1->GetBinContent(i);
	  err_ratio1[i-2]=ratio1[i-2]*(ULproj1->GetBinError(i)/ULproj1->GetBinContent(i)+Lproj1->GetBinError(i)/Lproj1->GetBinContent(i));
	  mee1[i-2]=0.01*i;
	}
    }

  TF1 *fit1 = new TF1("fit1","pol1",0,1);
  TGraphErrors *r1 = new TGraphErrors(58,mee1,ratio1,0,err_ratio1);
  r1->SetMarkerStyle(20);
  r1->SetMarkerSize(0.4);
  r1->Fit("fit1", "  ", "  ",0.35,0.6);
  r1->Draw("AP");


  c12->cd(2)->SetLogy();

  Double_t ratio2[60];
  Double_t err_ratio2[60]; 
  Double_t mee2[60];
 
  for(Int_t i=2; i<61; i++)
    {
      if(Lproj2->GetBinContent(i)!=0)
	{
	  ratio2[i-2] = ULproj2->GetBinContent(i)/Lproj2->GetBinContent(i);
	  err_ratio2[i-2]=ratio2[i-2]*(ULproj2->GetBinError(i)/ULproj2->GetBinContent(i)+Lproj2->GetBinError(i)/Lproj2->GetBinContent(i));
	  mee2[i-2]=0.01*i;
	}
    }

  TF1 *fit2 = new TF1("fit2","pol1",0,1);
  TGraphErrors *r2 = new TGraphErrors(58,mee2,ratio2,0,err_ratio2);
  r2->SetMarkerStyle(20);
  r2->SetMarkerSize(0.4);
  r2->Fit("fit2", "  ", "  ",0.35,0.6);
  r2->Draw("AP");

 
  c12->cd(3)->SetLogy();

  Double_t ratio3[60];
  Double_t err_ratio3[60]; 
  Double_t mee3[60];
 
  for(Int_t i=2; i<61; i++)
    {
      if(Lproj3->GetBinContent(i)!=0)
	{
	  ratio3[i-2] = ULproj3->GetBinContent(i)/Lproj3->GetBinContent(i);
	  err_ratio3[i-2]=ratio3[i-2]*(ULproj3->GetBinError(i)/ULproj3->GetBinContent(i)+Lproj3->GetBinError(i)/Lproj3->GetBinContent(i));
	  mee3[i-2]=0.01*i;
	}
    }

  TF1 *fit3 = new TF1("fit3","pol1",0,1);
  TGraphErrors *r3 = new TGraphErrors(58,mee3,ratio3,0,err_ratio3);
  r3->SetMarkerStyle(20);
  r3->SetMarkerSize(0.4);
  r3->Fit("fit3", "  ", "  ",0.35,0.6);
  r3->Draw("AP");


  c12->cd(4)->SetLogy();

  Double_t ratio4[60];
  Double_t err_ratio4[60]; 
  Double_t mee4[60];
 
  for(Int_t i=2; i<61; i++)
    {
      if(Lproj4->GetBinContent(i)!=0)
	{
	  ratio4[i-2] = ULproj4->GetBinContent(i)/Lproj4->GetBinContent(i);
	  err_ratio4[i-2]=ratio4[i-2]*(ULproj4->GetBinError(i)/ULproj4->GetBinContent(i)+Lproj4->GetBinError(i)/Lproj4->GetBinContent(i));
	  mee4[i-2]=0.01*i;
	}
    }

  TF1 *fit4 = new TF1("fit4","pol1",0,1);
  TGraphErrors *r4 = new TGraphErrors(58,mee4,ratio4,0,err_ratio4);
  r4->SetMarkerStyle(20);
  r4->SetMarkerSize(0.4);
  r4->Fit("fit4", "  ", "  ",0.35,0.6);
  r4->Draw("AP");


  c12->cd(5)->SetLogy();

  Double_t ratio5[60];
  Double_t err_ratio5[60]; 
  Double_t mee5[60];
 
  for(Int_t i=2; i<61; i++)
    {
      if(Lproj5->GetBinContent(i)!=0)
	{
	  ratio5[i-2] = ULproj5->GetBinContent(i)/Lproj5->GetBinContent(i);
	  err_ratio5[i-2]=ratio5[i-2]*(ULproj5->GetBinError(i)/ULproj5->GetBinContent(i)+Lproj5->GetBinError(i)/Lproj5->GetBinContent(i));
	  mee5[i-2]=0.01*i;
	}
    }

  TF1 *fit5 = new TF1("fit5","pol1",0,1);
  TGraphErrors *r5 = new TGraphErrors(58,mee5,ratio5,0,err_ratio5);
  r5->SetMarkerStyle(20);
  r5->SetMarkerSize(0.4);
  r5->Fit("fit5", "  ", "  ",0.35,0.6);
  r5->Draw("AP");


  c12->cd(6)->SetLogy();

  Double_t ratio6[60];
  Double_t err_ratio6[60]; 
  Double_t mee6[60];
 
  for(Int_t i=2; i<61; i++)
    {
      if(Lproj6->GetBinContent(i)!=0)
	{
	  ratio6[i-2] = ULproj6->GetBinContent(i)/Lproj6->GetBinContent(i);
	  err_ratio6[i-2]=ratio6[i-2]*(ULproj6->GetBinError(i)/ULproj6->GetBinContent(i)+Lproj6->GetBinError(i)/Lproj6->GetBinContent(i));
	  mee6[i-2]=0.01*i;
	}
    }

  TF1 *fit6 = new TF1("fit6","pol1",0,1);
  TGraphErrors *r6 = new TGraphErrors(58,mee6,ratio6,0,err_ratio6);
  r6->SetMarkerStyle(20);
  r6->SetMarkerSize(0.4);
  r6->Fit("fit6", "  ", "  ",0.35,0.6);
  r6->Draw("AP");  
  
  c12->cd(7)->SetLogy();

  Double_t ratio7[60];
  Double_t err_ratio7[60]; 
  Double_t mee7[60];
 
  for(Int_t i=2; i<61; i++)
    {
      if(Lproj7->GetBinContent(i)!=0)
	{
	  ratio7[i-2] = ULproj7->GetBinContent(i)/Lproj7->GetBinContent(i);
	  err_ratio7[i-2]=ratio7[i-2]*(ULproj7->GetBinError(i)/ULproj7->GetBinContent(i)+Lproj7->GetBinError(i)/Lproj7->GetBinContent(i));
	  mee7[i-2]=0.01*i;
	}
    }

  TF1 *fit7 = new TF1("fit7","pol1",0,1);
  TGraphErrors *r7 = new TGraphErrors(58,mee7,ratio7,0,err_ratio7);
  r7->SetMarkerStyle(20);
  r7->SetMarkerSize(0.4);
  r7->Fit("fit7", "  ", "  ",0.35,0.6);
  r7->Draw("AP");

  c12->cd(8)->SetLogy();

  Double_t ratio8[60];
  Double_t err_ratio8[60]; 
  Double_t mee8[60];
 
  for(Int_t i=2; i<61; i++)
    {
      if(Lproj8->GetBinContent(i)!=0)
	{
	  ratio8[i-2] = ULproj8->GetBinContent(i)/Lproj8->GetBinContent(i);
	  err_ratio8[i-2]=ratio8[i-2]*(ULproj8->GetBinError(i)/ULproj8->GetBinContent(i)+Lproj8->GetBinError(i)/Lproj8->GetBinContent(i));
	  mee8[i-2]=0.01*i;
	}
    }

  TF1 *fit8 = new TF1("fit8","pol1",0,1);
  TGraphErrors *r8 = new TGraphErrors(58,mee8,ratio8,0,err_ratio8);
  r8->SetMarkerStyle(20);
  r8->SetMarkerSize(0.4);
  r8->Fit("fit8", "  ", "  ",0.35,0.6);
  r8->Draw("AP");


  c12->cd(9)->SetLogy();

  Double_t ratio9[60];
  Double_t err_ratio9[60]; 
  Double_t mee9[60];
 
  for(Int_t i=2; i<61; i++)
    {
      if(Lproj9->GetBinContent(i)!=0)
	{
	  ratio9[i-2] = ULproj9->GetBinContent(i)/Lproj9->GetBinContent(i);
	  err_ratio9[i-2]=ratio9[i-2]*(ULproj9->GetBinError(i)/ULproj9->GetBinContent(i)+Lproj9->GetBinError(i)/Lproj9->GetBinContent(i));
	  mee9[i-2]=0.01*i;
	}
    }

  TF1 *fit9 = new TF1("fit9","pol1",0,1);
  TGraphErrors *r9 = new TGraphErrors(58,mee9,ratio9,0,err_ratio9);
  r9->SetMarkerStyle(20);
  r9->SetMarkerSize(0.4);
  r9->Fit("fit9", "  ", "  ",0.35,0.6);
  r9->Draw("AP");

  TCanvas *c13 = new TCanvas("c13","c13",10,10,1200,900);
  c13->Clear();
 
  Double_t ratio_pt[9];
  Double_t err_ratio_pt[9];

  ratio_pt[0] = fit1->GetParameter(0);
  err_ratio_pt[0] = fit1->GetParError(0);
  ratio_pt[1] = fit2->GetParameter(0);
  err_ratio_pt[1] = fit2->GetParError(0);
  ratio_pt[2] = fit3->GetParameter(0);
  err_ratio_pt[2] = fit3->GetParError(0);
  ratio_pt[3] = fit4->GetParameter(0);
  err_ratio_pt[3] = fit4->GetParError(0);
  ratio_pt[4] = fit5->GetParameter(0);
  err_ratio_pt[4] = fit5->GetParError(0);
  ratio_pt[5] = fit6->GetParameter(0);
  err_ratio_pt[5] = fit6->GetParError(0);
  ratio_pt[6] = fit7->GetParameter(0);
  err_ratio_pt[6] = fit7->GetParError(0);
  ratio_pt[7] = fit8->GetParameter(0);
  err_ratio_pt[7] = fit8->GetParError(0);
  ratio_pt[8] = fit9->GetParameter(0);
  err_ratio_pt[8] = fit9->GetParError(0);

  Double_t pt_bin[9];
  for(Int_t i=0; i<9; i++)
    pt_bin[i]=i+1;

  TGraphErrors *r_pt = new TGraphErrors(9,pt_bin,ratio_pt,0,err_ratio_pt);
  r_pt->SetMarkerStyle(20);
  r_pt->Fit("pol1");
  r_pt->Draw("AP");

  /*-------------------------------- K0L pt distributions from MC truth ------------------------------------------*/
  /*
  TCanvas *c9 = new TCanvas("c9","c9",10,10,1200,900);
  c9->Clear();
  c9->Divide(2,3);
  c9->cd(1);  
  l->FindObject("fHistelfromK0LPttrue")->Draw();
  c9->cd(2);
  l->FindObject("fHistelfromK0LPt")->Draw();
  c9->cd(3);
  l->FindObject("fHistPifromK0LPttrue")->Draw();
  c9->cd(4);
  l->FindObject("fHistPifromK0LPt")->Draw();
  c9->cd(5);
  TH2D *helPip = (TH2D*)spelPip->Projection(0,1);
  helPip->SetStats(0);
  helPip->SetTitle("P_{T} distribution for e^{-} #Pi^{+} from K^{0}_{L} ");
  helPip->Draw("lego2z");
  c9->cd(6);
  TH2D *hposPim = (TH2D*)spposPim->Projection(0,1);
  hposPim->SetStats(0);
  hposPim->SetTitle("P_{T} distribution for e^{+} #Pi^{-} from K^{0}_{L} ");
  hposPim->Draw("lego2z");

  /*------------------------------- K0L invariant mass distributions --------------------------------------------*/
  /*
  TCanvas *c10 = new TCanvas("c10","c10",10,10,1200,900);
  c10->Clear();
  c10->Divide(2,1);
  c10->cd(1);
  l->FindObject("fHistKe3Mass_me")->Draw();
  c10->cd(2);
  l->FindObject("fHistKe3Mass_mPi")->Draw();

  /*-----------------------------------------------------------*/
  /*
  TCanvas *c11 = new TCanvas("c11","c11",10,10,1200,900);
  c11->Clear();
  l->FindObject("HistMonitor")->Draw();
  */

  /*--------------------------- Real Data--------------------------------*/

  TCanvas *ppdata = new TCanvas("ppdata", "ppdata", 10,10,1200,900); 
  ppdata->SetLogy();

  TH1F *hULdata = (TH1F*)f2->Get("ULdata");
  TH1F *hLSdata = (TH1F*)f2->Get("LSdata");
  hULdata->Sumw2();
  hLSdata->Sumw2();

  TH1F *hData =(TH1F*)hULdata->Clone();
  hData->Add(hLSdata,-1);
  hData->SetStats(0);
  hData->SetMarkerSize(0.8);
  hData->SetMarkerColor(2);
  hData->SetLineColor(2);

   TF1 *dist_data = new TF1("dist_data","Pi02+Eta2+gamma2",0,1);
   dist_data->SetLineColor(2);
   dist_data->SetParameters(n_pi0_esd2,m_pi0_esd2,n_eta_esd2,m_eta_esd2,par0_esd2,par1_esd2, par2_esd2, par3_esd2); 

   dist_data->FixParameter(1,0.13957);//m_pi0_esd);
   dist_data->FixParameter(3,0.54785);//m_eta_esd);
   //dist_data->FixParameter(0,n_pi0_esd2);
   //dist_data->FixParameter(2,n_eta_esd2);
   // dist_data->FixParameter(6,par2_esd2);
   dist_data->FixParameter(5,0.001);//par1_esd2);
   //dist_data->FixParameter(7,par3_esd2);

  TFitResultPtr fitresultsData = hData->Fit("dist_data","SWLI");
   hData->Draw("E1E0");
   
   TF1 *Pi0Data = new TF1("Pi0Data","dalitz_Pi0",0.,1.);
   TF1 *EtaData = new TF1("EtaData","dalitz_Eta",0.,1.);
   TF1 *gammaData = new TF1("gammaData","fGaussExp",0.,1);

   Pi0Data->SetLineColor(1);
   EtaData->SetLineColor(4);
   gammaData->SetLineColor(3);
   Pi0Data->SetLineWidth(2);
   EtaData->SetLineWidth(2);
   gammaData->SetLineWidth(2);
   Pi0Data->SetLineStyle(10);
   EtaData->SetLineStyle(3);
   gammaData->SetLineStyle(9);

   Pi0Data->SetParameters(dist_data->GetParameter(0),dist_data->GetParameter(1));
   EtaData->SetParameters(dist_data->GetParameter(2),dist_data->GetParameter(3));  
   gammaData->SetParameters(dist_data->GetParameter(4),dist_data->GetParameter(5),dist_data->GetParameter(6),dist_data->GetParameter(7));

   Pi0Data->Draw("same");
   EtaData->Draw("same");
   gammaData->Draw("same");
   
   cout<<"NDF = " << dist_data->GetNDF()<<endl;   

  //Integrals and errors
  TMatrixDSym covData(8);
  covData = fitresultsData->GetCovarianceMatrix();    
  covData.Print();

  //gamma
  Double_t par_gammaData[4];
  TMatrixDSym cov_gammaData(4);

  for(Int_t i = 0; i<4; i++)
    {
      par_gammaData[i]=gammaData->GetParameter(i);
      for(Int_t j = 0; j< 4; j++)
	{
	  cov_gammaData(i,j) = covData(i+4,j+4);
	  if(i==j && i==3)
	    cov_gammaData(i,j)=0;
	}
    }

  //Pi0
  Double_t par_Pi0Data[2];
  TMatrixDSym cov_Pi0Data(2);

  for(Int_t i = 0; i<2; i++)
    {
      par_Pi0Data[i]=Pi0Data->GetParameter(i);
      for(Int_t j = 0; j< 2; j++)
	{
	  cov_Pi0Data(i,j) = covData(i,j);
	}
    }

  //Eta
  Double_t par_EtaData[2];
  TMatrixDSym cov_EtaData(2);

  for(Int_t i = 0; i<2; i++)
    {
      par_EtaData[i]=EtaData->GetParameter(i);
      for(Int_t j = 0; j< 2; j++)
	{
	  cov_EtaData(i,j) = covData(i+2,j+2);
	}
    }

  Double_t gamma_areaData = gammaData->Integral(0,0.6)*h4->GetNbinsX()/0.6;
  Double_t Pi0_areaData = Pi0Data->Integral(0,0.6)*h4->GetNbinsX()/0.6;
  Double_t Eta_areaData = EtaData->Integral(0,0.6)*h4->GetNbinsX()/0.6;
  Double_t gamma_error_areaData = gammaData->IntegralError(0,0.6,par_gammaData,cov_gammaData.GetMatrixArray())*h4->GetNbinsX()/0.6;
  Double_t Pi0_error_areaData = Pi0Data->IntegralError(0,0.6,par_Pi0Data,cov_Pi0Data.GetMatrixArray())*h4->GetNbinsX()/0.6;
  Double_t Eta_error_areaData = EtaData->IntegralError(0,0.6,par_EtaData,cov_EtaData.GetMatrixArray())*h4->GetNbinsX()/0.6;

  cout << "\n\nIntegrate area of gammas = " << gamma_areaData << " +/- " << gamma_error_areaData << endl;
  cout << "Integrate area of pions = " << Pi0_areaData << " +/- " << Pi0_error_areaData <<endl;
  cout << "Integrate area of etas = " << Eta_areaData << " +/- " << Eta_error_areaData <<endl;
  cout << "Pi0/Eta ratio = " << Pi0_areaData/Eta_areaData << " +/- " << Pi0_areaData/Eta_areaData*TMath::Sqrt(TMath::Power(Pi0_error_areaData/Pi0_areaData,2)+TMath::Power(Eta_error_areaData/Eta_areaData,2)) << endl;
  cout << "gamma/Pi0 ratio = " << gamma_areaData/Pi0_areaData <<" +/- " << gamma_areaData/Pi0_areaData*TMath::Sqrt(TMath::Power(gamma_error_areaData/gamma_areaData,2)+TMath::Power(Pi0_error_areaData/Pi0_areaData,2))<< endl;
  cout << "gamma/Eta ratio = " << gamma_areaData/Eta_areaData << " +/- " << gamma_areaData/Eta_areaData*TMath::Sqrt(TMath::Power(gamma_error_areaData/gamma_areaData,2)+TMath::Power(Eta_error_areaData/Eta_areaData,2))<< endl;
  cout << "gamma/(Pi0+Eta) ratio = " << gamma_areaData/(Eta_areaData+Pi0_areaData) << " +/- " << gamma_areaData/(Eta_areaData+Pi0_areaData)*TMath::Sqrt(TMath::Power(gamma_error_areaData/gamma_areaData,2)+TMath::Power(TMath::Sqrt(TMath::Power(Pi0_error_areaData,2)+TMath::Power(Eta_error_areaData,2))/(Pi0_areaData+Eta_areaData),2)) << "\n" << endl; 

  TLatex *latexData = new TLatex();
  latexData->SetTextSize(0.04);
  latexData->SetTextAlign(33);  //align at top
  latexData->SetTextFont(42);
  latexData->DrawLatex(0.15,200000,"pp, #sqrt{s} = 7 TeV,");
  latexData->DrawLatex(0.15,100000,"N_{e^{#pm}e^{#mp}} - N_{e^{#pm}e^{#pm}}");
  
  TLegend *legendData = new TLegend(0.65,0.55,0.89,0.89);
  legendData->SetBorderSize(0);
  legendData->SetFillColor(kWhite);
  legendData->SetTextSizePixels(1000);
  legendData->SetTextFont(42);
  legendData->SetTextSize(0.045);
  legendData->AddEntry(hData,"N_{e^{#pm}e^{#mp}}","lep");
  legendData->AddEntry(dist_esd2,"all sources ","l");
  legendData->AddEntry(gammaData,"#gamma #rightarrow e^{+} e^{-}","l");
  legendData->AddEntry(Pi0Data,"#pi^{0} #rightarrow e^{+} e^{-} #gamma","l");
  legendData->AddEntry(EtaData,"#eta #rightarrow e^{+} e^{-} #gamma","l");
  legendData->Draw();

  /*--------------- OUTPUT FILES --------------------------*/

  
  MC_dist->SaveAs("Inv_mass_MC.pdf");
  esd_dist->SaveAs("Inv_mass_ESD.pdf");
  LsUs->SaveAs("like-unlike.pdf");
  LsUs_norm->SaveAs("like-unlike_norm.pdf");
  esd_comb->SaveAs("esd_comb.pdf");
  IncAss->SaveAs("incl-ass.pdf");
  LSUS_ratio->SaveAs("USLS_ratio.pdf");
  ppdata->SaveAs("ppdata.pdf");
}



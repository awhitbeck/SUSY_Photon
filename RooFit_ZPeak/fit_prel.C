using namespace RooFit;

class model{
 public:
 
    RooBreitWigner* bw;

    RooRealVar* mean;
    RooRealVar* sigma;
    RooRealVar* cbsig;
    RooRealVar* n;
    RooRealVar* alpha;
    RooCBShape* cball;
    RooGaussian* gauss;

    RooRealVar* a;
    RooRealVar* a1;
    RooRealVar* a2;
    RooRealVar* a3;
    RooRealVar* c;
    RooBernstein* poly;
    RooExponential* exp;

    RooRealVar* bwcbsig;
    RooRealVar* b;

    int nfit_params;

    RooFFTConvPdf* sig_conv;

    RooAbsPdf* background_model,*signal_model;
    RooAddPdf* sum;

    enum sig_model {kBWGauss,kBWCrystBall,kNumSigModels};
    TString sig_names[kNumSigModels]{"kBWGauss","kBWCrystBall"};
    enum bkg_model {kPoly1,kPoly2,kPoly3,kPoly4,kExpPoly1,kNumBkgModels};
    TString bkg_names[kNumBkgModels]{"kPoly1","kPoly2","kPoly3","kPoly4","kExpPoly1"};

    model(RooRealVar* x,
          RooRealVar* bwmean,
          RooRealVar* bwsigma,
          TString tag = "",
          sig_model sig_model_=kBWGauss,
          bkg_model bkg_model_=kPoly3
          ){

        nfit_params=0;

        //////////////////////////////
        // inputs for signal models //
        //////////////////////////////
        //bwmean = new RooRealVar("bwmean_"+tag, "bwmean" , 0.0, 50.,150.);
        //bwsigma = new RooRealVar("bwsigma_"+tag, "bwsigma" , 2.0, 0.01, 40.0) ;
        bw = new RooBreitWigner("bw_"+tag,"bw", *x, *bwmean, *bwsigma);

        mean = new RooRealVar("mean_"+tag, "mean" , 0.0, -0.001, 0.001) ;
        mean->setConstant(true);
        sigma = new RooRealVar("sigma_"+tag, "sigma" , 2.0, 0.01, 40.0) ;
        n = new RooRealVar("n_"+tag,"n",0.1,0.0001,10000.);
        alpha = new RooRealVar("alpha_"+tag,"alpha",1.,-10.,10.);

        cball = new RooCBShape("cball_"+tag, "crystal ball", *x, *mean, *sigma, *alpha, *n);

        gauss = new RooGaussian("gauss_"+tag,"gauss",*x,*mean,*sigma);

        //////////////////////////////////
        // inputs for background models //
        //////////////////////////////////
        a = new RooRealVar("a_"+tag, "a",0.5,0., 1.0);
        a1 = new RooRealVar("a1_"+tag,"a1",0.5,0.,1.0);
        a2 = new RooRealVar("a2_"+tag,"a2",0.5,0.,1.0);
        a3 = new RooRealVar("a3_"+tag,"a3",0.5,0.,1.0);
        c = new RooRealVar("c_"+tag,"c",1.,-1000.,1000.);
        b = new RooRealVar("background_"+tag, "background yield", 10, 0, 10000000);

        switch(bkg_model_){
        case kPoly1 : 
            std::cout << "kPoly1" << std::endl;
            poly = new RooBernstein("poly_"+tag,"poly",*x,RooArgList(*a,*a1));
            background_model = poly;
            nfit_params+=3;
            break;
        case kPoly2 : 
            std::cout << "kPoly2" << std::endl;
            poly = new RooBernstein("poly_"+tag,"poly",*x,RooArgList(*a,*a1,*a2));
            background_model = poly;
            nfit_params+=4;
            break;
        case kPoly3 : 
            std::cout << "kPoly3" << std::endl;
            poly = new RooBernstein("poly_"+tag,"poly",*x,RooArgList(*a,*a1,*a2,*a3));
            background_model = poly;
            nfit_params+=5;
            break;
        case kExpPoly1 : 
            std::cout << "kExpPoly1" << std::endl;
            poly = new RooBernstein("poly_"+tag,"poly",*x,RooArgList(*a,*a1));
            exp = new RooExponential("exp_"+tag,"exp",*poly,*c);
            background_model = exp;
            nfit_params+=4;
            break;
        default : 
            std:: cout << "Error: this option has not been implemented" << std::endl;
            assert(0);
            break;
        }

        bwcbsig = new RooRealVar("bwcbsig_"+tag, "bwcbsignal", 10, 0, 10000000);

        switch(sig_model_){
        case kBWGauss : 
            std::cout << "kBWGauss" << std::endl;
            sig_conv = new RooFFTConvPdf("sig_conv_"+tag,"sig_conv",*x,*bw,*gauss);            
            signal_model = sig_conv;
            nfit_params+=3;
            break;
        case kBWCrystBall : 
            std::cout << "kBWCrystBall" << std::endl;
            sig_conv = new RooFFTConvPdf("sig_conv_"+tag,"sig_conv",*x,*bw,*cball);            
            signal_model = sig_conv;
            nfit_params+=5;
            break;
        default : 
            std:: cout << "Error: this option has not been implemented" << std::endl;
            assert(0);
            break;
        }
        
        sum = new RooAddPdf("sum_"+tag, "", RooArgList(*signal_model, *background_model), RooArgList(*bwcbsig, *b));

    }
};

void fit_prel(bool run_data = false,
              int qmult_low_cut = 0,
              int qmult_high_cut = 10,
              model::sig_model sig_model_=model::kBWCrystBall,
              model::bkg_model bkg_model_=model::kPoly3 ){

  gROOT->ProcessLine(".L ~/tdrstyle.C");
  gROOT->ProcessLine("setTDRStyle()");

  TString base_dir = "/eos/uscms/store/user/vhegde/GMSB_skims_ST_RA2b_TreesV12/MS_EleFakes/tnpTrees/";
  TString data_file_name = "Run2016_03Feb2017_SingleElectron_tnpTree.root";
  TString mc_dy_file_name = "MS_FR_DYJetsToLL_HT-100toInf_tnpTree.root";
  TString mc_wgjets_file_name = "runList_MS_FR_WGJets*.root";
  TString mc_wjets_file_name = "runList_MS_FR_WJetsToLNu*.root";
  TString mc_ttjets_file_name = "runList_MS_FR_TTJets*.root";
  TString mc_ttgjets_file_name = "runList_MS_FR_TTGJets*.root";
  TString tree_name = "tnpTree";

  TChain* t = new TChain(tree_name);
  if( run_data ) 
      t->Add(base_dir+data_file_name);
  else{
      t->Add(base_dir+mc_dy_file_name);
      t->Add(base_dir+mc_wgjets_file_name);
      t->Add(base_dir+mc_wjets_file_name);
      t->Add(base_dir+mc_ttjets_file_name);
      t->Add(base_dir+mc_ttgjets_file_name);
  }

  RooRealVar* x = new RooRealVar("ZMass", "m_{e,e/#gamma} [GeV]", 50., 130.);
  RooRealVar* w = new RooRealVar("EvtWt", "event weight",0);
  RooRealVar* reg = new RooRealVar("ProbeIsPhoton", "0: ee region, 1 egamma region",1);
  RooRealVar* qmult = new RooRealVar("QMultProbeJet", "Probe charge multiplicity",1);
  RooRealVar* isbkg = new RooRealVar("IsBGEvent","Background truth",1);

  RooRealVar* bwmean = new RooRealVar("bwmean", "bwmean" , 90.0, 20, 180);
  RooRealVar* bwsigma = new RooRealVar("bwsigma", "bwsigma" , 4.12, 0.1, 100.0);

  RooCategory sample("sample","sample");
  sample.defineType("ele");
  sample.defineType("pho");

  cout << "///////////////////" << endl;
  cout << "// getting data  //" << endl;
  cout << "///////////////////" << endl;
  RooDataSet* data_ele;
  RooDataSet* data_pho;
  RooDataSet* data_pho_background;
  RooDataSet* data_ele_background;

  char cut_string_ele[256],cut_string_pho[256];
  sprintf(cut_string_ele,"ProbeIsPhoton==0&&QMultProbeJet>=%i&&QMultProbeJet<=%i",qmult_low_cut,qmult_high_cut);
  sprintf(cut_string_pho,"ProbeIsPhoton==1&&QMultProbeJet>=%i&&QMultProbeJet<=%i",qmult_low_cut,qmult_high_cut);
  data_ele = new RooDataSet("data_ele","data_ele",RooArgSet(*x,*reg,*w,*qmult),Import(*t),Cut(cut_string_ele),WeightVar("EvtWt"));
  data_pho = new RooDataSet("data_pho","data_pho",RooArgSet(*x,*reg,*w,*qmult),Import(*t),Cut(cut_string_pho),WeightVar("EvtWt"));

  if( ! run_data ){
      sprintf(cut_string_ele,"ProbeIsPhoton==0&&QMultProbeJet>=%i&&QMultProbeJet<=%i&&IsBGEvent==1",qmult_low_cut,qmult_high_cut);
      sprintf(cut_string_pho,"ProbeIsPhoton==1&&QMultProbeJet>=%i&&QMultProbeJet<=%i&&IsBGEvent==1",qmult_low_cut,qmult_high_cut);
      data_ele_background = new RooDataSet("data_ele_background","data_ele_background",RooArgSet(*x,*reg,*w,*qmult,*isbkg),Import(*t),Cut(cut_string_ele),WeightVar("EvtWt"));
      data_pho_background = new RooDataSet("data_pho_background","data_pho_background",RooArgSet(*x,*reg,*w,*qmult,*isbkg),Import(*t),Cut(cut_string_pho),WeightVar("EvtWt"));
  }

  RooDataSet data("data","data",RooArgSet(*x,*reg,*w,*qmult),Index(sample),Import("ele",*data_ele),Import("pho",*data_pho));
  //RooDataHist* data_binned = data.binnedClone();
  cout << "////////////////" << endl;
  cout << "// got data  //" << endl;
  cout << "///////////////" << endl;

  model ele_pdf(x,bwmean,bwsigma,"ele",sig_model_,bkg_model_);
  model pho_pdf(x,bwmean,bwsigma,"pho",sig_model_,bkg_model_);
 
  RooSimultaneous simPdf("simPdf","simultaneous pdf",sample) ;
  simPdf.addPdf(*(ele_pdf.sum),"ele") ;
  simPdf.addPdf(*(pho_pdf.sum),"pho") ;
  
  RooFitResult* r;
  if( ! run_data ){
      r = simPdf.fitTo(data,RooFit::SumW2Error(true));
  }else{
      r = simPdf.fitTo(data);
  }

  TCanvas *c1=new TCanvas("c1","",1000,500);
  c1->Divide(2,1);

  cout << "drawing pho plot" << endl;
  TPad* pad_left = (TPad*)c1->cd(1);
  RooPlot *xframe_pho=x->frame(80);
  if( ! run_data ) data_pho_background->plotOn(xframe_pho,RooFit::MarkerStyle(4));
  data_pho->plotOn(xframe_pho);
  pho_pdf.sum->plotOn(xframe_pho, LineColor(kMagenta));
  pho_pdf.sum->plotOn(xframe_pho, RooFit::Components(*(pho_pdf.poly)), RooFit::LineStyle(kDashed));
  pho_pdf.sum->plotOn(xframe_pho, Components(RooArgSet(*(pho_pdf.poly))));//,DrawOption("F"), FillColor(kGreen));
  data_pho->plotOn(xframe_pho);
  xframe_pho->GetYaxis()->SetRangeUser(0.1,xframe_pho->GetMaximum()*2.);
  xframe_pho->Draw();

  cout << "draw ele plot" << endl;

  TPad* pad_right = (TPad*)c1->cd(2);
  RooPlot *xframe_ele=x->frame(80);
  if( ! run_data ) data_ele_background->plotOn(xframe_ele,RooFit::MarkerStyle(4));
  data_ele->plotOn(xframe_ele);
  ele_pdf.sum->plotOn(xframe_ele, LineColor(kMagenta));
  ele_pdf.sum->plotOn(xframe_ele, RooFit::Components(*(ele_pdf.poly)), RooFit::LineStyle(kDashed));
  ele_pdf.sum->plotOn(xframe_ele, Components(RooArgSet(*(ele_pdf.poly))));//,DrawOption("F"), FillColor(kGreen));
  data_ele->plotOn(xframe_ele);
  xframe_ele->GetYaxis()->SetRangeUser(0.1,xframe_ele->GetMaximum()*2.);
  xframe_ele->Draw();
  
  cout << "Extracting info from fits..." << endl;

  TH1F* h_ele_pdf = (TH1F*) ele_pdf.sum->createHistogram("ele_pdf_histo",*x) ;
  TH1F* h_ele_background_pdf = (TH1F*) ele_pdf.background_model->createHistogram("ele_background_pdf_histo",*x);
  TH1F* h_data_ele = (TH1F*) data_ele->createHistogram("data_ele_histo",*x);
  h_ele_pdf->Scale(h_data_ele->Integral()/h_ele_pdf->Integral());
  h_ele_background_pdf->Scale(h_data_ele->Integral()*ele_pdf.b->getVal()/(ele_pdf.b->getVal()+ele_pdf.bwcbsig->getVal()));
  TH1F* h_data_ele_background;
  if( !run_data ){
      h_data_ele_background = (TH1F*) data_ele_background->createHistogram("data_ele_background_histo",*x);
      h_data_ele_background->Print("all");
      h_ele_background_pdf->Print("all");
      cout << "Cross check ele background - MC: " << h_data_ele_background->Integral(38,63) << " model: " << h_ele_background_pdf->Integral(38,63) << endl;
  }

  TH1F* h_pho_pdf = (TH1F*) pho_pdf.sum->createHistogram("pho_pdf_histo",*x);
  TH1F* h_data_pho = (TH1F*) data_pho->createHistogram("data_pho_histo",*x);
  TH1F* h_pho_background_pdf = (TH1F*) pho_pdf.background_model->createHistogram("pho_background_pdf_histo",*x);  
  h_pho_pdf->Scale(h_data_pho->Integral()/h_pho_pdf->Integral());
  h_pho_background_pdf->Scale(h_data_pho->Integral()*pho_pdf.b->getVal()/(pho_pdf.b->getVal()+pho_pdf.bwcbsig->getVal()));
  TH1F* h_data_pho_background;
  if( !run_data ){
      h_data_pho_background = (TH1F*) data_pho_background->createHistogram("data_pho_background_histo",*x);
      h_data_pho_background->Print("all");
      h_pho_background_pdf->Print("all");
      cout << "Cross check pho background - MC: " << h_data_pho_background->Integral(38,63) << " model: " << h_pho_background_pdf->Integral(38,63) << endl;
  }

  pad_left->SetLogy(false);
  pad_right->SetLogy(false);
  char save_file_name[256];
  if( run_data ) 
      sprintf(save_file_name,"tag_probe_fit_qmult_%s_%s_%i_%i",ele_pdf.sig_names[sig_model_].Data(),ele_pdf.bkg_names[bkg_model_].Data(),qmult_low_cut,qmult_high_cut);
  else 
      sprintf(save_file_name,"tag_probe_MC_fit_qmult_%s_%s_%i_%i",ele_pdf.sig_names[sig_model_].Data(),ele_pdf.bkg_names[bkg_model_].Data(),qmult_low_cut,qmult_high_cut);
  c1->SaveAs(TString(save_file_name)+".png");
  c1->SaveAs(TString(save_file_name)+".pdf");

  pad_left->SetLogy(true);
  pad_right->SetLogy(true);
  if( run_data )
      sprintf(save_file_name,"tag_probe_fit_qmult_%s_%s_%i_%i_LogY",ele_pdf.sig_names[sig_model_].Data(),ele_pdf.bkg_names[bkg_model_].Data(),qmult_low_cut,qmult_high_cut);
  else
      sprintf(save_file_name,"tag_probe_MC_fit_qmult_%s_%s_%i_%i_LogY",ele_pdf.sig_names[sig_model_].Data(),ele_pdf.bkg_names[bkg_model_].Data(),qmult_low_cut,qmult_high_cut);
  c1->SaveAs(TString(save_file_name)+".png");
  c1->SaveAs(TString(save_file_name)+".pdf");

  cout <<  "//////////////" << endl;
  cout << "// SUMMARY  //" << endl;
  cout << "//////////////" << endl;
  
  x->setRange("SR",80,100);
  auto ele_yield_sr = ele_pdf.sig_conv->createIntegral(RooArgSet(*x),"SR");
  auto pho_yield_sr = pho_pdf.sig_conv->createIntegral(RooArgSet(*x),"SR");
  auto ele_yield = ele_pdf.sig_conv->createIntegral(RooArgSet(*x));
  auto pho_yield = pho_pdf.sig_conv->createIntegral(RooArgSet(*x));
  cout << "pho bwcbsig: " << pho_pdf.bwcbsig->getVal() << endl;
  cout << "pho_int: " << pho_yield->getVal() << " SR: " << pho_yield_sr->getVal() << endl;
  cout << "ele bwcbsig: " << ele_pdf.bwcbsig->getVal() << endl;
  cout << "ele_int: " << ele_yield->getVal() << " SR: " << ele_yield_sr->getVal() << endl;

  double fake_rate = (pho_yield_sr->getVal()/pho_yield->getVal())*pho_pdf.bwcbsig->getVal()/(ele_yield_sr->getVal()/ele_yield->getVal()*ele_pdf.bwcbsig->getVal()) ;
  cout << ele_pdf.sig_names[sig_model_] + " --- " + ele_pdf.bkg_names[bkg_model_] << endl;
  cout << "fake rate: " << fake_rate << " +/- " << sqrt(pow(pho_pdf.bwcbsig->getError()/pho_pdf.bwcbsig->getVal(),2)+pow(ele_pdf.bwcbsig->getError()/ele_pdf.bwcbsig->getVal(),2))*fake_rate << endl;

  double my_chisquared_pho = 0.;
  double my_chisquared_ele = 0.;
  for( int i = 1 ; i <= h_pho_pdf->GetNbinsX() ; i++ ){
      // cout << "i: " << i << endl;
      // cout << " observed: " << h_data_pho->GetBinContent(i) << " +/- " << ((h_data_pho->GetBinContent(i)==0)?1.8:h_data_pho->GetBinError(i)) << endl;
      // cout << " model: " << h_pho_pdf->GetBinContent(i) << endl;
      // cout << " delta-chi-squared: " << pow(h_pho_pdf->GetBinContent(i)-h_data_pho->GetBinContent(i),2)/pow((h_data_pho->GetBinContent(i)==0)?1.8:h_data_pho->GetBinError(i),2) << endl;
      my_chisquared_pho+=pow(h_pho_pdf->GetBinContent(i)-h_data_pho->GetBinContent(i),2)/pow((h_data_pho->GetBinContent(i)==0)?1.8:h_data_pho->GetBinError(i),2);
  }
  cout << "my chi-squared (pho): " << my_chisquared_pho << " ndof: " << h_pho_pdf->GetNbinsX()-1 << " and model params: " << pho_pdf.nfit_params << endl;

  for( int i = 1 ; i <= h_ele_pdf->GetNbinsX() ; i++ ){
      // cout << "i: " << i << endl;
      // cout << " observed: " << h_data_ele->GetBinContent(i) << " +/- " << ((h_data_ele->GetBinContent(i)==0)?1.8:h_data_ele->GetBinError(i)) << endl;
      // cout << " model: " << h_ele_pdf->GetBinContent(i) << endl;
      // cout << " delta-chi-squared: " << pow(h_ele_pdf->GetBinContent(i)-h_data_ele->GetBinContent(i),2)/pow((h_data_ele->GetBinContent(i)==0)?1.8:h_data_ele->GetBinError(i),2) << endl;
      my_chisquared_ele+=pow(h_ele_pdf->GetBinContent(i)-h_data_ele->GetBinContent(i),2)/pow((h_data_ele->GetBinContent(i)==0)?1.8:h_data_ele->GetBinError(i),2);
  }
  cout << "my chi-squared (ele): " << my_chisquared_ele << " ndof: " << h_ele_pdf->GetNbinsX()-1 << " and model params: " << ele_pdf.nfit_params << endl;

}

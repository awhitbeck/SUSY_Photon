using namespace RooFit;

class model{
 public:
 
    RooRealVar* bwsig;
    RooBreitWigner* bwgauss;

    RooRealVar* cbmean;
    RooRealVar* cbsigma;
    RooRealVar* cbsig;
    RooRealVar* n;
    RooRealVar* alpha;
    RooCBShape* cball;
    
    RooRealVar* a;
    RooRealVar* a1;
    RooRealVar* a2;
    RooBernstein* poly;
    RooRealVar* b;
    
    RooRealVar* bwcbsig;
    RooFFTConvPdf* bwcbconv;
    RooAddPdf* sum;
    
    model(RooRealVar* x,
          RooRealVar* bwmean,
          RooRealVar* bwsigma,
          TString tag = ""){

        bwsig = new RooRealVar("bwsig_"+tag, "signal", 10, 0, 1000000);
        bwgauss = new RooBreitWigner("bwgauss_"+tag,"bwgauss", *x, *bwmean, *bwsigma);

        cbmean = new RooRealVar("cbmean_"+tag, "cbmean" , 0.0, -0.001, 0.001) ;
        cbsigma = new RooRealVar("cbsigma_"+tag, "cbsigma" , 2.0, 1.0, 40.0) ;
        cbsig = new RooRealVar("cbsig_"+tag, "cbsignal", 10, 0, 1000000);
        n = new RooRealVar("n_"+tag,"", 5.1);
        alpha = new RooRealVar("alpha_"+tag,"", 1.3);
        cball = new RooCBShape("cball_"+tag, "crystal ball", *x, *cbmean, *cbsigma, *alpha, *n);

        a = new RooRealVar("a_"+tag, "a",0., 1.0);
        a1 = new RooRealVar("a1_"+tag,"a1",0.,1.0);
        a2 = new RooRealVar("a2_"+tag,"a2",0.,1.0);
        poly = new RooBernstein("poly_"+tag,"poly",*x,RooArgList(*a,*a1,*a2));
        b = new RooRealVar("background_"+tag, "background yield", 10, 0, 10000000);
        
        bwcbsig = new RooRealVar("bwcbsig_"+tag, "bwcbsignal", 10, 0, 10000000);
        bwcbconv = new RooFFTConvPdf("bwcbconv_"+tag,"bwcbconv",*x,*bwgauss,*cball);
        sum = new RooAddPdf("sum_"+tag, "", RooArgList(*bwcbconv, *poly), RooArgList(*bwcbsig, *b));
    }
};

void fit_prel(bool run_data = false,
              int qmult_low_cut = 0,
              int qmult_high_cut = 10){

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

  RooRealVar* x = new RooRealVar("ZMass", "zmass", 50, 130);
  RooRealVar* w = new RooRealVar("EvtWt", "event weight",0);
  RooRealVar* reg = new RooRealVar("ProbeIsPhoton", "0: ee region, 1 egamma region",1);
  RooRealVar* qmult = new RooRealVar("QMultProbeJet", "Probe charge multiplicity",1);

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
  char cut_string_ele[256],cut_string_pho[256];
  sprintf(cut_string_ele,"ProbeIsPhoton==0&&QMultProbeJet>=%i&&QMultProbeJet<=%i",qmult_low_cut,qmult_high_cut);
  sprintf(cut_string_pho,"ProbeIsPhoton==1&&QMultProbeJet>=%i&&QMultProbeJet<=%i",qmult_low_cut,qmult_high_cut);
  if( run_data ){
      data_ele = new RooDataSet("data_ele","data_ele",RooArgSet(*x,*reg,*w,*qmult),Import(*t),Cut(cut_string_ele));
      data_pho = new RooDataSet("data_pho","data_pho",RooArgSet(*x,*reg,*w,*qmult),Import(*t),Cut(cut_string_pho));
  }else{
      data_ele = new RooDataSet("data_ele","data_ele",RooArgSet(*x,*reg,*w,*qmult),Import(*t),Cut(cut_string_ele),WeightVar("EvtWt"));
      data_pho = new RooDataSet("data_pho","data_pho",RooArgSet(*x,*reg,*w,*qmult),Import(*t),Cut(cut_string_pho),WeightVar("EvtWt"));
  }
  RooDataSet data("data","data",RooArgSet(*x,*reg,*w,*qmult),Index(sample),Import("ele",*data_ele),Import("pho",*data_pho));
  cout << "////////////////" << endl;
  cout << "// got data  //" << endl;
  cout << "///////////////" << endl;

  model ele_pdf(x,bwmean,bwsigma,"ele");
  model pho_pdf(x,bwmean,bwsigma,"pho");
 
  RooSimultaneous simPdf("simPdf","simultaneous pdf",sample) ;
  simPdf.addPdf(*(ele_pdf.sum),"ele") ;
  simPdf.addPdf(*(pho_pdf.sum),"pho") ;
  
  auto r = simPdf.fitTo(data);

  TCanvas *c1=new TCanvas("c1","",1000,500);
  c1->Divide(2,1);

  TPad* pad_left = (TPad*)c1->cd(1);
  RooPlot *xframe_pho=x->frame();
  data_pho->plotOn(xframe_pho);
  pho_pdf.sum->plotOn(xframe_pho, LineColor(kMagenta));
  pho_pdf.sum->plotOn(xframe_pho, RooFit::Components(*(pho_pdf.poly)), RooFit::LineStyle(kDashed));
  pho_pdf.sum->plotOn(xframe_pho, Components(RooArgSet(*(pho_pdf.poly))),DrawOption("F"), FillColor(kGreen));
  data_pho->plotOn(xframe_pho);
  xframe_pho->Draw();

  TPad* pad_right = (TPad*)c1->cd(2);
  RooPlot *xframe_ele=x->frame();
  data_ele->plotOn(xframe_ele);
  ele_pdf.sum->plotOn(xframe_ele, LineColor(kMagenta));
  ele_pdf.sum->plotOn(xframe_ele, RooFit::Components(*(ele_pdf.poly)), RooFit::LineStyle(kDashed));
  ele_pdf.sum->plotOn(xframe_ele, Components(RooArgSet(*(ele_pdf.poly))),DrawOption("F"), FillColor(kGreen));
  data_ele->plotOn(xframe_ele);
  xframe_ele->Draw();

  pad_left->SetLogy(false);
  pad_right->SetLogy(false);
  char save_file_name[256];
  if( run_data ) 
      sprintf(save_file_name,"tag_probe_fit_qmult_%i_%i.png",qmult_low_cut,qmult_high_cut);
  else 
      sprintf(save_file_name,"tag_probe_MC_fit_qmult_%i_%i.png",qmult_low_cut,qmult_high_cut);
  c1->SaveAs(save_file_name);

  pad_left->SetLogy(true);
  pad_right->SetLogy(true);
  if( run_data )
      sprintf(save_file_name,"tag_probe_fit_qmult_%i_%i_LogY.png",qmult_low_cut,qmult_high_cut);
  else
      sprintf(save_file_name,"tag_probe_MC_fit_qmult_%i_%i_LogY.png",qmult_low_cut,qmult_high_cut);
  c1->SaveAs(save_file_name);

  cout <<  "//////////////" << endl;
  cout << "// SUMMARY  //" << endl;
  cout << "//////////////" << endl;
  
  x->setRange("SR",80,100);
  auto ele_yield_sr = ele_pdf.bwcbconv->createIntegral(RooArgSet(*x),"SR");
  auto pho_yield_sr = pho_pdf.bwcbconv->createIntegral(RooArgSet(*x),"SR");
  auto ele_yield = ele_pdf.bwcbconv->createIntegral(RooArgSet(*x));
  auto pho_yield = pho_pdf.bwcbconv->createIntegral(RooArgSet(*x));
  cout << "pho bwcbsig: " << pho_pdf.bwcbsig->getVal() << endl;
  cout << "pho_int: " << pho_yield->getVal() << " SR: " << pho_yield_sr->getVal() << endl;
  cout << "ele bwcbsig: " << ele_pdf.bwcbsig->getVal() << endl;
  cout << "ele_int: " << ele_yield->getVal() << " SR: " << ele_yield_sr->getVal() << endl;

  cout << "fake rate: " << (pho_yield_sr->getVal()/pho_yield->getVal())*pho_pdf.bwcbsig->getVal()/(ele_yield_sr->getVal()/ele_yield->getVal()*ele_pdf.bwcbsig->getVal()) << endl;

}


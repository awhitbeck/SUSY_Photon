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
        poly = new RooBernstein("poly_"+tag,"poly",*x,RooArgList(*a,*a1,*a2));//,a3));
        b = new RooRealVar("background_"+tag, "background yield", 0, 0, 10000000);
        
        bwcbsig = new RooRealVar("bwcbsig_"+tag, "bwcbsignal", 0, 0, 10000000);
        bwcbconv = new RooFFTConvPdf("bwcbconv_"+tag,"bwcbconv",*x,*bwgauss,*cball);
        sum = new RooAddPdf("sum_"+tag, "", RooArgList(*bwcbconv, *poly), RooArgList(*bwcbsig, *b));
    }
};

void fit_prel(TString filename="MS_FR_SigBG_ST200NoZWindow.root"){

  gROOT->ProcessLine(".L ~/tdrstyle.C");
  gROOT->ProcessLine("setTDRStyle()");

  TFile * hFile = new TFile(filename);
  TH1F * DataHist_Pho = (TH1F*) hFile->Get("ZMass_Pho");
  TH1F * DataHist_Ele = (TH1F*) hFile->Get("ZMass_Ele");
  
  RooRealVar* x = new RooRealVar("x", "Mass (GeV/c^{2})", 50, 130);
  RooRealVar* bwmean = new RooRealVar("bwmean", "bwmean" , 90.0, 20, 180);
  RooRealVar* bwsigma = new RooRealVar("bwsigma", "bwsigma" , 4.12, 0.1, 100.0);

  RooCategory sample("sample","sample");
  sample.defineType("ele");
  sample.defineType("pho");

  RooDataHist data_ele("data_ele","data_ele",*x,DataHist_Ele);
  RooDataHist data_pho("data_pho","data_pho",*x,DataHist_Pho);
  RooDataHist data("data","data",*x,Index(sample),Import("ele",data_ele),Import("pho",data_pho));

  model ele_pdf(x,bwmean,bwsigma,"ele");
  model pho_pdf(x,bwmean,bwsigma,"pho");
 
  RooSimultaneous simPdf("simPdf","simultaneous pdf",sample) ;
  simPdf.addPdf(*(ele_pdf.sum),"ele") ;
  simPdf.addPdf(*(pho_pdf.sum),"pho") ;
  
  simPdf.fitTo(data);

  TCanvas *c1=new TCanvas("c1","",1000,500);
  c1->Divide(2,1);

  TPad* pad_left = (TPad*)c1->cd(1);
  RooPlot *xframe_pho=x->frame();
  data_pho.plotOn(xframe_pho);
  pho_pdf.sum->plotOn(xframe_pho, LineColor(kMagenta));
  pho_pdf.sum->plotOn(xframe_pho, RooFit::Components(*(pho_pdf.poly)), RooFit::LineStyle(kDashed));
  pho_pdf.sum->plotOn(xframe_pho, Components(RooArgSet(*(pho_pdf.poly))),DrawOption("F"), FillColor(kGreen));
  data_pho.plotOn(xframe_pho);
  xframe_pho->Draw();

  TPad* pad_right = (TPad*)c1->cd(2);
  RooPlot *xframe_ele=x->frame();
  data_ele.plotOn(xframe_ele);
  ele_pdf.sum->plotOn(xframe_ele, LineColor(kMagenta));
  ele_pdf.sum->plotOn(xframe_ele, RooFit::Components(*(ele_pdf.poly)), RooFit::LineStyle(kDashed));
  ele_pdf.sum->plotOn(xframe_ele, Components(RooArgSet(*(ele_pdf.poly))),DrawOption("F"), FillColor(kGreen));
  data_ele.plotOn(xframe_ele);
  xframe_ele->Draw();

  pad_left->SetLogy(false);
  pad_right->SetLogy(false);
  c1->SaveAs("test.png");

  pad_left->SetLogy(true);
  pad_right->SetLogy(true);
  c1->SaveAs("test_LogY.png");
}


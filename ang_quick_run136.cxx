/*
 quick and fast way to check quasi simulation outcomes
 26/06/17 Yun
*/
#include "iostream.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TSystem.h"
#include "TLine.h"
#include "TStyle.h"
#include "snprintf.h"
#include<string>

#include <TCanvas.h>
#include <TMath.h>
#include <TGraph.h>
#include "TF1.h"
#include "TLine.h"

using namespace std;

int loaded=0;

counts(){

  gStyle->SetOptStat(0);
  gROOT->ProcessLine(".L madTools_david_V9_Y.C");
  
  TFile* input1 = new TFile("YGSI3.run315.PEG.root"); // Input sorted root file

  TH2D *plot1=(TH2D *)input1->Get("Hist2D/mrthcm13_mir");// Experiment  [counts]
  plot1->RebinX(10); // binning factor; 1000/10=100
  plot1->RebinY(30); // binning factor; 1800/30=60
  
  //---------------------------------------------------------------------------------------
  //comparison regions
  TCutG *cutg1 = new TCutG("rgate",5);
   cutg1->SetPoint(0,0.15,0);
   cutg1->SetPoint(1,0.5,0);
   cutg1->SetPoint(2,0.5,180);
   cutg1->SetPoint(3,0.15,180);
   cutg1->SetPoint(4,0.15,0);
  
  TCutG *cutg2 = new TCutG("mrgate",5); //gating on 0.1 < MR < 0.9 && 15 < thetaCM < 90 
   cutg2->SetPoint(0,0.152,10);
   cutg2->SetPoint(1,0.851,10);
   cutg2->SetPoint(2,0.851,170);
   cutg2->SetPoint(3,0.152,170);
   cutg2->SetPoint(4,0.152,10);

   TH1D *px1 =plot1->ProjectionX("px1",0,-1,"[mrgate]"); //Exp
   //TH1D *px1 =plot1->ProjectionX("px1",0,-1,""); //Exp
   TH1D *py1 =plot1->ProjectionY("py1",0,-1,"[rgate]"); //

 //-----------------------------------------------------------------------------------------
 //Ouput file figures name and setting options

  TString saveFilePath1 = "Error_corrected/"; // figures saved in captures directory
  TString fileName = "54Cr208Pb_e263_run241";
  TString format = "png"; // e.g pdf, eps, jpg etc...
  TString title = "^{54}Cr + ^{208}Pb     E_{cm} = 206.06 MeV";
  TString energy = "E_{cm} = 206.06 MeV";
  //-----------------------------------------------------------------------------------------
 //Plot conditions-
  plot1->SetMinimum(1);
  plot1->SetMaximum(1000); //HERE
  px1->SetMaximum(4500);
  py1->SetMaximum(4000);

 //-----------------------------------------------------------------------------------------

  TCanvas *c1 = new TCanvas("c1", title,1,77,1585,524);
  c1->Divide(3,1);
  madColour();

 c1->cd(1);
  c1_1->SetLogz();
  plot1->SetTitle("Experiment  " + energy + "  [counts]"); //HERE
  plot1->GetXaxis()->SetTitle("M_{R}");
  plot1->GetXaxis()->CenterTitle(true);
  plot1->GetXaxis()->SetTitleOffset(0.87);
  plot1->GetYaxis()->SetTitle("#theta_{cm} [deg]");
  plot1->GetYaxis()->SetTitleOffset(1.3);
  plot1->GetYaxis()->CenterTitle(true);
  plot1->Draw("colz");
  rgate->Draw("same");
  rgate->SetLineColor(kBlue);
  rgate->SetLineWidth(3);
  // mrgate->Draw("same");
  //mrgate->SetLineColor(kBlack);
  //mrgate->SetLineWidth(2);
  tickStyle2D();

  
  c1->cd(2);
  py1->SetLineColor(kBlue+1);//Red
  py1->SetLineWidth(2);
  py1->SetMarkerColor(kBlue+1);
  py1->SetMarkerSize(1);
  py1->GetYaxis()->SetTitleOffset(1.3);
  py1->SetMarkerStyle(kFullCircle);
  py1->SetTitle("Angular Distribution [counts]"); //HERE
  py1->GetXaxis()->SetTitle("#theta_{cm} [deg]");
  py1->GetYaxis()->SetTitle("counts");
  py1->GetXaxis()->CenterTitle(true);
  py1->GetYaxis()->CenterTitle(true);
  py1->GetYaxis()->SetTitleOffset(1.5);
  py1->Draw("EP ");
  tickStyle1D();
     
 TLegend *legy = new TLegend(0.5857927,0.7509875,0.8320683,0.8456221,NULL,"brNDC");
   legy->SetBorderSize(1);
   legy->SetLineColor(0);
   legy->SetLineStyle(1);
   legy->SetLineWidth(2);
   legy->SetFillColor(0);
   legy->SetFillStyle(1001);
   TLegendEntry *entry=legy->AddEntry(py1,"Experiment","p"); //HERE
   entry->SetLineColor(4);
   entry->SetLineStyle(1);
   entry->SetLineWidth(2);
   entry->SetTextFont(42);

   legy->Draw();

   c1->cd(3);

   px1->SetLineColor(kBlue+1);
   px1->SetLineWidth(2);
   px1->SetMarkerStyle(kFullCircle);
   px1->SetMarkerColor(kBlue+1);
   px1->SetMarkerSize(1.0);
   px1->SetTitle("Mass Distribution"); //HERE
   px1->GetXaxis()->SetTitle("M_{R}");
   px1->GetYaxis()->SetTitle("counts");
   px1->GetXaxis()->CenterTitle(true);
   px1->GetYaxis()->CenterTitle(true);
   px1->GetYaxis()->SetTitleOffset(1.5);
   px1->Draw("EP");
  
   tickStyle1D();

   TLegend *legx = new TLegend(0.5857927,0.7509875,0.8320683,0.8456221,NULL,"brNDC"); 
   legx->SetLineColor(0);
   legx->SetLineStyle(1);
   legx->SetLineWidth(2);
   legx->SetFillColor(0);
   legx->SetFillStyle(1001);
   TLegendEntry *entry=legx->AddEntry(px1,"Experiment","p");//HERE
   entry->SetLineColor(4);
   entry->SetLineStyle(1);
   entry->SetLineWidth(2);
   entry->SetTextFont(42);

   legx->Draw();

    c1->SaveAs(saveFilePath1 + "exp_counts_"+ fileName + "." + format); 
   //-----------------------------------------------------------------------------------------

}
//-------------------------------------------------------------------------------
//-------------------------------------------------------------------------------
angdist(){

  gStyle->SetOptStat(0);
  gROOT->ProcessLine(".L madTools_david_V9_Y.C");
  
  TFile* input1 = new TFile("YGSI3.run119.PEG.root"); // Input sorted root file
  TFile* input2 = new TFile("quasisim_out.root"); // Input sorted root file

  TH2D *plot1=(TH2D *)input1->Get("Xsec_pub/Exp_MAD");// Experiment  [counts] already binned MAD
  TH2D *plot2=(TH2D *)input1->Get("Xsec_group/MrThCMxsec_T_mir_g");// Experiment  [dsigdth] for ang dist
  TH2D *plot4=(TH2D *)input1->Get("Xsec_pub/MrThCMxsec_T_mir");// Experiment  [dsigdthdmr] for noramlised MAD
  TH2D *plot5=(TH2D *)input2->Get("quasi_mrth3");// Simulated MAD
  TH1D *plot6=(TH1D *)input2->Get("quasi_time");// Simulated time
  TH1D *plot7=(TH1D *)input2->Get("quasi_jdist");// Simulated J
  TH2D *plot8=(TH2D *)input2->Get("quasi_Lth");// Simulated J vs Th
//  plot4->RebinY(2);// binning factor; 60/30=22
  //---------------------------------------------------------------------------------------
  //comparison regions
//  TCutG *cutg1 = new TCutG("rgate",5);
//   cutg1->SetPoint(0,0.15,0);
//   cutg1->SetPoint(1,0.5,0);
//   cutg1->SetPoint(2,0.5,180);
//   cutg1->SetPoint(3,0.15,180);
//   cutg1->SetPoint(4,0.15,0);

   TCutG *cutg1 = new TCutG("rgate",5);
   cutg1->SetPoint(0,0.26,0);
   cutg1->SetPoint(1,0.50,0);
   cutg1->SetPoint(2,0.50,180);
   cutg1->SetPoint(3,0.26,180);
   cutg1->SetPoint(4,0.26,0);
  
//  TCutG *cutg2 = new TCutG("mrgate",5); //gating on 0.1 < MR < 0.9 && 15 < thetaCM < 90 
//   cutg2->SetPoint(0,0.152,10);
//   cutg2->SetPoint(1,0.851,10);
//   cutg2->SetPoint(2,0.851,170);
//   cutg2->SetPoint(3,0.152,170);
//   cutg2->SetPoint(4,0.152,10);

   TCutG *cutg2 = new TCutG("mrgate",5); //gating on 0.1 < MR < 0.9 && 90 < thetaCM < 160 
   cutg2->SetPoint(0,0.1,90);
   cutg2->SetPoint(1,0.9,90);
   cutg2->SetPoint(2,0.9,170);
   cutg2->SetPoint(3,0.1,170);
   cutg2->SetPoint(4,0.1,90);
   
  TH1D *px1 =plot1->ProjectionX("px1",0,-1,"[mrgate]"); //Experiment in counts
  TH1D *py1 =plot1->ProjectionY("py1",0,-1,"[rgate]"); //
  
  TH1D *px2 =plot2->ProjectionX("px2",0,-1,"[mrgate]"); //Experiment in mb/rad.
  TH1D *py2 =plot2->ProjectionY("py2",0,-1,"[rgate]"); //
  TH1D *px4 =plot4->ProjectionX("px4",0,-1,"[mrgate]"); //Experiment in mb/rad.
  TH1D *py4 =plot4->ProjectionY("py4",0,-1,"[rgate]"); //

  TH1D *px5 =plot5->ProjectionX("px5",0,-1,"[mrgate]"); //Simulation.
  TH1D *py5 =plot5->ProjectionY("py5",0,-1,"[rgate]"); //Simulation.

//-----------------------------------------------------------------------------------------
 //Computing errors between histograms

  //Double_t MADerror =plot5->KolmogorovTest(plot1,"M","N")
  //Double_t MRkolmogorov =px5->KolmogorovTest(px1,"M","N")
  //Double_t MRchi2 =px5->Chi2Test(px1,"WW","CHI2/NDF")
  
  //cout << "\n MAD error is" << MADerror << "\n";
  //cout << "\n Mass Ratio K-S error is" << MRkolmogorov << "\n";
  //cout << "\n Mass Ratio Chi2 error is" << MRchi2 << "\n";

 //-----------------------------------------------------------------------------------------
 //Ouput file figures name and setting options

  TString saveFilePath1 = "50Ti+238U/"; // figures saved in captures directory
  TString fileName = "50Ti238U_e258_run119";
  TString format = "png"; // e.g pdf, eps, jpg etc...
  TString title = "^{50}Ti + ^{238}U     E_{cm} = 212.14 MeV";

  filename1 = saveFilePath1 + fileName +"_expAngDist.dat"; //only half
  filename2 = saveFilePath1 + fileName +"_expMRDist.dat"; //only half
 
  ofstream fileout1(filename1);
  ofstream fileout2(filename2);

  //-----------------------------------------------------------------------------------------
 //Plot conditions
  double scalefactor = 1.0;//1.0
  double sim_scalefactor = 0.0007;// Simulation is compared with  publication version 
  
//  plot4->SetMinimum(1);
//  plot4->SetMaximum(100); //HERE
  
  px2->SetMinimum(0);//disdth
  px2->SetMaximum(560);//
  py2->SetMinimum(0);//
  py2->SetMaximum(480);//
  
 //-----------------------------------------------------------------------------------------

  int c_bins = plot1->GetYaxis()->GetNbins(); 
//    int c_bins = plot2->GetYaxis()->GetNbins(); 
  int x_bins = plot2->GetYaxis()->GetNbins(); //normased 

  int mr_bins = plot2->GetXaxis()->GetNbins(); //normased 
  
   cout << "c_bins = " << c_bins << "\t" << "x_bins = " << x_bins << endl; // should be the same bin size

  const int n=c_bins;
  const int m=mr_bins;
  
  double thCM,fis_yield, xsec,xsec_err;
  double mr,fis_yield_m,xsec_m, xsec_err_m; 
  
  //double scalefactor = 100.;

  //-----------------------error calculation -----------------------------------
  //This information is just for calculting error bars 

  double cal_mon_eff =0.889807;
  double cal_cube_eff =0.842332;
  double cal_monsum = 165512;

  double fis_mon_eff = 0.854057;
  double fis_cube_eff = 0.815727;
  double fis_monsum = 374128;
 
  //----------------------------------------------------------------------------
  
  TCanvas *c2 = new TCanvas("c2",title,1,77,1585,524);
  TCanvas *c3 = new TCanvas("c3","c3",675,700);//2D Hist. Exp. Angle without lables
  TCanvas *c4 = new TCanvas("c4","c4",675,700);//2D Hist. Exp. Mass ratio without lables
  TCanvas *c5 = new TCanvas("c5","c5",675,700);//2D Hist. Exp. MAD without lables
  TCanvas *c6 = new TCanvas("c6","c6",675,700);//2D Hist. Sim. MAD without lables
  TCanvas *c7 = new TCanvas("c7","c7",675,700);//1D Hist. Sim. Time without lables
  TCanvas *c8 = new TCanvas("c8","c8",675,700);//1D Hist. Sim. J without lables
  TCanvas *c9 = new TCanvas("c9","c9",675,700);//2D Hist. Sim. J vs th without lables


 c3->SetRightMargin(0.003);
 c3->SetLeftMargin(0.004);
 c3->SetTopMargin(0.005);
 c3->SetBottomMargin(0.005); 

 c4->SetRightMargin(0.003);
 c4->SetLeftMargin(0.004);
 c4->SetTopMargin(0.005);
 c4->SetBottomMargin(0.005);

 c5->SetRightMargin(0.003);
 c5->SetLeftMargin(0.004);
 c5->SetTopMargin(0.005);
 c5->SetBottomMargin(0.005);//0.005

 c6->SetRightMargin(0.003);
 c6->SetLeftMargin(0.004);
 c6->SetTopMargin(0.005);
 c6->SetBottomMargin(0.005);//0.005

 c7->SetRightMargin(0.003);
 c7->SetLeftMargin(0.004);
 c7->SetTopMargin(0.005);
 c7->SetBottomMargin(0.005);//0.005

 c8->SetRightMargin(0.003);
 c8->SetLeftMargin(0.004);
 c8->SetTopMargin(0.005);
 c8->SetBottomMargin(0.005);//0.005

 c9->SetRightMargin(0.003);
 c9->SetLeftMargin(0.004);
 c9->SetTopMargin(0.005);
 c9->SetBottomMargin(0.005);//0.005
 
  c2->Divide(3,1);
  madColour();

  
  TGraphErrors *gr_ang = new TGraphErrors(); //# points,x,y,dx,dy
  TGraphErrors *gr_mr = new TGraphErrors(); //# points,x,y,dx,dy

  /////////////////////////////// //Angular distribution ///////////////////////////
  for(int i=0; i<=n;i++){

    if(py2-> GetBinContent(i)>0 && (plot1->GetYaxis()->GetBinCenter(i))>0 ){

   thCM = plot1->GetYaxis()->GetBinCenter(i);
    fis_yield = (py1->GetBinContent(i))/fis_cube_eff;

  //   thCM = plot2->GetYaxis()->GetBinCenter(i);
   // fis_yield = py2->GetBinContent(i);
    xsec = scalefactor*(py2->GetBinContent(i));//dsigdth     //READ FROM GROUP FIGURE!!!!!!!!!!!!!!!
    xsec_err= xsec*sqrt(1./fis_yield + 1./cal_monsum + 1./fis_monsum);
        
    gr_ang->SetPoint(i,thCM,xsec);
    gr_ang->SetPointError(i,0,xsec_err);

    //cout << thCM << "\t\t " << fis_yield << "\t\t " << xsec << "\t\t " << xsec_err << endl;
    fileout1 << thCM << " \t" <<  xsec << " \t" << xsec_err <<endl;
   }
    //Simulation data
    // fileout2 << py3->GetBinCenter(i) << " \t" << py3->GetBinContent(i) << " \t"<< py3->GetBinError(i) << endl;

  }
 /////////////////////////////// //Mass distribution ///////////////////////////////

  for (int j=0;j<m;j++){

    if (px1->GetBinContent(j)>0){
    mr = plot1->GetXaxis()->GetBinCenter(j);
    fis_yield_m = (px1->GetBinContent(j))/fis_cube_eff;// counts
    xsec_m = scalefactor*(px2-> GetBinContent(j));//dsigdth
    xsec_err_m= xsec_m*sqrt(1./fis_yield_m + 1./cal_monsum + 1./fis_monsum);
      
    gr_mr->SetPoint(j,mr,xsec_m);
    gr_mr->SetPointError(j,0,xsec_err_m);
    }
    //cout << mr << "\t\t " << fis_yield_m << "\t\t " << xsec_m << "\t\t " << xsec_err_m << endl;
    fileout2 << mr << "\t\t " << "\t\t " << xsec_m << "\t\t " << xsec_err_m << endl;
  
  }
  //////////////////////////////////////////////////////////////////////////////////
  
  c2->cd(1);
  c2_1->SetLogz();
  //Publication figuare
 // plot4->SetTitle("Experiment [SAnorm - d^{2}#sigma/(d#theta dM_{R})]"); //HERE
  plot4->SetTitle(""); //HERE
  plot4->GetXaxis()->SetTitle("         M_{R}    [bin width 0.01]");
  plot4->GetXaxis()->CenterTitle(true);
  plot4->GetXaxis()->SetTitleOffset(0.87);
  plot4->GetYaxis()->SetTitle("         #theta_{cm} [deg]      [bin width 3.0]");
  plot4->GetYaxis()->SetTitleOffset(1.3);
  plot4->GetXaxis()->SetTitleOffset(1.3);
  plot4->GetYaxis()->CenterTitle(true);
plot4->GetXaxis()->SetTitleSize(0.04);
plot4->GetYaxis()->SetTitleSize(0.04);
//  plot4->Scale(scalefactor);
  plot4->Draw("colz");
  //rgate->Draw("same");
  //rgate->SetLineColor(kRed);
  //rgate->SetLineWidth(5);
  //mrgate->Draw("same");
  //mrgate->SetLineColor(kBlue);
  //mrgate->SetLineWidth(5);
  tickStyle2D();


 
 
 Int_t ci;   // for color index setting
 ci = TColor::GetColor("#dcdcf0");
 


  c2->cd(2);
  //Experimental angdist (dsigdth)
//  py2->SetTitle("Angular Distribution");
   py2->SetTitle("");
  py2->GetXaxis()->SetTitle("#theta_{CM} [deg]");//ejectile angle
  py2->GetYaxis()->SetTitle("d#sigma/d#theta [mb/rad]");//recoile angle
py2->GetXaxis()->SetTitleSize(0.04);
py2->GetYaxis()->SetTitleSize(0.04);
//   py2->GetYaxis()->SetTitle("d#sigma/d#Omega [mb/rad]");//recoile angle
  py2->GetYaxis()->SetTitleOffset(1.5);
  py2->GetYaxis()->CenterTitle(true);
  py2->GetXaxis()->CenterTitle(true);
  py2->SetFillColor(ci);
//  py2->Draw("P");
  py2->Draw("EBAR");
  gr_ang->SetLineColor(kBlue+1);
  gr_ang->SetLineWidth(2);
  gr_ang->SetMarkerColor(kBlue+1);
  gr_ang->SetMarkerStyle(21);
  gr_ang->SetMarkerSize(1.0);
  py5->SetMarkerSize(1.3);
  py5->SetMarkerStyle(kFullStar);
  py5->SetMarkerColor(kRed);
    py5->Scale(sim_scalefactor);
  gr_ang->Draw("P same");
  py5->Draw("P same");
//    gr_ang->Draw("P");
  tickStyle1D();
 
  TLegend *legy = new TLegend(0.5857927,0.7509875,0.8320683,0.8456221,NULL,"brNDC"); 
  legy->SetBorderSize(1);
  legy->SetLineColor(0);
  legy->SetLineStyle(1);
  legy->SetLineWidth(2);
  legy->SetFillColor(0);
  legy->SetFillStyle(1001);
  TLegendEntry *entry=legy->AddEntry(gr_ang,"Experiment","p"); //HERE
  entry->SetLineColor(4);
  entry->SetLineStyle(1);
  entry->SetLineWidth(2);
  entry->SetTextFont(42);
  
  legy->Draw();
 
  c2->cd(3);
   //Experimental angdist (dsigdth)

//   px2->SetTitle("Mass Distribution");
   px2->SetTitle("");
   px2->GetXaxis()->SetTitle("M_{R}");//ejectile angle
   px2->GetYaxis()->SetTitle("d#sigma/dM_{R} [mb]");//recoile angle
   px2->GetYaxis()->SetTitleOffset(1.5);
   px2->GetYaxis()->CenterTitle(true);
   px2->GetXaxis()->CenterTitle(true);
   px2->GetXaxis()->SetTitleSize(0.04);
   px2->GetYaxis()->SetTitleSize(0.04);
   px2->SetFillColor(ci);
//   px2->Draw("P");
   px2->Draw("EBAR");
   gr_mr->SetLineColor(kBlue+1);
   gr_mr->SetLineWidth(2);
   gr_mr->SetMarkerColor(kBlue+1);
   gr_mr->SetMarkerStyle(21);
   gr_mr->SetMarkerSize(1.0);
     px5->SetMarkerSize(1.3);
  px5->SetMarkerStyle(kFullStar);
  px5->SetMarkerColor(kRed);
  px5->Scale(sim_scalefactor);
   gr_mr->Draw("P same");
//   gr_mr->Draw("P");
   px5->Draw("P same");
   tickStyle1D();

   TLegend *legx = new TLegend(0.5857927,0.7509875,0.8320683,0.8456221,NULL,"brNDC");
   legx->SetLineColor(0);
   legx->SetLineStyle(1);
   legx->SetLineWidth(2);
   legx->SetFillColor(0);
   legx->SetFillStyle(1001);
   TLegendEntry *entry=legx->AddEntry(gr_mr,"Experiment","p");//HERE
   entry->SetLineColor(4);
   entry->SetLineStyle(1);
   entry->SetLineWidth(2);
   entry->SetTextFont(42);

   legx->Draw();
   
   c2->SaveAs(saveFilePath1 + "exp_norm_"+ fileName + "." + format);
   c3->cd();
   py2->Scale(scalefactor);
    py2->GetXaxis()->SetLimits(0,180);
    py2->GetYaxis()->SetRangeUser(0,100);
   py2->Draw("hist bar");
   gr_ang->Draw("P, same");
   py5->Draw("P, same");
  py2->GetXaxis()->SetNdivisions(-306);
//   tickStyle1D();
TLegend *legend1 = new TLegend(0.78,0.80,0.88,0.90);

legend1->AddEntry(gr_ang," expt", "p");//ANU
legend1->AddEntry(py5," Sim ", "p");//Dub08

legend1->SetBorderSize(0);

legend1->SetTextSize(0.05);
legend1->Draw("same");

   c3->SaveAs(saveFilePath1 + "angle_"+ fileName + "." + format);

     c4->cd();
   px2->Scale(scalefactor);
   px2->GetYaxis()->SetRangeUser(0,100);
   px2->Draw("hist bar");
   gr_mr->Draw("P, same");
   px5->Draw("P, same");


TLegend *legend2 = new TLegend(0.18,0.80,0.28,0.90);

legend2->AddEntry(gr_mr," expt", "p");//ANU
legend2->AddEntry(px5," Sim ", "p");//Dub08

legend2->SetBorderSize(0);

legend2->SetTextSize(0.05);
legend2->Draw("same");
   c4->SaveAs(saveFilePath1 + "mass_"+ fileName + "." + format);

      c5->cd();
      c5->SetLogz();
   plot4->SetMinimum(1);
  plot4->SetMaximum(100000); //HERE
//  gStyle->SetLabelSize(1, "Z");
//  gStyle->SetTitleFontSize(0.5);
     plot4->GetZaxis()->SetLabelSize(0.04);
    plot4->Draw("col");
    rgate->SetLineColor(kRed);
    rgate->SetLineWidth(5);
//    rgate->Draw("same");
   c5->SaveAs(saveFilePath1 + "MAD_"+ fileName + "." + format);

      c6->cd();
      c6->SetLogz();
   plot5->SetMinimum(1);
  plot5->SetMaximum(10000000); //HERE
//    plot5->Scale(sim_scalefactor);
//  gStyle->SetLabelSize(1, "Z");
//  gStyle->SetTitleFontSize(0.5);

     plot5->GetZaxis()->SetLabelSize(0.04);
    plot5->SetTitle("");

    plot5->Draw("col");
         tickStyle2D();
        rgate->SetLineColor(kRed);
    rgate->SetLineWidth(5);
//    rgate->Draw("same");
   c6->SaveAs(saveFilePath1 + "Sim_MAD_"+ fileName + "." + format);

  c7->cd();
 plot6->GetXaxis()->SetLimits(0,50e-21);
    plot6->SetTitle("");
    plot6->GetYaxis()->SetNdivisions(505);
  plot6->Draw();
   c7->SaveAs(saveFilePath1 + "Time_"+ fileName + "." + format);

  c8->cd();
    plot7->SetTitle("");
 plot7->Draw();
 plot7->GetXaxis()->SetNdivisions(505);
  plot7->GetYaxis()->SetNdivisions(505);
 plot7->GetXaxis()->SetLimits(0,150);
 c8->SaveAs(saveFilePath1 + "Angular_Momentum_"+ fileName + "." + format);

      c9->cd();
      c9->SetLogz();
   plot8->SetMinimum(1);
  plot8->SetMaximum(10000000); //HERE
//    plot5->Scale(sim_scalefactor);
//  gStyle->SetLabelSize(1, "Z");
//  gStyle->SetTitleFontSize(0.5);

     plot8->GetZaxis()->SetLabelSize(0.04);
    plot8->SetTitle("");

    plot8->Draw("col");
         tickStyle2D();
        rgate->SetLineColor(kRed);
    rgate->SetLineWidth(5);
//    rgate->Draw("same");
   c9->SaveAs(saveFilePath1 + "Sim_JTh_"+ fileName + "." + format);
  
}
   //----------------------------------------------------------------------------------
void colour_mad()
//void madColour()
{
  const Int_t NRGBs = 20;
  const Int_t NCont = 20;
  // Colour scheme to match "Hinde special" in dagui
    Double_t stops[NRGBs] = { 0.00, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 1.00};//simulation data
    Double_t red[NRGBs]   = { 0.05, 0.05, 0.15, 0.30, 0.75, 1.00, 1.00, 1.00, 1.00, 1.00, 0.85, 0.55, 0.50, 0.50, 0.70, 1.00, 1.00, 1.00, 0.95, 1.00};
    Double_t green[NRGBs] = { 0.10, 0.35, 0.45, 0.85, 1.00, 1.00, 0.95, 0.70, 0.35, 0.00, 0.00, 0.00, 0.10, 0.10, 0.10, 0.25, 0.40, 0.60, 0.90, 1.00};
    Double_t blue[NRGBs]  = { 0.70, 1.00, 1.00, 0.50, 0.30, 0.40, 0.20, 0.00, 0.00, 0.00, 0.10, 0.25, 0.60, 0.70, 0.90, 1.00, 1.00, 1.00, 0.95, 1.00};

  TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
  gStyle->SetNumberContours(NCont);
  gStyle->SetLineWidth(2);



}

#include<stdio.h>
#include<stddef.h>
#include<stdlib.h>
#include<TFile.h>
#include<TH2.h>
#include<TMath.h>
#include<TRandom.h>
#include<TRandom2.h>
#include<TRandom3.h>

#include"ldistro_TiU_212.h" // if using generated L distro from CC-full

void quasisim(void){

  gStyle->SetOptStat(0);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1); // use 2,3,4 to get axis label as well 

  // input stuff (can be read from input files) 
  int Z1 = 22;                // Projectile
  int A1 = 50;
  
  int Z2 = 92;                // Target 
  int A2 = 238;
  
  double ELAB = 258.00;          // MeV
  double ECM = (double) A2/(A1+A2) * ELAB;
  
  //int NMAX = 10000000;
  int NMAX =1800000;          // Total number of fission events
  int MEL = 0;               // Multiplier for number of Elasics, ETOT = MEL * NMAX 100
  int JMIN = 1;               // Needs to be JMIN > 0 
  int JMAX = 79;
  int THMIN = 15.0;             // Need to have THMIN > 0 (Rutherford -> infinity at 0)
  int THMAX = 165.0; 
  bool ccfullinput = true;    // read L distro from CC-full output file

  // elastics sigma
  double M0 = 0.015;          // sigma at 0 deg 0.008
  double M180 = 0.015;        // sigma at 180 deg 0.008
  
  // MR sigma
  double MRSIGMAMIN = 0.010;   // sigma at MRSTART 0.025
  double MRSIGMAMAX = 0.1;  // sigma at MR = 0.5
  double Scale_factor_CCFULL = 1.0;

  // fission times
  bool usegausexp = true;        // Use gausexp or gausgaus
  double xd1 = 1.0 ; //0.7            // Ratio of Landau/Gaus     
  double sh_frac = 0.78; //0.4, 0.65           // Fraction of total events that are shell dependent
  double time_offset = 0.00e-21;     //0 sec; t = time_offset + F(t) where F(t) is Landau() or Gaus() dep on xd1 above 
  //double mL = 5e-21;        // sec; most prob value - Landau   
  //double sL = 20e-22;       // sec; sigma - Landau   
  double mean = 7.300e-21;      // 6.8 sec; mean  - gausgaus/gausexp   
  double sigma = 2.800e-21;       //1.35 sec; sigma - gaus (gausgaus/gausexp)   
  double lifetime = 2.800e-21;    //1.9 sec; lifetime decay - exp/sigma (gausexp/gausgaus)   
  const double m_sh = 0.280;        // Mr; mean value - mass component that is shell dependent
  const double sig_sh = 0.05;     // Mr; sigma -  mass component that is shell dependent 0.05

//    double m_sh = 0.0;        // Mr; mean value - mass component that is shell dependent
//  double sig_sh = 0.0;     // Mr; sigma -  mass component that is shell dependent
  

  double mG = 0.5e-18;         // sec; mean value - Gaus   
  double sG = 0.3e-18;       // sec; sigma - Gaus   

//    cout << "----------" << endl;
//    cout << "ZP=" << Z1 << " AP=" << A1 << " ZT=" << Z2 << " AT=" << A2 << endl;
//    cout << "ELAB=" << ELAB << " ECM=" << ECM << " NMAX=" << NMAX << " MEL="<< MEL << endl;
//    if(ccfullinput) cout << "Using L-distro from CCFULL" << endl;
//    else cout << JMIN << " < L < " << JMAX << endl;
//    cout << "MO="<< M0 << " M180=" << M180 << " MRSMIN="<< MRSIGMAMIN << " MRSMAX=" << MRSIGMAMAX << endl;
   cout<<mean<<"\t"<<sigma<<"\t"<<lifetime<<"\t"<<xd1<<endl;

  //double timeAxisMax = time_offset + 2*mL + 3*sL;  // rough estimate of maximum time ... 
  double timeAxisMax = 1e-18;  // rough estimate of maximum time ... 
  if(xd1!=1.0) timeAxisMax += 2*mG + 3*sG;         // ... and add for Gaus ... 

  // The system parameters   
  int ATOT = A1+A2;                // A for Compound Nucleus
  int ZTOT = Z1+Z2;                // Z for Compound Nucleus
  
  double MR1IN = (double) A1/ATOT; // massratio particle 1
  double MR2IN = (double) A2/ATOT; // massratio particle 2

  double MRSTART = MR1IN;          // start value of MR 
  
  // constants and converters
  double const d2r = TMath::Pi()/180.0;  
  double const r2d = 180.0/TMath::Pi();
  
  // for mass excange 
  double const tau = 5.2e-21;           // sec - from Nucl Phys A440 (1985) 327-365 ... (p361) 5.2e-21

  // angle stuff 
  double const hbar = 6.58e-16;         // planckyboy        -   eV*s
  double const mn = 931.5e6;              // nucleon mass      -   eV/c^2
  double const c = 3e8;                 // light             -   m/s
  double const massn = mn/pow(c,2);     // nucleon mass      -   eV*s^2/m^2

  double const r0 = 1.2e-15;            // radii stuff       -   m
  double const r1 = r0*pow(A1,1./3.);   // radii nucleus 1   
  double const r2 = r0*pow(A2,1./3.);   // radii nucleus 2   

  double const x1 = (r1+r2)/(1+(double)A1/A2);  // displacement nucleus 1
  double const x2 = (r1+r2)/(1+(double)A2/A1);  // displacement nucleus 2

  //  double const test = (2./5.)*pow(r1,2);
  //  cout << "test " << test <<endl;

   //Aditya: Begin edit for deformed moment of inertia calculations.
  // lengths of axes for deformed nucleus
  double const r_minor = 6.79142e-15; // calculated from beta2 
  double const r_major = 8.7269e-15; 

  // from TDHF: CM distance between nuclei at neck density of 0.08fm-3
  double const r_axis = 2.4386e-14; 
  double const r_equator = 1.5542e-14;
  double const r_z = 1.5542e-14;

  //displacements for collisions with deformed nucleus using TDHF
  double const x1_axis_TDHF = r_axis/(1+(double)A1/A2); 
  double const x2_axis_TDHF = r_axis/(1+(double)A2/A1); 
  double const x1_equator_TDHF = r_equator/(1+(double)A1/A2); 
  double const x2_equator_TDHF = r_equator/(1+(double)A2/A1); 
  double const x1_z_TDHF = r_z/(1+(double)A1/A2); 
  double const x2_z_TDHF = r_z/(1+(double)A2/A1); 
  
  //displacements for collisions with deformed nucleus using gemoetric model
  double const x1_axis = (r1+r_major)/(1+(double)A1/A2);
  double const x2_axis = (r1+r_major)/(1+(double)A2/A1); 
  double const x1_equator = (r1+r_minor)/(1+(double)A1/A2);
  double const x2_equator = (r1+r_minor)/(1+(double)A2/A1);
  double const x1_z = (r1+r_minor)/(1+(double)A1/A2);
  double const x2_z = (r1+r_minor)/(1+(double)A2/A1);
  
  double inertia = 0.0; 
  //  double const inertia = 6.852e-35; 
  double const inertia_spherical = massn * (A1*(pow(x1,2)+2./5.*pow(r1,2)) + A2*(pow(x2,2)+2./5.*pow(r2,2)));  // momentum of inertia old(spherical) correct coding
  double const inertia_old = massn * (A1*(pow(x1,2)+2/5*pow(r1,2)) + A2*(pow(x2,2)+2/5*pow(r2,2))); //old (spherical) coding that didn't give correct value

 double inertia_tdhf = 1.711196e-34;
  // moments of inertia for collisions with deformed nucleus using TDHF cm distances
  double const inertia_axis = massn * (A1*(pow(x1_axis_TDHF,2)+2./5.*pow(r1,2)) + A2*(pow(x2_axis_TDHF,2)+1./5.*(pow(r_minor,2)+pow(r_major,2)))); 
  double const inertia_equator = massn * (A1*(pow(x1_equator_TDHF,2)+2./5.*pow(r1,2)) + A2*(pow(x2_equator_TDHF,2)+1./5.*(pow(r_minor,2)+pow(r_major,2))));
  double const inertia_z = massn * (A1*(pow(x1_z_TDHF,2)+2./5.*pow(r1,2)) + A2*(pow(x2_z_TDHF,2)+1./5.*(pow(r_minor,2)+pow(r_minor,2)))); 

  // moments of inertia for collisions with deformed nucleus using geometric cm distances
  double const inertia_axis_g = massn * (A1*(pow(x1_axis,2)+2./5.*pow(r1,2)) + A2*(pow(x2_axis,2)+1./5.*(pow(r_minor,2)+pow(r_major,2)))); 
  double const inertia_equator_g = massn * (A1*(pow(x1_equator,2)+2./5.*pow(r1,2)) + A2*(pow(x2_equator,2)+1./5.*(pow(r_minor,2)+pow(r_major,2))));
  double const inertia_z_g = massn * (A1*(pow(x1_z,2)+2./5.*pow(r1,2)) + A2*(pow(x2_z,2)+1./5.*(pow(r_minor,2)+pow(r_minor,2)))); 

  cout<<"old_spherical"<<"\t"<<"spherical"<<"\t"<<"TDHF_axis"<<"\t"<<"TDHF_equator"<<"\t"<<"TDHF_z"<<"\t"<<"\t"<<"geom_axis"<<"\t"<<"geom_equator"<<"\t"<<"geom_z"<<endl; 
  cout<<inertia_old<<"\t"<<inertia<<"\t"<<inertia_axis<<"\t"<<inertia_equator<<"\t"<<inertia_z<<"\t"<<inertia_axis_g<<"\t"<<inertia_equator_g<<"\t"<<inertia_z_g<<endl;
  
  //inertia = inertia_axis_g;
	inertia = inertia_tdhf*1.0;
	cout<<"inertia_tdhf"<<endl;
	cout<< inertia<<endl;
  //Aditya: End edit for deformed moment of inertia calculations.
  
  double const red_mass = massn*A1*A2/(A1+A2);
  double theta0 = 0;

  //opening a root file to store
  TFile *fg_sim=new TFile("quasisim_out.root","RECREATE");

  // Functions used
  // mass exchange with time (for quasis)
  TF1 *f1 = new TF1("f1","([0]-0.5)*exp(-x/[1])+0.5",0,1e-19);
  f1->SetParameter(0,MRSTART);
  f1->SetParameter(1,tau);

  // mass exchange for shell dependent component (for quasis)
  TF1 *f4 = new TF1("f4","([0]- 0.280)*exp(-x/[1])+ 0.280",0,1e-19); // for 266Sg is 0.21804
  f4->SetParameter(0,MRSTART);
  f4->SetParameter(1,tau);
  
  // rutherford scattering (for elastics, no azimuthal dependence ie dOmega = 2*Pi*sin(x)*dx)
  TF1 *f2 = new TF1("f2","sin(x/180*TMath::Pi())/(pow(sin(x/180*TMath::Pi()/2),4))",0,180);

  // time distro
  if(usegausexp) {
    TF1 *f3 = new TF1("f3",gausexp,0,100e-21,3);
    f3->SetParameters(mean, sigma, lifetime);
  }else{
    TF1 *f3 = new TF1("f3",gausgaus,0,100e-21,3);
    f3->SetParameters(mean, sigma, lifetime);
  }
  cout << "Random is " << f3->GetRandom() << "\n" << endl;
  // Randomize input
  TRandom3 *rand0 = new TRandom3(0);
  TRandom3 *rand1 = new TRandom3(0);
  TRandom3 *rand2 = new TRandom3(0);
  TRandom3 *rand3 = new TRandom3(0);
  TRandom3 *rand4 = new TRandom3(0);
  
  //  double Range = RangeMean;
  //  Energy = r3->Gaus(EnergyMean,EnergySigma); // rand

  double MR = 0;
  double MRSIGMA = 0;
  //  double MRSIGMAMAX = 0.1;
  double time = 0;
  double theta = 0;
  double massratio = 0;
  
  double phiCoulomb0 = 0;
  double phiCoulomb1 = 0;
  double wt = 0; 
  double thetaTot = 0; 

  double D = 1.44*Z1*Z2/ECM;
  double B1 = 0;
  double D1 = 0;
  double TKE = 0;
  double B2 = 0;
  double D2 = 0;

  // calculate how much each J need to be normalized with to reproduce total number of events
  float cntTot = 0;
  if(ccfullinput){
    cntTot = SUMCOL2; // from generated headerfile
    JMIN = (int) COL1[1];
    JMAX = (int) COL1[NROW-1];
  }else{
    for(int J=JMIN;J<=JMAX;J++){
      cntTot += 2*J+1;                // total sum
    }
  }

  int cntNorm = (int) NMAX/cntTot;  // normalization factor

  //defining histograms
  //  TH1F *quasi_mr = new TH1F("quasi_mr","MassRatio;MassRatio;Counts",100,0,1);
    TH1F *quasi_time = new TH1F("quasi_time","Time;Time [s];Counts",1000,-1e-21,100e-21); 
    TH1F *quasi_jdist = new TH1F("quasi_jdist","J-distribution;J;Counts",JMAX-JMIN+1,JMIN-0.5,JMAX+0.5);
  TH2F *quasi_mrth1 = new TH2F("quasi_mrth1","MassRatio - Theta;MassRatio;#Theta [deg]",100,0,1,60,0,180);
  TH2F *quasi_mrth2 = new TH2F("quasi_mrth2","MassRatio - Theta;MassRatio;#Theta [deg]",100,0,1,60,0,180);
  TH2F *quasi_mrth3 = new TH2F("quasi_mrth3","MassRatio - Theta;MassRatio;#Theta [deg]",100,0,1,60,0,180);
  TH2F *quasi_Lth = new TH2F("quasi_Lth","Angular Momentum - Theta;AngularMomentum;#Theta [deg]",150,0,150,60,0,180);
  quasi_mrth1->SetOption("COLSCATZ");
  quasi_mrth2->SetOption("COLSCATZ");
  quasi_mrth3->SetOption("COLSCATZ");
  quasi_Lth->SetOption("COLSCATZ");


  // elastic stuff
  double w, MREL, MRELSIG;
  for(int i=0;i<=180;i++){
    if(i<THMIN || i>THMAX) continue;

    w = f2->Eval(i)/f2->Integral(THMIN,THMAX);
    MRELSIG = (M180-M0)/180.*i + M0; 
    cout << "Elastics, Theta = " << i << " out of " << THMAX << "\r" << flush;

    for(int j=0;j<MEL*NMAX*w;j++){
      MREL = rand0->Gaus(MR1IN,MRELSIG);
      quasi_mrth1->Fill(MREL,i);
      quasi_mrth1->Fill(1-MREL,180-i);
      quasi_mrth2->Fill(MREL,i);
      quasi_mrth2->Fill(1-MREL,180-i);
      quasi_mrth3->Fill(MREL,i);
      quasi_mrth3->Fill(1-MREL,180-i);
    }
  }

  cout << endl; 

  int halfturns = 0;
  for(int J=JMIN;J<=JMAX;J++){
 //   J = J*sqrt(Scale_factor_CCFULL);
    // first approx with theta0 - simple method
    //theta0 = asin(sqrt(J*(J+1)/(2*red_mass*Ecm*1e6))*hbar / (r1+r2));
    
    // more 'serious' approach following D. Hindes recipe  
    //B1 = 6.32*J*ATOT/(A1*A2*1.385*sqrt(ECM/A1));       // fm using A1 to be projectile nucleus
    B1 = 6.32*J*ATOT/(A1*A2*1.385*sqrt(ELAB/A1));       // fm using A1 to be projectile nucleus
    D1 = D/2*(1+sqrt(1+pow((2*B1/D),2)));              // fm 
    
    cout << "Fissions, J = " << J << " out of " << JMAX << "\r" << flush;
    //int nevt = 1000;
 
    int nevt = 0;
    if(ccfullinput){
      nevt = cntNorm*sqrt(Scale_factor_CCFULL)*COL2[J];
    }else{
      nevt = cntNorm*(2*J+1);
    }
    
    for(int i=0;i<nevt;i++){
     
      // Use flat time distro    
      // for(int i=0;i<5000;i++){
      //    time = i * 1e-23;
      // Use a distro for time - 1000 events
 
      // for(int i=0;i<nevt;i++){
      //if(i>=0  && i<(xd1*nevt)) time = time_offset+rand1->Landau(mL,sL);
      if(i>=0  && i<(xd1*nevt)) time = (time_offset+f3->GetRandom())*(exp(-0.0095*J));
      if(i>=(xd1*nevt) && i<nevt) time = (time_offset+rand2->Gaus(mG,sG))*(exp(-0.0095*J));
      
 //     if(time<0) continue;
       quasi_time->Fill(time);
      quasi_jdist->Fill(J);

      // adding third mass component that is dependent on shell closures
      // this is to simulate events that may be trapped in a 'potential pocket' 
      // a fixed number of events whos masses are dependent on shells are added to MR which is then randomised
      // and filled as massratio in the histograms
      if(i>=0  && i<(xd1*nevt) && i<(sh_frac*nevt) )
	{
	  MR = f4->Eval(time); // for mass evolving to 208Pb
	   //MR = f1->Eval(time);
	   if(MR<m_sh)
	    {
	      //MR = f1->Eval(time);
	        MR= rand4->Gaus(m_sh,sig_sh);
	    }
	  //else
	  //  {
	  //    MR= rand4->Gaus(m_sh,sig_sh); // for mass evolving independent of time after reaching m_sh with mass width sig_sh
	  //  }
	}
      else
	{
	  MR = f1->Eval(time);
	}
      // uncomment this next line to go back to old code with no shell stuff
      //      MR = f1->Eval(time);
      // end old code massratio determination

      //MRSIGMA = MRSIGMAMAX/(0.5-MR2IN)*MR - MRSIGMAMAX/(0.5-MR2IN)*MR2IN;
      MRSIGMA = (MRSIGMAMAX-MRSIGMAMIN)/(0.5-MRSTART)*MR + MRSIGMAMIN - (MRSIGMAMAX-MRSIGMAMIN)/(0.5-MRSTART)*MRSTART;

      // more 'serious' approach with theta0 continued  
      TKE = 0.785*ZTOT*ZTOT*MR*(1-MR) / (pow(ATOT,(1./3.))*(pow(MR,(1./3.))+pow((1-MR),(1./3.))));
      D2 = 1.834*pow(ATOT,(1./3.))*(pow(MR,(1./3.))+pow((1-MR),(1./3.)));                                     
      //      B2 = 5*B1*A2/ATOT*sqrt(ECM*A1) / (7*sqrt(TKE*ATOT)*MR*(1-MR)*(pow(((1-MR)/MR),0.5)+pow((MR/(1-MR)),0.5)));
      B2 = 5*B1*A2/ATOT*sqrt(ELAB*A1) / (7*sqrt(TKE*ATOT)*MR*(1-MR)*(pow(((1-MR)/MR),0.5)+pow((MR/(1-MR)),0.5)));
      
      // previous approach with halfturns
      //theta0 = atan(D1/(2*B1)) + atan(D2/(2*B2));
      
      // new approach using different theta depending on value of wt
      phiCoulomb0 = 0.5*TMath::Pi() - atan(D1/(2*B1));
      phiCoulomb1 = 0.5*TMath::Pi() - atan(D2/(2*B2));
      wt = Scale_factor_CCFULL*sqrt(J*(J+1))*hbar/inertia*time;
      thetaTot = phiCoulomb0 + wt + phiCoulomb1;
      thetaTot = fmod(thetaTot,2*TMath::Pi());    // use fmod for modulus with float
      
      if(thetaTot <= TMath::Pi()) theta = r2d * (TMath::Pi() - thetaTot); 
      if(thetaTot > TMath::Pi()) theta = r2d * (thetaTot - TMath::Pi());  
      
      //      if(i%100==0) cout << MR << "\t" << MRSIGMA << "\t" << MR2IN << endl;
      //theta = r2d * (theta0 + sqrt(J*(J+1))*hbar/inertia*time);    // simple model 
      //theta = r2d * (theta0 - sqrt(J*(J+1))*hbar/inertia*time);    // first attempt (strange when over 180 deg)
      
      //      cout << time << "\t" << r2d * thetaTot << "\t" << theta << endl; 
      
      //if(i%1000==0) cout << MR << "\t" << MRSIGMA << "\t" << theta0 << "\t"  << theta*d2r << "\t" << time ;
      

      //       halfturns = theta/180;
      //       theta = fmod(theta,180);               // fmod for modulus with float
      //       //   if (theta<0) theta = - theta;
      //       cout << "\t" << halfturns << "\t" << theta << endl;
      
      //       if(halfturns%2) theta = 180 - theta;
      
      // 	if(i==0 && cnt==0) cout << "\n" << MR1IN << "\t" << MR2IN << endl;
      //if(i==0 ) cout << "MR\t\tTKE\tD\tB1\t\tD1\tB2\t\tD2\ttheta" << endl;
      //if((i==0 || i==10 || i==20 || i==30 || i==40 || i==50) ) cout << MR << "\t" << TKE << "\t" << D << "\t" << B1 << "\t" << D1 << "\t" << B2 << "\t" << D2 << "\t" << theta << endl;
      

      if(theta<THMIN || theta>THMAX) continue;

      // HISTO FILLING STUFF FROM HERE 
      
      //  fill only for J=10
      //  if(J==10){
      //   quasi_mrth1->Fill(MR,theta);
      //   quasi_mrth1->Fill(1-MR,180-theta);
      //  }
      
      quasi_mrth1->Fill(MR,theta);
      quasi_mrth1->Fill(1-MR,180-theta);
      
      // randomize massratio for particular theta
      //    for(int j=0;j<10;j++){
      massratio = rand3->Gaus(MR,MRSIGMA); // radn
      
      // quasi_mr->Fill(massratio);
      
      
      // fill only for J=10
      if(J==10){
	quasi_mrth2->Fill(massratio,theta);
	quasi_mrth2->Fill(1-massratio,180-theta);
      }
      
      // fill sum of all
      quasi_mrth3->Fill(massratio,theta);
      quasi_mrth3->Fill(1-massratio,180-theta);
      quasi_Lth->Fill(J,theta);
//      quasi_Lth->Fill(J,180-theta);
      // }
    }
  }
  
  cout << "\nDONE!\n";

  TExec *ex1 = new TExec("ex1","colour_mad()");  
//  TExec *ex1 = new TExec("ex1","madColour()");
    
  // draw histo in separate canvases
  TCanvas *can[10];
  
//   can[1]=new TCanvas();
//   can[1]->SetWindowSize(400,400);
//   can[1]->SetWindowPosition(10,10);
//   can[1]->SetLogz(1);
//   quasi_mrth1->Draw();
  
//   can[2]=new TCanvas();
//   can[2]->SetWindowSize(400,400);
//   can[2]->SetWindowPosition(410,10);
//   can[2]->SetLogz(1);
//   quasi_mrth2->Draw();
  
  can[3]=new TCanvas();
  can[3]->SetWindowSize(600,600);
  can[3]->SetWindowPosition(10,10);
  can[3]->SetLogz(1);
  quasi_mrth3->SetMaximum(200000);

  //  gStyle->SetPalette (colour_mad(20,20));
  quasi_mrth3->GetListOfFunctions()->Add(ex1);
//  quasi_mrth3->Draw();
  quasi_mrth3->Draw("col");
  tickStyle2D();

//  quasi_Lth->Draw("col");
  tickStyle2D();
//   can[4]=new TCanvas();
//   can[4]->SetWindowSize(400,400);
//   can[4]->SetWindowPosition(10,410);
//   quasi_mr->Draw();
  
  can[5]=new TCanvas();
  can[5]->SetWindowSize(400,400);
  can[5]->SetWindowPosition(10,410);
  quasi_time->Draw();
  
  
//   can[6]=new TCanvas();
//   can[6]->SetWindowSize(400,400);
//   can[6]->SetWindowPosition(410,410);
//   can[6]->Divide(2,2);
//   can[6]->cd(1);
//   can[6]->GetPad(1)->SetLogz(1);
//   quasi_mrth1->Draw();
//   can[6]->cd(2);
//   can[6]->GetPad(2)->SetLogz(1);
//   quasi_mrth2->Draw();
//   can[6]->cd(3);
//   can[6]->GetPad(3)->SetLogz(1);
//   quasi_mrth3->Draw();
//   can[6]->cd(4);
//   quasi_time->Draw();

    can[7]=new TCanvas();
   can[7]->SetWindowSize(400,400);
   can[7]->SetWindowPosition(810,410);
   quasi_jdist->Draw();

  
  fg_sim->Write();
  // fg_sim->Close();
    
  fg_sim->IsOpen();
}



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
    // old colour scheme
    // Double_t stops[NRGBs] = { 0.00, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 1.00};
    // Double_t red[NRGBs]   = { 0.00, 0.00, 0.00, 0.00, 0.50, 1.00, 1.00, 1.00, 1.00, 1.00, 0.70, 0.60, 0.40, 0.40, 0.70, 1.00, 0.50, 0.70, 0.90, 1.00};
    //  Double_t green[NRGBs] = { 0.00, 0.00, 0.00, 0.50, 0.80, 1.00, 0.85, 0.70, 0.40, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.60, 0.80, 0.90, 1.00};
    //  Double_t blue[NRGBs]  = { 0.50, 0.70, 1.00, 0.60, 0.10, 0.00, 0.10, 0.00, 0.00, 0.00, 0.40, 0.60, 0.40, 0.70, 0.90, 1.00, 1.00, 1.00, 0.90, 1.00};
  TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
  gStyle->SetNumberContours(NCont);
  gStyle->SetLineWidth(2);
  //old scheme
  //const Int_t NRGBs = 20;
  //const Int_t NCont = 20;
  //Double_t stops[NRGBs] = { 0.00, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55,
  //			    0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 1.00};
  //Double_t red[NRGBs]   = { 0.00, 0.00, 0.00, 0.00, 0.50, 1.00, 1.00, 1.00, 1.00, 1.00, 0.70, 0.60,
  //			    0.40, 0.40, 0.70, 1.00, 0.50, 0.70, 0.90, 1.00};
  //Double_t green[NRGBs] = { 0.00, 0.00, 0.00, 0.50, 0.80, 1.00, 0.85, 0.70, 0.40, 0.00, 0.00, 0.00,
  //			    0.00, 0.00, 0.00, 0.00, 0.60, 0.80, 0.90, 1.00};
  //Double_t blue[NRGBs]  = { 0.50, 0.70, 1.00, 0.60, 0.10, 0.00, 0.10, 0.00, 0.00, 0.00, 0.40, 0.60,
  //			    0.40, 0.70, 0.90, 1.00, 1.00, 1.00, 0.90, 1.00};
  //TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
  //gStyle->SetNumberContours(NCont);






// Double_t stops[NRGBs] = { 0.00, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 1.00};//expt data
// Double_t red[NRGBs]   = {0,0,0,0,0.25,0.55,1,1,1,1,0.86,0.67,0.47,0.47,0.525,0.57,0.62,0.75,0.825,0.95};
// Double_t green[NRGBs] = {0,0,0.5,0.75,0.9,1,1,0.75,0.5,0.35,0.19,0.12,0,0.1,0.27,0.43,0.55,0.72,0.815,0.95};
// Double_t blue[NRGBs]  = {0.5,0.85,0.5,0.3,0.1,0,0,0,0,0,0,0,0.18,0.48,0.59,0.7,0.81,0.88,0.94,0.99};


}

void tickStyle2D(){

   TH2 *h2;
   TObject *obj;
   TIter next(gPad->GetListOfPrimitives());
   while ((obj = next())){
      if(obj->IsA() == TH2C::Class() || 
	 obj->IsA() == TH2S::Class() || 
	 obj->IsA() == TH2I::Class() || 
	 obj->IsA() == TH2F::Class() || 
	 obj->IsA() == TH2D::Class()){
	
	gStyle->SetLineWidth(2);
	
	h2 = (TH2*) obj;
	h2->SetNdivisions(-304,"Y");
	h2->SetNdivisions(-205,"X");
	h2->SetTickLength(0.06,"XY");
	//h2->GetXaxis()->CenterTitle(true);
	//h2->GetYaxis()->CenterTitle(true);

      }
   }
}

void colour_def(){
  gStyle->SetPalette(1);
  gStyle->SetNumberContours(20);
}

void colour_def_highres(){
  gStyle->SetPalette(1);
  gStyle->SetNumberContours(100);
}

double gausgaus(double *x, double *par){
  float xx=x[0];
  double f = 0;
  if(xx > 0 && xx < par[0]){
 //   f = 1./sqrt(2*TMath::Pi()*par[1]) * exp(-0.5*pow((xx-par[0])/par[1],2));
   f = exp(-0.5*pow((xx-par[0])/par[1],2));
  }
  else{
//    f = 1./sqrt(2*TMath::Pi()*par[1]) * exp(-0.5*pow((xx-par[0])/par[2],2));
    f = exp(-0.5*pow((xx-par[0])/par[2],2));
  }
  return f;
}

double gausexp(double *x, double *par){
  float xx=x[0];
  double f = 0;
  if(xx > 0 && xx < par[0]){
    //f = 1./sqrt(2*TMath::Pi()*par[1]) * exp(-0.5*pow((xx-par[0])/par[1],2));
    f = exp(-0.5*pow((xx-par[0])/par[1],2));
  }
  else{
    //f = 1./sqrt(2*TMath::Pi()*par[1]) * exp(par[2]*(par[0]-xx));
    f = exp(1./par[2]*(par[0]-xx));
  }
  return f;
}


void myfunc(){
TCanvas *c4 = new TCanvas("c1");
 c4->cd();
//   TF1 *f = new TF1("f",gausexp,0,100,3);
//   f->SetParameters(3,0.3,0.03);
  TF1 *f = new TF1("f",gausexp,0,100e-21,3);
  f->SetParameters(15.6e-21,1.2e-21,0.50e-21);
  f->Draw();
}

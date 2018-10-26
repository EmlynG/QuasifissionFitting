//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
//
//	Fitting code for quasisim_tip.cxx variables
//	Author: Emlyn Graham
//	Instructions: Components to change located under ********* in the code.
//	Options to change include:
//		1. Experiment: Global input files above the fitting function need to
//		   be changed
//		2. Variables: Variable starting values, step size, min and max value
//		   should be defined (min max optional, MINUIT is faster without
//		   these but may find non-physical values). Variables may also be fixed
//		   at a specific value if the user does not wish to fit them.
//
//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------


#include "iostream.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TSystem.h"
#include "TLine.h"
#include "TStyle.h"
#include "snprintf.h"
#include <string>
#include <typeinfo.h>

#include <TCanvas.h>
#include <TMath.h>
#include <TGraph.h>
#include "TF1.h"
#include "TLine.h"
#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <TFile.h>
#include <TH2.h>
#include <TMath.h>
#include <TRandom.h>
#include <TRandom2.h>
#include <TRandom3.h>
#include <TMinuit.h>

#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "TError.h"
#include <iostream>
#include "TCollection.h"

#include"ldistro_TiU_212.h" // if using generated L distro from CC-full



//-----------------------------------------------------------------------------------------
//	Needed to explicitly define sqrt in the MINUIT subroutine
//-----------------------------------------------------------------------------------------

using namespace std;

using std::sqrt;

inline int
sqrt( int in )
{
    return static_cast<int>( std::sqrt( static_cast<double>( in ) ) );
}

//-----------------------------------------------------------------------------------------
//	Function to set parameter details in Minuit fitting routine
//-----------------------------------------------------------------------------------------

void SetParameter(int number, const std::string& name, double startingval, double stepsize, double minval, double maxval, int ierflg)
{
	gMinuit->mnparm(number, name, startingval, stepsize, minval, maxval, ierflg);
}

//-----------------------------------------------------------------------------------------
//	Distribution functions
//-----------------------------------------------------------------------------------------

double gausgaus(double *x, double *par){
  double xx=x[0];
  double f = 0;
  if(xx > 0 && xx < par[0]){
   f = exp(-0.5*pow((xx-par[0])/par[1],2));
  }
  else{
    f = exp(-0.5*pow((xx-par[0])/par[2],2));
  }
  return f;
}

double gausexp(double *x, double *par){
  double xx=x[0];
  double f = 0.;
  if(xx > 0 && xx < par[0]){
    f = exp(-0.5*pow((xx-par[0])/par[1],2));
  }
  else{
    f = exp(1./par[2]*(par[0]-xx));
  }
  return f;
}


//-----------------------------------------------------------------------------------------
//	Colour scheme for MAD representation
//-----------------------------------------------------------------------------------------

void colour_mad()
{
  const Int_t NRGBs = 20;
  const Int_t NCont = 20;
  // Colour scheme to match "Hinde special" in dagui
    Double_t stops[NRGBs] = { 0.00, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45
				, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90
				, 1.00};//simulation data
    Double_t red[NRGBs]   = { 0.05, 0.05, 0.15, 0.30, 0.75, 1.00, 1.00, 1.00, 1.00, 1.00
				, 0.85, 0.55, 0.50, 0.50, 0.70, 1.00, 1.00, 1.00, 0.95
				, 1.00};
    Double_t green[NRGBs] = { 0.10, 0.35, 0.45, 0.85, 1.00, 1.00, 0.95, 0.70, 0.35, 0.00
				, 0.00, 0.00, 0.10, 0.10, 0.10, 0.25, 0.40, 0.60, 0.90
				, 1.00};
    Double_t blue[NRGBs]  = { 0.70, 1.00, 1.00, 0.50, 0.30, 0.40, 0.20, 0.00, 0.00, 0.00
				, 0.10, 0.25, 0.60, 0.70, 0.90, 1.00, 1.00, 1.00, 0.95
				, 1.00};
  TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
  gStyle->SetNumberContours(NCont);
  gStyle->SetLineWidth(2);
}


//-----------------------------------------------------------------------------------------
//	Function to measure fit between simulated and experimental MADs
//-----------------------------------------------------------------------------------------
//	Arguments:
//		avg	-	the average for the set expexp or expgaus distribution
//		INCLUDE FURTHER ARGUMENTS HERE
//-----------------------------------------------------------------------------------------

double madfitfcn(double avg, double sig, double shellfrac){
	//cout << "Average is " << avg << endl;
	//cout << "Sigma is " << sig << endl;
	//cout << "Shell fraction is " << shellfrac << endl;
	
	//------------------------------------------------------------------------------
	//	SET OF CHECKS FOR FITTING
	//------------------------------------------------------------------------------

	// double inputgaus[3] = {1.0, 2.0, 3.0};
	// double inputmean[] = {5.};
	// cout << "Output var sum: " << inputgaus[0]+inputgaus[1]+inputgaus[2] << endl;
	// cout << "Output mean: " << inputmean[0] << endl;
	// double gausout = gausexp(inputmean, inputgaus);
	// cout << "Output gaus: "<< gausout << endl;
	// CURRENT ISSUE THAT I'VE NOTICED IS THAT THE gausexp cannot be called, either return 
	// error or  
	// just 0 when it is used to populate a TH1F and then get random value
	// cout << "GausExp " << gausexp(5., inputgaus) << "n\" << endl;
	// cout << "Input mean is: " << *avg << endl;

	//------------------------------------------------------------------------------
	
	//gROOT->SetBatch(kTRUE);
	gROOT->SetBatch(0);
	gStyle->SetOptStat(0);
	gStyle->SetPadTickX(1);
	gStyle->SetPadTickY(1); // use 2,3,4 to get axis label as well 
	
	double scalefactor = 1.0;//1.0
  	double sim_scalefactor = 0.0007;// Simulation is compared with  publication version 

	//------------------------------------------------------------------------------
	//	Input stuff (can be read from input files) 
	//------------------------------------------------------------------------------

	int Z1 = 22;                // Projectile
	int A1 = 50;

	int Z2 = 92;                // Target 
	int A2 = 238;
  
	double ELAB = 258.00;          // MeV
	double ECM = (double) A2/(A1+A2) * ELAB;
	//int NMAX =1800000;          // Total number of fission events
	int NMAX =180000;          // Total number of fission events
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
	double sh_frac = shellfrac;
	//double sh_frac = 0.78; //0.4, 0.65      // Fraction of total events that are shell 
						// dependent
	double time_offset = 0.00e-21;  //0 sec; t = time_offset + F(t) where F(t) is Landau() 
					// or Gaus() dep on xd1 above 
	// Landau	
	//double mL = 5e-21;        // sec; most prob value - Landau   
	//double sL = 20e-22;       // sec; sigma - Landau   
	
	// Set mean
	//double mean = 7.300e-21;      // 6.8 sec; mean  - gausgaus/gausexp 
	double mean = avg;		// Mean controlled by fitting variable  
	//double sigma = 2.800e-21;       //1.35 sec; sigma - gaus (gausgaus/gausexp)   
	double sigma = sig;		// Sigma controlled by fitting variable
	double lifetime = 2.800e-21;    //1.9 sec; lifetime decay - exp/sigma (gausexp/gausgaus)   
	const double m_sh = 0.280;        // Mr; mean value - mass component that is shell dependent
	const double sig_sh = 0.05;     // Mr; sigma -  mass component that is shell dependent 0.05

	// double m_sh = 0.0;        // Mr; mean value - mass component that is shell dependent
	// double sig_sh = 0.0;     // Mr; sigma -  mass component that is shell dependent

	double mG = 0.5e-18;         // sec; mean value - Gaus   
	double sG = 0.3e-18;       // sec; sigma - Gaus   

	// double timeAxisMax = time_offset + 2*mL + 3*sL;  // rough estimate of maximum time ... 
	double timeAxisMax = 1e-18;  // rough estimate of maximum time ... 
	if(xd1!=1.0) timeAxisMax += 2*mG + 3*sG;         // ... and add for Gaus ... 

	//------------------------------------------------------------------------------
	//	The system parameters   
	//------------------------------------------------------------------------------

	int ATOT = A1+A2;                // A for Compound Nucleus
	int ZTOT = Z1+Z2;                // Z for Compound Nucleus
  
	double MR1IN = (double) A1/ATOT; // massratio particle 1
	double MR2IN = (double) A2/ATOT; // massratio particle 2

	double MRSTART = MR1IN;          // start value of MR 
  
	// constants and converters
	double const d2r = TMath::Pi()/180.0;  
	double const r2d = 180.0/TMath::Pi();
  
	// for mass excange 
	double const tau = 5.2e-21;      // sec - from Nucl Phys A440 (1985) 327-365 ... 
					 // (p361)	-	5.2e-21
	// angle stuff 
	double const hbar = 6.58e-16;         // planckyboy        -   eV*s
	double const mn = 931.5e6;            // nucleon mass      -   eV/c^2
	double const c = 3e8;                 // light             -   m/s
	double const massn = mn/pow(c,2);     // nucleon mass      -   eV*s^2/m^2

	double const r0 = 1.2e-15;            // radii stuff       -   m
	double const r1 = r0*pow(A1,1./3.);   // radii nucleus 1   
	double const r2 = r0*pow(A2,1./3.);   // radii nucleus 2   

	double const x1 = (r1+r2)/(1+(double)A1/A2);  // displacement nucleus 1
	double const x2 = (r1+r2)/(1+(double)A2/A1);  // displacement nucleus 2

	// Aditya: Begin edit for deformed moment of inertia calculations.
	// lengths of axes for deformed nucleus
	double const r_minor = 6.79142e-15; // calculated from beta2 
	double const r_major = 8.7269e-15; 

	// from TDHF: CM distance between nuclei at neck density of 0.08fm-3
	double const r_axis = 2.4386e-14; 
	double const r_equator = 1.5542e-14;
	double const r_z = 1.5542e-14;

	// displacements for collisions with deformed nucleus using TDHF
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

	// momentum of inertia old (spherical) correct coding	
	double const inertia_spherical = massn * (A1*(pow(x1,2)+2./5.*pow(r1,2))
					  + A2*(pow(x2,2)+2./5.*pow(r2,2)));
	
	//old (spherical) coding that didn't give correct value
	double const inertia_old = massn * (A1*(pow(x1,2)+2/5*pow(r1,2))
				    + A2*(pow(x2,2)+2/5*pow(r2,2)));


	double inertia_tdhf = 1.711196e-34;
	// moments of inertia for collisions with deformed nucleus using TDHF cm distances
	double const inertia_axis = massn * (A1*(pow(x1_axis_TDHF,2)+2./5.*pow(r1,2))
				    + A2*(pow(x2_axis_TDHF,2)+1./5.*(pow(r_minor,2)
				    + pow(r_major,2))));
 
	double const inertia_equator = massn * (A1*(pow(x1_equator_TDHF,2)+2./5.*pow(r1,2))
					 + A2*(pow(x2_equator_TDHF,2)+1./5.*(pow(r_minor,2)
					 + pow(r_major,2))));

	double const inertia_z = massn * (A1*(pow(x1_z_TDHF,2)+2./5.*pow(r1,2))
				 + A2*(pow(x2_z_TDHF,2)+1./5.*(pow(r_minor,2)
				 + pow(r_minor,2)))); 

	// Moments of inertia for collisions with deformed nucleus using geometric cm distances
	double const inertia_axis_g = massn * (A1*(pow(x1_axis,2)+2./5.*pow(r1,2))
					+ A2*(pow(x2_axis,2)+1./5.*(pow(r_minor,2)
					+ pow(r_major,2)))); 

	double const inertia_equator_g = massn * (A1*(pow(x1_equator,2)+2./5.*pow(r1,2))
					 + A2*(pow(x2_equator,2)+1./5.*(pow(r_minor,2)
					 + pow(r_major,2))));

	double const inertia_z_g = massn * (A1*(pow(x1_z,2)+2./5.*pow(r1,2))
					+ A2*(pow(x2_z,2)+1./5.*(pow(r_minor,2)
					+ pow(r_minor,2)))); 

	// Choice of MoI

	// inertia = inertia_axis_g;
	inertia = inertia_tdhf*1.0;

	// Aditya: End edit for deformed moment of inertia calculations.
  
	double const red_mass = massn*A1*A2/(A1+A2);
	double theta0 = 0;

	//------------------------------------------------------------------------------
	//	Functions used
	//------------------------------------------------------------------------------

	// Mass exchange with time (for quasis)
	TF1 *f1 = new TF1("f1","([0]-0.5)*exp(-x/[1])+0.5",0,1e-19);
	f1->SetParameter(0,MRSTART);
	f1->SetParameter(1,tau);

	// Mass exchange for shell dependent component (for quasis)
	TF1 *f4 = new TF1("f4","([0]- 0.280)*exp(-x/[1])+ 0.280",0,1e-19); // for 266Sg is 0.21804
	f4->SetParameter(0,MRSTART);
	f4->SetParameter(1,tau);
  
	// Rutherford scattering (for elastics, no azimuthal dependence ie dOmega = 2*Pi*sin(x)*dx)
	TF1 *f2 = new TF1("f2","sin(x/180*TMath::Pi())/(pow(sin(x/180*TMath::Pi()/2),4))",0,180);

	// Time distro
	if(usegausexp) {
			TF1 *f3 = new TF1("f3",gausexp,0,100e-21,3);
			f3->SetParameters(mean, sigma, lifetime);
			f3->SetParNames("mean", "sigma", "lifetime");
	}
	else{
		TF1 *f3 = new TF1("f3",gausgaus,0,100e-21,3);
		f3->SetParameters(mean, sigma, lifetime);
		f3->SetParNames("mean", "sigma", "lifetime");
	}
	
	//cout << "Random val is " << f3->GetRandom() << endl;

	// Randomize input
	TRandom3 *rand0 = new TRandom3(0);
	TRandom3 *rand1 = new TRandom3(0);
	TRandom3 *rand2 = new TRandom3(0);
	TRandom3 *rand3 = new TRandom3(0);
	TRandom3 *rand4 = new TRandom3(0);
  
	// Double Range = RangeMean;
	// Energy = r3->Gaus(EnergyMean,EnergySigma); // rand

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

	//---------------------------------------------------------------------------------
	//	Calculate how much each J need to be normalized with to reproduce total 
	//	number of events
	//---------------------------------------------------------------------------------
	float cntTot = 0;
	if(ccfullinput){
			cntTot = SUMCOL2; // from generated headerfile
			JMIN = (int) COL1[1];
			JMAX = (int) COL1[NROW-1];
	}
	else{
		for(int J=JMIN;J<=JMAX;J++){
		cntTot += 2*J+1;                // total sum
	}
	}

	int cntNorm = (int) NMAX/cntTot;  // normalization factor

	// Defining histograms
	// TH1F *quasi_mr = new TH1F("quasi_mr","MassRatio;MassRatio;Counts",100,0,1);
	TH1F *quasi_time = new TH1F("quasi_time","Time;Time [s];Counts",1000,-1e-21,100e-21); 
	TH1F *quasi_jdist = new TH1F("quasi_jdist","J-distribution;J;Counts",JMAX-JMIN+1
				,JMIN-0.5,JMAX+0.5);

	TH2F *quasi_mrth1 = new TH2F("quasi_mrth1","MassRatio - Theta;MassRatio;#Theta [deg]"
					,100,0,1,60,0,180);

	TH2F *quasi_mrth2 = new TH2F("quasi_mrth2","MassRatio - Theta;MassRatio;#Theta [deg]"
					,100,0,1,60,0,180);

	TH2F *quasi_mrth3 = new TH2F("quasi_mrth3","MassRatio - Theta;MassRatio;#Theta [deg]"
					,100,0,1,60,0,180);
	
	TH2F *quasi_Lth = new TH2F("quasi_Lth","Angular Momentum - Theta;AngularMomentum;#Theta [deg]"
					,150,0,150,60,0,180);
	
	quasi_mrth1->SetOption("COLSCATZ");
	quasi_mrth2->SetOption("COLSCATZ");
	quasi_mrth3->SetOption("COLSCATZ");
	quasi_Lth->SetOption("COLSCATZ");


	// Elastic stuff
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
		// J = J*sqrt(Scale_factor_CCFULL);

		// First approx with theta0 - simple method

		// theta0 = asin(sqrt(J*(J+1)/(2*red_mass*Ecm*1e6))*hbar / (r1+r2));
    
		// more 'serious' approach following D. Hindes recipe  
		// B1 = 6.32*J*ATOT/(A1*A2*1.385*sqrt(ECM/A1)); 	// fm using A1 to be 
									// projectile nucleus
		B1 = 6.32*J*ATOT/(A1*A2*1.385*sqrt(ELAB/A1));   // fm using A1 to be 
								// projectile nucleus
		D1 = D/2*(1+sqrt(1+pow((2*B1/D),2)));           // fm 
    
		cout << "Fissions, J = " << J << " out of " << JMAX << "\r" << flush;
		// int nevt = 1000;
 
		int nevt = 0;
		if(ccfullinput){
			nevt = cntNorm*sqrt(Scale_factor_CCFULL)*COL2[J];
		}
		else{
			nevt = cntNorm*(2*J+1);
		}
    
		for(int i=0;i<nevt;i++){
     
			// Use flat time distro    
			
			// for(int i=0;i<5000;i++){
				// time = i * 1e-23;
				
				// Use a distro for time - 1000 events
 
				// for(int i=0;i<nevt;i++){
					// if(i>=0  && i<(xd1*nevt)) time = time_offset+rand1->
					//					Landau(mL,sL);
      					if(i>=0  && i<(xd1*nevt)) time =
					 	(time_offset+f3->GetRandom())*(exp(-0.0095*J));
					if(i>=(xd1*nevt) && i<nevt) time =
						(time_offset+rand2->Gaus(mG,sG))*(exp(-0.0095*J));
      
					// if(time<0) continue;
					quasi_time->Fill(time);
					quasi_jdist->Fill(J);

		// adding third mass component that is dependent on shell closures
		// this is to simulate events that may be trapped in a 'potential pocket' 
		// a fixed number of events whos masses are dependent on shells are added
		// to MR which is then randomised and filled as massratio in the histograms
					if(i>=0  && i<(xd1*nevt) && i<(sh_frac*nevt) ){
						MR = f4->Eval(time); // for mass evolving to 208Pb
					//	MR = f1->Eval(time);
						if(MR<m_sh){
					//		MR = f1->Eval(time);
							MR = rand4->Gaus(m_sh,sig_sh);
						}
					//	else{
					//		MR= rand4->Gaus(m_sh,sig_sh);
						// for mass evolving independent of time
						// after reaching m_sh with mass width sig_sh
					//	}
					}
					else{
						MR = f1->Eval(time);
					}
		}
					// Uncomment this next line to go back to old code with no 
					// shell stuff
					// MR = f1->Eval(time);
					// end old code massratio determination

		// MRSIGMA = MRSIGMAMAX/(0.5-MR2IN)*MR - MRSIGMAMAX/(0.5-MR2IN)*MR2IN;
		MRSIGMA = (MRSIGMAMAX-MRSIGMAMIN)/(0.5-MRSTART)*MR + MRSIGMAMIN
				- (MRSIGMAMAX-MRSIGMAMIN)/(0.5-MRSTART)*MRSTART;

		// more 'serious' approach with theta0 continued  
		TKE = 0.785*ZTOT*ZTOT*MR*(1-MR) / (pow(ATOT,(1./3.))*(pow(MR,(1./3.))
			+ pow((1-MR),(1./3.))));
		D2 = 1.834*pow(ATOT,(1./3.))*(pow(MR,(1./3.))
			+ pow((1-MR),(1./3.)));                                     
		// B2 = 5*B1*A2/ATOT*sqrt(ECM*A1) / (7*sqrt(TKE*ATOT)*MR*(1-MR)
		//	 *(pow(((1-MR)/MR),0.5) + pow((MR/(1-MR)),0.5)));
		B2 = 5*B1*A2/ATOT*sqrt(ELAB*A1) / (7*sqrt(TKE*ATOT)*MR*(1-MR)*(pow(((1-MR)/MR),0.5)
			+ pow((MR/(1-MR)),0.5)));
      
		// Previous approach with halfturns
		// theta0 = atan(D1/(2*B1)) + atan(D2/(2*B2));
      
		// New approach using different theta depending on value of wt
		phiCoulomb0 = 0.5*TMath::Pi() - atan(D1/(2*B1));
		phiCoulomb1 = 0.5*TMath::Pi() - atan(D2/(2*B2));
		wt = Scale_factor_CCFULL*sqrt(J*(J+1))*hbar/inertia*time;
		thetaTot = phiCoulomb0 + wt + phiCoulomb1;
		thetaTot = fmod(thetaTot,2*TMath::Pi());    // use fmod for modulus with float
      
		if(thetaTot <= TMath::Pi()) theta = r2d * (TMath::Pi() - thetaTot); 
		if(thetaTot > TMath::Pi()) theta = r2d * (thetaTot - TMath::Pi());  
      
		// if(i%100==0) cout << MR << "\t" << MRSIGMA << "\t" << MR2IN << endl;
		// theta = r2d * (theta0 + sqrt(J*(J+1))*hbar/inertia*time);    // simple model 
		// theta = r2d * (theta0 - sqrt(J*(J+1))*hbar/inertia*time);    // first attempt
								// (strange when over 180 deg)


		// halfturns = theta/180;
		// theta = fmod(theta,180);               // fmod for modulus with float
		// if (theta<0) theta = - theta;
		// if(halfturns%2) theta = 180 - theta;
	      
		// if(i==0 && cnt==0) cout << "\n" << MR1IN << "\t" << MR2IN << endl;
		// if(i==0 ) cout << "MR\t\tTKE\tD\tB1\t\tD1\tB2\t\tD2\ttheta" << endl;
		// if((i==0 || i==10 || i==20 || i==30 || i==40 || i==50) )
		//	cout << MR << "\t" << TKE << "\t" << D << "\t" << B1 <<
		//	"\t" << D1 << "\t" << B2 << "\t" << D2 << "\t" << theta << endl;
	      

		if(theta<THMIN || theta>THMAX) continue;

		//---------------------------------------------------------------------------------
		//	HISTO FILLING STUFF FROM HERE 
		//---------------------------------------------------------------------------------

		// Fill only for J=10
		// if(J==10){
		// quasi_mrth1->Fill(MR,theta);
		// quasi_mrth1->Fill(1-MR,180-theta);
		// }
	      
		quasi_mrth1->Fill(MR,theta);
		quasi_mrth1->Fill(1-MR,180-theta);
	      
		// Randomise massratio for particular theta
		// for(int j=0;j<10;j++){
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
		// quasi_Lth->Fill(J,180-theta);
		// }
		// }
	}
  
	cout << "\nDONE!\n";

	TExec *ex1 = new TExec("ex1","colour_mad()");
	quasi_mrth3->SetMaximum(200000);
	quasi_mrth3->GetListOfFunctions()->Add(ex1);


	gROOT->FindObject("plot1");
	gROOT->FindObject("plot2");
	gROOT->FindObject("plot4");
	TH2D *plot5 = quasi_mrth3->Clone();
  
	gStyle->SetOptStat(0);
	gROOT->ProcessLine(".L madTools_david_V9_Y.C");
	TString title = "^{50}Ti + ^{238}U     E_{cm} = 212.14 MeV";
	TString saveFilePath1 = "50Ti+238U/"; // figures saved in captures directory
	TString fileName = "50Ti238U_e258_run119";
	TString format = "png"; // e.g pdf, eps, jpg etc...
	filename1 = saveFilePath1 + fileName +"_expAngDist.dat"; //only half
	filename2 = saveFilePath1 + fileName +"_expMRDist.dat"; //only half
 
	ofstream fileout1(filename1);
	ofstream fileout2(filename2);

	//---------------------------------------------------------------------------------	
	//	Comparison regions
	//---------------------------------------------------------------------------------	

	TCutG *cutg1 = new TCutG("rgate",5);
	cutg1->SetPoint(0,0.26,0);
	cutg1->SetPoint(1,0.50,0);
	cutg1->SetPoint(2,0.50,180);
	cutg1->SetPoint(3,0.26,180);
	cutg1->SetPoint(4,0.26,0);

	TCutG *cutg2 = new TCutG("mrgate",5); //gating on 0.1 < MR < 0.9 && 90 < thetaCM < 160 
	cutg2->SetPoint(0,0.1,90);
	cutg2->SetPoint(1,0.9,90);
	cutg2->SetPoint(2,0.9,170);
	cutg2->SetPoint(3,0.1,170);
	cutg2->SetPoint(4,0.1,90);
   
	TH1D *px1 =plot1->ProjectionX("px1",0,-1,"[mrgate]"); // Experiment in counts
	TH1D *py1 =plot1->ProjectionY("py1",0,-1,"[rgate]");  //
	
	TH1D *px2 =plot2->ProjectionX("px2",0,-1,"[mrgate]"); // Experiment in mb/rad.
	TH1D *py2 =plot2->ProjectionY("py2",0,-1,"[rgate]");  //
	TH1D *px4 =plot4->ProjectionX("px4",0,-1,"[mrgate]"); // Experiment in mb/rad.
	TH1D *py4 =plot4->ProjectionY("py4",0,-1,"[rgate]");  //

	TH1D *px5 =plot5->ProjectionX("px5",0,-1,"[mrgate]"); // Simulation.
	TH1D *py5 =plot5->ProjectionY("py5",0,-1,"[rgate]");  // Simulation.


	px2->SetMinimum(0);//disdth
	px2->SetMaximum(560);//
	py2->SetMinimum(0);//
	py2->SetMaximum(480);//


	//py5->Scale(sim_scalefactor);
	//px5->Scale(sim_scalefactor);
	//py2->Scale(scalefactor);
	//px2->Scale(scalefactor);
	//plot5->Scale(sim_scalefactor*sim_scalefactor);
	//plot1->Scale(scalefactor*scalefactor);
	//py5->Scale(1.);
	//px5->Scale(1.);
	//py2->Scale(1.);
	//px2->Scale(1.);
	//plot5->Scale(1.);
	//plot1->Scale(1.);
	

	int c_bins = plot1->GetYaxis()->GetNbins(); 
	int x_bins = plot2->GetYaxis()->GetNbins(); //normased 
	int mr_bins = plot2->GetXaxis()->GetNbins(); //normased 
	const int n=c_bins;
	const int m=mr_bins;
  
	double thCM,fis_yield, xsec,xsec_err;
	double mr,fis_yield_m,xsec_m, xsec_err_m; 

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

			xsec = scalefactor*(py2->GetBinContent(i));//dsigdth     //READ FROM GROUP FIGURE!!!!!!!!!!!!!!!
			xsec_err= xsec*sqrt(1./fis_yield + 1./cal_monsum + 1./fis_monsum);
        
			gr_ang->SetPoint(i,thCM,xsec);
			gr_ang->SetPointError(i,0,xsec_err);

			fileout1 << thCM << " \t" <<  xsec << " \t" << xsec_err <<endl;
		}
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
	plot4->Draw("colz");
	rgate->Draw("same");
	rgate->SetLineColor(kBlue);
	rgate->SetLineWidth(3);
	tickStyle2D();

 
	Int_t ci;   // for color index setting
	ci = TColor::GetColor("#dcdcf0");
 


	c2->cd(2);
	//Experimental angdist (dsigdth)
	//  py2->SetTitle("Angular Distribution");
	py2->SetTitle("");
	py2->GetXaxis()->SetTitle("#theta_{CM} [deg]");//ejectile angle
	py2->GetYaxis()->SetTitle("d#sigma/d#theta [mb/rad]");//recoile angle
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
   
	//c2->SaveAs(saveFilePath1 + "exp_norm_"+ fileName + "." + format);
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

	//c3->SaveAs(saveFilePath1 + "angle_"+ fileName + "." + format);
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
	//c4->SaveAs(saveFilePath1 + "mass_"+ fileName + "." + format);


	// Scale simulation to match experimental counts

	//---------------------------------------------------------------------------------	
	//	Computing errors between histograms
	//---------------------------------------------------------------------------------	

	// Define various measures of fit	
	// NEED TO FIGURE OUT WHY ONLY THE KStot MEASURE WORKS PROPERLY
	//Double_t MADerror = plot5->KolmogorovTest(plot1,"M WW");
	Double_t ChiNDFError = plot5->Chi2Test(plot1, "UU CHI2/NDF");
	Double_t ChiError = plot5->Chi2Test(plot1, "UU CHI2");
	Double_t Xchi2ndf = px5->Chi2Test(px1,"UU CHI2/NDF");
	Double_t Ychi2ndf = py5->Chi2Test(py1,"UU CHI2/NDF");
	Double_t Xchi2 = px5->Chi2Test(px1,"UU CHI2");
	Double_t Ychi2 = py5->Chi2Test(py1,"UU CHI2");
	//Double_t XKS = px5->KolmogorovTest(px1,"M WW");
	//Double_t YKS = py5->KolmogorovTest(py1,"M WW");
	//Double_t XAD = px5->AndersonDarlingTest(px1);
	//Double_t YAD = py5->AndersonDarlingTest(py1);

	//TH1* px1c = px1->GetCumulative();
	//TH1* px5c = px5->GetCumulative();
	//scale hint1 to the pad coordinates
	
	
	
	// draw histogram together with its cumulative distribution
	TCanvas* canva1 = new TCanvas;
	canva->Divide(1,2);
	canva->cd(1);
	//px1->Scale(1/px1->Integral(), "width");
	px1->Draw();
	px1->SetLineColor(kBlue);
	//Float_t rightmax = 1.1*px5->GetMaximum();
	//Float_t scaleplot = gPad->GetUymax()/rightmax;
	//px5->Scale(scaleplot);
	//px5->Scale(1/px5->Integral(), "width");
	px5->Draw("same");
	px5->SetLineColor(kRed);	
	canva1->Update();
	canva->cd(2);
	TCanvas* canva2 = new TCanvas;
	//px1c->Scale(1/px1c->Integral(), "width");
	px1c->Draw();
	px1c->SetLineColor(kBlue);
	//Float_t rightmax2 = 1.1*px5c->GetMaximum();
	//Float_t scaleplot2 = gPad->GetUymax()/rightmax2;
	//px5c->Scale(scaleplot2);
	//px5c->Scale(1/px5c->Integral(), "width");
	px5c->Draw("same");
	px5c->SetLineColor(kRed);
	canva2->Update();
	// Delete histograms to prevent memory leaks when fitting	
	delete quasi_time;
	delete quasi_jdist;
	delete quasi_mrth1;
	delete quasi_mrth2;
	delete quasi_mrth3;
	delete quasi_Lth;
	delete plot5;

	// Set measures of fit and return one	
	//Double_t KStot = (XKS+YKS)*0.5;
	Double_t CHI2tot = (Xchi2+Ychi2)*0.5;
	Double_t CHI2NDFtot = (Xchi2ndf+Ychi2ndf)*0.5;
	//Double_t ADtot = (XAD+YAD)*0.5;

	//cout << "KS MAD error is: " << MADerror << endl;
	cout << "Chi2 MAD error is: " << ChiError << endl;
	cout << "Chi2/NDF MAD error is: " << ChiNDFError << endl;
	//cout << "KS MR error is: " << XKS << endl;	
	//cout << "KS proj error is: " << KStot << endl;
	cout << "Chi2 proj error is: " << CHI2tot << endl;
	cout << "Chi2NDF proj error is: " << CHI2NDFtot << endl;
	//cout << "AD proj error is: " << ADtot << endl;

	Double_t objective = ChiError;
	return objective;
}




//-----------------------------------------------------------------------------------------
//	Define function for MINUIT to fit
//-----------------------------------------------------------------------------------------
//	Arguments automatically dealt with by MINUIT, par values are inputs
//	for each loop of the minimiser
//-----------------------------------------------------------------------------------------


void minuitFCN(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag){
	double mean = par[0];
	double sigma = par[1];
	double sh_frac = par[2];
	cout << "Mean is " << mean << endl;
	cout << "Sigma is " << sigma << endl;
	cout << "Shell fraction is " << sh_frac << endl;
	f = madfitfcn(mean, sigma, sh_frac);
        cout << "\n Objective value is " << f << "\n" << endl;
}

//******************************************************************************************

//-----------------------------------------------------------------------------------------
//	Open experiment files globally
//-----------------------------------------------------------------------------------------

TFile* input1 = new TFile("YGSI3.run119.PEG.root"); // Input sorted root file

TH2D *plot1=(TH2D *)input1->Get("Xsec_pub/Exp_MAD");
						// Experiment  [counts] already binned MAD
TH2D *plot2=(TH2D *)input1->Get("Xsec_group/MrThCMxsec_T_mir_g");
						// Experiment  [dsigdth] for ang dist
TH2D *plot4=(TH2D *)input1->Get("Xsec_pub/MrThCMxsec_T_mir");
						// Experiment  [dsigdthdmr] for noramlised MAD


//-----------------------------------------------------------------------------------------
//	Minimisation routine
//-----------------------------------------------------------------------------------------

void minimiseMAD(){	
	
	//---------------------------------------------------------------------------------	
	//	Initialise Minuit object, set objective function for fitting, and
	//	initialise arguments
	//---------------------------------------------------------------------------------

	// Minuit object	
	TMinuit *gMinuit = new TMinuit(1);

	// Set function for minimisation
	gMinuit -> SetFCN(minuitFCN);

	// Initialise arguments
	bool converged;
   	Int_t ierflg = 0;
	Double_t arglist[10];
	
	// Set print level (1 - verbose, 0 - normal, -1 - quiet (suppress warnings))
	arglist[0] = 1;
	gMinuit->mnexcm("SET PRINT", arglist, 1, ierflg);

	// Define error as Chi2 error
	// (in future change to 0.5 for log-likelihood measure)
	arglist[0] = 1;
	gMinuit->mnexcm("SET ERR", arglist, 1, ierflg);
	

	//---------------------------------------------------------------------------------
	//	Parameter setting
	//---------------------------------------------------------------------------------

	// Set all parameters
	// Inputs:
	//	- Index (int): Start at 0 increment for additional parameters
	//	- Name (str): Name given to parameter, use something relevant to variable
	//	- Starting value (double): Fitting routine starts at this, set near where
	//	  the assumed fitted value is
	//	- Step value (double): Increment used by the routine, vary depending on
	//	  time available and confidence in starting guess
	//	- Lower bound (double): Minimum value the routine checks. Set to 0.
	//	  to remove bound, only change if fitting fails
	//	- Upper bound (double): Maximum value the routine checks. Set to 0.
	//	  to remove bound, only change if fitting fails

	SetParameter(0, "meanval", 7.3e-21, 1.000e-23, 1.000e-23, 3.000e-20, ierflg);
	SetParameter(1, "sigma", 2.800e-21, 1.000e-23, 1.000e-21, 1.000e-20, ierflg);
	SetParameter(2, "shellfraction", 0.78, 0.001, 0.5, 1.0, ierflg);

	// Fix parameters which don't need to be fitted (at starting value)
	// Input: Index (int) of parameter as set above
	gMinuit->FixParameter(2);
	gMinuit->FixParameter(1);

	
	
	//---------------------------------------------------------------------------------
	//	Minimisation step
	//---------------------------------------------------------------------------------
	
	// Initial function run
	gMinuit->mnexcm("CALL FCN", arglist ,1,ierflg);
	
	// Max number of iterations
	arglist[0]=50000;
	arglist[1]=1.;
	
	// Uses simple SIMPLEX approach to find good starting guess
	//gMinuit->mnexcm("SIMPLEX", arglist ,1 ,ierflg);
	//gMinuit->mnexcm("SCAN", arglist ,1 ,ierflg);
	// Continues search with more robust MIGRAD
   	//gMinuit->mnexcm("MIGRAD", arglist ,2 ,ierflg);	
	gMinuit->mnimpr();	
	converged = true;
	// Check MIGRAD convergence	
	//if (ierflg == 4) {

	//	converged = false;

	//}
	// If converged=false, finish fitting routine, print error
	// (try rerunning with higher iteration number)
	//if (!converged){
	//	cout << "MAD fit routine didn't converge" << endl;
	//	cout << "Check SIMPLEX results" << endl;
		
	//}
	// If converged=true, run MINOS for better estimate of minimum and errors (optional)
	//else{
	// Calculate Hesse matrix
	//gMinuit->mnexcm("HESSE", arglist, 2, ierflg);	
	// Improve error accuracy
	//gMinuit->mnexcm("MINOS", arglist ,0,ierflg);
	//}
	gROOT->SetBatch(0);
	gMinuit->SetGraphicsMode(kTRUE);
   	//gMinuit->mnscan();
	//gMinuit->mnplot();
	//gMinuit->mncont(0,0,10);
	//TGraph *graph = (TGraph*)gMinuit -> Contour(0,1);
	//TGraph->Draw();
	//gMinuit->mncomd("scan 0 100 0.5 1.0",ierflg);
	//TGraph *gr = (TGraph*)gMinuit->GetPlot();
	//gr->SetMarkerColor(1);
	//gr->SetMarkerStyle(20);
	//gr->Draw("AP");

	gMinuit->mnexcm("EXIT", arglist, 1, ierflg);
	
	//---------------------------------------------------------------------------------
	//	Print results
	//---------------------------------------------------------------------------------
	Double_t amin, edm, errdef;
	Int_t nvpar, nparx, icstat;
	gMinuit-> mnstat(amin,edm,errdef,nvpar,nparx,icstat);
	gMinuit->mnprin(3,amin);

}

	




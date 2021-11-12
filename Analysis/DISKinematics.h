//! D. Britzger, MPP, 2021
#ifndef __DISKinematics_H
#define __DISKinematics_H

#include <utility>
#include <vector>
#include <TLorentzVector.h>

// --------------------------------------------------------------------- //
//!
//! Calculate the DIS kinematic variables from three input observables
//!
//! All formulae assume massless particles.
//! For exact (within numerical precision) agreement among all definitions 
//! one needs to have energy-momentum conservation in the transverse plane
//! and no longitudinal energy-loss (e.g. ISR photon).
//!
// --------------------------------------------------------------------- //

class DISKinematics
{
   
public:

   DISKinematics(){};
   ~DISKinematics(){};

   void SetBeamElec  (const TLorentzVector& elec0 ) { fE0=elec0.E(); }
   void SetBeamElecEn(double Elec0En ) { fE0 = Elec0En; }         //!< Set electron beam energy
   void SetBeamProtEn(double Prot0En ) { fEp = Prot0En; }         //!< Set proton beam energy (only needed for x)
   void SetElecEn    (double ElecEn  ) { fElEn  = ElecEn;  }      //!< Set scattered electron energy
   void SetElecTh    (double ElecTh  ) { fElTh  = ElecTh;  }      //!< Set scattered electron polar angle theta
   void SetHFSSigma  (double Sigma  )  { fSigma = Sigma;  }       //!< Set Sigma of HFS, so E-Pz of all HFS particles
   //void SetHFSgamma  (double gamma  )  { fgamma = gamma;  }       //!< Set gamma of HFS, so polar angle of HFS 
   void SetHFSPt     (double pT  )     { fT     = pT;     }       //!< Set transverse momenta of HFS
   double GetHFSgamma () const { return 2.*atan(fSigma/fT);};

   
protected:     
   double fE0 = 0;      //! beam energy
   double fT  = 0;      //! pT of HFS
   //double fgamma = 0;   //! gamma of HFS
   double fSigma = 0;   //! Sigma of HFS
   double fElEn = 0;    //! Electron energy
   double fElTh = 0;    //! Electron angle
   double fEp = 920.;    //! Electron angle

public:
   double GetQ2_E0_E_theta() const {
      return 4*fE0*fElEn*pow(cos(fElTh/2.),2);
   }

   double GetQ2_E0_E_Sigma() const {
      return 4*fE0*fElEn-4*fE0*fE0+2.*fE0*fSigma;
   }

   // double GetQ2_E0_E_gamma() const {
   //    double theta1 = asin(fT/fElEn);
   //    double theta2 = acos(fT/fElEn)+M_PI/2.;
   //    cout<<"real: "<<fElTh<<"\tth: "<<theta1<<"\tth2: "<<theta2<<endl;
   //    cout<<"Q2  :             "<<"\tth: "<< 4*fE0*fElEn*pow(cos(theta1/2.),2)<<"\tth2: "<<4*fE0*fElEn*pow(cos(theta2/2.),2)<<endl;
   //    cout<<"Q2  :             "
   //        <<"\tth1: "<< 2*fE0*(fElEn + sqrt( fElEn*fElEn - fT*fT) )
   //        <<"\tth2: "<< 2*fE0*(fElEn - sqrt( fElEn*fElEn - fT*fT) )<<endl;
   //    cout<<"cos(th1): "<<cos(theta1)<<"\tcos2: "<<cos(theta2)<<endl;
   //    cout<<"Q2  :             "
   //        <<"\tth1: "<< 2*fE0*(fElEn + fElEn*cos(theta1) )
   //        <<"\tth2: "<< 2*fE0*(fElEn + fElEn*cos(theta2) )<<endl;
   //    if ( fElTh > M_PI/2.) 
   //       return 0;
   //    return 4*fE0*fElEn*pow(cos(theta2/2.),2);
   // }

   pair<double,double> GetQ2_E0_E_T() const {
      double ElmT = fElEn*fElEn - fT*fT;
      if ( ElmT < 0  ) ElmT = 0;// because of numerical limiatations, like 1.0001
      double sqrtE2T2 = sqrt(ElmT);
      double Q2a = 2*fE0*(fElEn - sqrtE2T2 ); // == fT*fT/(1-GetY_E0_En_T().first);
      double Q2b = 2*fE0*(fElEn + sqrtE2T2 ); // == fT*fT/(1-GetY_E0_En_T().second)
      return {Q2a,Q2b};
   }

   double GetQ2_E0_theta_Sigma() const {
      return 2.*fE0*(2.*fE0-fSigma)/pow(tan(fElTh/2.),2);
   }

   double GetQ2_E0_theta_gamma() const {
      //! double angle method
      return 4*fE0*fE0 / tan(fElTh/2.) / (tan(GetHFSgamma()/2.)+tan(fElTh/2.));
   }

   double GetQ2_E0_Sigma_gamma() const {
      return fT*fT/ (1 - GetYh());
   }

   double GetQ2_E_theta_Sigma() const {
      return fElEn*fElEn*pow(sin(fElTh),2)/(1-GetYs());
   }

   pair<double,double> GetQ2_E_Sigma_T() const {
      auto yy = GetY_E_Sigma_T();
      return {fT*fT/(1-yy.first), fT*fT/(1-yy.second)};
   }

   double GetQ2_E_theta_gamma() const {
      double tanth2 = tan(fElTh/2.);
      double tangam = tan(GetHFSgamma()/2.);
      //double gamma = GetHFSgamma();
      //double Q2_IDA = fElEn*fElEn * tanth2 * ( tangam + tanth2 )  / (1./tanth2 + tanth2); // bassler,bernardi 94 (not correct?)
      double Q2_IDA = fElEn*fElEn * sin(fElTh)*(1.+cos(fElTh)) * ( tangam + tanth2 );
      // double Q2_IDA2 = fElEn*fElEn * pow(sin(fElTh),2) / ( 1- tangam/(tangam+tan(fElTh/2.)) );
      // double Q2_IDA3 = fElEn*fElEn * ( tangam + tanth2 )/tanth2 * pow(sin(fElTh),2);
      // double Q2_yDA  = fElEn*fElEn * pow(sin(fElTh),2) / ( 1. - GetYda());
      return Q2_IDA ;
   }
   
   double GetQ2_theta_Sigma_gamma() const {
      return fT*fT/(1-GetYda());
   }



   // double GetY_En_Sigma_gamma() const {
   //    //double tantheta = asin(fT/fElEn);
   //    double TovE = fT/fElEn;
   //    if ( TovE > 1.) TovE=1; // because of numerical limiatations, like 1.00001
   //    double theta = asin(TovE);
   //    double theta2 = M_PI - asin(TovE);

   //    double th1 = M_PI/2. - (asin(TovE) - M_PI/2.);
   //    double th2 = M_PI/2. + (asin(TovE) - M_PI/2.);
   //    double tangam = tan(GetHFSgamma()/2.);
   //    double y0 = tangam/(tangam+tan(th2/2.));
   //    double y1 = tangam/(tangam+tan(th1/2.));
   //    double tanth1 = fT/(fElEn + sqrt(fElEn*fElEn - fT*fT ));
   //    double tanth2 = fT/(fElEn - sqrt(fElEn*fElEn - fT*fT ));
   //    double y2 = tangam/(tangam+tanth1);
   //    double y3 = tangam/(tangam+tanth2);
   //    double y4 = fSigma/(fSigma + (fT*tanth1) );
   //    double y5 = fSigma/(fSigma + (fT*tanth2) );

   //    return tangam/(tangam+tan(theta/2.));
   // }

   // ---------------------------------------------------------------------- //


   pair<double,double> GetY_E_Sigma_T() const {
      double ElmT = fElEn*fElEn - fT*fT;
      if ( ElmT < 0  ) ElmT = 0;// because of numerical limiatations, like 1.0001
      double sqrtE2T2 = sqrt(ElmT);
      double y1 = fSigma/(fSigma+fElEn+sqrtE2T2);
      double y2 = fSigma/(fSigma+fElEn-sqrtE2T2);
      return {y1,y2};
   }


   pair<double,double> GetY_E0_E_T() const {
      double ElmT = fElEn*fElEn - fT*fT;
      if ( ElmT < 0  ) ElmT = 0;// because of numerical limiatations, like 1.0001
      double y1 = 1. - (fElEn+sqrt(ElmT) )/(2.*fE0);
      double y2 = 1. - (fElEn-sqrt(ElmT) )/(2.*fE0);
      return {y1,y2};
   }

   double GetYh() const {
      return fSigma/(2.*fE0);
   }

   double GetYe() const {
      return 1. - fElEn/fE0*pow(sin(fElTh/2.),2);
   }

   double GetYeSigma() const {
      return 2.*fE0*fSigma/pow(fSigma+fElEn*(1.-cos(fElTh)),2);
   }

   double GetYs() const {
      return fSigma/(fSigma + fElEn*(1.-cos(fElTh)));
   }

   double GetYda() const {
      double tangam = tan(GetHFSgamma()/2.);
      return tangam/(tangam+tan(fElTh/2.));
   }


   // ---------------------------------------------------------------------- //
   double GetX_E0_E_theta() {
      //! Electron method
      return fE0*fElEn*pow(cos(fElTh/2.),2) / (fE0 - fElEn*pow(sin(fElTh/2.),2)) / fEp;
   }

   double GetX_E0_E_Sigma() {
      return GetQ2_E0_E_Sigma()/2./fSigma / fEp;
   }

   double GetX_E0_theta_Sigma() {
      return GetQ2_E0_theta_Sigma()/2./fSigma / fEp;;
   }

   double GetX_E0_theta_gamma() { 
      //! double-angle method
      return GetQ2_E0_theta_gamma() / ( 4.*fE0 * GetYda() )  / fEp;
   }

   double GetX_E0_Sigma_T() {
      //! Hadron method (JB)
      return fE0 * fT*fT / (2.*fE0*fSigma - fSigma*fSigma)  / fEp;
   }

   double GetX_E_theta_Sigma() {
      return fElEn * pow(cos(fElTh/2.),2) / GetYs()  / fEp;
   }

   double GetX_E_theta_gamma() {
      //! IDA Method
      double tanth2 = tan(fElTh/2.);
      double tangam = tan(GetHFSgamma()/2.);
      //return fElEn * (1./tangam + 1./tanth2) / (1./tanth2 - tanth2) / fEp;
      return fElEn * pow(cos(fElTh/2.),2) * (tangam + tanth2) / (tangam) / fEp;
   }

   double GetX_theta_Sigma_gamma() {
      return GetQ2_theta_Sigma_gamma() /2. / fSigma / fEp;
   }

   pair<double,double> GetX_E_Sigma_T() {
      return { GetQ2_E_Sigma_T().first /2./fSigma / fEp,
               GetQ2_E_Sigma_T().second/2./fSigma / fEp};
   } 

   pair<double,double> GetX_E0_E_T() {
      return {
         GetQ2_E0_E_T().first /2./fSigma / fEp,
            GetQ2_E0_E_T().second/2./fSigma / fEp};
   }
   

};

#endif

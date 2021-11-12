//! D. Britzger, MPP, 2021
#ifndef __H1BestKinematics_H
#define __H1BestKinematics_H

#include <utility>
#include <vector>
#include <TLorentzVector.h>

#include "H1Mods/H1PartMC.h"
#include "H1Calculator/H1Calculator.h"
#include "H1PhysUtils/H1BoostedJets.h"
#include "H1Calculator/H1CalcHad.h"
#include "H1Calculator/H1CalcHad.h"
#include "DISKinematics.h"
#include "H1Pointers/H1ShortPtr.h"
#include "H1Pointers/H1IntPtr.h"
// --------------------------------------------------------------------- //
//!
//! H1BestKinematics
//!
// --------------------------------------------------------------------- //

class H1BestKinematics
{
   
public:

   H1BestKinematics(){};
   ~H1BestKinematics(){};

   void CalculateAll(bool force = false) {
      if ( fThisEvNumber == gH1Calc->GetEventNumber() && !force ) 
         return;
      
      fThisEvNumber = gH1Calc->GetEventNumber();
      CalcGen();
      CalcRec();
   }

   const TLorentzVector& GetGenElec() {
      CalculateAll();
      return fGenElec;
   }

   const TLorentzVector& GetGenPhoton() {
      CalculateAll();
      return fGenPhoton;
   }
   
   const TLorentzVector& GetGenHFS() {
      CalculateAll();
      return fGenHFS;
   }

   const vector<H1PartMC*>& GetGenHFSarray() {
      CalculateAll();
      return fHadronarray;
   }

   double GetGenEmpz() {
      CalculateAll();
      return (fGenElec.E() - fGenElec.Pz()) + (fGenHFS.E()-fGenHFS.Pz());
   }

   DISKinematics GetGenKine() {
      CalculateAll();
      return fGenKine;
   }


   void SetUseGenIDA(bool UseIDA) {
      // use IDA method for gen-level, or
      bUseIDA = UseIDA;
      CalculateAll(true);
   }

   void SetUseGenISigma(bool UseISigma) {
      // use ISimga method for gen-level
      bUseIDA = !UseISigma;
      CalculateAll(true);
   }

   double GetGenQ2() {
      CalculateAll();
      return fGenQ2;
   }

   double GetGenY() {
      CalculateAll();
      return fGenY;
   }

   double GetGenX() {
      CalculateAll();
      return fGenX;
   }


   double GetRecEmpz() { 
      CalculateAll();
      return (fRecElec.E() - fRecElec.Pz() ) + (fRecHFS.E() - fRecHFS.Pz());
   }

   const TLorentzVector& GetRecElec() {
      CalculateAll();
      return fRecElec;
   }

   const TLorentzVector& GetRecHFS() {
      CalculateAll();
      return fRecHFS;
   }

   const vector<TLorentzVector>& GetRecParticleFourvectors() {
      CalculateAll();
      return fParticleFourvectors;
   }
   
   DISKinematics GetRecKinematics() {
      CalculateAll();
      return fRecKine;
   }

   //return 
   //vector<H1PartCand*>    fParticlearray;

protected:
   int fThisEvNumber=-1;

   // -------------------------------------------------- 
   // gen-level
   // -------------------------------------------------- 
   TLorentzVector       fElecBeam;
   TLorentzVector       fGenElec;
   TLorentzVector       fGenPhoton;
   TLorentzVector       fGenHFS;
   vector<H1PartMC*>    fHadronarray;
   bool                 fElecComb = false;
   bool bUseIDA = false;

   DISKinematics fGenKine;

   // preferred Gen-Kinematics
   double fGenQ2 = 0;
   double fGenY  = 0;
   double fGenX  = 0;
   // -------------------------------------------------- 


   // -------------------------------------------------- 
   // rec-level
   // -------------------------------------------------- 
   TLorentzVector         fRecElec;
   TLorentzVector         fRecHFS;
   vector<H1PartCand*>    fParticlearray;
   vector<TLorentzVector> fParticleFourvectors;
   // preferred Gen-Kinematics
   DISKinematics fRecKine;
   double fRecQ2 = 0;
   double fRecY  = 0;
   double fRecX  = 0;
   // -------------------------------------------------- 

protected:
   void CalcRec() {
      //! calculate all detector-level variables
      // bug!
      ////////////////////fRecElec       = gH1Calc->Elec()->GetFirstElectronGen();
      fRecElec       = gH1Calc->Elec()->GetFirstElectron();
      fParticlearray.clear();
      fParticleFourvectors.clear();
      auto tmparray = to_vector<H1PartCand*>(H1BoostedJets::Instance()->GetHFSArray());
      fRecHFS.SetPxPyPzE(0,0,0,0);
      
      // find photon
      H1PartCand* photonCand =NULL;
      double photrminval = 10000.;
      const double cut_photEmin     = 3.0; // minimum photon energy
      const double cut_photThetamin = 0.8; // minimum polar angle
      for ( H1PartCand* part : tmparray ) {
         if ( part->IsPhoton() && part->GetTheta()> cut_photThetamin && part->GetE()>cut_photEmin) {
            double rmin = min(fRecElec.Angle(part->GetFourVector().Vect()),  (M_PI - part->GetTheta()) );
            if ( photonCand == NULL  ||  rmin<photrminval )  {
               photonCand  = part;
               photrminval = rmin;
            }
         }
      }

      // treat photon
      if ( photonCand != NULL ) {
         double angleBeam = M_PI - photonCand->GetTheta();
         double angleElec = fRecElec.Angle(photonCand->GetFourVector().Vect());
         if ( angleBeam < M_PI/3. && angleBeam<angleElec ) { // possibly ISR
            if ( angleBeam > 0.8 ) // it is not an ISR
               photonCand = NULL;
         }
         else { // possibly FSR
            if ( angleElec < 1.2 ) { // merge it with elec
               fRecElec = fRecElec + photonCand->GetFourVector(); 
            }
            else 
               photonCand = NULL;
         }
      }

      // fill fourvectors
      for ( H1PartCand* part : tmparray ) {
         // is it possibly a photon ?
         if ( part == photonCand ) {
            continue;
         } // end photon

         TLorentzVector part_v4 = part->GetFourVector();

         // is it too close to electron
         if ( fRecElec.DeltaR( part_v4 ) < 0.2 ) {
            continue;
         }

         if ( part_v4.M() > 1.1 ) {  // check for reasonable mass assumption
            //cout<<"large mass!  M="<<part_v4.M()<<endl;
            //part->Print();
            static const double m_pion = 0.13957;
            part_v4 = part->GetFourVector(m_pion);
         }

         fParticlearray.     push_back(part); // keep PartCand
         fParticleFourvectors.push_back(part_v4); // keep PartCand
         fRecHFS += part_v4;
      }

      fRecKine.SetBeamElecEn( fElecBeam.E() ); // 27.6                                                                  
      fRecKine.SetElecEn( fRecElec.E() );
      fRecKine.SetElecTh( fRecElec.Theta() );
      fRecKine.SetHFSPt( fRecHFS.Pt() );
      fRecKine.SetHFSSigma( fRecHFS.E() - fRecHFS.Pz() );

   }

   void CalcGen() {
      //! calculate all generator-level variables
      
      // beam
      static H1FloatPtr ElecBeamEnergy("EBeamE");
      fElecBeam.SetPxPyPzE(0,0,-(*ElecBeamEnergy),(*ElecBeamEnergy));
      fElecComb = false;

      // photon first !
      fGenPhoton = MakeGen_Photon();

      // electron !
      bool IsIgnore = false;
      fGenElec   = MakeGen_UncombElectron();
      if (fGenPhoton.E() != 0 ) {
         double angleEscat = fGenPhoton.Angle(fGenElec.Vect());
         double angleBeam  = M_PI - fGenPhoton.Theta();
         if ( !(angleBeam < M_PI/3. && angleBeam<angleEscat) ) { // FSR
            fGenElec = fGenElec+fGenPhoton; // recombine gamma+elec
            fElecComb = true;
         }
         else 
            IsIgnore = true;
      }
      
      // HFS
      //const TLorentzVector& HFS  = gH1Calc->Had()->GetHadGen();
      fGenHFS = gH1Calc->Had()->GetHadGen(); 
      H1BoostedJets::Instance()->UseGenElecGammaCombined(fElecComb);
      H1BoostedJets::Instance()->ExcludeGenGamma(true); // must be called AFTER 'UseGenElecGammaCombined'
      fHadronarray = to_vector<H1PartMC*>(H1BoostedJets::Instance()->GetHadronArray());
      TLorentzVector HFS2 = H1BoostedJets::Instance()->FillHadronArray();
      //cout<<"HFS ratio: iscomb: "<<fElecComb<<"\tisIgnore: "<<IsIgnore <<"\tHFSratio"<<fGenHFS.E() /HFS2.E()<<endl;

      // kinematics
      fGenKine.SetBeamElecEn( (*ElecBeamEnergy) ); // 27.6                                                                  
      fGenKine.SetElecEn( fGenElec.E() );
      fGenKine.SetElecTh( fGenElec.Theta() );
      fGenKine.SetHFSPt( fGenHFS.Pt() );
      fGenKine.SetHFSSigma( fGenHFS.E() - fGenHFS.Pz() );
      
      // 'truth' kinematics
      if ( bUseIDA ) {
         // use IDA method
         fGenQ2= fGenKine.GetQ2_E_theta_gamma();
         fGenY = fGenKine.GetYda();
         fGenX = fGenKine.GetX_E_theta_gamma();
      }
      else { // use ISigma method 
         fGenQ2= fGenKine.GetQ2_E_theta_Sigma();
         fGenY = fGenKine.GetYs();
         fGenX = fGenKine.GetX_E_theta_Sigma();
      }
   }


   TLorentzVector MakeGen_UncombElectron() const {
      //! get the generated (uncombined) electron  
      static H1FloatPtr EPtr("GenEnElecUncombined");
      static H1FloatPtr ThPtr("GenThElecUncombined");
      static H1FloatPtr PhPtr("GenPhElecUncombined");
      Float_t elecE = (*EPtr);
      Float_t elecTh = (*ThPtr);
      Float_t elecPh = (*PhPtr);
      Float_t elecP  = elecE; // ignore the mass                                                                            

      Float_t elecPx = elecP*TMath::Sin(elecTh)*TMath::Cos(elecPh);
      Float_t elecPy = elecP*TMath::Sin(elecTh)*TMath::Sin(elecPh);
      Float_t elecPz = elecP*TMath::Cos(elecTh);
      return TLorentzVector(elecPx, elecPy, elecPz, elecE);
   }
   
   TLorentzVector MakeGen_Photon() const {
      static H1ShortPtr GenRadType("GenRad");
      if ((*GenRadType)==0) // is there no radiation?
         return TLorentzVector(0,0,0,0);

      // get the photon
      static H1FloatPtr EnPhPtr("GenEnPhoton");
      static H1FloatPtr ThPhPtr("GenThPhoton");
      static H1FloatPtr PhPhPtr("GenPhPhoton");
      Float_t phEn = (*EnPhPtr);
      Float_t phTh = (*ThPhPtr);
      Float_t phPh = (*PhPhPtr);

      Float_t phPx = phEn*TMath::Sin(phTh)*TMath::Cos(phPh);
      Float_t phPy = phEn*TMath::Sin(phTh)*TMath::Sin(phPh);
      Float_t phPz = phEn*TMath::Cos(phTh);
      return TLorentzVector(phPx, phPy, phPz, phEn);
   }


   
private:
   //! convert a TObjArray into a std::vector, using proper type_cast
   template<class TObj, class TArr>
   std::vector<TObj> to_vector(TArr* array) {
      std::vector<TObj> ret(array->GetEntries());
      for ( int i = 0 ; i<array->GetEntries() ; i++ )
         ret[i] = static_cast<TObj>(array->At(i));
      return ret;
   }


};

#endif

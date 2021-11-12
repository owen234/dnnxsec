// (c) MPI 2020
#include "AnalysisEventShapes.h"
 
// // C++ includes
// #include <map>
// #include <set>
#include <vector>

// Root includes
// #include <TH1.h>
// #include <TH2.h>
 
// // H1 includes
#include "H1PhysUtils/H1BoostedJets.h"
#include "H1HadronicCalibration/H1HadronicCalibration.h"

#include "H1Calculator/H1CalcGenericInterface.h"
using namespace H1CalcGenericInterface;

#include "H1Calculator/H1Calculator.h"
#include "H1Calculator/H1CalcTrig.h"
#include "H1Calculator/H1CalcWeight.h"
#include "H1Calculator/H1CalcVertex.h"
#include "H1Calculator/H1CalcEvent.h"
#include "H1Calculator/H1CalcKine.h"
#include "H1Calculator/H1CalcElec.h"
#include "H1Calculator/H1CalcFs.h"
#include "H1Calculator/H1CalcHad.h"
#include "H1Calculator/H1CalcTrack.h"
#include "H1Calculator/H1CalcSystematic.h"
#include "H1Mods/H1PartCandArrayPtr.h"

// #include "H1PhysUtils/H1MakeKine.h"
#include "EventshapeTools.h"
#include "JetTools.h"         // JetsAtHighQ2
#include "H1Mods/H1PartSelTrack.h"
#include "H1JetFinder/H1EventShape.h"
//#include "UsefulTools.h"


// fjcontrib
#include "fastjet/PseudoJet.hh"
#include "fastjet/contrib/Centauro.hh"


//--- owen
#include "H1Mods/H1PartMCArrayPtr.h"
#include "H1Mods/H1PartMC.h"

#include "H1Mods/H1PartCandArrayPtr.h"
#include "H1Mods/H1PartCand.h"

#include "H1Mods/H1PartSelTrackArrayPtr.h"
#include "H1Mods/H1PartSelTrack.h"

#include "H1Skeleton/H1Tree.h"


float calc_dr( double phi1, double phi2, double eta1, double eta2 ) {
   float deta = fabs( eta1 - eta2 ) ;
   float dphi = phi1 - phi2 ;
   if ( dphi < -3.14159265 ) dphi += 2*3.14159265 ;
   if ( dphi >  3.14159265 ) dphi -= 2*3.14159265 ;
   return sqrt( deta*deta + dphi*dphi ) ;
}


using namespace std;
using namespace fastjet;

// _______________________________________________________ //
//! Constructor
AnalysisEventShapes::AnalysisEventShapes(TString chain) : AnalysisBase(chain) {
//--- owen
   do_dump = false ;
   //do_dump = true ;
   if ( do_dump) {
      printf("\n\n Opening dump file: dump-file.txt\n\n") ;
      dump_file = fopen( "dump-file.txt", "w" ) ;
   }
}

// _______________________________________________________ //
AnalysisEventShapes::~AnalysisEventShapes() {
//--- owen
   if ( do_dump ) {
      printf("\n\n closing dump-file.txt\n\n") ;
      fclose( dump_file ) ;
   }
}



// _______________________________________________________ //
//!
//!  AnalysisEventShapes::DoInitialSettings()
//!
//!  This function is called by the main program
//!  for the very first event once
//!
void AnalysisEventShapes::DoInitialSettings() {
   // nothing todo
   // ... initialize reweightings etc...
}



// _______________________________________________________ //
//!
//!  AnalysisEventShapes::DoReset()
//!
//!  This function is called by the main program
//!  at the beginning of the event loop
//!
//!  Reset all members.
//!  Note: Also AnalysisBase::DoReset() is called
//!
void AnalysisEventShapes::DoReset() {
   // -- reset event quantites
   fGen     = CrossSectionQuantities();
   fRec     = CrossSectionQuantities();
   fTreeVar = TreeVariables();

}



// _______________________________________________________ //
//!
//!  AnalysisEventShapes::DoAnalysisCutsGen()
//!
//!  This function is called by the main program
//!  at the beginning of the event loop
//!
//!  Define analysis specific generator level cuts
//!
bool AnalysisEventShapes::DoAnalysisCutsGen() {
   // -- reset event quantites
   fAnalysisCutsGen = true;

   //Exclude MC Events to avoid double counting
   // Q2<4 is covered by the Pythia photo production
   // Q2>60 is covered by Django and Rapgap


   if ( IsBkgMC && fChainName=="DjBkg" ){
      if( gH1Calc->Kine()->GetQ2Gen() > 60 || gH1Calc->Kine()->GetQ2Gen() < 4) {
         fAnalysisCutsGen = false;
      }
   }

   // set fGen.IsGood
   fGen.IsGood  =  fAnalysisCutsGen && fBasicCutsGen;
   return fAnalysisCutsGen;
}

// _______________________________________________________ //
//!
//!  AnalysisEventShapes::DoAnalysisCutsRec()
//!
//!  This function is called by the main program
//!  at the beginning of the event loop
//!
//!  Define analysis specific detector level cuts
//!
bool AnalysisEventShapes::DoAnalysisCutsRec() {
   // -- reset event quantites
   fAnalysisCutsRec = true;
   
   // set fRec.IsGood
   fRec.IsGood  =  fAnalysisCutsRec && fBasicCutsRec;

   return fAnalysisCutsRec;
}



// _______________________________________________________ //
//!
//!  AnalysisEventShapes::DoCrossSectionObservablesGen()
//!
//!  This function is called by the main program
//!  during the event loop
//!
//!  Set observables needed to make the cross section
//!  histograms. Store them in (GenLevQuantities) fGen
//!
void AnalysisEventShapes::DoCrossSectionObservablesGen() {



   //--- owen

   bool has_isr ;
   bool has_fsr ;


   int se_pi = -1 ;
   int rad_gam_pi = -1 ;

   int beam_p_pi = -1 ;
   int beam_e_pi = -1 ;


   has_isr = false ;
   has_fsr = false ;

   float beam_p_e ;
   float beam_e_e ;


   static H1PartMCArrayPtr partMc ;
   if ( do_dump ) {
      fprintf( dump_file, "\n\n ============ Number of MC particles: %d\n", partMc.GetEntries() ) ;
      for ( int pi = 0; pi < partMc.GetEntries(); pi++ ) {
          H1PartMC* pmc = partMc[pi] ;
          fprintf( dump_file, "  %4d :  PDG %8d  M1 %4d M1 %4d D1 %4d D2 %4d  stat %4d  Px %9.5f Py %9.5f Pz %9.5f\n",
               pi, pmc->GetPDG(), pmc->GetMother1(), pmc->GetMother2(), pmc->GetDaughter1(), pmc->GetDaughter2(),
               pmc->GetStatus(),
               pmc->GetPx(), pmc->GetPy(), pmc->GetPz() ) ;
      } // pi
   }

   for ( int pi = 0; pi < partMc.GetEntries(); pi++ ) {

       H1PartMC* pmc = partMc[pi] ;

       int pid = pmc->GetPDG() ;
       int status = pmc->GetStatus() ;
       int mom1_pi = pmc->GetMother1() ;

       if ( status == 201 && abs(pid) == 11 ) {  //-- usually (always?) pi == 0
          beam_e_e = pmc->GetE() ;
          beam_e_pi = pi ;
          continue ;
       }
       if ( status == 201 && pid == 2212 ) {  //-- usually (always?) pi == 1
          beam_p_e = pmc->GetE() ;
          beam_p_pi = pi ;
          continue ;
       }
       if ( status == 0 && abs(pid) == 11 && se_pi < 0 ) { //-- usually (always?) pi == 3
          se_pi = pi ;
          continue ;
       }
       if ( status == 202 && pid == 22 && rad_gam_pi < 0 ) { //-- usually (always?) pi == 4, if present
          rad_gam_pi = pi ;
          if ( mom1_pi == beam_e_pi ) {
             has_isr = true ;
          }
          if ( mom1_pi == se_pi    ) {
             has_fsr = true ;
          }
          continue ;
       }

   } // pi

   if ( do_dump ) {
      fprintf( dump_file, "  se_pi %4d  rad_gam_pi %4d  beam_e_pi %4d beam_p_pi %4d\n",
        se_pi,   rad_gam_pi,   beam_e_pi, beam_p_pi ) ;
      if ( has_isr ) fprintf( dump_file, "  Has ISR.\n") ;
      if ( has_fsr ) fprintf( dump_file, "  Has FSR.\n") ;
      if ( rad_gam_pi >= 0 ) {
         float dr_se = calc_dr( partMc[rad_gam_pi]->GetPhi(), partMc[se_pi]->GetPhi(),  partMc[rad_gam_pi]->GetEta(), partMc[se_pi]->GetEta() ) ;
         float cos_angle_beam = ( partMc[rad_gam_pi]->GetPx() * partMc[beam_e_pi]->GetPx()
                            + partMc[rad_gam_pi]->GetPy() * partMc[beam_e_pi]->GetPy()
                            + partMc[rad_gam_pi]->GetPz() * partMc[beam_e_pi]->GetPz()) / ( partMc[rad_gam_pi]->GetP() * partMc[beam_e_pi]->GetP() ) ;
         float angle_beam = acos( cos_angle_beam ) ;
         fprintf( dump_file, "   radiated photon, dr_se = %9.5f  angle_beam = %9.5f\n", dr_se, angle_beam ) ;
      }
   }

   if ( beam_p_pi < 0 ) {
      ////printf("\n\n *** No beam proton in gen particle list!\n\n") ;
      return ;
   }
   if ( beam_e_pi < 0 ) {
      ////printf("\n\n *** No beam electron in gen particle list!\n\n") ;
      return ;
   }
   if ( se_pi < 0 ) {
      ////printf("\n\n *** No scattered electron in gen particle list!\n\n") ;
      return ;
   }

   TLorentzVector tlv_p( partMc[beam_p_pi]->GetPx(),  partMc[beam_p_pi]->GetPy(),  partMc[beam_p_pi]->GetPz(),  partMc[beam_p_pi]->GetE() ) ;
   TLorentzVector tlv_l( partMc[beam_e_pi]->GetPx(),  partMc[beam_e_pi]->GetPy(),  partMc[beam_e_pi]->GetPz(),  partMc[beam_e_pi]->GetE() ) ;
   TLorentzVector tlv_lprime( partMc[se_pi]->GetPx(),  partMc[se_pi]->GetPy(),  partMc[se_pi]->GetPz(),  partMc[se_pi]->GetE() ) ;
   TLorentzVector tlv_rad_pho(0.,0.,0.,0.) ;
   if ( has_isr && rad_gam_pi >= 0 ) {
      TLorentzVector tlv_pho( partMc[rad_gam_pi]->GetPx(),  partMc[rad_gam_pi]->GetPy(),  partMc[rad_gam_pi]->GetPz(),  partMc[rad_gam_pi]->GetE() ) ;
      tlv_rad_pho = tlv_pho ;
      tlv_l = tlv_l - tlv_pho ; // subtract gen radiated photon from gen beam electron for ISR to get post-radiation beam electron 4-vector.
   }
   if ( has_fsr && rad_gam_pi >= 0 ) {
      TLorentzVector tlv_pho( partMc[rad_gam_pi]->GetPx(),  partMc[rad_gam_pi]->GetPy(),  partMc[rad_gam_pi]->GetPz(),  partMc[rad_gam_pi]->GetE() ) ;
      tlv_rad_pho = tlv_pho ;
      tlv_lprime = tlv_lprime + tlv_pho ; // add gen radiated photon to gen scattered electron for FSR to get pre-radiation scattered electron 4-vector.
   }

   TLorentzVector tlv_q = tlv_l - tlv_lprime ;

   float from_tlv_gen_Q2  ;
   float from_tlv_gen_y  ;
   float from_tlv_gen_x  ;

   from_tlv_gen_Q2 = -1 * tlv_q * tlv_q ;
   from_tlv_gen_y = ( tlv_p * tlv_q ) / ( tlv_p * tlv_l ) ;
   from_tlv_gen_x = -1 * ( tlv_q * tlv_q ) / ( 2. * tlv_p * tlv_q ) ;



   static EventshapeTools ESTools;
   // event weight
   // fGen.wgt      = gH1Calc->Weight()->GetWeightGen(); // gH1Calc->GetFloatVariable(Weight_Weight);
   //fGen.BoostToBreit = ESTools.BoostToBreitFrame(fGen.Q2, fGen.Y, fGen.X, ScatElecGen.Phi());
   //fGen.BoostToLab = ESTools.BoostToLabFrame(fGen.BoostToBreit);
   //mini-tree                                                                                                                                                                                            
  

   //  TLorentzVector virtualphoton = ebeam - ScatElec;

   // --------------- general kinematics -------------------                                                                                                                                                
   TLorentzVector ebeam, pbeam;
   H1BoostedJets::Instance()->BeamVectors(&pbeam, &ebeam);

   fTreeVar.event_weight = gH1Calc->Weight()->GetWeightGen();
   TLorentzVector ScatElecGen = gH1Calc->Elec()->GetFirstElectronGen();
   TLorentzVector virtual_photon = ebeam - ScatElecGen;

   vector<H1PartMC*> hadronarray =
     to_vector<H1PartMC*>(H1BoostedJets::Instance()->GetHadronArray());
   const TLorentzVector& HFS  = gH1Calc->Fs()->GetAllFsLessElectronGen();

   fTreeVar.gen_event_Q2 =  gH1Calc->Kine()->GetQ2Gen();
   fTreeVar.gen_event_y  =  gH1Calc->Kine()->GetYGen();
   fTreeVar.gene_px = ScatElecGen.Px();
   fTreeVar.gene_py = ScatElecGen.Py();
   fTreeVar.gene_pz = ScatElecGen.Pz();
   fTreeVar.gene_eta = ScatElecGen.Eta();
   fTreeVar.genHFS_px = HFS.Px();
   fTreeVar.genHFS_py = HFS.Py();
   fTreeVar.genHFS_pz = HFS.Pz();
   fTreeVar.genHFS_E = HFS.E();
   fTreeVar.genHFS_eta = HFS.Eta();

   fTreeVar.gen_e_e = ScatElecGen.E() ;
   fTreeVar.gen_e_pz = ScatElecGen.Pz() ;
   fTreeVar.gen_e_pt = ScatElecGen.Pt() ;
   fTreeVar.gen_e_phi = ScatElecGen.Phi() ;
   fTreeVar.gen_e_eta = ScatElecGen.Eta() ;

   fTreeVar.gen_hfs_e = HFS.E() ;
   fTreeVar.gen_hfs_pz = HFS.Pz() ;
   fTreeVar.gen_hfs_pt = HFS.Pt() ;
   fTreeVar.gen_hfs_phi = HFS.Phi() ;
   fTreeVar.gen_hfs_eta = HFS.Eta() ;

   float dphi = ScatElecGen.Phi() - HFS.Phi() ;
   if ( dphi < -3.14159265 ) dphi += 2 * 3.14159265 ;
   if ( dphi >  3.14159265 ) dphi -= 2 * 3.14159265 ;
   if ( dphi < 0 ) dphi += 2 * 3.14159265 ;  //-- center this at pi.
   fTreeVar.gen_dphi = dphi ;




   fGen.Q2 = gH1Calc->Kine()->GetQ2Gen();

   fTreeVar.from_tlv_gen_Q2 = from_tlv_gen_Q2 ;
   fTreeVar.from_tlv_gen_x = from_tlv_gen_x ;
   fTreeVar.from_tlv_gen_y = from_tlv_gen_y ;

   fTreeVar.has_isr = has_isr ;
   fTreeVar.has_fsr = has_fsr ;

   fTreeVar.beam_e_e = beam_e_e ;
   fTreeVar.beam_p_e = beam_p_e ;

}



// _______________________________________________________ //
//!
//!  AnalysisEventShapes::DoCrossSectionObservablesRec()
//!
//!  This function is called by the main program
//!  during the event loop
//!
//!  Set observables needed to make the cross section
//!  histograms. Store them in (RecLevQuantities) fRec
//!
void AnalysisEventShapes::DoCrossSectionObservablesRec() {

   static EventshapeTools ESTools;
   // event weight
   fRec.wgt      = gH1Calc->Weight()->GetWeight();
   // initial and final electron and quark
   TLorentzVector ebeam, pbeam;
   H1BoostedJets::Instance()->BeamVectors(&pbeam, &ebeam);
   TLorentzVector ScatElec = gH1Calc->Elec()->GetFirstElectron();
   TLorentzVector virtual_photon = ebeam - ScatElec;



   vector<H1PartCand*> particlearray = to_vector<H1PartCand*>(H1BoostedJets::Instance()->GetHFSArray());

   //minitree

   Float_t hadptda = gH1Calc->Fs()->GetHadPtDa();
   TLorentzVector HFS = gH1Calc->Fs()->GetAllFsLessElectron();


   fTreeVar.event_Q2_e = gH1Calc->Kine()->GetQ2e();
   fTreeVar.event_y_e  = gH1Calc->Kine()->GetYe();
   fTreeVar.event_Q2_sigma = gH1Calc->Kine()->GetQ2s();
   fTreeVar.event_y_sigma  = gH1Calc->Kine()->GetYs();
   fTreeVar.event_Q2_da = gH1Calc->Kine()->GetQ2da();
   fTreeVar.event_y_da  = gH1Calc->Kine()->GetYda();
   fTreeVar.event_Q2_esigma = gH1Calc->Kine()->GetQ2es();
   fTreeVar.event_y_esigma  = gH1Calc->Kine()->GetYes();

   fTreeVar.event_y_h = gH1Calc->Kine()->GetYh();
   fTreeVar.event_Q2_h = gH1Calc->Kine()->GetQ2h();

   fTreeVar.Empz     = gH1Calc->Fs()->GetEmpz();
   fTreeVar.e_px = ScatElec.Px();
   fTreeVar.e_py = ScatElec.Py();
   fTreeVar.e_pz = ScatElec.Pz();
   fTreeVar.e_eta = ScatElec.Eta();
   fTreeVar.HFS_px = HFS.Px();
   fTreeVar.HFS_py = HFS.Py();
   fTreeVar.HFS_pz = HFS.Pz();
   fTreeVar.HFS_E = HFS.E();
   fTreeVar.HFS_eta = HFS.Eta();

   fTreeVar.ptmiss =  gH1Calc->Fs()->GetPtMiss();
   fTreeVar.pth    =  gH1Calc->Fs()->GetPtCalo();
   fTreeVar.vertex_z = gH1Calc->Vertex()->GetZ();
   fTreeVar.ptratio_da = HFS.Pt()/hadptda;
   fTreeVar.ptratio_ele = HFS.Pt()/ScatElec.Pt();

   fTreeVar.obs_e_e = ScatElec.E() ;
   fTreeVar.obs_e_pz = ScatElec.Pz() ;
   fTreeVar.obs_e_pt = ScatElec.Pt() ;
   fTreeVar.obs_e_phi = ScatElec.Phi() ;
   fTreeVar.obs_e_eta = ScatElec.Eta() ;

   fTreeVar.obs_hfs_e = HFS.E() ;
   fTreeVar.obs_hfs_pz = HFS.Pz() ;
   fTreeVar.obs_hfs_pt = HFS.Pt() ;
   fTreeVar.obs_hfs_phi = HFS.Phi() ;
   fTreeVar.obs_hfs_eta = HFS.Eta() ;

   float dphi = ScatElec.Phi() - HFS.Phi() ;
   if ( dphi < -3.14159265 ) dphi += 2 * 3.14159265 ;
   if ( dphi >  3.14159265 ) dphi -= 2 * 3.14159265 ;
   if ( dphi < 0 ) dphi += 2 * 3.14159265 ;  //-- center this at pi.
   fTreeVar.obs_dphi = dphi ;


   //-- make sure beam energies are set for data.
   if ( !(gH1Tree ->IsMC()) ) {
      fTreeVar.beam_e_e = ebeam.E() ;
      fTreeVar.beam_p_e = pbeam.E() ;
   }

   int im ;

   float e ;
   float e0 ;
   float theta ;
   float Sigma ;
   float T ;
   float tan_gamma_over_2 ;
   float ep ;

   e = fTreeVar.obs_e_e ;
   e0 = fTreeVar.beam_e_e ;
   theta = atan2( fTreeVar.obs_e_pt, fTreeVar.obs_e_pz ) ;
   Sigma = fTreeVar.obs_hfs_e - fTreeVar.obs_hfs_pz  ;
   T = fTreeVar.obs_hfs_pt ;
   tan_gamma_over_2 = Sigma / T ;
   ep = fTreeVar.beam_p_e ;



  //-- Method  0 : electron
   im = 0 ;

   fTreeVar.obs_y[im] = 1. - (e/e0) * pow( sin(theta/2.), 2 ) ;

   fTreeVar.obs_Q2[im] = 4. * e0 * e * pow( cos(theta/2.), 2 ) ;

   fTreeVar.obs_x[im] = ( e0 * e * pow( cos(theta/2.), 2 ) ) / ( ep * ( e0 - e * pow( sin(theta/2), 2 ) ) ) ;




  //-- Method  1 : unpub1
   im = 1 ;

   fTreeVar.obs_y[im] = Sigma / ( 2. * e0 ) ; // yh

   fTreeVar.obs_Q2[im] = 4. * e0 * e - 4. * e0 * e0 + 2 * e0 * Sigma ;

   fTreeVar.obs_x[im] = fTreeVar.obs_Q2[im] / ( 2. * ep * Sigma ) ;


  //-- Method  2 : unpub2
   im = 2 ;

   fTreeVar.obs_y[im] = Sigma / ( 2. * e0 ) ; // yh

   fTreeVar.obs_Q2[im] = 2. * e0 * ( 2. * e0 - Sigma ) / pow( tan(theta/2.), 2. ) ;

   fTreeVar.obs_x[im] = fTreeVar.obs_Q2[im] / ( 2. * ep * Sigma ) ;




  //-- Method  3 : DA
   im = 3 ;

   fTreeVar.obs_y[im] = tan_gamma_over_2 / ( tan_gamma_over_2 + tan(theta/2.) ) ;

   fTreeVar.obs_Q2[im] = 4. * e0 * e0 * (1./tan(theta/2.)) / ( tan_gamma_over_2 + tan(theta/2.) ) ;

   fTreeVar.obs_x[im] = fTreeVar.obs_Q2[im] / ( 4. * ep * e0 * fTreeVar.obs_y[im] ) ;




  //-- Method  4 : h
   im = 4 ;

   fTreeVar.obs_y[im] = Sigma / ( 2. * e0 ) ;

   fTreeVar.obs_Q2[im] = T * T / ( 1. - fTreeVar.obs_y[im] ) ;

   fTreeVar.obs_x[im] = fTreeVar.obs_Q2[im] / ( 2. * ep * Sigma ) ;




  //-- Method  5 : ISigma
   im = 5 ;


   fTreeVar.obs_y[im] = Sigma / ( Sigma + e * ( 1. - cos(theta) ) ) ;

   fTreeVar.obs_Q2[im] = e * e * pow( sin(theta), 2 ) / ( 1. - fTreeVar.obs_y[im] ) ;

   fTreeVar.obs_x[im] = e * pow( cos(theta/2.), 2 ) / ( ep * fTreeVar.obs_y[im] ) ;




  //-- Method  6 : IDA
   im = 6 ;

   fTreeVar.obs_y[im] = tan_gamma_over_2 / ( tan_gamma_over_2 + tan(theta/2.) ) ; // yDA

   //fTreeVar.obs_Q2[im] = e * e * tan(theta/2.) * ( tan_gamma_over_2 + tan(theta/2.) ) / ( 1./tan(theta/2.) + tan(theta/2.) ) ; // *** this is wrong somehow
   fTreeVar.obs_Q2[im] = 4. * e * e * (1./tan(theta/2.)) * ( tan_gamma_over_2 + tan(theta/2.) ) / pow( (1./tan(theta/2.) + tan(theta/2.)), 2 ) ;

   fTreeVar.obs_x[im] = e * pow( cos(theta/2.), 2. ) / ( ep * fTreeVar.obs_y[im] ) ;


  //-- Method  7 : unpub3
   im = 7 ;

   fTreeVar.obs_y[im] = tan_gamma_over_2 / ( tan_gamma_over_2 + tan(theta/2.) ) ; // yDA

   fTreeVar.obs_Q2[im] = pow( Sigma / tan_gamma_over_2, 2 ) / ( 1. - fTreeVar.obs_y[im] ) ;

   fTreeVar.obs_x[im] = e * pow( cos(theta/2.), 2. ) / ( ep * fTreeVar.obs_y[im] ) ;





  //-- Method  8 : eSigma
   im = 8 ;

   fTreeVar.obs_y[im] = 2. * e0 * Sigma / pow( ( Sigma + e*(1.-cos(theta)) ) , 2 ) ;

   fTreeVar.obs_Q2[im] = 4. * e0 * e * pow( cos(theta/2.), 2 ) ; // Q2e

   fTreeVar.obs_x[im] = ( e * pow( cos(theta/2), 2 ) + (e * e/(2.*Sigma)) * pow( sin(theta), 2 ) ) / ep ; // *** fixed typo in note







   //---- Find and save the scattered electron track.

   static H1PartSelTrackArrayPtr partTrk ;
   float trk_min_dr = 999. ;
   float trk_min_ti = -1 ;
   for ( int ti=0; ti<partTrk.GetEntries(); ti++ ) {
      H1PartSelTrack* trk = partTrk[ti] ;
      float dr = calc_dr( fTreeVar.obs_e_phi,  trk->GetPhi(),   fTreeVar.obs_e_eta, trk->GetEta() ) ;
      if ( dr < trk_min_dr ) {
         trk_min_dr = dr ;
         trk_min_ti = ti ;
      }
   } // ti
   if ( trk_min_ti >= 0 && trk_min_dr < 0.05 ) {
      H1PartSelTrack* trk = partTrk[trk_min_ti] ;
      fTreeVar.obs_e_trk_e = trk->GetE() ;
      fTreeVar.obs_e_trk_eta = trk->GetEta() ;
      fTreeVar.obs_e_trk_phi = trk->GetPhi() ;
   } else {
      fTreeVar.obs_e_trk_e = -1. ;
      fTreeVar.obs_e_trk_eta = 5. ;
      fTreeVar.obs_e_trk_phi = 5. ;
   }


   if ( do_dump ) {
      fprintf( dump_file, "\n\n Number of SelTrack cands : %d\n", partTrk.GetEntries() ) ;
      for ( int ti=0; ti<partTrk.GetEntries(); ti++ ) {
         H1PartSelTrack* trk = partTrk[ti] ;
         fprintf( dump_file, "   track %3d : E %9.3f  Pt %9.3f  Pz %9.3f  Eta %9.5f  Phi %9.5f",
           ti, trk->GetE(), trk->GetPt(), trk->GetPz(), trk->GetEta(), trk->GetPhi() ) ;
         if ( ti == trk_min_ti ) { fprintf( dump_file, " *** scattered electron track, dr = %9.5f", trk_min_dr ) ; }
         fprintf( dump_file, "\n" ) ;
      } // ti
   }




  //======== Stuff for radiated photons below here.


   static H1PartCandArrayPtr partCand ;
   if ( do_dump ) {
      fprintf( dump_file, "\n\n ------ Number of reconstructed particle candidates: %d\n", partCand.GetEntries() ) ;
      for ( int pi = 0; pi < partCand.GetEntries(); pi++ ) {
         H1PartCand* pcp = partCand[pi] ;
         fprintf( dump_file, "  %3d :  EMfrac %9.3f  LArfrac %9.3f  EmHadSepClassifier %9.3f  E %9.3f  Pt %9.3f  Pz %9.3f  Eta %9.5f  Phi %9.5f  IsEm %d IsCluster %d IsSelTrack %d\n",
            pi, pcp->GetEMFraction(), pcp->GetLArFraction(), pcp->GetEmHadSepClassifier(),
            pcp->GetE(), pcp->GetPt(), pcp->GetPz(), pcp->GetEta(), pcp->GetPhi(),
            pcp->IsEm(), pcp->IsCluster(), pcp->IsSelTrack() ) ;

         float cand_phi = pcp->GetPhi() ;
         float cand_eta = pcp->GetEta() ;

         float min_dr(999.) ;
         int   min_dr_mcpi(-1) ;
         static H1PartMCArrayPtr partMc ;
         for ( int mcpi = 0; mcpi < partMc.GetEntries(); mcpi++ ) {
            H1PartMC* pmc = partMc[mcpi] ;
            if ( !( pmc->GetStatus() == 0 || pmc->GetStatus() == 202 )  ) continue ;
            float mc_phi = pmc->GetPhi() ;
            float mc_eta = pmc->GetEta() ;
            float dr = calc_dr( cand_phi, mc_phi, cand_eta, mc_eta ) ;
            if ( dr < min_dr ) {
               min_dr = dr ;
               min_dr_mcpi = mcpi ;
            }
         } // mcpi
         if ( min_dr_mcpi >= 0 && min_dr < 0.05 ) {
            H1PartMC* pmc = partMc[min_dr_mcpi] ;
            fprintf( dump_file, "                         %4d     MC match  dr = %9.5f  PID = %8d  E %9.3f  Pt %9.3f  Pz %9.3f  Eta %9.5f  Phi %9.5f",
               min_dr_mcpi, min_dr, pmc->GetPDG(), pmc->GetE(), pmc->GetPt(), pmc->GetPz(), pmc->GetEta(), pmc->GetPhi() ) ;
            if ( abs(pmc->GetPDG()) == 11 && min_dr_mcpi < 4 ) fprintf( dump_file, " +++ scattered electron" ) ;
            if ( pmc->GetPDG() == 22 && pmc->GetStatus() == 202 ) fprintf( dump_file, "  *** ISR/FSR") ;
            if ( pmc->GetPDG() == 22 ) fprintf( dump_file, " photon ***" ) ;
            fprintf( dump_file, "\n") ;
         }
      } // pi
   }

   float pho_min_eta_val = 9. ;
   int   pho_min_eta_pi = -1 ;
   float se_40_em_sum(0.) ;
   int   se_40_em_n(0) ;

   for ( int pi = 0; pi < partCand.GetEntries(); pi++ ) {
      H1PartCand* pcp = partCand[pi] ;
      //-- include the scattered electron in the sums.
      if ( pcp->IsSelTrack() ) {
         float dr = calc_dr( fTreeVar.obs_e_phi, pcp->GetPhi(),   fTreeVar.obs_e_eta, pcp->GetEta() ) ;
         if ( dr < 0.05 ) {
            se_40_em_sum += pcp->GetE() ;
            se_40_em_n ++ ;
         }
      }
      if ( pcp->IsSelTrack() ) continue ;
      if ( pcp->GetEMFraction() < 1. ) continue ;
      if ( pcp->GetEta() < pho_min_eta_val ) {
         pho_min_eta_val = pcp->GetEta() ;
         pho_min_eta_pi = pi ;
      }
      float dr = calc_dr( fTreeVar.obs_e_phi, pcp->GetPhi(),   fTreeVar.obs_e_eta, pcp->GetEta() ) ;
      if ( dr < 0.4 ) {
         se_40_em_sum += pcp->GetE() ;
         se_40_em_n ++ ;
      }
   } // pi
   fTreeVar.tower_sum_40 = se_40_em_sum ;
   fTreeVar.n_towers_40 = se_40_em_n ;
   if ( pho_min_eta_pi >= 0 ) {
      fTreeVar.phi_pho_closest_to_ebeam = partCand[pho_min_eta_pi]->GetPhi() ;
      fTreeVar.eta_pho_closest_to_ebeam = partCand[pho_min_eta_pi]->GetEta() ;
      fTreeVar.e_pho_closest_to_ebeam = partCand[pho_min_eta_pi]->GetE() ;
   } else {
      fTreeVar.e_pho_closest_to_ebeam = -1. ;
      fTreeVar.eta_pho_closest_to_ebeam = 9. ;
      fTreeVar.phi_pho_closest_to_ebeam = 0. ;
   }


}//end main function



// _______________________________________________________ //
//!
//!  AnalysisEventShapes::DoControlPlotsGen()
//!
//!  This function is called by the main program
//!  during the event loop
//!
void AnalysisEventShapes::DoControlPlotsGen() {
   
  return;
}



// _______________________________________________________ //
//!
//!  AnalysisEventShapes::DoControlPlotsRec()
//!
//!  This function is called by the main program
//!  during the event loop
//!
void AnalysisEventShapes::DoControlPlotsRec() {
   
   static EventshapeTools ESTools;

   // --- event weight
   //double wgt = fRec.wgt;

   return;

}



// _______________________________________________________ //
//!
//!  AnalysisEventShapes::DoControlPlotsGenRec()
//!
//!  This function is called by the main program
//!  during the event loop
//!
void AnalysisEventShapes::DoControlPlotsGenRec() {
  return;
}



// _______________________________________________________ //
//!
//!  AnalysisEventShapes::DoCrossSectionsGenRec()
//!
//!  This function is called by the main program
//!  during the event loop
//!
//!  Note: Fill histograms, but make USE ONLY of
//!  observables stored previously in CrossSectionQuantities
//! 
void AnalysisEventShapes::DoCrossSectionsGenRec() {
   

   //Determine acceptance and purity
   //Create 3D Histos with bins in Q2, X and tau_zQ
   
}



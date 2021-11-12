#ifndef __H2020HISTMANAGER_H
#define __H2020HISTMANAGER_H

#include <map>
#include <unordered_map>
#include <iostream>
#include <string>
#include <vector>

#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include "TDirectory.h"

// _________________________________________________________ //
//!
//! H2020HistManager
//!
//! Small helper class to store histograms
//! and to fill them easily in a 1-line call
//!
// _________________________________________________________ //
class H2020HistManager { 
   
public:

   //! Constructors
   H2020HistManager(const std::string& HMname, const std::string& dirname="");
   ~H2020HistManager();

   //! Return name
   const std::string& GetName() const { return fHMname; }
   const std::string& GetDirname() const { return fDirname; }

   //! Write all histograms into current TDirectory
   void Write();

   // 1D
   template<class TH>
   TH* Get(const std::string& histname, const std::string& title, int nbinsx, double xlow, double xup, int id=9999);
   template<class TH>
   TH* Get(const std::string& histname, const std::string& title, const std::vector<double>& xbins, int id=9999);
   // 2D
   template<class TH>
   TH* Get(const std::string& histname, const std::string& title, int nbinsx, double xlow, double xup, int nbinsy, double ylow, double yup, int id=9999);
   template<class TH>
   TH* Get(const std::string& histname, const std::string& title, const std::vector<double>& xbins, const std::vector<double>& ybins, int id=9999);

   // 3D
   template<class TH>
      TH* Get(const std::string& histname, const std::string& title, int nbinsx, double xlow, double xup, int nbinsy, double ylow, double yup, int nbinsz, double zlow, double zup, int id=9999);
   template<class TH>
      TH* Get(const std::string& histname, const std::string& title, const std::vector<double>& xbins, const std::vector<double>& ybins, const std::vector<double>& zbins,int id=9999);

   // compatibility with h1oo
   template<class TH>
   TH* Get(const std::string& histname, int nbins, const double* xbins, const std::string& xtitle, const std::string& ytitle);
   template<class TH>
   TH* Get(const std::string& histname, int notused, const std::vector<double>& xbins, const std::string& xtitle, const std::string& ytitle);
   template<class TH>
   TH* Get(const std::string& histname, int nbinsx, double xlow, double xup, const std::string& xtitle, const std::string& ytitle);


protected:   
   std::string fHMname;   //!< Name of this HistManager to identify it by a string
   std::string fDirname;  //!< TDirectory to store the histograms (default: dirname==HMname)

   //std::unordered_map<const std::string*,std::unordered_map<int, std::unordered_map<char, std::map<std::string,TH1*> > > > fHistmap;
   //std::unordered_map<const std::string*, std::unordered_map<int, std::map<std::string,TH1*> > > fHistmap;
   //std::unordered_map<const std::string*, std::unordered_map<int, TH1* > > fHistmap;
   std::map<const std::string*, std::map<int, TH1* > > fHistmap;
   //std::map<const std::string*, TH1* > fHistmap;

public:
   static std::vector<double> MakeLogBinning(int n_bins, double xmin, double xmax)  { //!< helper function for compatibility with h1oo
         // return pointer to an array of bins in logscale
         // usage:  Double_t* q2bins = MakeLogBinning(50,100,10000.);
         //        ... TH1F("q2Log","q2 of the event",50,q2bins);
      std::vector<double> binning(n_bins+1);
      double delta_x_log = (std::log(xmax)-std::log(xmin))/n_bins;
      binning[0]=xmin;
      for(int i=1;i<=n_bins;++i) binning[i] = std::exp(std::log(binning[i-1])+ delta_x_log);
      return binning;
   }

};


// // ________________________________________ TH1 ______________________________________________ //
//! return histogram pointer from fHistmap;
//! Make new histogram if not existent.
template <class TH>
TH* H2020HistManager::Get (const std::string& histname, const std::string& title, int nbinsx, double xlow, double xup, int id) {
   TH1* hptr = fHistmap[&histname][id];
   //TH1* hptr = fHistmap[&histname][id][histname];
   if ( hptr == nullptr ) {
      std::string newhistname = histname;
      if ( id!=9999 ) newhistname += "_"+std::to_string(id);
      printf("New histogram: %s  %s\n",newhistname.c_str(), title.c_str());
      // fHistmap[&histname][id][back][histname] = new TH(histname.c_str(),title.c_str(),nbinsx,xlow,xup);
      // hptr = fHistmap[&histname][id][back][histname];
      // fHistmap[&histname][id][histname] = new TH(histname.c_str(),title.c_str(),nbinsx,xlow,xup);
      // hptr = fHistmap[&histname][id][histname];
      if ( histname.find("_lx") != std::string::npos ) { // log-x histogram
         std::vector<double> lbins(nbinsx+1);
         lbins[0] = xlow;
         double dlogx = (log(xup)-log(xlow))/nbinsx;
         for ( int i=1 ; i<=nbinsx ; ++i ) lbins[i] = std::exp(std::log(lbins[i-1]) + dlogx);
         fHistmap[&histname][id] = new TH(histname.c_str(),title.c_str(),nbinsx,&lbins[0]);
      }
      else {
         fHistmap[&histname][id] = new TH(histname.c_str(),title.c_str(),nbinsx,xlow,xup);
      }
      hptr = fHistmap[&histname][id];
   }
   return ((TH*)hptr);
}

//! return histogram pointer from fHistmap;
//! Make new histogram if not existent.
template <class TH>
TH* H2020HistManager::Get (const std::string& histname, const std::string& title, const std::vector<double>& xbins, int id) {
   TH1* hptr = fHistmap[&histname][id];
   if ( hptr == nullptr ) {
      std::string newhistname = histname;
      if ( id!=9999 ) newhistname += "_"+std::to_string(id);
      printf("New histogram: %s  %s\n",newhistname.c_str(), title.c_str());
      fHistmap[&histname][id] = new TH(histname.c_str(),title.c_str(),xbins.size()-1,&xbins[0]);
      hptr = fHistmap[&histname][id];
   }
   return ((TH*)hptr);
}


// ---------------- compatibility with H1Analysis --------------------- //
template <class TH>
TH* H2020HistManager::Get (const std::string& histname, int notused, const std::vector<double>& xbins, const std::string& xtitle, const std::string& ytitle){
   return Get<TH>(histname,int(xbins.size()-1),&xbins[0],xtitle,ytitle);
}

template <class TH>
TH* H2020HistManager::Get (const std::string& histname, int nbins, const double* xbins, const std::string& xtitle, const std::string& ytitle){
   const int id = 9999;
   TH1* hptr = fHistmap[&histname][id];
   if ( hptr == nullptr ) {
      printf("New histogram: %s  %s\n",histname.c_str(), histname.c_str());
      fHistmap[&histname][id] = new TH(histname.c_str(),histname.c_str(),nbins,xbins);
      hptr = fHistmap[&histname][id];
      hptr->GetXaxis()->SetTitle(xtitle.c_str());
      hptr->GetYaxis()->SetTitle(ytitle.c_str());
   }
   return ((TH*)hptr);
}

template <class TH>
TH* H2020HistManager::Get (const std::string& histname, int nbinsx, double xlow, double xup, const std::string& xtitle, const std::string& ytitle) {
   const int id = 9999;
   TH1* hptr = fHistmap[&histname][id];
   if ( hptr == nullptr ) {
      printf("New histogram: %s  %s\n", histname.c_str(), histname.c_str());
      fHistmap[&histname][id] = new TH(histname.c_str(),histname.c_str(),nbinsx,xlow,xup);
      hptr = fHistmap[&histname][id];
      hptr->GetXaxis()->SetTitle(xtitle.c_str());
      hptr->GetYaxis()->SetTitle(ytitle.c_str());
   }
   return ((TH*)hptr);
}


// ________________________________________ TH2 ______________________________________________ //
//! return histogram pointer from fHistmap;
//! Make new histogram if not existent.
template <class TH>
TH* H2020HistManager::Get (const std::string& histname, const std::string& title, int nbinsx, double xlow, double xup, int nbinsy, double ylow, double yup, int id) {
   TH1* hptr = fHistmap[&histname][id];
   if ( hptr == nullptr ) {
      std::string newhistname = histname;
      if ( id!=9999 ) newhistname += "_"+std::to_string(id);
      printf("New histogram: %s  %s\n",newhistname.c_str(), title.c_str());
      fHistmap[&histname][id] = new TH(histname.c_str(),title.c_str(),nbinsx,xlow,xup, nbinsy, ylow, yup);
      hptr = fHistmap[&histname][id];
   }
   return ((TH*)hptr);
}


//! return histogram pointer from fHistmap;
//! Make new histogram if not existent.
template <class TH>
TH* H2020HistManager::Get (const std::string& histname, const std::string& title, const std::vector<double>& xbins, const std::vector<double>& ybins, int id) {
   TH1* hptr = fHistmap[&histname][id];
   if ( hptr == nullptr ) {
      std::string newhistname = histname;
      if ( id!=9999 ) newhistname += "_"+std::to_string(id);
      printf("New histogram: %s  %s\n",newhistname.c_str(), title.c_str());
      fHistmap[&histname][id] = new TH(histname.c_str(),title.c_str(),xbins.size()-1,&xbins[0],ybins.size()-1,&ybins[0]);
      hptr = fHistmap[&histname][id];
   }
   return ((TH*)hptr);
}

// ________________________________________ TH3 ______________________________________________ //
//! return histogram pointer from fHistmap;
//! Make new histogram if not existent.
template <class TH>
TH* H2020HistManager::Get (const std::string& histname, const std::string& title, int nbinsx, double xlow, double xup, int nbinsy, double ylow, double yup, int nbinsz, double zlow, double zup, int id) {
   TH1* hptr = fHistmap[&histname][id];
   if ( hptr == nullptr ) {
      std::string newhistname = histname;
      if ( id!=9999 ) newhistname += "_"+std::to_string(id);
      printf("New histogram: %s  %s\n",newhistname.c_str(), title.c_str());
      fHistmap[&histname][id] = new TH(histname.c_str(),title.c_str(),nbinsx,xlow,xup, nbinsy, ylow, yup, nbinsz, zlow, zup);
      hptr = fHistmap[&histname][id];
   }
   return ((TH*)hptr);
}

//! return histogram pointer from fHistmap;
//! Make new histogram if not existent.
template <class TH>
TH* H2020HistManager::Get (const std::string& histname, const std::string& title, const std::vector<double>& xbins, const std::vector<double>& ybins, const std::vector<double>& zbins, int id) {
   TH1* hptr = fHistmap[&histname][id];
   if ( hptr == nullptr ) {
      std::string newhistname = histname;
      if ( id!=9999 ) newhistname += "_"+std::to_string(id);
      printf("New histogram: %s  %s\n",newhistname.c_str(), title.c_str());
      fHistmap[&histname][id] = new TH(histname.c_str(),title.c_str(),xbins.size()-1,&xbins[0],ybins.size()-1,&ybins[0],zbins.size()-1,&zbins[0]);
      hptr = fHistmap[&histname][id];
   }
   return ((TH*)hptr);
}




// _________________________________________________________ //
//!
//! HistMaster
//!
//! Singleton class to manage all histogram managers
//!
// _________________________________________________________ //
class HistMaster {
public:
   static HistMaster* Instance () {
      static HistMaster fInstance;
      return &fInstance;
   }
   ~HistMaster() { }

   //! Write a single histogram manager
   void WriteHistManager(TDirectory* dir, const std::string& HMname, const std::string& dirname) {
      if ( fHMs.count(HMname) != 0 ) {
         TDirectory* tmp = gDirectory;
         if ( dirname != "" ) {
            dir->mkdir(dirname.c_str(),dirname.c_str(),true);
            dir->cd(dirname.c_str());
         }  
         fHMs.at(HMname).Write();
         tmp->cd(); // to back
      }
      else
         std::cout<<"WARNING in HistMaster::NewHistManager! Cannot write. HistogramManager with name '"<<HMname<<"' not existent."<<std::endl;
   }
   //! write all histogram managers
   void WriteAll(TDirectory* dir) { 
      for ( auto& [HMname, HM] : fHMs ) {
         auto dirname = HM.GetDirname();
         WriteHistManager(dir,HMname,dirname);
      }
   }

   //! Get a histogram manager
   //! Instantiate a new one, of not yet existent
   H2020HistManager& GetHistManager(const std::string& HMname, const std::string& dirname = "") {
      if ( fHMs.count(HMname) == 0 ) 
         fHMs.insert({HMname, H2020HistManager(HMname,dirname)});
      // todo: possibly add a check for different dirnames here (slow?!)
      return fHMs.at(HMname);
   }
   
private:
   static HistMaster* fInstance;
   HistMaster() {}
   std::map<std::string,H2020HistManager> fHMs;
};



#endif

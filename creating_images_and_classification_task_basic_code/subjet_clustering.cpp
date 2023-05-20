
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <utility>
#include <vector>
#include <math.h>

#include <TFile.h>

#include "TROOT.h"
#include "TSystem.h"
#include "TApplication.h"

#include "TString.h"

#include "TH2.h"
#include "THStack.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TClonesArray.h"
#include "TLorentzVector.h"

#include "TCanvas.h"

#include "modules/Delphes.h"
#include "classes/DelphesClasses.h"
#include "classes/DelphesFactory.h"

#include "ExRootAnalysis/ExRootTreeReader.h"
#include "ExRootAnalysis/ExRootTreeWriter.h"
#include "ExRootAnalysis/ExRootTreeBranch.h"
#include "ExRootAnalysis/ExRootResult.h"
#include "ExRootAnalysis/ExRootUtilities.h"

#include "fastjet/PseudoJet.hh"
#include "fastjet/JetDefinition.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/contrib/Nsubjettiness.hh"
#include "fastjet/contrib/Njettiness.hh"
#include "fastjet/contrib/NjettinessPlugin.hh"
//#include "fastjet/contrib/WinnerTakeAllRecombiner.hh"

#include <fastjet/Selector.hh>
#include "fastjet/tools/Filter.hh"
#include "fastjet/tools/Pruner.hh"
#include "fastjet/tools/MassDropTagger.hh"

#include "fastjet/contrib/Nsubjettiness.hh"
#include "fastjet/contrib/Njettiness.hh"
#include "fastjet/contrib/NjettinessPlugin.hh"
//#include "fastjet/contrib/WinnerTakeAllRecombiner.hh"

#include <sstream>
#include "fastjet/contrib/EnergyCorrelator.hh"

using namespace std;
using namespace fastjet;
using namespace fastjet::contrib;


#include <fastjet/tools/JHTopTagger.hh>

//.................................................................................
//......... Random number generator .......

#include <chrono>

/* double rangen()
{
    std::mt19937_64 rng;
    // initialize the random number generator with time-dependent seed
    uint64_t timeSeed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
    std::seed_seq ss{uint32_t(timeSeed & 0xffffffff), uint32_t(timeSeed>>32)};
    rng.seed(ss);
    // initialize a uniform distribution between 0 and 1
    std::uniform_real_distribution<double> unif(0, 1);
    // ready to generate random numbers
     const int nSimulations = 1;
     double currentRandomNumber;
    for (int i = 0; i < nSimulations; i++)
    {
        currentRandomNumber = unif(rng);
   //     std::cout << currentRandomNumber << std::endl;
    }
   // return 0;
     return currentRandomNumber;
}
*///.................................................................................



int main()
//main loop begins
{     vector<PseudoJet> inputList, inputList2, inputList3;
    vector<PseudoJet> constituents;
    vector<PseudoJet> outputList,outputList2,outputList3;
    PseudoJet jetphoton,jethadron,jettrack; //track, emcal, hcal pseudojets which is pushed into inputlist
    PseudoJet jetparticle,jettowerEM,jettowerHAD; //track, emcal, hcal pseudojets which is pushed into inputlist
    
    TFile file1("/Users/apple/Delphes-3.5.0/ppaa.root","READ");
    
    gSystem->Load("libDelphes");
    
    
    // Create chain of root trees
    TChain chain("Delphes");
    chain.Add("/Users/apple/Delphes-3.5.0/ppaa.root");
    //Create object of class ExRootTreeReader
    ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
    Long64_t numberOfEntries = treeReader->GetEntries();
    //    TClonesArray *branchphoton = treeReader->UseBranch("photonisolation");//ECAL output
    TClonesArray *branchelectron = treeReader->UseBranch("Electron");//ECAL output
    TClonesArray *branchmuon = treeReader->UseBranch("Muon");//ECAL output
    TClonesArray *branchphoton = treeReader->UseBranch("Photon");
    
    
    TClonesArray *branchtower = treeReader->UseBranch("Tower");
    TClonesArray *branchtrack = treeReader->UseBranch("Track");
    TClonesArray *branchParticle = treeReader->UseBranch("Particle");
    TClonesArray *branchJet = treeReader->UseBranch("Jet");
    TClonesArray *branchGenParticle = treeReader->UseBranch("GenParticle");
    TClonesArray *branchMissingET = treeReader->UseBranch("MissingET");
    
    TClonesArray *branchEFlowTrack = treeReader->UseBranch("EFlowTrack");
    TClonesArray *branchEFlowTower = treeReader->UseBranch("EFlowPhoton");
    TClonesArray *branchEFlowNhadron = treeReader->UseBranch("EFlowHadron");
    
    // TFile* fileout=new TFile("diphotonout_CA.root","recreate");
    TFile* fileout=new TFile("ppaasjout.root","recreate");
    
    
    TTree *tree = new TTree("TreeB","Background");
    TH1F* h1 = new TH1F("h1", "histogram", 10000, 0.0, 1000);
    
    //jet and event variables
    
    Float_t tau1=0;
    Float_t higgspt=0;
    Float_t top1pt=0;
    Float_t top2pt=0;
    Float_t gluonpt=0;
    Float_t mass=0;
    Float_t edge=0;
    Float_t edge1=0;
    Float_t missingET=0;
    Float_t edge_tprime=0;
    Float_t edge_hprime=0;
    Float_t transversemass_hprime =0;
    Float_t transversemass_tprime=0;
    Float_t transversemass_tprime_real_max=0;
    Float_t transversemass_tprime_real_min=0;
    Float_t invmassmax=0;
    Float_t invmassmin=0;
    Float_t invmass1=0;
    Float_t phileading=0;
    
    Float_t invmass2=0;
    
    Float_t invmass3=0;
    
    Float_t invmasstot=0;
    Float_t thetaJ0=0;
    Float_t thetaJ1=0;
    Float_t thetaJ2=0;
    Float_t edge_real=0;
    Float_t edge_real1=0;
    Float_t min_invmass=0;
    
    Float_t mass_total=0;
    Float_t edge_real_hprime=0;
    
    Float_t tprimept=0;
    
    Float_t tau21_beta1=1;
    Float_t tau31_beta1=1;
    Float_t tau32_beta1=1;
    
    //subleading jet subjettiness
    Float_t tau1s=0;
    Float_t tau21_beta1s=1;
    Float_t tau31_beta1s=1;
    Float_t tau32_beta1s=1;
    
    
    Float_t Area=0;
    Float_t thetaJ=0;
    Float_t trackno=0;
    Float_t LambdaJ=0;
    Float_t rho=0;
    Float_t epsJ=0;
    Float_t pionpT=0;
    
    Float_t pTbyM=0;
    Float_t Rdelgamgam=100;
    Float_t Rdelphoton=0;
    
    
    Float_t ecf0=0;
    Float_t ecf1=0;
    Float_t ecf2=0;
    Float_t ecf3=0;
    Float_t ecf4=0;
    Float_t ecf5=0;
    
    Float_t ecfr0=1;
    Float_t ecfr1=1;
    Float_t ecfr2=1;
    Float_t ecfr3=1;
    Float_t ecfr4=1;
    
    
    Float_t etag1 =0;
    Float_t phig1 = 0;
    Float_t etag2 = 0;
    Float_t phig2 =0;
    
    Float_t ecfdr1=1;
    Float_t ecfdr2=1;
    Float_t ecfdr3=1;
    Float_t ecfdr4=1;
    Float_t jetno=0;
    Float_t jetno1=0;
    Float_t  pt3=0;
    Float_t  pt3a=0;
    
    Float_t deltaR1=0;
    Float_t deltaR2=0;
    Float_t deltaR3=0;
    
    
    Float_t deltaR = 0;
    
    
    //photon variables
    Float_t diphotonmin=0;
    Float_t photonnumber=0;
    Float_t photonnumber1=0;
    Float_t photonpT0=0;
    Float_t photonpT1=0;
    Float_t photonpT2=0;
    Float_t mass1=0;
    Float_t mass2=0;
    Float_t mass3=0;
    Float_t bjet=0;
    Float_t bmeson=0;
    
    Float_t etaphoton1=0;
    Float_t etaphoton2=0;
    Float_t deltaeta1=0;
    Float_t deltaeta2=0;
    
    Float_t jetm1=0;
    Float_t jetm2=0;
    Float_t deltaetasubjet=0;
    Float_t ptsubjetplot=0;
    Float_t costheta=0;
    Float_t leadingjet=0;
    Float_t jet_number =0;
    Float_t higgsjet=0;
    
    Float_t diphotonminsubjet=0;
    
    Float_t transverse_momentum=0;
    tree->Branch("transverse_momentum", &transverse_momentum, "transverse_momentum/F");
    Float_t numtracks=0;
    tree->Branch("numtracks", &numtracks, "numtracks/F");
  
    ofstream variables7;
    ///  variables7.open("/Users/abhishekiyer/Dropbox/chargedvsneutral/HLfiles/zpmu14_4.txt",  std::ios::app);
    variables7.open("/Users/apple/Documents/softwares/iyer/ppaasjout.txt",  std::ios::app);
    
    
    cout<<numberOfEntries<<"\n";
    
    TLorentzVector momentum;
    TObject *object;
    Track *track;
    Tower *tower;
    GenParticle *particle1;
    int count1=0;
    int count2=0;
    int count3=0;
    int count4=0;
    int count=0;
    int count4a=0;
    int eventcount=0;
    int eventcount1=0;
    
    
    
    /// To Collect and Expport some data in csv form
    fstream fout;
    /*   fout.open("bbbargamma30ptout.csv", ios::out | ios::app); ///
     /// Adding Columns to data.csv
     fout << "invmassphoton" << ","
     << "Einv" <<
     
     "\n"; */
    ///
    Float_t  totalnumjets =0;
    
    
    //event loop begins
    for(Int_t entry=0; entry < numberOfEntries; ++entry)
    { treeReader->ReadEntry(entry);
        
        Float_t pxtot=0;
        Float_t pytot=0;
       // int tracksize = 0;
        Float_t ptl1=0;
        Float_t ptl2=0;
        Float_t phil1=0;
        Float_t phil2=0;
        Float_t etal1=0;
        Float_t etal2=0;
        
        Float_t ljeta = 0;
        Float_t ljphi =0;
        Float_t sjeta = 0;
        Float_t sjphi =0;
        
        Float_t coneta = 0;
        Float_t conphi = 0;
        Float_t conEnergy = 0;
        Float_t transverse_energy = 0;
        
        Float_t sjconeta = 0;
        Float_t sjconphi = 0;
        Float_t sjconEnergy = 0;
        Float_t sjtransverse_energy = 0;
        
        
        Float_t numjets=0;
        
        double ptvector [100];
        double phivector[100];
        double etavector[100];
        
        Float_t pxelectron=0;
        Float_t pyelectron=0;
        Float_t pzelectron=0;
        Float_t Eelectron=0;
        Float_t Ptelectron=0;
        
        Float_t pxmuon=0;
        Float_t pymuon=0;
        Float_t pzmuon=0;
        Float_t Emuon=0;
        Float_t Ptmuon=0;
        
        Float_t pnx=0;
        Float_t pny=0;
        Float_t pnz=0;
        Float_t pne=0;
        
        Float_t pxtop2=0;
        Float_t pytop2=0;
        Float_t pztop2=0;
        Float_t Etop2=0;
        Float_t etatop2=0;
        
        Float_t pxtop1=0;
        Float_t pytop1=0;
        Float_t pztop1=0;
        Float_t Etop1=0;
        
        Float_t pxgluon=0;
        Float_t pygluon=0;
        Float_t pzgluon=0;
        Float_t Egluon=0;
        Float_t zprimemass=0;
        
        int k=0;
        Float_t pxtau[20]={};
        Float_t pytau[20]={} ;
        Float_t pztau[20]={};
        Float_t energytau[20]={};
        Float_t etatau[20]={};
        Float_t phitau[20]={};
        Float_t pttau[20]={};
        int counttaujet=0;
        
        
        Float_t missET=0;
       
        //   numjets = branchJet->GetEntriesFast();
        //  Float_t numparticles = branchParticle->GetEntriesFast();
        //   cout<<entry<<"\t"<<numjets<<"\n";
        
   
        
        for(int ii=0; ii < branchtrack -> GetEntriesFast(); ++ii)
                {
                    Track *track2=(Track*) branchtrack->At(ii);
            
                    Float_t pt1 = track2->PT;
                    Float_t eta1 = track2->Eta;
                    Float_t phi1 = track2->Phi;
                    
                    
                    //defining 4 momentum variables
                    Float_t eT1 = sqrt(pow(pt1,2));
                    Float_t pe1 = pt1*cosh(eta1);
                    Float_t pz1 = pt1*sinh(eta1);
                    Float_t px1 = pt1*cos(phi1);
                    Float_t py1 = pt1*sin(phi1);
                    
                    //         re-scale the track four momentum
                  //  if (pt1 >=.1)
                    if (pt1 >=2)
                    {  //tracks = tracks+1;
                        
                        Float_t eps = pow(10,5);
                        
                        jettrack = PseudoJet(px1/eps, py1/eps, pz1/eps, pe1/eps);
                        jettrack.set_user_index(0);
                        inputList.push_back(jettrack);
                    }
                }//track branch ends here
                
                //   cout<<branchtrack->GetEntriesFast()<<"\n";
                //-----------------------------------------------------------------------------------------
                
                //Pusing tower elements into input list
                
                for(int i=0; i < branchtower -> GetEntriesFast(); ++i)
                {
                    Tower *tower2=(Tower*) branchtower->At(i);
                    
                    Float_t eta1 = tower2->Eta;
                    Float_t phi1 = tower2->Phi;
                    Float_t  Ehad = tower2->Ehad;
                    Float_t Eem = tower2->Eem;
                   
                    
                    
                    //defining 4 momentum variables for ecal
                    Float_t eTecal = Eem/cosh(eta1);
                    Float_t peecal = Eem;
                    Float_t pzecal = eTecal*sinh(eta1);
                    Float_t pxecal =  eTecal*cos(phi1);
                    Float_t pyecal =  eTecal*sin(phi1);

                        if(eTecal>0.1)
                        {
                            jetphoton = PseudoJet(pxecal, pyecal, pzecal, peecal);
                            jetphoton.set_user_index(1);
                            inputList.push_back(jetphoton);
                            
                            
                        }
                        //defining 4 momentum variables for hcal
                        Float_t eThcal = Ehad/cosh(eta1);
                        Float_t pehcal = Ehad;
                        Float_t pzhcal = eThcal*sinh(eta1);
                        Float_t pxhcal =  eThcal*cos(phi1);
                        Float_t pyhcal =  eThcal*sin(phi1);
                        
                        
                        if(eThcal>0.5)
                        {
                            jethadron = PseudoJet(pxhcal, pyhcal, pzhcal, pehcal);
                            jethadron.set_user_index(2);
                            inputList.push_back(jethadron);
                         
                        }
                   
                }//hadron branch ends here
        
        
        for (unsigned int jj=0; jj<branchParticle -> GetEntriesFast(); ++jj)
        {
            //GenParticle *genparticle=(GenParticle*) branchParticle->At(jj);
            
            
            particle1 =(GenParticle*) branchParticle->At(jj);
            
          //  cout << "checkiing if gen particle is working" <<"\n";
   
   
   
       if ((particle1->PID==9000005)&& particle1->Status==22)
       //  if ((particle1->PID==5) && particle1->Status==21)
        {
            cout << "checkiing if gen particle is working" <<"\n";
       Float_t  pxd1=particle1->Px;
       
       Float_t  pyd1=particle1->Py;
       
       Float_t   pzd1=particle1->Pz;
       
       Float_t  ped1=particle1->E;
       
       //       phivector[jj]= genparticle->Phi;
       
   // transverse_momentum = particle1->PT;
       
       Float_t Einvd=ped1;
       Float_t pxinvd=pxd1;
       Float_t pyinvd=pyd1;
       Float_t pzinvd=pzd1;
       
       
        }
   }
   
        double R = 0.4;
                JetDefinition jet_def1(antikt_algorithm, R);
               //non area clustering
                ClusterSequence sequence1(inputList,jet_def1);
        
        outputList.clear();
        outputList = sorted_by_pt(sequence1.inclusive_jets(2.0));
        
            //  cout << "output list size is " << outputList.size() << "\n" ;
        
        numjets = outputList.size();
           
           //  cout << " pt pseudo-rapidity phi Energy" << endl;
             
             for (unsigned i = 0; i < outputList.size(); i++) {
            //     cout << "jet " << i << ": "<< outputList[0].perp() << "\n ";
            //     << outputList1[i].pseudorapidity() << " " << outputList1[i].phi() << " "<< outputList1[i].E()<< endl;
           //    Float_t  etag2 = outputList[0].eta();
          //   Float_t phig2 = outputList[0].phi();
              //   variables7<<setprecision(3)<<outputList[0].phi()<< "\t" << outputList[0].eta() << "\n";
                    
                // etag2 = 12;
                // phig2 = 23;
                 transverse_momentum = outputList[0].perp();
               
                 for (unsigned j = 0; j < constituents.size(); j++) {
                    // cout << " constituent " << j << "’s pt: "<< constituents[j].perp() << endl;
                    // for(unsigned t = 0; t < constituents.size(); t++)
                  //   {
                  //      if(jettrack.set_user_index == 0)
                  //      {cout << "something" << "\n" ;}
                     //}
                   //
                 
                    
                 }
               
             }
        
     //   Float_t deltaR = sqrt(abs(pow((phig1-phig2),2)+pow((etag1-etag2),2)));
     //   cout << "delta R value is " << deltaR << "\n";
        int tracksize = 0;
        
      if (outputList.size() > 0)
      {
          
          ljeta = outputList[0].eta();
          ljphi = outputList[0].phi();
          //  cout << "ljeta and ljphi " << ljeta << "\t" << ljphi <<  "\n" ;
          
          vector<PseudoJet> constituents = outputList[0].constituents();
          //    cout << "constituents size is " << constituents.size()<< "\n";
          
          
          
          vector<PseudoJet>  subjet;
          int Nprefilter=5;
          double R1=1.0;
          //   JetDefinition jet_def_sub(cambridge_algorithm,R1);
          JetDefinition jet_def_sub(kt_algorithm,R1);
          
          //non area clustering
          ClusterSequence clust_seq(constituents, jet_def_sub);
          subjet.clear();
          subjet = sorted_by_pt(clust_seq.inclusive_jets(2.0));
          
          for (unsigned mm = 0; mm < subjet.size(); mm++) {
              
              transverse_momentum = subjet[0].perp();
              
          }
          
          
          if(subjet.size()>0)
          {
              
              sjeta = subjet[0].eta();
              sjphi = subjet[0].phi();
              
              vector<PseudoJet> sjconstituents = subjet[0].constituents();
              cout<<entry<<"\t"<<sjconstituents.size()<<"\n";
              
              for (int nn = 0; nn < sjconstituents.size(); nn++)
              {
                  sjconeta = sjconstituents[nn].eta();
                  sjconphi = sjconstituents[nn].phi();
                  sjconEnergy = sjconstituents[nn].E();
                  sjtransverse_energy = sjconEnergy/cosh(sjconeta);
                  
                  if( sjconstituents[nn].user_index() == 0)
                  {
                      tracksize = tracksize + 1;
                      //  cout << entry << "checking  " << tracksize << "\n";
                      
                  }
                  Float_t deltaetasj =0;
                  Float_t deltaphisj = 0;
                  deltaetasj = sjeta - sjconeta;
                  deltaphisj = sjphi - sjconphi;
                  cout << "deltaetaj and deltaphij and energy" << "\t" << deltaetasj << "\t" << deltaphisj <<"\t" << sjtransverse_energy <<  "\n" ;
                  //  cout<<entry<<"\t"<<constituents.size()<<"\n";
                  variables7<<setprecision(3)<<deltaetasj<< "\t" << deltaphisj<< "\t" << sjtransverse_energy << "\n";
              }
              
              
          }
      }
      
          numtracks = tracksize;
   
      
        
             inputList.clear();

   
     
           
        // }// branch jet loop ends here inside event loop
        tree ->Fill();
       
    } //event loop ends here
 //   cout<<"Number of di-`electron' events"<<"\t"<<eventcount<<"\n";
    tree->Write();
    delete fileout;
    variables7.close();
    return 0;
    
} // main loop ends here

        
       

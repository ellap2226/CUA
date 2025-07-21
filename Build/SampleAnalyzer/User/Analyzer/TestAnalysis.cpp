#include "SampleAnalyzer/User/Analyzer/TestAnalysis.h"
#include <TLorentzVector.h>
#include "TVector3.h"
#include <TFile.h>
#include <cmath>

using namespace MA5;
using namespace std;

// -----------------------------------------------------------------------------
// Initialize
// function called one time at the beginning of the analysis
// -----------------------------------------------------------------------------
bool TestAnalysis::Initialize(const MA5::Configuration& cfg, const std::map<std::string,std::string>& parameters)
{
  cout << "BEGIN Initialization" << endl;
  // initialize variables, histos
  cout << "END   Initialization" << endl;
  //Histogram to fill
  ptaamtotal = new TH1F("ptaamtotal","p_{T}^{aa}/m_{total}",100,0.,1.);
  ptbbmtotal = new TH1F("ptbbmtotal","p_{T}^{bb}/m_{total}",100,0.,1.);
  ptamaal = new TH1F("ptamaal","p_{T}^{a}/m_{aa} leading",100,0.,2.);
  ptamaasl = new TH1F("ptamaasl","p_{T}^{a}/m_{aa} subleading",100,0.,2.);
  ptbmbbl = new TH1F("ptbmbbl","p_{T}^{b}/m_{bb} leading",100,0.,2.);
  ptbmbbsl = new TH1F("ptbmbbsl","p_{T}^{b}/m_{bb} subleading",100,0.,2.);
  minAngularDist = new TH1F("minAngularDist","Minimum Angular Distance (photon&quark)",100,0.,3.);
  cosThetaY = new TH1F("cosThetaY","costheta^{*} Y",100,0.,3.5);
  cosThetaH = new TH1F("cosThetaH","costheta^{*} H",100,0.,3.5);
  return true;
}

// -----------------------------------------------------------------------------
// Finalize
// function called one time at the end of the analysis
// -----------------------------------------------------------------------------
void TestAnalysis::Finalize(const SampleFormat& summary, const std::vector<SampleFormat>& files)
{
  cout << "BEGIN Finalization" << endl;
  TFile* outFile = new TFile("first_trial.root","recreate");
  ptaamtotal->Write();
  ptbbmtotal->Write();
  ptamaal->Write();
  ptamaasl->Write();
  ptbmbbl->Write();
  ptbmbbsl->Write();
  minAngularDist->Write();
  cosThetaY->Write();
  cosThetaH->Write();
  outFile->Close();
  // saving histos
  cout << "END   Finalization" << endl;
}

// -----------------------------------------------------------------------------
// Execute
// function called each time one event is read
// -----------------------------------------------------------------------------
bool TestAnalysis::Execute(SampleFormat& sample, const EventFormat& event)
{
  // ***************************************************************************
  // Example of analysis with generated particles
  // Concerned samples : LHE/STDHEP/HEPMC
  // ****************************************************************************
  if (event.mc()!=0)
  {
    cout << "*---------------NEW EVENT-------------------*" << endl;
    
    //PT
    TLorentzVector pho1(0,0,0,0);
    TLorentzVector pho2(0,0,0,0);
    TLorentzVector b(0,0,0,0);
    TLorentzVector ab(0,0,0,0);
    TLorentzVector x(0,0,0,0);
    TLorentzVector y(0,0,0,0);
    TLorentzVector h(0,0,0,0);

    //Total Mass
    //TLorentzVector totalMass(0,0,0,0);

    // Initial state
    for (MAuint32 i=0;i<event.mc()->particles().size();i++)
    {
      const MCParticleFormat& part = event.mc()->particles()[i];

      cout << "----------------------------------" << endl;
      cout << "MC particle" << endl;
      cout << "----------------------------------" << endl;

      // display index particle
      cout << "index=" << i+1;

      // display the status code
      cout << " Status Code=" << part.statuscode() << endl;
      if (PHYSICS->Id->IsInitialState(part)) cout << " (Initial state) ";
      else if (PHYSICS->Id->IsFinalState(part)) cout << " (Final state) ";
      else cout << " (Intermediate state) ";
      cout << endl;

    // Find and set x particle
    if(part.pdgid()==45){
	      x.SetPxPyPzE(part.px(),
		    part.py(),
		    part.pz(),
		    part.e());
    } else if(part.pdgid()==35){
        y.SetPxPyPzE(part.px(),
		    part.py(),
		    part.pz(),
		    part.e());
    } else if(part.pdgid()==25){
        h.SetPxPyPzE(part.px(),
		    part.py(),
		    part.pz(),
		    part.e());
    }

      if(!PHYSICS->Id->IsFinalState(part))continue;

      // pdgid
      cout << "pdg id=" << part.pdgid() << endl;

      // Set momentums for pho1, pho2, b, ab, x
    if(part.pdgid()==22 && pho1.E()<1e-8){
	      pho1.SetPxPyPzE(part.px(),
		    part.py(),
		    part.pz(),
		    part.e());
    } else if(part.pdgid()==22 && pho2.E()<1e-8){
        pho2.SetPxPyPzE(part.px(),
		    part.py(),
		    part.pz(),
		    part.e());
    } else if(part.pdgid()==5 && b.E()<1e-8){
	      b.SetPxPyPzE(part.px(),
		    part.py(),
		    part.pz(),
		    part.e());
    } else if(part.pdgid()==-5 && ab.E()<1e-8){
        ab.SetPxPyPzE(part.px(),
		    part.py(),
		    part.pz(),
		    part.e());
    }

      // display kinematics information
      cout << "px=" << part.px()
                << " py=" << part.py()
                << " pz=" << part.pz()
                << " e="  << part.e()
                << " m="  << part.m() << endl;
      cout << "pt=" << part.pt() 
                << " eta=" << part.eta() 
                << " phi=" << part.phi() << endl;


      // display particle mother id
      if (part.mothers().empty()) 
      {
        cout << "particle with no mother." << endl;
      }
      else
      {
        std::cout << "particle coming from the decay of";
        for(MAuint32 j=0;j<part.mothers().size();j++)
        {
          const MCParticleFormat* mother = part.mothers()[j];
          cout << " " << mother->pdgid();
        }
        std::cout << "." << endl;
      }

    }
    
    cout << "----------------PT(aa)/m_total------------------" << endl;

    // Total Transverse Momentum
    double pTa = (pho1 + pho2).Pt();

    // Total Mass
    double totalMass = pho1.M()+pho2.M()+b.M()+ab.M();
    double fourBodyMass = (pho1+pho2+b+ab).M();
    
    // Final Number
    double finalValueaat = pTa/fourBodyMass;
    
    // Print outs
    cout << "Pho1 Mass = " << pho1.M() << endl;
    cout << "Pho2 Mass = " << pho2.M() << endl;
    cout << "b Mass = " << b.M() << endl;
    cout << "ab Mass = " << ab.M() << endl;
    cout << "Total Mass = " << totalMass << endl;
    cout << "Four body mass = " << fourBodyMass << endl;
    cout << "PT Photons = " << pTa << endl;
    cout << "Final Value = " << finalValueaat << endl;

    ptaamtotal->Fill(finalValueaat);

    cout << endl;
    cout << "----------------PT(bb)/m_total------------------" << endl;

    // Total Transverse Momentum
    double pTb = (b + ab).Pt();
    
    // Final Number
    double finalValuebbt = pTb/fourBodyMass;
    
    // Print outs
    cout << "Pho1 Mass = " << pho1.M() << endl;
    cout << "Pho2 Mass = " << pho2.M() << endl;
    cout << "b Mass = " << b.M() << endl;
    cout << "ab Mass = " << ab.M() << endl;
    cout << "Total Mass = " << totalMass << endl;
    cout << "Four body mass = " << fourBodyMass << endl;
    cout << "PT quarks = " << pTb << endl;
    cout << "Final Value = " << finalValuebbt << endl;
    
    ptbbmtotal->Fill(finalValuebbt);

    cout << "----------------PT(a)/m_aa leading------------------" << endl;
    
    // Photon Mass
    double photonM = (pho1 + pho2).M();

    // Find leading photon
    double leadingPho;
    double subleadingPho;

    if(pho1.Pt()>pho2.Pt()){
      leadingPho = pho1.Pt();
      subleadingPho = pho2.Pt();
    } else {
      leadingPho = pho2.Pt();
      subleadingPho = pho1.Pt();
    }

    // Final Number
    double finalValueaaal = leadingPho/photonM;
    
    // Print outs
    cout << "Photon Mass = " << photonM << endl;
    cout << "Final Value = " << finalValueaaal << endl;
    
    ptamaal->Fill(finalValueaaal);

    cout << "----------------PT(a)/m_aa subleading------------------" << endl;
    
    // Final Number
    double finalValueaaasl = subleadingPho/photonM;
    
    // Print outs
    cout << "Photon Mass = " << photonM << endl;
    cout << "Final Value = " << finalValueaaasl << endl;
    
    ptamaasl->Fill(finalValueaaasl);

    cout << "----------------PT(b)/m_bb leading------------------" << endl;

    // Photon Mass
    double qM = (b+ab).M();

    // Find leading quark
    double leadingb;
    double subleadingb;

    if(b.Pt()>ab.Pt()){
      leadingb = b.Pt();
      subleadingb = ab.Pt();
    } else {
      leadingb = ab.Pt();
      subleadingb = b.Pt();
    }


    // Final Number
    double finalValuebbbl = leadingb/qM;
    
    // Print outs
    cout << "Quark Mass = " << qM << endl;
    cout << "Final Value = " << finalValueaaasl << endl;
    
    ptbmbbl->Fill(finalValuebbbl);

    cout << "----------------PT(b)/m_bb subleading------------------" << endl;

    // Final Number
    double finalValuebbbsl = subleadingb/qM;
    
    // Print outs
    cout << "Quark Mass = " << qM << endl;
    cout << "Final Value = " << finalValueaaasl << endl;
    
    ptbmbbsl->Fill(finalValuebbbsl);

    cout << "----------------Minimum Angular Distance (photon & quark)------------------" << endl;

    double comboArray [4] = {pho1.DeltaR(b), pho1.DeltaR(ab), pho2.DeltaR(b), pho2.DeltaR(ab)};

    double minDist = comboArray[0];

    for(int combo = 1; combo < 4; combo++){
      if(minDist>comboArray[combo]){
        minDist = comboArray[combo];
      }
    }
    
    // Print Outs
    /*
    cout << "MinDist" << minDist << endl;
    cout << "Full Array" << endl;
    for (int i = 0; i < 4; i ++){
      cout << " " << comboArray[i] << " ";
    }
      */

    minAngularDist->Fill(minDist);
///*
    cout << "----------------Cos Theta------------------" << endl;

    TVector3 boostToRest = -x.BoostVector();

    // Copy of original y and h
    TLorentzVector yInRest = y;
    TLorentzVector hInRest = h;

    // Boost both Y and H to rest
    yInRest.Boost(boostToRest);
    hInRest.Boost(boostToRest);

    // Get direction of Y and H (in Rest) and X
    TVector3 yDirection = yInRest.Vect();
    TVector3 hDirection = hInRest.Vect();
    TVector3 xDirection = x.Vect();

    double costhetaY = yDirection.Angle(xDirection);
    double costhetaH = hDirection.Angle(xDirection);

    cout << "cos theta Y and X " << costhetaY << endl;
    cout << "cos theta H and X " << costhetaH << endl;

    /*
    cout << "y and h in x rest frame" << endl;
    cout << endl;

    cout << "y in rest = " << endl;
    yInRest.Print();
    cout << endl;
    cout << "y original = " << endl;
    y.Print();
    cout << endl;
    cout << "h in rest = " << endl;
    hInRest.Print();
    cout << endl;
    cout << "h original = " << endl;
    h.Print();
    cout << endl;
    */

    cosThetaY->Fill(costhetaY);
    cosThetaH->Fill(costhetaH);

//*/
    cout << endl;
    cout << endl;

  return true;
  }

  return true;
}

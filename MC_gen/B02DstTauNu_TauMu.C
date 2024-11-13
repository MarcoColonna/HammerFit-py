
/***********************************************************************
* Copyright 1998-2020 CERN for the benefit of the EvtGen authors       *
*                                                                      *
* This file is part of EvtGen.                                         *
*                                                                      *
* EvtGen is free software: you can redistribute it and/or modify       *
* it under the terms of the GNU General Public License as published by *
* the Free Software Foundation, either version 3 of the License, or    *
* (at your option) any later version.                                  *
*                                                                      *
* EvtGen is distributed in the hope that it will be useful,            *
* but WITHOUT ANY WARRANTY; without even the implied warranty of       *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the        *
* GNU General Public License for more details.                         *
*                                                                      *
* You should have received a copy of the GNU General Public License    *
* along with EvtGen.  If not, see <https://www.gnu.org/licenses/>.     *
***********************************************************************/

//
//  Sample test program for running EvtGen
//

#define EVTGEN_HEPMC3 
#define EVTGEN_EXTERNAL

#include "EvtGen/EvtGen.hh"

#include "EvtGenBase/EvtAbsRadCorr.hh"
#include "EvtGenBase/EvtDecayBase.hh"
#include "EvtGenBase/EvtHepMCEvent.hh"
#include "EvtGenBase/EvtMTRandomEngine.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtParticleFactory.hh"
#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtRandom.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenBase/EvtSimpleRandomEngine.hh"

#include "TH1D.h"
#include "TFile.h"
#include "TRandom.h"
#include "TRandom3.h"
#include "TTree.h"
#include "TLorentzVector.h"

#ifdef EVTGEN_EXTERNAL
#include "EvtGenExternal/EvtExternalGenList.hh"
#include "EvtGenExternal/EvtExternalGenFactory.hh"
#endif

#include <iostream>
#include <list>
#include <string>

// Random smearing
TRandom3 smear, kill;

struct event {
    std::vector <double> p_B0, p_Dst, p_Dst_D0, p_D0_pi, p_D0_K, p_Dst_pi, p_tau, p_tau_mu, p_tau_nu1, p_tau_nu2, p_nu;
    double B0_PE, Dst_PE, Dst_D0_PE, D0_pi_PE, D0_K_PE, Dst_pi_PE, tau_PE, tau_mu_PE, tau_nu1_PE, tau_nu2_PE, nu_PE;
    double B0_PX, Dst_PX, Dst_D0_PX, D0_pi_PX, D0_K_PX, Dst_pi_PX, tau_PX, tau_mu_PX, tau_nu1_PX, tau_nu2_PX, nu_PX;
    double B0_PY, Dst_PY, Dst_D0_PY, D0_pi_PY, D0_K_PY, Dst_pi_PY, tau_PY, tau_mu_PY, tau_nu1_PY, tau_nu2_PY, nu_PY;
    double B0_PZ, Dst_PZ, Dst_D0_PZ, D0_pi_PZ, D0_K_PZ, Dst_pi_PZ, tau_PZ, tau_mu_PZ, tau_nu1_PZ, tau_nu2_PZ, nu_PZ;
    int B0_ID, Dst_ID, Dst_D0_ID, D0_pi_ID, D0_K_ID, Dst_pi_ID, tau_ID, tau_mu_ID, tau_nu1_ID, tau_nu2_ID, nu_ID;
    //std::vector <double> p_B0_tag, p_Dst_tag, p_Dst_D0_tag, p_Dst_pi_tag, p_D0_K_tag, p_D0_pi_tag, p_tau_tag, p_nu_tag;
    // Spectral moments we want
   // double El, MX, q2;
    // Smeared versions
   // double El_reco, MX_reco, q2_reco;
    // alternative Smeared versions
   // double El_reco2, MX_reco2, q2_reco2;    
    // Lepton in barrel
   // bool inAcceptance;
    double mmiss2,q2,El;
};

TLorentzVector SmearTLV(TLorentzVector in, double res );

TLorentzVector LoopOverDaughters(EvtParticle *part);

bool KeepDaughter(EvtParticle *part);
TLorentzVector SmearDaughter(EvtParticle *part);
TLorentzVector AddBeamBackground();


void DumpResults(std::vector <event> v_evt);

int main( int /*argc*/, char** /*argv*/ )
{
    EvtParticle* parent( 0 );
    // Define the random number generator
    EvtRandomEngine* eng = 0;

    eng = new EvtMTRandomEngine();

    EvtRandom::setRandomEngine( eng );

    //    std::string fname = "tauola";
    std::string fname = "B02DstTauNu_TauMu";
    // Set the Tauola external generator
    //EvtAbsExternalGen* eng = 0;

    //EvtExternalGenFactory* externalGenerators = EvtExternalGenFactory::getInstance();

    // We are using TAUOLA physics codes in the decay.dec file(s).
    //bool convertPhysCode(true);
    //externalGenerators->defineTauolaGenerator(convertPhysCode);

    //eng = EvtExternalGenFactory::getInstance()->getGenerator(EvtExternalGenFactory::TauolaGenId);

    EvtAbsRadCorr* radCorrEngine = 0;
    std::list<EvtDecayBase*> extraModels;

    #ifdef EVTGEN_EXTERNAL
        
        EvtExternalGenList genList(true, "", "gamma", true);
        extraModels = genList.getListOfModels();
    
    #endif

    //Initialize the generator - read in the decay table and particle properties
    EvtGen myGenerator( "./"+fname+".dec", "./evt.pdl", eng, radCorrEngine,
                        &extraModels );
    //If I wanted a user decay file, I would read it in now.
    myGenerator.readUDecay("./"+fname+".dec");
    static EvtId PS = EvtPDL::getId( std::string( "B0" ) );

    int nEvents( 10000 );
    //int nEvents( 200000 );
    //    int nEvents( 1000000 );
    
    std::vector <event> v_events;
    
    // Loop to create nEvents, starting from an Upsilon(4S)
    int i;
    
    for ( i = 0; i < nEvents; i++ ) {
        if(i%10000 == 0) std::cout << "Event number " << i << std::endl;
        
        // Set up the parent particle
        EvtVector4R pInit( EvtPDL::getMass( PS ), 0.0, 0.0, 0.0 );
        parent = EvtParticleFactory::particleFactory( PS, pInit );
        //parent->setVectorSpinDensity(); << why is this problematic?
        
        // Generate the event
        myGenerator.generateDecay( parent );
        
        int nDaug = parent->getNDaug();
        // Event Information
        event ev;

        // Loop over the daughters of PS4
        //TLorentzVector tl_B0(parent->getP4Lab().get(1), parent->getP4Lab().get(2), parent->getP4Lab().get(3), parent->getP4Lab().get(0));
        std::vector<TLorentzVector> v_tl_D;
        std::vector<std::vector<TLorentzVector>> vv_tl_DD;
        
        std::vector<TLorentzVector> v_tl_DDD;

        for ( int iDaug = 0; iDaug < nDaug; iDaug++ ) {
            std::vector<TLorentzVector> v_tl_DD_temp;
            int nDDaug = parent->getDaug(iDaug)->getNDaug();
            for(int jDaug = 0; jDaug<nDDaug; jDaug++){
                if(iDaug!=0 || jDaug!=0){
                    TLorentzVector tl_DD_temp(parent->getDaug(iDaug)->getDaug(jDaug)->getP4Lab().get(1),
                                              parent->getDaug(iDaug)->getDaug(jDaug)->getP4Lab().get(2),
                                              parent->getDaug(iDaug)->getDaug(jDaug)->getP4Lab().get(3),
                                              parent->getDaug(iDaug)->getDaug(jDaug)->getP4Lab().get(0));
                    v_tl_DD_temp.push_back(tl_DD_temp);
                }
                else{
                    TLorentzVector tl_DDD_1_temp(parent->getDaug(iDaug)->getDaug(jDaug)->getDaug(0)->getP4Lab().get(1),
                                                 parent->getDaug(iDaug)->getDaug(jDaug)->getDaug(0)->getP4Lab().get(2),
                                                 parent->getDaug(iDaug)->getDaug(jDaug)->getDaug(0)->getP4Lab().get(3),
                                                 parent->getDaug(iDaug)->getDaug(jDaug)->getDaug(0)->getP4Lab().get(0));

                    TLorentzVector tl_DDD_2_temp(parent->getDaug(iDaug)->getDaug(jDaug)->getDaug(1)->getP4Lab().get(1),
                                                 parent->getDaug(iDaug)->getDaug(jDaug)->getDaug(1)->getP4Lab().get(2),
                                                 parent->getDaug(iDaug)->getDaug(jDaug)->getDaug(1)->getP4Lab().get(3),
                                                 parent->getDaug(iDaug)->getDaug(jDaug)->getDaug(1)->getP4Lab().get(0));
                    
                    TLorentzVector tl_DD_temp = tl_DDD_1_temp + tl_DDD_2_temp;
                    
                    v_tl_DD_temp.push_back(tl_DD_temp);

                    v_tl_DDD.push_back(tl_DDD_1_temp);
                    v_tl_DDD.push_back(tl_DDD_2_temp);
                }
            }
            vv_tl_DD.push_back(v_tl_DD_temp);
            if(iDaug == 0){
                TLorentzVector tl_D_temp = v_tl_DD_temp[0] + v_tl_DD_temp[1];
                v_tl_D.push_back(tl_D_temp);
            }
            if(iDaug == 1){ 
                TLorentzVector tl_D_temp = v_tl_DD_temp[0] + v_tl_DD_temp[1] + v_tl_DD_temp[2];
                v_tl_D.push_back(tl_D_temp);
            }
            if(iDaug == 2){
                TLorentzVector tl_D_temp(parent->getDaug(iDaug)->getP4Lab().get(1), 
                                         parent->getDaug(iDaug)->getP4Lab().get(2), 
                                         parent->getDaug(iDaug)->getP4Lab().get(3), 
                                         parent->getDaug(iDaug)->getP4Lab().get(0));
                v_tl_D.push_back(tl_D_temp);
            }
        }    

        TLorentzVector tl_B0 = v_tl_D[0] + v_tl_D[1] + v_tl_D[2];

        ev.p_B0.push_back(tl_B0.Px());
        ev.p_B0.push_back(tl_B0.Py());
        ev.p_B0.push_back(tl_B0.Pz());
        ev.p_B0.push_back(tl_B0.E());
        ev.p_B0.push_back(tl_B0.M());

        ev.B0_PE = tl_B0.E();
        ev.B0_PX = tl_B0.Px();
        ev.B0_PY = tl_B0.Py();
        ev.B0_PZ = tl_B0.Pz();

        ev.B0_ID = parent->getPDGId();

        // B0->DstTauNu

        ev.p_Dst.push_back(v_tl_D[0].Px());
        ev.p_Dst.push_back(v_tl_D[0].Py());
        ev.p_Dst.push_back(v_tl_D[0].Pz());
        ev.p_Dst.push_back(v_tl_D[0].E());
        ev.p_Dst.push_back(v_tl_D[0].M());

        ev.Dst_PE = v_tl_D[0].E();
        ev.Dst_PX = v_tl_D[0].Px();
        ev.Dst_PY = v_tl_D[0].Py();
        ev.Dst_PZ = v_tl_D[0].Pz();

        ev.Dst_ID = parent->getDaug(0)->getPDGId();

        ev.p_tau.push_back(v_tl_D[1].Px());
        ev.p_tau.push_back(v_tl_D[1].Py());
        ev.p_tau.push_back(v_tl_D[1].Pz());
        ev.p_tau.push_back(v_tl_D[1].E());
        ev.p_tau.push_back(v_tl_D[1].M());

        ev.tau_PE = v_tl_D[1].E();
        ev.tau_PX = v_tl_D[1].Px();
        ev.tau_PY = v_tl_D[1].Py();
        ev.tau_PZ = v_tl_D[1].Pz();

        ev.tau_ID = parent->getDaug(1)->getPDGId();

        ev.p_nu.push_back(v_tl_D[2].Px());
        ev.p_nu.push_back(v_tl_D[2].Py());
        ev.p_nu.push_back(v_tl_D[2].Pz());
        ev.p_nu.push_back(v_tl_D[2].E());
        ev.p_nu.push_back(v_tl_D[2].M());

        ev.nu_PE = v_tl_D[2].E();
        ev.nu_PX = v_tl_D[2].Px();
        ev.nu_PY = v_tl_D[2].Py();
        ev.nu_PZ = v_tl_D[2].Pz();

        ev.nu_ID = parent->getDaug(2)->getPDGId();

        // Dst->D0pi

        ev.p_Dst_D0.push_back(vv_tl_DD[0][0].Px());
        ev.p_Dst_D0.push_back(vv_tl_DD[0][0].Py());
        ev.p_Dst_D0.push_back(vv_tl_DD[0][0].Pz());
        ev.p_Dst_D0.push_back(vv_tl_DD[0][0].E());
        ev.p_Dst_D0.push_back(vv_tl_DD[0][0].M());

        ev.Dst_D0_PE = vv_tl_DD[0][0].E();
        ev.Dst_D0_PX = vv_tl_DD[0][0].Px();
        ev.Dst_D0_PY = vv_tl_DD[0][0].Py();
        ev.Dst_D0_PZ = vv_tl_DD[0][0].Pz();

        ev.Dst_D0_ID = parent->getDaug(0)->getDaug(0)->getPDGId();

        ev.p_Dst_pi.push_back(vv_tl_DD[0][1].Px());
        ev.p_Dst_pi.push_back(vv_tl_DD[0][1].Py());
        ev.p_Dst_pi.push_back(vv_tl_DD[0][1].Pz());
        ev.p_Dst_pi.push_back(vv_tl_DD[0][1].E());
        ev.p_Dst_pi.push_back(vv_tl_DD[0][1].M());

        ev.Dst_pi_PE = vv_tl_DD[0][1].E();
        ev.Dst_pi_PX = vv_tl_DD[0][1].Px();
        ev.Dst_pi_PY = vv_tl_DD[0][1].Py();
        ev.Dst_pi_PZ = vv_tl_DD[0][1].Pz();

        ev.Dst_pi_ID = parent->getDaug(0)->getDaug(1)->getPDGId();

        //D0->Kpi

        ev.p_D0_K.push_back(v_tl_DDD[0].Px());
        ev.p_D0_K.push_back(v_tl_DDD[0].Py());
        ev.p_D0_K.push_back(v_tl_DDD[0].Pz());
        ev.p_D0_K.push_back(v_tl_DDD[0].E());
        ev.p_D0_K.push_back(v_tl_DDD[0].M());

        ev.D0_K_PE = v_tl_DDD[0].E();
        ev.D0_K_PX = v_tl_DDD[0].Px();
        ev.D0_K_PY = v_tl_DDD[0].Py();
        ev.D0_K_PZ = v_tl_DDD[0].Pz();
        ev.D0_K_ID = parent->getDaug(0)->getDaug(0)->getDaug(0)->getPDGId();

        ev.p_D0_pi.push_back(v_tl_DDD[1].Px());
        ev.p_D0_pi.push_back(v_tl_DDD[1].Py());
        ev.p_D0_pi.push_back(v_tl_DDD[1].Pz());
        ev.p_D0_pi.push_back(v_tl_DDD[1].E());
        ev.p_D0_pi.push_back(v_tl_DDD[1].M());

        ev.D0_pi_PE = v_tl_DDD[1].E();
        ev.D0_pi_PX = v_tl_DDD[1].Px();
        ev.D0_pi_PY = v_tl_DDD[1].Py();
        ev.D0_pi_PZ = v_tl_DDD[1].Pz();
        ev.D0_pi_ID = parent->getDaug(0)->getDaug(0)->getDaug(1)->getPDGId();

        //Tau->3PiNu

        ev.p_tau_mu.push_back(vv_tl_DD[1][0].Px());
        ev.p_tau_mu.push_back(vv_tl_DD[1][0].Py());
        ev.p_tau_mu.push_back(vv_tl_DD[1][0].Pz());
        ev.p_tau_mu.push_back(vv_tl_DD[1][0].E());
        ev.p_tau_mu.push_back(vv_tl_DD[1][0].M());
        
        ev.tau_mu_PE = vv_tl_DD[1][0].E();
        ev.tau_mu_PX = vv_tl_DD[1][0].Px();
        ev.tau_mu_PY = vv_tl_DD[1][0].Py();
        ev.tau_mu_PZ = vv_tl_DD[1][0].Pz();
        ev.tau_mu_ID = parent->getDaug(1)->getDaug(0)->getPDGId();

        ev.p_tau_nu1.push_back(vv_tl_DD[1][1].Px());
        ev.p_tau_nu1.push_back(vv_tl_DD[1][1].Py());
        ev.p_tau_nu1.push_back(vv_tl_DD[1][1].Pz());
        ev.p_tau_nu1.push_back(vv_tl_DD[1][1].E());
        ev.p_tau_nu1.push_back(vv_tl_DD[1][1].M());

        ev.tau_nu1_PE = vv_tl_DD[1][1].E();
        ev.tau_nu1_PX = vv_tl_DD[1][1].Px();
        ev.tau_nu1_PY = vv_tl_DD[1][1].Py();
        ev.tau_nu1_PZ = vv_tl_DD[1][1].Pz();
        ev.tau_nu1_ID = parent->getDaug(1)->getDaug(1)->getPDGId();

        ev.p_tau_nu2.push_back(vv_tl_DD[1][2].Px());
        ev.p_tau_nu2.push_back(vv_tl_DD[1][2].Py());
        ev.p_tau_nu2.push_back(vv_tl_DD[1][2].Pz());
        ev.p_tau_nu2.push_back(vv_tl_DD[1][2].E());
        ev.p_tau_nu2.push_back(vv_tl_DD[1][2].M());
        
        ev.tau_nu2_PE = vv_tl_DD[1][2].E();
        ev.tau_nu2_PX = vv_tl_DD[1][2].Px();
        ev.tau_nu2_PY = vv_tl_DD[1][2].Py();
        ev.tau_nu2_PZ = vv_tl_DD[1][2].Pz();
        ev.tau_nu2_ID = parent->getDaug(1)->getDaug(2)->getPDGId();

        ev.q2 = (tl_B0 - v_tl_D[0]).M2();
        ev.mmiss2 = (v_tl_D[2]+vv_tl_DD[1][2]+vv_tl_DD[1][1]).M2();

        TVector3 boost0=tl_B0.BoostVector();
        TLorentzVector lep_syst = v_tl_D[1] + v_tl_D[2];
        
        lep_syst.Boost(-boost0);
        vv_tl_DD[1][0].Boost(-boost0);

        TVector3 boost1 = lep_syst.BoostVector();
        vv_tl_DD[1][0].Boost(-boost1);

        ev.El = vv_tl_DD[1][0].E();

        v_events.push_back(ev);        
        parent->deleteTree();
        
        
    }
    
    // Store simulated events
    DumpResults(v_events);
    
    delete eng;
    
    return 0;
}



void DumpResults(std::vector <event> v_evt) {

    TFile f("../root/B02DstTauNu_TauMu.root","RECREATE");    
    TString TreeName = "Decay";

    //TFile f("test.root","UPDATE");
    //TString TreeName = "evtgen";

    TTree *t_evts = new TTree(Form("%sTree",TreeName.Data()),Form("%sTree",TreeName.Data()),-99); t_evts->SetDirectory(0);

    event ev;

    t_evts->Branch("B0_PE", &ev.B0_PE);
    t_evts->Branch("B0_PX", &ev.B0_PX);
    t_evts->Branch("B0_PY", &ev.B0_PY);
    t_evts->Branch("B0_PZ", &ev.B0_PZ);

    t_evts->Branch("Dst_PE", &ev.Dst_PE);
    t_evts->Branch("Dst_PX", &ev.Dst_PX);
    t_evts->Branch("Dst_PY", &ev.Dst_PY);
    t_evts->Branch("Dst_PZ", &ev.Dst_PZ);

    t_evts->Branch("Dst_D0_PE", &ev.Dst_D0_PE);
    t_evts->Branch("Dst_D0_PX", &ev.Dst_D0_PX);
    t_evts->Branch("Dst_D0_PY", &ev.Dst_D0_PY);
    t_evts->Branch("Dst_D0_PZ", &ev.Dst_D0_PZ);

    t_evts->Branch("D0_pi_PE", &ev.D0_pi_PE);
    t_evts->Branch("D0_pi_PX", &ev.D0_pi_PX);
    t_evts->Branch("D0_pi_PY", &ev.D0_pi_PY);
    t_evts->Branch("D0_pi_PZ", &ev.D0_pi_PZ);

    t_evts->Branch("D0_K_PE", &ev.D0_K_PE);
    t_evts->Branch("D0_K_PX", &ev.D0_K_PX);
    t_evts->Branch("D0_K_PY", &ev.D0_K_PY);
    t_evts->Branch("D0_K_PZ", &ev.D0_K_PZ);

    t_evts->Branch("Dst_pi_PE", &ev.Dst_pi_PE);
    t_evts->Branch("Dst_pi_PX", &ev.Dst_pi_PX);
    t_evts->Branch("Dst_pi_PY", &ev.Dst_pi_PY);
    t_evts->Branch("Dst_pi_PZ", &ev.Dst_pi_PZ);

    t_evts->Branch("tau_PE", &ev.tau_PE);
    t_evts->Branch("tau_PX", &ev.tau_PX);
    t_evts->Branch("tau_PY", &ev.tau_PY);
    t_evts->Branch("tau_PZ", &ev.tau_PZ);

    t_evts->Branch("tau_mu_PE", &ev.tau_mu_PE);
    t_evts->Branch("tau_mu_PX", &ev.tau_mu_PX);
    t_evts->Branch("tau_mu_PY", &ev.tau_mu_PY);
    t_evts->Branch("tau_mu_PZ", &ev.tau_mu_PZ);

    t_evts->Branch("tau_nu1_PE", &ev.tau_nu1_PE);
    t_evts->Branch("tau_nu1_PX", &ev.tau_nu1_PX);
    t_evts->Branch("tau_nu1_PY", &ev.tau_nu1_PY);
    t_evts->Branch("tau_nu1_PZ", &ev.tau_nu1_PZ);

    t_evts->Branch("tau_nu2_PE", &ev.tau_nu2_PE);
    t_evts->Branch("tau_nu2_PX", &ev.tau_nu2_PX);
    t_evts->Branch("tau_nu2_PY", &ev.tau_nu2_PY);
    t_evts->Branch("tau_nu2_PZ", &ev.tau_nu2_PZ);


    t_evts->Branch("nu_PE", &ev.nu_PE);
    t_evts->Branch("nu_PX", &ev.nu_PX);
    t_evts->Branch("nu_PY", &ev.nu_PY);
    t_evts->Branch("nu_PZ", &ev.nu_PZ);

    t_evts->Branch("B0_ID", &ev.B0_ID);

    t_evts->Branch("Dst_ID", &ev.Dst_ID);
    t_evts->Branch("Dst_D0_ID", &ev.Dst_D0_ID);
    t_evts->Branch("Dst_pi_ID", &ev.Dst_pi_ID);
    
    t_evts->Branch("D0_K_ID", &ev.D0_K_ID);
    t_evts->Branch("D0_pi_ID", &ev.D0_pi_ID);
    
    t_evts->Branch("tau_ID", &ev.tau_ID);

    t_evts->Branch("tau_mu_ID", &ev.tau_mu_ID);
    t_evts->Branch("tau_nu1_ID", &ev.tau_nu1_ID);
    t_evts->Branch("tau_nu2_ID", &ev.tau_nu2_ID);

    t_evts->Branch("nu_ID", &ev.nu_ID);

    t_evts->Branch("mmiss2", &ev.mmiss2);
    t_evts->Branch("q2", &ev.q2);
    t_evts->Branch("El", &ev.El);


   for(auto ev0 : v_evt) {
    
     // Resize covariance matrix
     ev = ev0;
     // Now let's fill the tree   
     t_evts->Fill();
   }

   t_evts->Write();   
   f.Write(); f.Close();

   std::cout << " Done writing " << "test.root" << std::endl;
    
   t_evts->Delete();

} 

TLorentzVector SmearTLV(TLorentzVector in, double res ) {
    
    TLorentzVector out;
    
    out.SetPx(  smear.Gaus( in.Px(), in.Px() * res )  );
    out.SetPy(  smear.Gaus( in.Py(), in.Py() * res )  );
    out.SetPz(  smear.Gaus( in.Pz(), in.Pz() * res )  );
    out.SetE( sqrt( pow(out.Px(),2) + pow(out.Py(),2) + pow(out.Pz(),2) + pow(in.M(),2) ) );                    
    
    return out;
    
}

TLorentzVector AddBeamBackground() {
 
    
    TLorentzVector p_beambkg;
    
    double phi = 2*TMath::Pi()*smear.Rndm();
    double theta = acos( smear.Rndm()*2 -1 );
    
    double gamma = 50;
    
    double Eph = -1./gamma * log(1 - smear.Rndm() );
    
    // With gamma of 50 ca. 70% of all photons will have an energy 
    //  below 50 MeV
    
    if (Eph > 0.025) {
        
        double Ex = Eph * sin(theta) * cos(phi), 
               Ey = Eph * sin(theta) * sin(phi), 
               Ez = Eph * cos(theta);
        
        p_beambkg.SetPxPyPzE( Ex, Ey, Ez, Eph);
        
    }
        
    return p_beambkg;

}


TLorentzVector LoopOverDaughters(EvtParticle *part) {
    
    // Let's first check if we have daughters
    int nDaug = part->getNDaug();
    // No? then return the parent four vector and we are done
    
    if (nDaug == 0)
        if (KeepDaughter(part))
            return SmearDaughter(part);
        //return TLorentzVector( part->getP4Lab().get(1), part->getP4Lab().get(2), part->getP4Lab().get(3), part->getP4Lab().get(0) );

    // If we have daughters, loop over them
    TLorentzVector tl_X;    
    
    
    for ( int iDaug = 0; iDaug < nDaug; iDaug++ ) {
            
         EvtParticle* daug = part->getDaug( iDaug );
                
         // Do the dauthers have daughters? Then let's call this iteratively
         if (daug->getNDaug() != 0) {
             if (KeepDaughter(daug))
                tl_X += LoopOverDaughters(daug);
         } else {  // No? Then we can sum their four-momentum up
             //TLorentzVector tl_Xd ( daug->getP4Lab().get(1), daug->getP4Lab().get(2), daug->getP4Lab().get(3), daug->getP4Lab().get(0) );
             //tl_X += tl_Xd;
             if (KeepDaughter(daug))
                tl_X += SmearDaughter(daug);
         }
             
             
    }
    
    // Let's add in up to 7 Beam background photons    
    tl_X += AddBeamBackground();
    tl_X += AddBeamBackground();
    tl_X += AddBeamBackground();
    tl_X += AddBeamBackground();
    tl_X += AddBeamBackground();
    tl_X += AddBeamBackground();
    tl_X += AddBeamBackground();
    
    return tl_X;
    
}

bool KeepDaughter(EvtParticle *part) {
    
    int pdg = abs(part->getPDGId());
    TLorentzVector tl_part ( part->getP4Lab().get(1), part->getP4Lab().get(2), part->getP4Lab().get(3), part->getP4Lab().get(0) );
    
    // Killing probability
    double killprob = kill.Rndm();
    
    //killprob = 0;
    
    if (pdg == 111) { // pi0
        
        if ( killprob < 0.50) return true;
            else return false;
        
    } else if (pdg == 22) { // Photon
        
        if ( killprob < sqrt(0.5)) return true;
            else return false;        
        
    } else if (pdg == 130 || pdg == 311) { // Klong
                
        if ( killprob < 0.50) return true;
            else return false;     
        
    } else if (pdg == 311) { // K0
                
        if ( killprob < 0.25) return true; // 50% of K0 are Klong, with an efficiency of 50%
            else return false;     
                
    } else if (pdg == 310) { // Kshort
        
        if ( killprob < 0.80) return true;
            else return false;          
        
    } else { // Kaon or Pion
        
        // Rough estimage of trk efficiency from 1507.03453
        if (tl_part.Pt() < 0.1 && killprob < 0.8)
            return true;
        else if (tl_part.Pt() > 0.1 && tl_part.Pt() < 0.25 && killprob < 0.9)
            return true;
        else if ( killprob < 0.98)
            return true;
        else 
            return false;
        
    }
    
}

 

TLorentzVector SmearDaughter(EvtParticle *part) {
    
    int pdg = abs(part->getPDGId());
    TLorentzVector tl_part ( part->getP4Lab().get(1), part->getP4Lab().get(2), part->getP4Lab().get(3), part->getP4Lab().get(0) );

    if (pdg == 111) { // pi0
        
        return SmearTLV(tl_part, 0.03);        
        
    } else if (pdg == 22) { // Photon
        
        return SmearTLV(tl_part, 0.015);
        
    } else if (pdg == 130) { // Klong
        
        return SmearTLV(tl_part, 0.30);
        
    } else if (pdg == 311) { // K0
        
        if (kill.Rndm() > 0.5) // Decide if this is a K-short or a K-long
            return SmearTLV(tl_part, 0.01);                    
        else
            return SmearTLV(tl_part, 0.30);        
        
    } else if (pdg == 310) { // Kshort
        
        return SmearTLV(tl_part, 0.01);
        
    } else { // Kaon or Pion
        
        return SmearTLV(tl_part, 0.005);
        
    }
    
}

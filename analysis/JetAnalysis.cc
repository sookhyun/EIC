/**
 * @file JetAnalysis.cc
 * @author Sookhyun Lee
 * @date May 2020
 */

#include "JetAnalysis.h"

#include <fastjet/ClusterSequence.hh>
#include <fastjet/Selector.hh>
#include <fastjet/contrib/SoftDrop.hh>
#include <fastjet/PseudoJet.hh>

#include "EventShape.h"
#include "Opts.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <map>
#include <stdio.h>

using namespace std;

JetAnalysis::JetAnalysis():
         infile_mctrue(0)
	,infile_mcreco(0)
	,outfile(0)
	,mctree(0)
	,outtree(0)
	,tevent(0)
	,sevent(0)
	,m_qjet_mode(true)
	,nfiles(0)
	,m_smear_mode(0)
	,verbosity(0)
{
        printf("Hello world!\n");
}

int JetAnalysis::open(std::string fname = "test.root"){

        std::size_t truefound = fname.find("true");
        std::size_t recofound = fname.find("reco");
        if (truefound!=std::string::npos) {
                TFile* infile_mctrue = new TFile(fname.c_str());
                if (!infile_mctrue) cout<<"MC True File not found!"<<endl;
		mctree=(TTree*) infile_mctrue->Get("EICTree");
		++nfiles;
        }
	if (recofound!=std::string::npos){
		infile_mcreco = new TFile(fname.c_str());
		if (!infile_mcreco) cout<<"MC Reco File not found!"<<endl;
		if(mctree){
		mctree->AddFriend("Smeared",fname.c_str());
		}else{
		infile_mcreco->GetObject("Smeared",mctree);
		}
		++nfiles;
	}

        return 1;
}

int JetAnalysis::init(){
	// Event level variables
	eventnum=0;
	outfile = new TFile("tau1.root","recreate");
	outtree = new TTree("T","tau1 tree");
	outtree->Branch("event_number",&eventnum,"event_number/i");
	outtree->Branch("true_tau1_a",&true_tau1_a,"true_tau1_a/F");
	outtree->Branch("true_tau1_b",&true_tau1_b,"true_tau1_b/F");
	outtree->Branch("true_tau1_c",&true_tau1_c,"true_tau1_c/F");
        outtree->Branch("true_rapgap_a",&true_rapgap_a,"true_rapgap_a/F");
        outtree->Branch("true_rapgap_b",&true_rapgap_b,"true_rapgap_b/F");
        outtree->Branch("true_rapgap_c",&true_rapgap_c,"true_rapgap_c/F");
        outtree->Branch("reco_tau1_a",&reco_tau1_a,"reco_tau1_a/F");
        outtree->Branch("reco_tau1_b",&reco_tau1_b,"reco_tau1_b/F");
        outtree->Branch("reco_tau1_c",&reco_tau1_c,"reco_tau1_c/F");
        outtree->Branch("reco_rapgap_a",&reco_rapgap_a,"reco_rapgap_a/F");
        outtree->Branch("reco_rapgap_b",&reco_rapgap_b,"reco_rapgap_b/F");
        outtree->Branch("reco_rapgap_c",&reco_rapgap_c,"reco_rapgap_c/F");

        outtree->Branch("true_x",&true_x,"true_x/F");
        outtree->Branch("true_Q2",&true_Q2,"true_Q2/F");
	outtree->Branch("true_e_y",&true_e_y,"true_e_y/F");
	outtree->Branch("true_e_pt",&true_e_pt,"true_e_pt/F");
	outtree->Branch("true_e_eta",&true_e_eta,"true_e_eta/F");
        outtree->Branch("true_e_y_cm",&true_e_y_cm,"true_e_y_cm/F");
        outtree->Branch("true_e_pt_cm",&true_e_pt_cm,"true_e_pt_cm/F");
        outtree->Branch("true_e_eta_cm",&true_e_eta_cm,"true_e_eta_cm/F");
        outtree->Branch("true_e_y_br",&true_e_y_br,"true_e_y_br/F");
        outtree->Branch("true_e_pt_br",&true_e_pt_br,"true_e_pt_br/F");
        outtree->Branch("true_e_eta_br",&true_e_eta_br,"true_e_eta_br/F");
        outtree->Branch("true_jet_y",&true_jet_y,"true_jet_y/F");
        outtree->Branch("true_jet_y_br",&true_jet_y_br,"true_jet_y_br/F");
        outtree->Branch("true_jet_y_cm",&true_jet_y_cm,"true_jet_y_cm/F");
        outtree->Branch("true_jet_pt",&true_jet_pt,"true_jet_pt/F");
        outtree->Branch("true_jet_pt_br",&true_jet_pt_br,"true_jet_pt_br/F");
        outtree->Branch("true_jet_pt_cm",&true_jet_pt_cm,"true_jet_pt_cm/F");
        outtree->Branch("true_jet_eta",&true_jet_eta,"true_jet_eta/F");
        outtree->Branch("true_jet_eta_br",&true_jet_eta_br,"true_jet_eta_br/F");
        outtree->Branch("true_jet_eta_cm",&true_jet_eta_cm,"true_jet_eta_cm/F");
        outtree->Branch("true_jet_p",&true_jet_p,"true_jet_p/F");
        outtree->Branch("true_jet_p_br",&true_jet_p_br,"true_jet_p_br/F");
        outtree->Branch("true_jet_p_cm",&true_jet_p_cm,"true_jet_p_cm/F");
        outtree->Branch("true_jet_e",&true_jet_e,"true_jet_e/F");
        outtree->Branch("true_jet_e_br",&true_jet_e_br,"true_jet_e_br/F");
        outtree->Branch("true_jet_e_cm",&true_jet_e_cm,"true_jet_e_cm/F");
        outtree->Branch("true_jet_m",&true_jet_m,"true_jet_m/F");
        outtree->Branch("true_jet_m_br",&true_jet_m_br,"true_jet_m_br/F");
        outtree->Branch("true_jet_m_cm",&true_jet_m_cm,"true_jet_m_cm/F");

        outtree->Branch("reco_x",&reco_x,"reco_x/F");
        outtree->Branch("reco_Q2",&reco_Q2,"reco_Q2/F");
        outtree->Branch("reco_e_y",&reco_e_y,"reco_e_y/F");
        outtree->Branch("reco_e_pt",&reco_e_pt,"reco_e_pt/F");
        outtree->Branch("reco_e_eta",&reco_e_eta,"reco_e_eta/F");
        outtree->Branch("reco_e_y_cm",&reco_e_y_cm,"reco_e_y_cm/F");
        outtree->Branch("reco_e_pt_cm",&reco_e_pt_cm,"reco_e_pt_cm/F");
        outtree->Branch("reco_e_eta_cm",&reco_e_eta_cm,"reco_e_eta_cm/F");
        outtree->Branch("reco_e_y_br",&reco_e_y_br,"reco_e_y_br/F");
        outtree->Branch("reco_e_pt_br",&reco_e_pt_br,"reco_e_pt_br/F");
        outtree->Branch("reco_e_eta_br",&reco_e_eta_br,"reco_e_eta_br/F");
        outtree->Branch("reco_jet_y",&reco_jet_y,"reco_jet_y/F");
        outtree->Branch("reco_jet_y_br",&reco_jet_y_br,"reco_jet_y_br/F");
        outtree->Branch("reco_jet_y_cm",&reco_jet_y_cm,"reco_jet_y_cm/F");
        outtree->Branch("reco_jet_pt",&reco_jet_pt,"reco_jet_pt/F");
        outtree->Branch("reco_jet_pt_br",&reco_jet_pt_br,"reco_jet_pt_br/F");
        outtree->Branch("reco_jet_pt_cm",&reco_jet_pt_cm,"reco_jet_pt_cm/F");
        outtree->Branch("reco_jet_eta",&reco_jet_eta,"reco_jet_eta/F");
        outtree->Branch("reco_jet_eta_br",&reco_jet_eta_br,"reco_jet_eta_br/F");
        outtree->Branch("reco_jet_eta_cm",&reco_jet_eta_cm,"reco_jet_eta_cm/F");
        outtree->Branch("reco_jet_p",&reco_jet_p,"reco_jet_p/F");
        outtree->Branch("reco_jet_p_br",&reco_jet_p_br,"reco_jet_p_br/F");
        outtree->Branch("reco_jet_p_cm",&reco_jet_p_cm,"reco_jet_p_cm/F");
        outtree->Branch("reco_jet_e",&reco_jet_e,"reco_jet_e/F");
        outtree->Branch("reco_jet_e_br",&reco_jet_e_br,"reco_jet_e_br/F");
        outtree->Branch("reco_jet_e_cm",&reco_jet_e_cm,"reco_jet_e_cm/F");
        outtree->Branch("reco_jet_m",&reco_jet_m,"reco_jet_m/F");
        outtree->Branch("reco_jet_m_br",&reco_jet_m_br,"reco_jet_m_br/F");
        outtree->Branch("reco_jet_m_cm",&reco_jet_m_cm,"reco_jet_m_cm/F");

        outtree->Branch("true_sf_a_lb",&true_sf_a_lb,"true_sf_a_lb/F");
        outtree->Branch("true_sf_a_br",&true_sf_a_br,"true_sf_a_br/F");
        outtree->Branch("true_sf_a_cm",&true_sf_a_cm,"true_sf_a_cm/F");
        outtree->Branch("true_sf_b_lb",&true_sf_b_lb,"true_sf_b_lb/F");
        outtree->Branch("true_sf_b_br",&true_sf_b_br,"true_sf_b_br/F");
        outtree->Branch("true_sf_b_cm",&true_sf_b_cm,"true_sf_b_cm/F");
        outtree->Branch("true_sf_c_lb",&true_sf_c_lb,"true_sf_c_lb/F");
        outtree->Branch("true_sf_c_br",&true_sf_c_br,"true_sf_c_br/F");
        outtree->Branch("true_sf_c_cm",&true_sf_c_cm,"true_sf_c_cm/F");

        outtree->Branch("reco_sf_a_lb",&reco_sf_a_lb,"reco_sf_a_lb/F");
        outtree->Branch("reco_sf_a_br",&reco_sf_a_br,"reco_sf_a_br/F");
        outtree->Branch("reco_sf_a_cm",&reco_sf_a_cm,"reco_sf_a_cm/F");
        outtree->Branch("reco_sf_b_lb",&reco_sf_b_lb,"reco_sf_b_lb/F");
        outtree->Branch("reco_sf_b_br",&reco_sf_b_br,"reco_sf_b_br/F");
        outtree->Branch("reco_sf_b_cm",&reco_sf_b_cm,"reco_sf_b_cm/F");
        outtree->Branch("reco_sf_c_lb",&reco_sf_c_lb,"reco_sf_c_lb/F");
        outtree->Branch("reco_sf_c_br",&reco_sf_c_br,"reco_sf_c_br/F");
        outtree->Branch("reco_sf_c_cm",&reco_sf_c_cm,"reco_sf_c_cm/F");

        outtree->Branch("true_qf_a",&true_qf_a,"true_qf_a/F");
        outtree->Branch("true_qf_b",&true_qf_b,"true_qf_b/F");
        outtree->Branch("true_qf_c",&true_qf_c,"true_qf_c/F");

	outtree->Branch("true_qj_a",true_qj_a,"true_qj_a[2]/F");
	outtree->Branch("true_qj_b",true_qj_b,"true_qj_b[2]/F");
	outtree->Branch("true_qj_c",true_qj_c,"true_qj_c[2]/F");
        outtree->Branch("true_qj_a_br",true_qj_a_br,"true_qj_a_br[2]/F");
        outtree->Branch("true_qj_b_br",true_qj_b_br,"true_qj_b_br[2]/F");
        outtree->Branch("true_qj_c_br",true_qj_c_br,"true_qj_c_br[2]/F");
        outtree->Branch("true_qj_a_cm",true_qj_a_cm,"true_qj_a_cm[2]/F");
        outtree->Branch("true_qj_b_cm",true_qj_b_cm,"true_qj_b_cm[2]/F");
        outtree->Branch("true_qj_c_cm",true_qj_c_cm,"true_qj_c_cm[2]/F");

	outtree->Branch("true_qb_a",true_qb_a,"true_qb_a[2]/F");
	outtree->Branch("true_qb_b",true_qb_b,"true_qb_b[2]/F");
	outtree->Branch("true_qb_c",true_qb_c,"true_qb_c[2]/F");
        outtree->Branch("true_qb_a_br",true_qb_a_br,"true_qb_a_br[2]/F");
        outtree->Branch("true_qb_b_br",true_qb_b_br,"true_qb_b_br[2]/F");
        outtree->Branch("true_qb_c_br",true_qb_c_br,"true_qb_c_br[2]/F");
        outtree->Branch("true_qb_a_cm",true_qb_a_cm,"true_qb_a_cm[2]/F");
        outtree->Branch("true_qb_b_cm",true_qb_b_cm,"true_qb_b_cm[2]/F");
        outtree->Branch("true_qb_c_cm",true_qb_c_cm,"true_qb_c_cm[2]/F");

	outtree->Branch("true_pjq_a_n",&true_pjq_a_n,"true_pjq_a_n/I");
	outtree->Branch("true_pjq_a_id",true_pjq_a_id,"true_pjq_a_id[true_pjq_a_n]/I");
	outtree->Branch("true_pjq_a_dot_lab",true_pjq_a_dot_lab,"true_pjq_a_dot_lab[true_pjq_a_n]/F");
        outtree->Branch("true_pjq_a_eta_lab",true_pjq_a_eta_lab,"true_pjq_a_eta_lab[true_pjq_a_n]/F");
        outtree->Branch("true_pjq_a_pt_lab",true_pjq_a_pt_lab,"true_pjq_a_pt_lab[true_pjq_a_n]/F");
        outtree->Branch("true_pjq_a_p_lab",true_pjq_a_p_lab,"true_pjq_a_p_lab[true_pjq_a_n]/F");
        outtree->Branch("true_pjq_a_dot_br",true_pjq_a_dot_br,"true_pjq_a_dot_br[true_pjq_a_n]/F");
        outtree->Branch("true_pjq_a_eta_br",true_pjq_a_eta_br,"true_pjq_a_eta_br[true_pjq_a_n]/F");
        outtree->Branch("true_pjq_a_pt_br",true_pjq_a_pt_br,"true_pjq_a_pt_br[true_pjq_a_n]/F");
        outtree->Branch("true_pjq_a_p_br",true_pjq_a_p_br,"true_pjq_a_p_br[true_pjq_a_n]/F");
        outtree->Branch("true_pjq_a_dot_cm",true_pjq_a_dot_cm,"true_pjq_a_dot_cm[true_pjq_a_n]/F");
        outtree->Branch("true_pjq_a_eta_cm",true_pjq_a_eta_cm,"true_pjq_a_eta_cm[true_pjq_a_n]/F");
        outtree->Branch("true_pjq_a_pt_cm",true_pjq_a_pt_cm,"true_pjq_a_pt_cm[true_pjq_a_n]/F");
        outtree->Branch("true_pjq_a_p_cm",true_pjq_a_p_cm,"true_pjq_a_p_cm[true_pjq_a_n]/F");

	outtree->Branch("true_pjq_b_n",&true_pjq_b_n,"true_pjq_b_n/I");
	outtree->Branch("true_pjq_b_dot_br",&true_pjq_b_dot_br,"true_pjq_b_dot_br[true_pjq_b_n]/F");
	outtree->Branch("true_pjq_b_eta_br",&true_pjq_b_eta_br,"true_pjq_b_eta_br[true_pjq_b_n]/F");
	outtree->Branch("true_pjq_b_pt_br",&true_pjq_b_pt_br,"true_pjq_b_pt_br[true_pjq_b_n]/F");
	outtree->Branch("true_pjq_b_p_br",&true_pjq_b_p_br,"true_pjq_b_p_br[true_pjq_b_n]/F");
        outtree->Branch("true_pjq_b_dot_lab",&true_pjq_b_dot_lab,"true_pjq_b_dot_lab[true_pjq_b_n]/F");
        outtree->Branch("true_pjq_b_eta_lab",&true_pjq_b_eta_lab,"true_pjq_b_eta_lab[true_pjq_b_n]/F");
        outtree->Branch("true_pjq_b_pt_lab",&true_pjq_b_pt_lab,"true_pjq_b_pt_lab[true_pjq_b_n]/F");
        outtree->Branch("true_pjq_b_p_lab",&true_pjq_b_p_lab,"true_pjq_b_p_lab[true_pjq_b_n]/F");
        outtree->Branch("true_pjq_b_dot_cm",&true_pjq_b_dot_cm,"true_pjq_b_dot_cm[true_pjq_b_n]/F");
        outtree->Branch("true_pjq_b_eta_cm",&true_pjq_b_eta_cm,"true_pjq_b_eta_cm[true_pjq_b_n]/F");
        outtree->Branch("true_pjq_b_pt_cm",&true_pjq_b_pt_cm,"true_pjq_b_pt_cm[true_pjq_b_n]/F");
        outtree->Branch("true_pjq_b_p_cm",&true_pjq_b_p_cm,"true_pjq_b_p_cm[true_pjq_b_n]/F");

	outtree->Branch("true_pjq_c_n",&true_pjq_c_n,"true_pjq_c_n/I");
	outtree->Branch("true_pjq_c_dot_cm",&true_pjq_c_dot_cm,"true_pjq_c_dot_cm[true_pjq_c_n]/F");
	outtree->Branch("true_pjq_c_eta_cm",&true_pjq_c_eta_cm,"true_pjq_c_eta_cm[true_pjq_c_n]/F");
	outtree->Branch("true_pjq_c_pt_cm",&true_pjq_c_pt_cm,"true_pjq_c_pt_cm[true_pjq_c_n]/F");
	outtree->Branch("true_pjq_c_p_cm",&true_pjq_c_p_cm,"true_pjq_c_p_cm[true_pjq_c_n]/F");
        outtree->Branch("true_pjq_c_dot_br",&true_pjq_c_dot_br,"true_pjq_c_dot_br[true_pjq_c_n]/F");
        outtree->Branch("true_pjq_c_eta_br",&true_pjq_c_eta_br,"true_pjq_c_eta_br[true_pjq_c_n]/F");
        outtree->Branch("true_pjq_c_pt_br",&true_pjq_c_pt_br,"true_pjq_c_pt_br[true_pjq_c_n]/F");
        outtree->Branch("true_pjq_c_p_br",&true_pjq_c_p_br,"true_pjq_c_p_br[true_pjq_c_n]/F");
        outtree->Branch("true_pjq_c_dot_lab",&true_pjq_c_dot_lab,"true_pjq_c_dot_lab[true_pjq_c_n]/F");
        outtree->Branch("true_pjq_c_eta_lab",&true_pjq_c_eta_lab,"true_pjq_c_eta_lab[true_pjq_c_n]/F");
        outtree->Branch("true_pjq_c_pt_lab",&true_pjq_c_pt_lab,"true_pjq_c_pt_lab[true_pjq_c_n]/F");
        outtree->Branch("true_pjq_c_p_lab",&true_pjq_c_p_lab,"true_pjq_c_p_lab[true_pjq_c_n]/F");


	outtree->Branch("true_pjb_a_n",&true_pjb_a_n,"true_pjb_a_n/I");
	outtree->Branch("true_pjb_a_id",true_pjb_a_id,"true_pjb_a_id[true_pjb_a_n]/I");
	outtree->Branch("true_pjb_a_dot_lab",true_pjb_a_dot_lab,"true_pjb_a_dot_lab[true_pjb_a_n]/F");
        outtree->Branch("true_pjb_a_eta_lab",true_pjb_a_eta_lab,"true_pjb_a_eta_lab[true_pjb_a_n]/F");
        outtree->Branch("true_pjb_a_pt_lab",true_pjb_a_pt_lab,"true_pjb_a_pt_lab[true_pjb_a_n]/F");
        outtree->Branch("true_pjb_a_p_lab",true_pjb_a_p_lab,"true_pjb_a_p_lab[true_pjb_a_n]/F");
        outtree->Branch("true_pjb_a_dot_br",true_pjb_a_dot_br,"true_pjb_a_dot_br[true_pjb_a_n]/F");
        outtree->Branch("true_pjb_a_eta_br",true_pjb_a_eta_br,"true_pjb_a_eta_br[true_pjb_a_n]/F");
        outtree->Branch("true_pjb_a_pt_br",true_pjb_a_pt_br,"true_pjb_a_pt_br[true_pjb_a_n]/F");
        outtree->Branch("true_pjb_a_p_br",true_pjb_a_p_br,"true_pjb_a_p_br[true_pjb_a_n]/F");
        outtree->Branch("true_pjb_a_dot_cm",true_pjb_a_dot_cm,"true_pjb_a_dot_cm[true_pjb_a_n]/F");
        outtree->Branch("true_pjb_a_eta_cm",true_pjb_a_eta_cm,"true_pjb_a_eta_cm[true_pjb_a_n]/F");
        outtree->Branch("true_pjb_a_pt_cm",true_pjb_a_pt_cm,"true_pjb_a_pt_cm[true_pjb_a_n]/F");
        outtree->Branch("true_pjb_a_p_cm",true_pjb_a_p_cm,"true_pjb_a_p_cm[true_pjb_a_n]/F");

        outtree->Branch("true_pjb_b_n",&true_pjb_b_n,"true_pjb_b_n/I");
        outtree->Branch("true_pjb_b_dot_lab",&true_pjb_b_dot_lab,"true_pjb_b_dot_lab[true_pjb_b_n]/F");
        outtree->Branch("true_pjb_b_eta_lab",&true_pjb_b_eta_lab,"true_pjb_b_eta_lab[true_pjb_b_n]/F");
        outtree->Branch("true_pjb_b_pt_lab",&true_pjb_b_pt_lab,"true_pjb_b_pt_lab[true_pjb_b_n]/F");
        outtree->Branch("true_pjb_b_p_lab",&true_pjb_b_p_lab,"true_pjb_b_p_lab[true_pjb_b_n]/F");
        outtree->Branch("true_pjb_b_dot_br",&true_pjb_b_dot_br,"true_pjb_b_dot_br[true_pjb_b_n]/F");
        outtree->Branch("true_pjb_b_eta_br",&true_pjb_b_eta_br,"true_pjb_b_eta_br[true_pjb_b_n]/F");
        outtree->Branch("true_pjb_b_pt_br",&true_pjb_b_pt_br,"true_pjb_b_pt_br[true_pjb_b_n]/F");
        outtree->Branch("true_pjb_b_p_br",&true_pjb_b_p_br,"true_pjb_b_p_br[true_pjb_b_n]/F");
        outtree->Branch("true_pjb_b_dot_cm",&true_pjb_b_dot_cm,"true_pjb_b_dot_cm[true_pjb_b_n]/F");
        outtree->Branch("true_pjb_b_eta_cm",&true_pjb_b_eta_cm,"true_pjb_b_eta_cm[true_pjb_b_n]/F");
        outtree->Branch("true_pjb_b_pt_cm",&true_pjb_b_pt_cm,"true_pjb_b_pt_cm[true_pjb_b_n]/F");
        outtree->Branch("true_pjb_b_p_cm",&true_pjb_b_p_cm,"true_pjb_b_p_cm[true_pjb_b_n]/F");

        outtree->Branch("true_pjb_c_n",&true_pjb_c_n,"true_pjb_c_n/I");
        outtree->Branch("true_pjb_c_dot_cm",&true_pjb_c_dot_cm,"true_pjb_c_dot_cm[true_pjb_c_n]/F");
        outtree->Branch("true_pjb_c_eta_cm",&true_pjb_c_eta_cm,"true_pjb_c_eta_cm[true_pjb_c_n]/F");
        outtree->Branch("true_pjb_c_pt_cm",&true_pjb_c_pt_cm,"true_pjb_c_pt_cm[true_pjb_c_n]/F");
        outtree->Branch("true_pjb_c_p_cm",&true_pjb_c_p_cm,"true_pjb_c_p_cm[true_pjb_c_n]/F");
        outtree->Branch("true_pjb_c_dot_br",&true_pjb_c_dot_br,"true_pjb_c_dot_br[true_pjb_c_n]/F");
        outtree->Branch("true_pjb_c_eta_br",&true_pjb_c_eta_br,"true_pjb_c_eta_br[true_pjb_c_n]/F");
        outtree->Branch("true_pjb_c_pt_br",&true_pjb_c_pt_br,"true_pjb_c_pt_br[true_pjb_c_n]/F");
        outtree->Branch("true_pjb_c_p_br",&true_pjb_c_p_br,"true_pjb_c_p_br[true_pjb_c_n]/F");
        outtree->Branch("true_pjb_c_dot_lab",&true_pjb_c_dot_lab,"true_pjb_c_dot_lab[true_pjb_c_n]/F");
        outtree->Branch("true_pjb_c_eta_lab",&true_pjb_c_eta_lab,"true_pjb_c_eta_lab[true_pjb_c_n]/F");
        outtree->Branch("true_pjb_c_pt_lab",&true_pjb_c_pt_lab,"true_pjb_c_pt_lab[true_pjb_c_n]/F");
        outtree->Branch("true_pjb_c_p_lab",&true_pjb_c_p_lab,"true_pjb_c_p_lab[true_pjb_c_n]/F");


	// arrays for a only, b and c to be implemented later
        outtree->Branch("reco_qj_a",reco_qj_a,"reco_qj_a[2]/F");
        outtree->Branch("reco_qj_b",reco_qj_b,"reco_qj_b[2]/F");
        outtree->Branch("reco_qj_c",reco_qj_c,"reco_qj_c[2]/F");
        outtree->Branch("reco_qj_a_br",reco_qj_a_br,"reco_qj_a_br[2]/F");
        outtree->Branch("reco_qj_b_br",reco_qj_b_br,"reco_qj_b_br[2]/F");
        outtree->Branch("reco_qj_c_br",reco_qj_c_br,"reco_qj_c_br[2]/F");
        outtree->Branch("reco_qj_a_cm",reco_qj_a_cm,"reco_qj_a_cm[2]/F");
        outtree->Branch("reco_qj_b_cm",reco_qj_b_cm,"reco_qj_b_cm[2]/F");
        outtree->Branch("reco_qj_c_cm",reco_qj_c_cm,"reco_qj_c_cm[2]/F");

        outtree->Branch("reco_qb_a",reco_qb_a,"reco_qb_a[2]/F");
        outtree->Branch("reco_qb_b",reco_qb_b,"reco_qb_b[2]/F");
        outtree->Branch("reco_qb_c",reco_qb_c,"reco_qb_c[2]/F");
        outtree->Branch("reco_qb_a_br",reco_qb_a_br,"reco_qb_a_br[2]/F");
        outtree->Branch("reco_qb_b_br",reco_qb_b_br,"reco_qb_b_br[2]/F");
        outtree->Branch("reco_qb_c_br",reco_qb_c_br,"reco_qb_c_br[2]/F");
        outtree->Branch("reco_qb_a_cm",reco_qb_a_cm,"reco_qb_a_cm[2]/F");
        outtree->Branch("reco_qb_b_cm",reco_qb_b_cm,"reco_qb_b_cm[2]/F");
        outtree->Branch("reco_qb_c_cm",reco_qb_c_cm,"reco_qb_c_cm[2]/F");

        outtree->Branch("reco_pjq_a_n",&reco_pjq_a_n,"reco_pjq_a_n/I");
	outtree->Branch("reco_pjq_a_id",reco_pjq_a_id,"reco_pjq_a_id[reco_pjq_a_n]/I");
        outtree->Branch("reco_pjq_a_dot_lab",reco_pjq_a_dot_lab,"reco_pjq_a_dot_lab[reco_pjq_a_n]/F");
        outtree->Branch("reco_pjq_a_eta_lab",reco_pjq_a_eta_lab,"reco_pjq_a_eta_lab[reco_pjq_a_n]/F");
        outtree->Branch("reco_pjq_a_pt_lab",reco_pjq_a_pt_lab,"reco_pjq_a_pt_lab[reco_pjq_a_n]/F");
        outtree->Branch("reco_pjq_a_p_lab",reco_pjq_a_p_lab,"reco_pjq_a_p_lab[reco_pjq_a_n]/F");
        outtree->Branch("reco_pjq_a_dot_br",reco_pjq_a_dot_br,"reco_pjq_a_dot_br[reco_pjq_a_n]/F");
        outtree->Branch("reco_pjq_a_eta_br",reco_pjq_a_eta_br,"reco_pjq_a_eta_br[reco_pjq_a_n]/F");
        outtree->Branch("reco_pjq_a_pt_br",reco_pjq_a_pt_br,"reco_pjq_a_pt_br[reco_pjq_a_n]/F");
        outtree->Branch("reco_pjq_a_p_br",reco_pjq_a_p_br,"reco_pjq_a_p_br[reco_pjq_a_n]/F");
        outtree->Branch("reco_pjq_a_dot_cm",reco_pjq_a_dot_cm,"reco_pjq_a_dot_cm[reco_pjq_a_n]/F");
        outtree->Branch("reco_pjq_a_eta_cm",reco_pjq_a_eta_cm,"reco_pjq_a_eta_cm[reco_pjq_a_n]/F");
        outtree->Branch("reco_pjq_a_pt_cm",reco_pjq_a_pt_cm,"reco_pjq_a_pt_cm[reco_pjq_a_n]/F");
        outtree->Branch("reco_pjq_a_p_cm",reco_pjq_a_p_cm,"reco_pjq_a_p_cm[reco_pjq_a_n]/F");

        outtree->Branch("reco_pjb_a_n",&reco_pjb_a_n,"reco_pjb_a_n/I");
	outtree->Branch("reco_pjb_a_id",&reco_pjb_a_id,"reco_pjb_a_id[reco_pjb_a_n]/I");
        outtree->Branch("reco_pjb_a_dot_lab",reco_pjb_a_dot_lab,"reco_pjb_a_dot_lab[reco_pjb_a_n]/F");
        outtree->Branch("reco_pjb_a_eta_lab",reco_pjb_a_eta_lab,"reco_pjb_a_eta_lab[reco_pjb_a_n]/F");
        outtree->Branch("reco_pjb_a_pt_lab",reco_pjb_a_pt_lab,"reco_pjb_a_pt_lab[reco_pjb_a_n]/F");
        outtree->Branch("reco_pjb_a_p_lab",reco_pjb_a_p_lab,"reco_pjb_a_p_lab[reco_pjb_a_n]/F");
        outtree->Branch("reco_pjb_a_dot_br",reco_pjb_a_dot_br,"reco_pjb_a_dot_br[reco_pjb_a_n]/F");
        outtree->Branch("reco_pjb_a_eta_br",reco_pjb_a_eta_br,"reco_pjb_a_eta_br[reco_pjb_a_n]/F");
        outtree->Branch("reco_pjb_a_pt_br",reco_pjb_a_pt_br,"reco_pjb_a_pt_br[reco_pjb_a_n]/F");
        outtree->Branch("reco_pjb_a_p_br",reco_pjb_a_p_br,"reco_pjb_a_p_br[reco_pjb_a_n]/F");
        outtree->Branch("reco_pjb_a_dot_cm",reco_pjb_a_dot_cm,"reco_pjb_a_dot_cm[reco_pjb_a_n]/F");
        outtree->Branch("reco_pjb_a_eta_cm",reco_pjb_a_eta_cm,"reco_pjb_a_eta_cm[reco_pjb_a_n]/F");
        outtree->Branch("reco_pjb_a_pt_cm",reco_pjb_a_pt_cm,"reco_pjb_a_pt_cm[reco_pjb_a_n]/F");
        outtree->Branch("reco_pjb_a_p_cm",reco_pjb_a_p_cm,"reco_pjb_a_p_cm[reco_pjb_a_n]/F");


	// Particle level variables
	// pt_cm pt_br pt 
	// y_jb
	return 1;
}

int JetAnalysis::save_results(){

	outfile->cd();
	outtree->Write();
	outfile->Close();
	return 1;
}

int JetAnalysis::analyze(){
	set_branches_mctree();

	unsigned int totentries= mctree->GetEntries();
	cout<<"Total entries "<<totentries<<endl;
	int tau_ab=0;
	int tau_bc=0;
	TLorentzVector qj_a; TLorentzVector qj_b; TLorentzVector qj_c;
	TLorentzVector qb_a; TLorentzVector qb_b; TLorentzVector qb_c;
//	htjpt= new TH1D("htjpt","htjpt",200,0.,50.);
//	htj_yphi= new TH2D("htjyphi","",200, -3.1415,3.1415,100,2.,5.);
	Frame fra = Frame::Lab;
	Frame frb = Frame::Breit;
	Frame frc = Frame::CM;
	int nq=0; int nb=0;
	for(unsigned int ii=0; ii<totentries; ++ii)
	{
		mctree->GetEntry(ii);
                EventShape es(*tevent,*sevent);
		es.set_verbosity(verbosity);
		cout<<"JetAnalysis::anlyze m_qjet_mode "<<m_qjet_mode<<endl;
		es.set_xQ2_method(m_Q2_method);
		es.set_qjet_mode(m_qjet_mode);
		es.set_smear_mode(m_smear_mode);
		es.set_tracking_mode(m_tracking_mode);
		es.set_pt_cut(m_pt_cut);
		for(int f=0; f<3; ++f){
			switch(f)	
			{
			case 0:
				// Truth
				true_tau1_a =es.get_tau1(TauType::A, EType::MCTrue);

                                true_qf_a = es.get_tau1_q_fraction(TauType::A);

				true_pjq_a_n=es.get_pjq(TauType::A, fra, Var::Dot, EType::MCTrue).size();
                                true_pjb_a_n=es.get_pjb(TauType::A, fra, Var::Dot, EType::MCTrue).size();

				nq += true_pjq_a_n;
				nb += true_pjb_a_n;

				fra = Frame::Lab;
				vector_to_array(true_pjq_a_id, es.get_pjq_id(TauType::A));
				vector_to_array(true_pjb_a_id, es.get_pjb_id(TauType::A));
				vector_to_array(true_pjq_a_dot_lab,es.get_pjq(TauType::A, fra, Var::Dot, EType::MCTrue));
				vector_to_array(true_pjb_a_dot_lab,es.get_pjb(TauType::A, fra, Var::Dot, EType::MCTrue));
				vector_to_array(true_pjq_a_eta_lab,es.get_pjq(TauType::A, fra, Var::Eta, EType::MCTrue));
				vector_to_array(true_pjb_a_eta_lab,es.get_pjb(TauType::A, fra, Var::Eta, EType::MCTrue));
				vector_to_array(true_pjq_a_pt_lab,es.get_pjq(TauType::A, fra, Var::Pt, EType::MCTrue));
				vector_to_array(true_pjb_a_pt_lab,es.get_pjb(TauType::A, fra, Var::Pt, EType::MCTrue));
				vector_to_array(true_pjq_a_p_lab,es.get_pjq(TauType::A, fra, Var::P, EType::MCTrue));
				vector_to_array(true_pjb_a_p_lab,es.get_pjb(TauType::A, fra, Var::P, EType::MCTrue));
				qj_a = es.get_qj(TauType::A, fra);		
				true_qj_a[0]= qj_a.Eta();
                                qb_a = es.get_qb(TauType::A, fra);
                                true_qb_a[0]= qb_a.Eta();
				true_sf_a_lb = es.get_suppression_factor(TauType::A, fra);

				fra =Frame::Breit;
                                vector_to_array(true_pjq_a_dot_br,es.get_pjq(TauType::A, fra, Var::Dot, EType::MCTrue));
                                vector_to_array(true_pjb_a_dot_br,es.get_pjb(TauType::A, fra, Var::Dot, EType::MCTrue));
                                vector_to_array(true_pjq_a_eta_br,es.get_pjq(TauType::A, fra, Var::Eta, EType::MCTrue));
                                vector_to_array(true_pjb_a_eta_br,es.get_pjb(TauType::A, fra, Var::Eta, EType::MCTrue));
                                vector_to_array(true_pjq_a_pt_br,es.get_pjq(TauType::A, fra, Var::Pt, EType::MCTrue));
                                vector_to_array(true_pjb_a_pt_br,es.get_pjb(TauType::A, fra, Var::Pt, EType::MCTrue));
                                vector_to_array(true_pjq_a_p_br,es.get_pjq(TauType::A, fra, Var::P, EType::MCTrue));
                                vector_to_array(true_pjb_a_p_br,es.get_pjb(TauType::A, fra, Var::P, EType::MCTrue));
                                qj_a = es.get_qj(TauType::A, fra);
                                true_qj_a_br[0]= qj_a.Eta();
                                qb_a = es.get_qb(TauType::A, fra);
                                true_qb_a_br[0]= qb_a.Eta();
                                true_sf_a_br = es.get_suppression_factor(TauType::A, fra);

				fra = Frame::CM;
                                vector_to_array(true_pjq_a_dot_cm,es.get_pjq(TauType::A, fra, Var::Dot, EType::MCTrue));
                                vector_to_array(true_pjb_a_dot_cm,es.get_pjb(TauType::A, fra, Var::Dot, EType::MCTrue));
                                vector_to_array(true_pjq_a_eta_cm,es.get_pjq(TauType::A, fra, Var::Eta, EType::MCTrue));
                                vector_to_array(true_pjb_a_eta_cm,es.get_pjb(TauType::A, fra, Var::Eta, EType::MCTrue));
                                vector_to_array(true_pjq_a_pt_cm,es.get_pjq(TauType::A, fra, Var::Pt, EType::MCTrue));
                                vector_to_array(true_pjb_a_pt_cm,es.get_pjb(TauType::A, fra, Var::Pt, EType::MCTrue));
                                vector_to_array(true_pjq_a_p_cm,es.get_pjq(TauType::A, fra, Var::P, EType::MCTrue));
                                vector_to_array(true_pjb_a_p_cm,es.get_pjb(TauType::A, fra, Var::P, EType::MCTrue));
                                qj_a = es.get_qj(TauType::A, fra);
                                true_qj_a_cm[0]= qj_a.Eta();
                                qb_a = es.get_qb(TauType::A, fra);
                                true_qb_a_cm[0]= qb_a.Eta();
                                true_sf_a_cm = es.get_suppression_factor(TauType::A, fra);

				// Reco
                                reco_tau1_a =es.get_tau1(TauType::A, EType::MCReco);
                                reco_pjq_a_n=es.get_pjq(TauType::A, fra, Var::Dot, EType::MCReco).size();
                                reco_pjb_a_n=es.get_pjb(TauType::A, fra, Var::Dot, EType::MCReco).size();

                                fra = Frame::Lab;
                                vector_to_array(reco_pjq_a_id, es.get_pjq_id(TauType::A));
                                vector_to_array(reco_pjb_a_id, es.get_pjb_id(TauType::A));
                                vector_to_array(reco_pjq_a_dot_lab,es.get_pjq(TauType::A, fra, Var::Dot, EType::MCReco));
                                vector_to_array(reco_pjb_a_dot_lab,es.get_pjb(TauType::A, fra, Var::Dot, EType::MCReco));
                                vector_to_array(reco_pjq_a_eta_lab,es.get_pjq(TauType::A, fra, Var::Eta, EType::MCReco));
                                vector_to_array(reco_pjb_a_eta_lab,es.get_pjb(TauType::A, fra, Var::Eta, EType::MCReco));
                                vector_to_array(reco_pjq_a_pt_lab,es.get_pjq(TauType::A, fra, Var::Pt, EType::MCReco));
                                vector_to_array(reco_pjb_a_pt_lab,es.get_pjb(TauType::A, fra, Var::Pt, EType::MCReco));
                                vector_to_array(reco_pjq_a_p_lab,es.get_pjq(TauType::A, fra, Var::P, EType::MCReco));
                                vector_to_array(reco_pjb_a_p_lab,es.get_pjb(TauType::A, fra, Var::P, EType::MCReco));
                                qj_a = es.get_qj(TauType::A, fra);
                                reco_qj_a[0]= qj_a.Eta();
                                qb_a = es.get_qb(TauType::A, fra);
                                reco_qb_a[0]= qb_a.Eta();
                                reco_sf_a_lb = es.get_suppression_factor(TauType::A, fra);


                                if(verbosity) cout<<"tau1 Aligned "<<true_tau1_a<<endl;
			break;
			case 1: 
				// Truth
				true_tau1_b =es.get_tau1(TauType::B, EType::MCTrue);

                                true_qf_b = es.get_tau1_q_fraction(TauType::B);
  
                                true_pjq_b_n=es.get_pjq(TauType::B, frb, Var::Dot, EType::MCTrue).size();
                                true_pjb_b_n=es.get_pjb(TauType::B, frb, Var::Dot, EType::MCTrue).size();

				frb = Frame::Breit;
                                vector_to_array(true_pjq_b_dot_br,es.get_pjq(TauType::B, frb, Var::Dot, EType::MCTrue));
                                vector_to_array(true_pjb_b_dot_br,es.get_pjb(TauType::B, frb, Var::Dot, EType::MCTrue));
                                vector_to_array(true_pjq_b_eta_br,es.get_pjq(TauType::B, frb, Var::Eta, EType::MCTrue));
                                vector_to_array(true_pjb_b_eta_br,es.get_pjb(TauType::B, frb, Var::Eta, EType::MCTrue));
                                vector_to_array(true_pjq_b_pt_br,es.get_pjq(TauType::B, frb, Var::Pt, EType::MCTrue));
                                vector_to_array(true_pjb_b_pt_br,es.get_pjb(TauType::B, frb, Var::Pt, EType::MCTrue));
                                vector_to_array(true_pjq_b_p_br,es.get_pjq(TauType::B, frb, Var::P, EType::MCTrue));
                                vector_to_array(true_pjb_b_p_br,es.get_pjb(TauType::B, frb, Var::P, EType::MCTrue));
                                qj_b = es.get_qj(TauType::B, frb);
                                true_qj_b_br[0]= qj_b.Eta();
                                qb_b = es.get_qb(TauType::B, frb);
                                true_qb_b_br[0]= qb_b.Eta();
                                true_sf_b_br = es.get_suppression_factor(TauType::B, fra);

                                frb= Frame::Lab;
                                vector_to_array(true_pjq_b_dot_lab,es.get_pjq(TauType::B, frb, Var::Dot, EType::MCTrue));
                                vector_to_array(true_pjb_b_dot_lab,es.get_pjb(TauType::B, frb, Var::Dot, EType::MCTrue));
                                vector_to_array(true_pjq_b_eta_lab,es.get_pjq(TauType::B, frb, Var::Eta, EType::MCTrue));
                                vector_to_array(true_pjb_b_eta_lab,es.get_pjb(TauType::B, frb, Var::Eta, EType::MCTrue));
                                vector_to_array(true_pjq_b_pt_lab,es.get_pjq(TauType::B, frb, Var::Pt, EType::MCTrue));
                                vector_to_array(true_pjb_b_pt_lab,es.get_pjb(TauType::B, frb, Var::Pt, EType::MCTrue));
                                vector_to_array(true_pjq_b_p_lab,es.get_pjq(TauType::B, frb, Var::P, EType::MCTrue));
                                vector_to_array(true_pjb_b_p_lab,es.get_pjb(TauType::B, frb, Var::P, EType::MCTrue));
                                qj_b = es.get_qj(TauType::B, frb);
                                true_qj_b[0]= qj_b.Eta();
                                qb_b = es.get_qb(TauType::B, frb);
                                true_qb_b[0]= qb_b.Eta();
                                true_sf_b_lb = es.get_suppression_factor(TauType::B, fra);

                                frb= Frame::CM;
                                vector_to_array(true_pjq_b_dot_cm,es.get_pjq(TauType::B, frb, Var::Dot, EType::MCTrue));
                                vector_to_array(true_pjb_b_dot_cm,es.get_pjb(TauType::B, frb, Var::Dot, EType::MCTrue));
                                vector_to_array(true_pjq_b_eta_cm,es.get_pjq(TauType::B, frb, Var::Eta, EType::MCTrue));
                                vector_to_array(true_pjb_b_eta_cm,es.get_pjb(TauType::B, frb, Var::Eta, EType::MCTrue));
                                vector_to_array(true_pjq_b_pt_cm,es.get_pjq(TauType::B, frb, Var::Pt, EType::MCTrue));
                                vector_to_array(true_pjb_b_pt_cm,es.get_pjb(TauType::B, frb, Var::Pt, EType::MCTrue));
                                vector_to_array(true_pjq_b_p_cm,es.get_pjq(TauType::B, frb, Var::P, EType::MCTrue));
                                vector_to_array(true_pjb_b_p_cm,es.get_pjb(TauType::B, frb, Var::P, EType::MCTrue));
                                qj_b = es.get_qj(TauType::B, frb);
                                true_qj_b_cm[0]= qj_b.Eta();
                                qb_b = es.get_qb(TauType::B, frb);
                                true_qb_b_cm[0]= qb_b.Eta();
                                true_sf_b_cm = es.get_suppression_factor(TauType::B, fra);

				// Reco
				reco_tau1_b =es.get_tau1(TauType::B, EType::MCReco);

				if(verbosity) cout<<"tau1 Breit "<<true_tau1_b<<" "<<reco_tau1_b<<endl;
			break;
			case 2:
				// Truth
                                true_tau1_c =es.get_tau1(TauType::C, EType::MCTrue);

                                true_qf_c = es.get_tau1_q_fraction(TauType::C);

                                true_pjq_c_n=es.get_pjq(TauType::C, frc, Var::Dot, EType::MCTrue).size();
                                true_pjb_c_n=es.get_pjb(TauType::C, frc, Var::Dot, EType::MCTrue).size();

				frc = Frame::CM;
                                vector_to_array(true_pjq_c_dot_cm,es.get_pjq(TauType::C, frc, Var::Dot, EType::MCTrue));
                                vector_to_array(true_pjb_c_dot_cm,es.get_pjb(TauType::C, frc, Var::Dot, EType::MCTrue));
                                vector_to_array(true_pjq_c_eta_cm,es.get_pjq(TauType::C, frc, Var::Eta, EType::MCTrue));
                                vector_to_array(true_pjb_c_eta_cm,es.get_pjb(TauType::C, frc, Var::Eta, EType::MCTrue));
                                vector_to_array(true_pjq_c_pt_cm,es.get_pjq(TauType::C, frc, Var::Pt, EType::MCTrue));
                                vector_to_array(true_pjb_c_pt_cm,es.get_pjb(TauType::C, frc, Var::Pt, EType::MCTrue));
                                vector_to_array(true_pjq_c_p_cm,es.get_pjq(TauType::C, frc, Var::P, EType::MCTrue));
                                vector_to_array(true_pjb_c_p_cm,es.get_pjb(TauType::C, frc, Var::P, EType::MCTrue));
                                qj_c = es.get_qj(TauType::C, frc);
                                true_qj_c_cm[0]= qj_c.Eta();
                                qb_c = es.get_qb(TauType::C, frc);
                                true_qb_c_cm[0]= qb_c.Eta();
                                true_sf_c_cm = es.get_suppression_factor(TauType::C, fra);

				frc = Frame::Breit;
                                vector_to_array(true_pjq_c_dot_br,es.get_pjq(TauType::C, frc, Var::Dot, EType::MCTrue));
                                vector_to_array(true_pjb_c_dot_br,es.get_pjb(TauType::C, frc, Var::Dot, EType::MCTrue));
                                vector_to_array(true_pjq_c_eta_br,es.get_pjq(TauType::C, frc, Var::Eta, EType::MCTrue));
                                vector_to_array(true_pjb_c_eta_br,es.get_pjb(TauType::C, frc, Var::Eta, EType::MCTrue));
                                vector_to_array(true_pjq_c_pt_br,es.get_pjq(TauType::C, frc, Var::Pt, EType::MCTrue));
                                vector_to_array(true_pjb_c_pt_br,es.get_pjb(TauType::C, frc, Var::Pt, EType::MCTrue));
                                vector_to_array(true_pjq_c_p_br,es.get_pjq(TauType::C, frc, Var::P, EType::MCTrue));
                                vector_to_array(true_pjb_c_p_br,es.get_pjb(TauType::C, frc, Var::P, EType::MCTrue));
                                qj_c = es.get_qj(TauType::C, frc);
                                true_qj_c_br[0]= qj_c.Eta();
                                qb_c = es.get_qb(TauType::C, frc);
                                true_qb_c_br[0]= qb_c.Eta();
                                true_sf_c_br = es.get_suppression_factor(TauType::C, fra);

				frc = Frame::Lab;
                                vector_to_array(true_pjq_c_dot_lab,es.get_pjq(TauType::C, frc, Var::Dot, EType::MCTrue));
                                vector_to_array(true_pjb_c_dot_lab,es.get_pjb(TauType::C, frc, Var::Dot, EType::MCTrue));
                                vector_to_array(true_pjq_c_eta_lab,es.get_pjq(TauType::C, frc, Var::Eta, EType::MCTrue));
                                vector_to_array(true_pjb_c_eta_lab,es.get_pjb(TauType::C, frc, Var::Eta, EType::MCTrue));
                                vector_to_array(true_pjq_c_pt_lab,es.get_pjq(TauType::C, frc, Var::Pt, EType::MCTrue));
                                vector_to_array(true_pjb_c_pt_lab,es.get_pjb(TauType::C, frc, Var::Pt, EType::MCTrue));
                                vector_to_array(true_pjq_c_p_lab,es.get_pjq(TauType::C, frc, Var::P, EType::MCTrue));
                                vector_to_array(true_pjb_c_p_lab,es.get_pjb(TauType::C, frc, Var::P, EType::MCTrue));
                                qj_c = es.get_qj(TauType::C, frc);
                                true_qj_c[0]= qj_c.Eta();
                                qb_c = es.get_qb(TauType::C, frc);
                                true_qb_c[0]= qb_c.Eta();
                                true_sf_c_lb = es.get_suppression_factor(TauType::C, fra);

				// Reco
				reco_tau1_c =es.get_tau1(TauType::C, EType::MCReco);

                                if(verbosity) cout<<"tau1 CM "<< true_tau1_c<<endl;
			break;

			default: 
			break;
			}//case
		}//three tau's

		true_x = es.get_parton_x(EType::MCTrue);
		true_Q2 = es.get_Q2(EType::MCTrue);
                reco_x = es.get_parton_x(EType::MCReco);
                reco_Q2 = es.get_Q2(EType::MCReco);		

		es.set_scattered_lepton(Frame::Lab,EType::MCTrue);
		true_e_y = es.get_e_y(Frame::Lab, EType::MCTrue);
		true_e_pt = es.get_e_pt(Frame::Lab, EType::MCTrue);
		true_e_eta = es.get_e_eta(Frame::Lab, EType::MCTrue);
		es.set_scattered_lepton(Frame::Breit,EType::MCTrue);
		true_e_y_br = es.get_e_y(Frame::Breit, EType::MCTrue);
		true_e_pt_br = es.get_e_pt(Frame::Breit, EType::MCTrue);
		true_e_eta_br = es.get_e_eta(Frame::Breit, EType::MCTrue);
		es.set_scattered_lepton(Frame::CM, EType::MCTrue);
		true_e_y_cm = es.get_e_y(Frame::CM, EType::MCTrue);
		true_e_pt_cm = es.get_e_pt(Frame::CM, EType::MCTrue);
		true_e_eta_cm = es.get_e_eta(Frame::CM, EType::MCTrue);

                es.set_scattered_lepton(Frame::Lab,EType::MCReco);
                reco_e_y = es.get_e_y(Frame::Lab, EType::MCReco);
                reco_e_pt = es.get_e_pt(Frame::Lab, EType::MCReco);
                reco_e_eta = es.get_e_eta(Frame::Lab, EType::MCReco);
                es.set_scattered_lepton(Frame::Breit,EType::MCReco);
                reco_e_y_br = es.get_e_y(Frame::Breit, EType::MCReco);
                reco_e_pt_br = es.get_e_pt(Frame::Breit, EType::MCReco);
                reco_e_eta_br = es.get_e_eta(Frame::Breit, EType::MCReco);
                es.set_scattered_lepton(Frame::CM, EType::MCReco);
                reco_e_y_cm = es.get_e_y(Frame::CM, EType::MCReco);
                reco_e_pt_cm = es.get_e_pt(Frame::CM, EType::MCReco);
                reco_e_eta_cm = es.get_e_eta(Frame::CM, EType::MCReco);

		true_jet_y = es.get_jet_y(Frame::Lab, EType::MCTrue);
		true_jet_y_br = es.get_jet_y(Frame::Breit, EType::MCTrue);
		true_jet_y_cm = es.get_jet_y(Frame::CM, EType::MCTrue);
		true_jet_pt = es.get_jet_pt(Frame::Lab, EType::MCTrue);
		true_jet_pt_br = es.get_jet_pt(Frame::Breit, EType::MCTrue);
 		true_jet_pt_cm = es.get_jet_pt(Frame::CM, EType::MCTrue);
		true_jet_eta = es.get_jet_eta(Frame::Lab, EType::MCTrue);
		true_jet_eta_br = es.get_jet_eta(Frame::Breit, EType::MCTrue);
		true_jet_eta_cm = es.get_jet_eta(Frame::CM, EType::MCTrue);
		true_jet_p = es.get_jet_p(Frame::Lab,EType::MCTrue);
		true_jet_p_br = es.get_jet_p(Frame::Breit,EType::MCTrue);
		true_jet_p_cm = es.get_jet_p(Frame::CM,EType::MCTrue);
		true_jet_e = es.get_jet_e(Frame::Lab,EType::MCTrue);
		true_jet_e_br = es.get_jet_e(Frame::Breit,EType::MCTrue);
		true_jet_e_cm = es.get_jet_e(Frame::CM,EType::MCTrue);
                true_jet_m = es.get_jet_mass(Frame::Lab,EType::MCTrue);
                true_jet_m_br = es.get_jet_mass(Frame::Breit,EType::MCTrue);
                true_jet_m_cm = es.get_jet_mass(Frame::CM,EType::MCTrue);
		
                reco_jet_y = es.get_jet_y(Frame::Lab, EType::MCReco);
                reco_jet_y_br = es.get_jet_y(Frame::Breit, EType::MCReco);
                reco_jet_y_cm = es.get_jet_y(Frame::CM, EType::MCReco);
                reco_jet_pt = es.get_jet_pt(Frame::Lab, EType::MCReco);
                reco_jet_pt_br = es.get_jet_pt(Frame::Breit, EType::MCReco);
                reco_jet_pt_cm = es.get_jet_pt(Frame::CM, EType::MCReco);
                reco_jet_eta = es.get_jet_eta(Frame::Lab, EType::MCReco);
                reco_jet_eta_br = es.get_jet_eta(Frame::Breit, EType::MCReco);
                reco_jet_eta_cm = es.get_jet_eta(Frame::CM, EType::MCReco);
                reco_jet_p = es.get_jet_p(Frame::Lab,EType::MCReco);
                reco_jet_p_br = es.get_jet_p(Frame::Breit,EType::MCReco);
                reco_jet_p_cm = es.get_jet_p(Frame::CM,EType::MCReco);
                reco_jet_e = es.get_jet_e(Frame::Lab,EType::MCReco);
                reco_jet_e_br = es.get_jet_e(Frame::Breit,EType::MCReco);
                reco_jet_e_cm = es.get_jet_e(Frame::CM,EType::MCReco);
                reco_jet_m = es.get_jet_mass(Frame::Lab,EType::MCReco);
                reco_jet_m_br = es.get_jet_mass(Frame::Breit,EType::MCReco);
                reco_jet_m_cm = es.get_jet_mass(Frame::CM,EType::MCReco);

		if(verbosity>5)cout<<"x "<<true_x<<" e_y "<< true_e_y<<endl;

		if(true_tau1_a>true_tau1_b) ++tau_ab;
		if(true_tau1_b>true_tau1_c) ++tau_bc;
                outtree->Fill();
//		cout<<"pjq a dot "<<true_pjq_a_dot[0]<<endl;
//		float deltar = sqrt(pow(tj_phi-rj_phi,2)+pow(tj_y-rj_eta,2));
//		htj_yphi->Fill(tj_phi,tj_y);
	
		if(m_do_histos)
		{
		
		}
	}
	cout<<"nq "<<nq <<" nb "<<nb<<endl;
	cout<<"tau a>b occurences "<<tau_ab<<endl;
	cout<<"tau b>c occurences "<<tau_bc<<endl; 
	
/*
	bool draw=true;
	if(draw){
	TCanvas* cpt=new TCanvas();
	cpt->cd(); htjpt->Draw();//htjpt->Draw();
	}
*/
	++eventnum;
	return 1;
}

int JetAnalysis::set_branches_mctree(){
	mctree->SetBranchAddress("event", &tevent);
	mctree->SetBranchAddress("eventS",&sevent);
	return 1;
}

void JetAnalysis::vector_to_array(float* arr, const std::vector<float>& vec)
{
	int size = vec.size();	
	//cout<<"size "<<size<<endl;
	for (int i=0;i<size; ++i)
	{arr[i]=vec[i];
	//cout<<arr[i]<<endl;
	}
}

void JetAnalysis::vector_to_array(int* arr, const std::vector<int>& vec)
{
        int size = vec.size();
        //cout<<"size "<<size<<endl;
        for (int i=0;i<size; ++i)
        {arr[i]=vec[i];
        //cout<<arr[i]<<endl;
        }
}



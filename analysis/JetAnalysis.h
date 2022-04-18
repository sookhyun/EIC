/**
 * @file JetAnalysis.h
 * @author Sookhyun Lee
 * @date May 2020
 */

#ifndef __JETANALYSIS_H_
#define __JETANALYSIS_H_

#include <eicsmear/erhic/EventPythia.h>
#include <eicsmear/smear/EventS.h>


#include <string>
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h>

class TFile;
class TH1D;
class TH2D;
class TGraph2DErrors;

static const int NPART_MAX = 100;

class JetAnalysis
{
  public:
    JetAnalysis();
    virtual ~JetAnalysis() {}

    int open(std::string fname);
    int init();
    int analyze();
    int save_results();    
    void set_qjet_mode(bool do_min){m_qjet_mode=do_min;}
    void set_xQ2_method(int method){m_Q2_method=method;}
    void set_smear_mode(int mode) {m_smear_mode=mode;}
    void set_tracking_mode(int mode){m_tracking_mode=mode;}
    void set_pt_cut(float ptcut){m_pt_cut=ptcut;}
    void set_verbosity(int _verb){verbosity=_verb;}
    void do_histos(bool _do){m_do_histos = _do;}


  private:
    void vector_to_array(float* arr, const std::vector<float>& vec);
    void vector_to_array(int* arr, const std::vector<int>& vec);
    int set_branches_mctree();

    unsigned int eventnum;
    float true_tau1_a;
    float true_tau1_b;
    float true_tau1_c;
    float true_rapgap_a;
    float true_rapgap_b;
    float true_rapgap_c;

    float reco_tau1_a;
    float reco_tau1_b;
    float reco_tau1_c;
    float reco_rapgap_a;
    float reco_rapgap_b;
    float reco_rapgap_c;

    float true_x;
    float true_Q2;
    float true_e_y;	float true_e_y_br;	float true_e_y_cm;
    float true_e_pt;   float true_e_pt_br;    float true_e_pt_cm;
    float true_e_eta;  float true_e_eta_br;   float true_e_eta_cm;
    float true_jet_y;  float true_jet_y_br;   float true_jet_y_cm;
    float true_jet_pt; float true_jet_pt_br;  float true_jet_pt_cm;
    float true_jet_eta;float true_jet_eta_br; float true_jet_eta_cm;
    float true_jet_p;  float true_jet_p_br;   float true_jet_p_cm;
    float true_jet_e;  float true_jet_e_br;   float true_jet_e_cm;
    float true_jet_m;  float true_jet_m_br;   float true_jet_m_cm;

    float reco_x;
    float reco_Q2;
    float reco_e_y;     float reco_e_y_br;      float reco_e_y_cm;
    float reco_e_pt;   float reco_e_pt_br;    float reco_e_pt_cm;
    float reco_e_eta;  float reco_e_eta_br;   float reco_e_eta_cm;
    float reco_jet_y;  float reco_jet_y_br;   float reco_jet_y_cm;
    float reco_jet_pt; float reco_jet_pt_br;  float reco_jet_pt_cm;
    float reco_jet_eta;float reco_jet_eta_br; float reco_jet_eta_cm;
    float reco_jet_p;  float reco_jet_p_br;   float reco_jet_p_cm;
    float reco_jet_e;  float reco_jet_e_br;   float reco_jet_e_cm;
    float reco_jet_m;  float reco_jet_m_br;   float reco_jet_m_cm;

    float true_sf_a_lb; float true_sf_a_br; float true_sf_a_cm;
    float true_sf_b_lb; float true_sf_b_br; float true_sf_b_cm;
    float true_sf_c_lb; float true_sf_c_br; float true_sf_c_cm;
    float reco_sf_a_lb; float reco_sf_a_br; float reco_sf_a_cm;
    float reco_sf_b_lb; float reco_sf_b_br; float reco_sf_b_cm;
    float reco_sf_c_lb; float reco_sf_c_br; float reco_sf_c_cm;

    float true_qf_a; float true_qf_b; float true_qf_c;
    // a: lab frame, b: breit frame, c: CM frame
    float true_qj_a[2]; float true_qj_b[2]; float true_qj_c[2];
    float true_qj_a_br[2]; float true_qj_b_br[2]; float true_qj_c_br[2];
    float true_qj_a_cm[2]; float true_qj_b_cm[2]; float true_qj_c_cm[2];
    float true_qb_a[2]; float true_qb_b[2]; float true_qb_c[2];
    float true_qb_a_br[2]; float true_qb_b_br[2]; float true_qb_c_br[2];
    float true_qb_a_cm[2]; float true_qb_b_cm[2]; float true_qb_c_cm[2];
    int true_pjq_a_n; int true_pjq_b_n; int true_pjq_c_n;
    int true_pjb_a_n; int true_pjb_b_n; int true_pjb_c_n;
    int true_pjq_a_id[NPART_MAX]; int true_pjb_a_id[NPART_MAX];
    float true_pjq_a_dot_lab[NPART_MAX]; float true_pjb_a_dot_lab[NPART_MAX];
    float true_pjq_a_eta_lab[NPART_MAX]; float true_pjb_a_eta_lab[NPART_MAX];
    float true_pjq_a_pt_lab[NPART_MAX]; float true_pjb_a_pt_lab[NPART_MAX];
    float true_pjq_a_p_lab[NPART_MAX]; float true_pjb_a_p_lab[NPART_MAX];
    float true_pjq_a_dot_br[NPART_MAX]; float true_pjb_a_dot_br[NPART_MAX];
    float true_pjq_a_eta_br[NPART_MAX]; float true_pjb_a_eta_br[NPART_MAX];
    float true_pjq_a_pt_br[NPART_MAX]; float true_pjb_a_pt_br[NPART_MAX];
    float true_pjq_a_p_br[NPART_MAX]; float true_pjb_a_p_br[NPART_MAX];
    float true_pjq_a_dot_cm[NPART_MAX]; float true_pjb_a_dot_cm[NPART_MAX];
    float true_pjq_a_eta_cm[NPART_MAX]; float true_pjb_a_eta_cm[NPART_MAX];
    float true_pjq_a_pt_cm[NPART_MAX]; float true_pjb_a_pt_cm[NPART_MAX];
    float true_pjq_a_p_cm[NPART_MAX]; float true_pjb_a_p_cm[NPART_MAX];

    float true_pjq_b_dot_br[NPART_MAX]; float true_pjb_b_dot_br[NPART_MAX];
    float true_pjq_b_eta_br[NPART_MAX]; float true_pjb_b_eta_br[NPART_MAX];
    float true_pjq_b_pt_br[NPART_MAX]; float true_pjb_b_pt_br[NPART_MAX];
    float true_pjq_b_p_br[NPART_MAX]; float true_pjb_b_p_br[NPART_MAX];
    float true_pjq_b_dot_lab[NPART_MAX]; float true_pjb_b_dot_lab[NPART_MAX];
    float true_pjq_b_eta_lab[NPART_MAX]; float true_pjb_b_eta_lab[NPART_MAX];
    float true_pjq_b_pt_lab[NPART_MAX]; float true_pjb_b_pt_lab[NPART_MAX];
    float true_pjq_b_p_lab[NPART_MAX]; float true_pjb_b_p_lab[NPART_MAX];
    float true_pjq_b_dot_cm[NPART_MAX]; float true_pjb_b_dot_cm[NPART_MAX];
    float true_pjq_b_eta_cm[NPART_MAX]; float true_pjb_b_eta_cm[NPART_MAX];
    float true_pjq_b_pt_cm[NPART_MAX]; float true_pjb_b_pt_cm[NPART_MAX];
    float true_pjq_b_p_cm[NPART_MAX]; float true_pjb_b_p_cm[NPART_MAX];

    float true_pjq_c_dot_cm[NPART_MAX]; float true_pjb_c_dot_cm[NPART_MAX];
    float true_pjq_c_eta_cm[NPART_MAX]; float true_pjb_c_eta_cm[NPART_MAX];
    float true_pjq_c_pt_cm[NPART_MAX]; float true_pjb_c_pt_cm[NPART_MAX];
    float true_pjq_c_p_cm[NPART_MAX]; float true_pjb_c_p_cm[NPART_MAX];
    float true_pjq_c_dot_lab[NPART_MAX]; float true_pjb_c_dot_lab[NPART_MAX];
    float true_pjq_c_eta_lab[NPART_MAX]; float true_pjb_c_eta_lab[NPART_MAX];
    float true_pjq_c_pt_lab[NPART_MAX]; float true_pjb_c_pt_lab[NPART_MAX];
    float true_pjq_c_p_lab[NPART_MAX]; float true_pjb_c_p_lab[NPART_MAX];
    float true_pjq_c_dot_br[NPART_MAX]; float true_pjb_c_dot_br[NPART_MAX];
    float true_pjq_c_eta_br[NPART_MAX]; float true_pjb_c_eta_br[NPART_MAX];
    float true_pjq_c_pt_br[NPART_MAX]; float true_pjb_c_pt_br[NPART_MAX];
    float true_pjq_c_p_br[NPART_MAX]; float true_pjb_c_p_br[NPART_MAX];

    // arrays for a only, b and c to be implemented laer 
    float reco_qj_a[2]; float reco_qj_b[2]; float reco_qj_c[2];
    float reco_qj_a_br[2]; float reco_qj_b_br[2]; float reco_qj_c_br[2];
    float reco_qj_a_cm[2]; float reco_qj_b_cm[2]; float reco_qj_c_cm[2];
    float reco_qb_a[2]; float reco_qb_b[2]; float reco_qb_c[2];
    float reco_qb_a_br[2]; float reco_qb_b_br[2]; float reco_qb_c_br[2];
    float reco_qb_a_cm[2]; float reco_qb_b_cm[2]; float reco_qb_c_cm[2];
    int reco_pjq_a_n; int reco_pjq_b_n; int reco_pjq_c_n;
    int reco_pjb_a_n; int reco_pjb_b_n; int reco_pjb_c_n;
    int reco_pjq_a_id[NPART_MAX]; int reco_pjb_a_id[NPART_MAX];
    float reco_pjq_a_dot_lab[NPART_MAX]; float reco_pjb_a_dot_lab[NPART_MAX];
    float reco_pjq_a_eta_lab[NPART_MAX]; float reco_pjb_a_eta_lab[NPART_MAX];
    float reco_pjq_a_pt_lab[NPART_MAX]; float reco_pjb_a_pt_lab[NPART_MAX];
    float reco_pjq_a_p_lab[NPART_MAX]; float reco_pjb_a_p_lab[NPART_MAX];
    float reco_pjq_a_dot_br[NPART_MAX]; float reco_pjb_a_dot_br[NPART_MAX];
    float reco_pjq_a_eta_br[NPART_MAX]; float reco_pjb_a_eta_br[NPART_MAX];
    float reco_pjq_a_pt_br[NPART_MAX]; float reco_pjb_a_pt_br[NPART_MAX];
    float reco_pjq_a_p_br[NPART_MAX]; float reco_pjb_a_p_br[NPART_MAX];
    float reco_pjq_a_dot_cm[NPART_MAX]; float reco_pjb_a_dot_cm[NPART_MAX];
    float reco_pjq_a_eta_cm[NPART_MAX]; float reco_pjb_a_eta_cm[NPART_MAX];
    float reco_pjq_a_pt_cm[NPART_MAX]; float reco_pjb_a_pt_cm[NPART_MAX];
    float reco_pjq_a_p_cm[NPART_MAX]; float reco_pjb_a_p_cm[NPART_MAX];


    bool matched;
    TFile* infile_mctrue;
    TFile* infile_mcreco;
    TFile* outfile;
    TTree* mctree;
    TTree* outtree;
    TH1D* htjpt;
    TH2D* htj_yphi;

    erhic::EventPythia* tevent;
    Smear::Event* sevent;

    float m_pt_cut;
    bool m_do_histos;
    bool m_qjet_mode;
    int nfiles;
    int m_Q2_method;
    int m_smear_mode;
    int m_tracking_mode;
    int verbosity;

};

#endif

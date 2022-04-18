/**
 * @file EventShape.h
 * @author Sookhyun Lee
 * @date May 2020
 */

#ifndef _EVENTSHAPE_H_
#define _EVENTSHAPE_H_

#include <eicsmear/erhic/EventPythia.h>
#include <eicsmear/smear/EventS.h>
#include <iostream>
#include <map>
#include <stdio.h>

#include <TLorentzVector.h>
#include <TRotation.h>
#include <TRandom2.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TF1.h>
#include "Opts.h"

class EventShape {
  public:
    EventShape(){}
    EventShape(erhic::EventPythia &tevent)
       : m_tevent(&tevent)
	,m_pt_cut(0.)
	,m_smear_mode(0)
	,m_tracking_mode(0)
	,verbosity(0)       
    {
    ftreff = new TF1("ftreff","[0]/(exp(-((x-[2])/[1]))+1.)",0.,20.0);
//    ftreff->SetParameters(0.99,0.3,1.4);
    rnd = new TRandom2();
    rnd->SetSeed();
    }

    EventShape(erhic::EventPythia &tevent, Smear::Event &sevent)
       : m_tevent(&tevent)
	,m_sevent(&sevent)
	,m_pt_cut(0.)
	,m_smear_mode(0)
	,m_tracking_mode(0)
	,verbosity(0)
    {
    ftreff = new TF1("ftreff","[0]/(exp(-((x-[2])/[1]))+1.)",0.,20.0);
//    ftreff->SetParameters(0.99,0.3,1.4);
    rnd = new TRandom2();
    rnd->SetSeed();
    }


    double get_tau1(TauType ttype, EType etype);
    double get_parton_x(EType etype);
    double get_Q2(EType etype);
    double get_e_y(Frame f, EType etype);
    double get_e_pt(Frame f, EType etype);
    double get_e_eta(Frame f, EType etype);
    double get_jet_y(Frame f, EType etype);
    double get_jet_pt(Frame f, EType etype);
    double get_jet_eta(Frame f, EType etype);
    double get_jet_p(Frame f, EType etype);
    double get_jet_e(Frame f, EType etype);
    double get_jet_mass(Frame f, EType etype);
    /// This block has to be called immediately after calling get_tau1
    double get_suppression_factor(TauType ttype, Frame f);
    double get_tau1_q_fraction(TauType ttype);
    TLorentzVector get_qj(TauType ttype, Frame f);
    TLorentzVector get_qb(TauType ttype, Frame f);
    std::vector<float> get_pjq(TauType ttype, Frame f, Var var, EType etype);
    std::vector<float> get_pjb(TauType ttype, Frame f, Var var, EType etype);  
    std::vector<int> get_pjq_id(TauType ttype);
    std::vector<int> get_pjb_id(TauType ttype);

    void set_qjet_mode(bool do_min){m_do_minimization=do_min;}
    void set_pt_cut(float ptcut){m_pt_cut=ptcut;}
    // smear_mode
    // 0:true p4, 
    // 1:default reco p3, true E, (effect of tracking)
    // 2:default reco p3, E form pion mass 
    void set_smear_mode(int smear_mode){m_smear_mode = smear_mode;}
    void set_tracking_mode(int tracking_mode);
    void set_xQ2_method(int method){m_xQ2_method = method;}
    void set_scattered_lepton(Frame f, EType etype);
    void set_verbosity(int _verb){verbosity=_verb;}

  private:
    double get_tau1_true(TauType ttype);
    double get_tau1_reco(TauType ttype);
    double compute_tau1(EType etype);
    double compute_tau1_jet(EType etype);
    void set_suppression_factor(TauType ttype);
    void set_passed();
    void set_q_jet(TauType ttype, EType etype);
    void set_q_beam(TauType ttype, EType etype);
    void find_qj(EType etype);
    void set_yes();
    void set_jet(EType etype);
    void set_particles(EType etype);
    void set_pj(TauType ttype);
    void clear_pj();

    TLorentzVector get_scattered_lepton(Frame f, EType etype);
    TLorentzVector get_jet(Frame f, EType etype);    
    TLorentzVector get_smeared_chargedH(unsigned int it, const TLorentzVector& p4);
    TLorentzVector get_smeared_EM(unsigned int it, const TLorentzVector& p4);
    TLorentzVector get_smeared_calo_only(unsigned int it, const TLorentzVector& p4);

    void set_kinematics(EType etype);
    void set_boson(const TLorentzVector& pboson)
	{m_boson=pboson;}
    void set_beam_hadron(const TLorentzVector& phadron)
	{m_beamhadron=phadron;}
    void set_beam_lepton(const TLorentzVector& plepton)
	{m_beamlepton = plepton;}
    void set_scattered_hadron(const TLorentzVector& pscathadron)
	{m_scathadron = pscathadron;}
    void set_scattered_lepton(const TLorentzVector& pscatlepton)
	{m_scatlepton =pscatlepton;}
    // set boost and rotatoin vector
    void set_parton_x(EType etype);
    void set_Q2(EType etype);
    void set_boost(Frame f, EType etype);
    void set_boost_lab2CM(const TLorentzVector& , const TLorentzVector& );
    void set_boost_lab2Breit(const TLorentzVector& , const TLorentzVector& );
    bool pass_cuts(unsigned int it, bool is_true_part);    
    bool is_pid_3sig(double p, double eta);

    // transport particles 
    TLorentzVector boost_particle_lab2CM(const TLorentzVector& pini);
    TLorentzVector boost_particle_lab2Breit(const TLorentzVector& pini);
    TLorentzVector boost_particle_Breit2lab(const TLorentzVector& pini);
    TLorentzVector boost_particle_CM2lab(const TLorentzVector& pini);

  protected:

    TF1* ftreff;
    TRandom2* rnd;
    TVector3 m_q_jet;
    TVector3 m_b1_cm;
    TVector3 m_b1_br;
    TVector3 m_b2_br;
    TRotation m_r1_cm;
    TRotation m_r1_br;
    // True info from generator
    TLorentzVector m_beamhadron;
    TLorentzVector m_beamlepton;
    TLorentzVector m_scathadron;
    TLorentzVector m_scatlepton; 
    TLorentzVector m_boson;

    TLorentzVector m_qj;
    TLorentzVector m_qb;

    // q axes in lab frame
    TLorentzVector m_qj_a;
    TLorentzVector m_qb_a;
    TLorentzVector m_qj_b;
    TLorentzVector m_qb_b;
    TLorentzVector m_qj_c;
    TLorentzVector m_qb_c;
    

    // Jet
    TLorentzVector m_true_jet_lb;
    TLorentzVector m_true_jet_br;
    TLorentzVector m_true_jet_cm;
    TLorentzVector m_reco_jet_lb;
    TLorentzVector m_reco_jet_br;
    TLorentzVector m_reco_jet_cm;

    TLorentzVector m_true_e_lb;
    TLorentzVector m_true_e_br;
    TLorentzVector m_true_e_cm;
    TLorentzVector m_reco_e_lb;
    TLorentzVector m_reco_e_br;
    TLorentzVector m_reco_e_cm;
    
    /// all particles in an event
    std::map<int, TLorentzVector> m_particles;
    std::map<int, bool> m_passed;

    /// particles belonging to quark jet
    std::vector<int> m_pjq_a_id;
    std::vector<TLorentzVector> m_pjq_a_lab;
    std::vector<TLorentzVector> m_pjq_a_breit;
    std::vector<TLorentzVector> m_pjq_a_cm;
    std::vector<TLorentzVector> m_pjq_b_lab;
    std::vector<TLorentzVector> m_pjq_b_breit;
    std::vector<TLorentzVector> m_pjq_b_cm;
    std::vector<TLorentzVector> m_pjq_c_lab;
    std::vector<TLorentzVector> m_pjq_c_breit;
    std::vector<TLorentzVector> m_pjq_c_cm;

    /// particles belonging to beam jet   
    std::vector<int> m_pjb_a_id;  
    std::vector<TLorentzVector> m_pjb_a_lab;
    std::vector<TLorentzVector> m_pjb_a_breit;
    std::vector<TLorentzVector> m_pjb_a_cm;
    std::vector<TLorentzVector> m_pjb_b_lab;
    std::vector<TLorentzVector> m_pjb_b_breit;
    std::vector<TLorentzVector> m_pjb_b_cm;
    std::vector<TLorentzVector> m_pjb_c_lab;
    std::vector<TLorentzVector> m_pjb_c_breit;
    std::vector<TLorentzVector> m_pjb_c_cm;

    double m_sf_a_lb; double m_sf_a_br; double m_sf_a_cm;
    double m_sf_b_lb; double m_sf_b_br; double m_sf_b_cm;
    double m_sf_c_lb; double m_sf_c_br; double m_sf_c_cm;
    double m_qf_a;  double m_qf_b; double m_qf_c;

    int m_npjq_a;
    int m_npjb_a;
    int m_npjq_b;
    int m_npjb_b;
    int m_npjq_c;
    int m_npjb_c;
    
    erhic::EventPythia *m_tevent;  
    Smear::Event *m_sevent;

    double m_true_parton_x;  double m_reco_parton_x;   
    double m_true_Q2; 	     double m_reco_Q2;

    float m_pt_cut;
    bool m_save_histos;
    bool m_do_minimization;
    int m_smear_mode;
    int m_tracking_mode;
    int m_xQ2_method;
    int verbosity;
    

};
#endif


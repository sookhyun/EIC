/**
 * @file EventShape.cc
 * @author Sookhyun Lee
 * @date May 2020
 */

#include "EventShape.h"
#include "TMinuit.h"

using namespace std;
std::map<int, TLorentzVector> m_mommap;
std::map<int, bool> m_yes;
double threshold = 1e-5;

double Term(double *x, double *p)
{
  const double y1 = p[0];//costheta
  const double y2 = p[1];//phi
//  double cosphi = 

//  double z=-y1/sqrt(1.-y1*y1);
//  if((y1*y1)>1.) z=-y1/1e-09;
  TLorentzVector pi(x[0],x[1],x[2],x[3]);
//  TLorentzVector dcostheta(z*cos(y2),z*sin(y2), 1.,0.);
  TLorentzVector dcostheta(sqrt(1.-y1*y1)*cos(y2),sqrt(1.-y1*y1)*sin(y2),y1, 1.);
  TLorentzVector dcosphi(0.,0.,0.,0.);
  double val = dcostheta.Dot(pi) + dcosphi.Dot(pi);

  return val;
}

void CostFn(Int_t &npar, Double_t *gin, Double_t &f, double *p, Int_t iflag) 
{

  double term;
  double cost = 0.;
  double x[4];// particle 4-vector
//cout<<"m_mommap size "<<m_mommap.size()<<endl;
//cout<<" theta "<<p[0]<<" phi "<<p[1]<<endl;
        for (std::map<int,TLorentzVector>::iterator iter=m_mommap.begin(); iter!=m_mommap.end(); ++iter)
	{
        TLorentzVector p4 = iter->second;
//        double qbdotp =m_qb.Dot(p4); double qjdotp = m_qj.Dot(p4);
//        double term = (qbdotp < qjdotp)? qbdotp : qjdotp ;
//        if(verbosity>7) cout<<"term " <<term<<endl;
//        tau1 += (2.*term/m_Q2);
	 // skip if particle belongs to beam jet
	//cout<<"m_yes "<<m_yes[iter->first]<<endl;
	 if(!m_yes[iter->first]) continue;
	 x[0] = p4.Px(); x[1]= p4.Py(); x[2]=p4.Pz(); x[3]=p4.E();
	//cout<<" px "<<x[0]<<" py "<<x[1]<< " pz "<<x[3]<<endl;
 	 term =Term(x,p);
    	 cost += term;
	//cout<<" term "<<term<<endl;
        }
  
  f = cost;
//  cout<<"     cost "<<cost<<endl;
}

double EventShape::get_tau1(TauType ttype, EType etype)
{
	double tau1;
	set_kinematics(etype);
	set_particles(etype); // in lab frame for set_pj & set_suppression_factor
	set_passed(); // detector coverage
	if(etype == MCTrue) tau1 = get_tau1_true(ttype);
	else if (etype == MCReco) tau1 = get_tau1_reco(ttype);
	return tau1;
}

double EventShape::get_tau1_true(TauType ttype)
{
	clear_pj();
	m_mommap.clear();
	m_yes.clear();
	// set boost & rotation matrix
		set_boost(Frame::Breit, EType::MCTrue);
		set_boost(Frame::CM, EType::MCTrue);

        // boost particles
	double tau1=0;
	unsigned int nentries=m_tevent->GetNTracks();
	for(unsigned int it =0; it< nentries;++it)
	{
		if(it<3) continue;
		const Particle *part=m_tevent->GetTrack(it);
		if(part->GetStatus() !=1) continue;
		if(part->GetE() == m_scatlepton.E()) continue;
		//if(!pass_cuts(it, true)) continue;
		TLorentzVector partp4 = part->PxPyPzE();
		
		// boost particles
		switch(ttype)
		{
			case (TauType::A):
			partp4 = boost_particle_lab2Breit(partp4);
			break;
			case (TauType::B):
			partp4 = boost_particle_lab2Breit(partp4);
			break;
			case (TauType::C):
			partp4 = boost_particle_lab2CM(partp4);
			break;
			default:
			break;
		}
		// add it to the particle momentum containter 
		m_mommap.insert(make_pair((int) it,partp4));
		m_yes.insert(make_pair((int) it,false));
	}

	set_q_beam(ttype, EType::MCTrue);
	set_q_jet(ttype, EType::MCTrue);
	set_yes();
	set_pj(ttype); // call get_pjq()/get_pjb() right after get_tau1_true() 
	tau1 = compute_tau1(EType::MCTrue);	
	if(ttype==TauType::A) set_jet(EType::MCTrue);
	set_suppression_factor(ttype); // call get_suppression_factor() right after get_tau1_true()
        return tau1;
}


double EventShape::get_tau1_reco(TauType ttype)
{
        clear_pj();
        m_mommap.clear();
        m_yes.clear();
        // set boost & rotation matrix
                set_boost(Frame::Breit, EType::MCReco);
                set_boost(Frame::CM, EType::MCReco);

        // boost particles
        double tau1=0;
        unsigned int nentries=m_sevent->GetNTracks();

        for(unsigned int it =0; it< nentries;++it)
        {
                if(it<3) continue;

                const Smear::ParticleMCS *part = m_sevent->GetTrack(it);
                if (!part) continue;
		if(part->GetStatus() !=1) continue;
                double pE = part->GetE(); double pP = part->GetP();
                if((fabs(pE) < threshold) && (fabs(pP) < threshold)) continue;
                if(part->GetE() == m_scatlepton.E()) continue;
                if(!pass_cuts(it, true)) continue;
                TLorentzVector partp4 = part->PxPyPzE();

                double mass = m_tevent->GetTrack(it)->GetM();

                if(m_smear_mode==0)
                {
                partp4.SetPxPyPzE(m_tevent->GetTrack(it)->GetPx(),
                              m_tevent->GetTrack(it)->GetPy(),
                              m_tevent->GetTrack(it)->GetPz(),
                              m_tevent->GetTrack(it)->GetE());
                }
                else
                {
			int pid = std::abs(m_tevent->GetTrack(it)->GetPdgCode());
                	if(fabs(pE) < threshold)
                	{//Tracking info only
			 double energy=0.;
			 if(m_smear_mode==3 || m_smear_mode==4)
			 { // 3 sig separation
                //partp4.SetPxPyPzE(m_tevent->GetTrack(it)->GetPx(),
                //              m_tevent->GetTrack(it)->GetPy(),
                //              m_tevent->GetTrack(it)->GetPz(),
                //              m_tevent->GetTrack(it)->GetE());

			  mass=0.139;
			  double p3 = m_tevent->GetTrack(it)->GetP();
			  double eta = m_tevent->GetTrack(it)->GetEta();
			  bool is_pkpi = ((pid==211) || (pid==321) || (pid ==2212));
			  bool is_3sig = is_pid_3sig(p3, eta);
			  //if(is_3sig && is_pkpi){ mass = (pid==321) ? 0.4937 : 0.9383;}
			  if(verbosity>7)cout<<"track only pid  " <<pid<<endl;

			  if (is_3sig && is_pkpi)
			  {
			  	 mass= m_tevent->GetTrack(it)->GetM();
				 energy = sqrt(pow(p3,2) + pow(mass,2));
			  }
			  
			//energy = m_tevent->GetTrack(it)->GetE();
			 }
			 if (m_smear_mode==2){ mass =0.139; energy = sqrt(pow(pP,2) + pow(mass,2));}
			 if (m_smear_mode==1){ //energy = m_tevent->GetTrack(it)->GetE();}
					       mass = m_tevent->GetTrack(it)->GetM(); 
					       energy = sqrt(pow(pP,2)+pow(mass,2));}// ensure time-like 4 vector
	                 partp4.SetE(energy);

			 if (m_tracking_mode>0) 
			 {
			 double pt = m_tevent->GetTrack(it)->GetPt();
			 double eff = ftreff->Eval(pt);
			 if (eff < rnd->Uniform(0.,1.)) continue; // no track, no calo, skip
			 }

	                }//end track only
	                else if(fabs(pP) < threshold)
	                {//Calo info only

			if(verbosity>7) cout<<"calo only pid  " <<pid<<endl;
			partp4 = get_smeared_calo_only(it, partp4);

	                }// end calo only
			else if(fabs(pP) >= threshold && fabs(pE) >= threshold)
			{  // both tracking and calo 	
                        if(pid!=11)
			{
			if(verbosity>7)cout<<"chargedH pid  " <<pid<<endl;
                        partp4 = get_smeared_chargedH(it,partp4 );
                        }else
                        partp4 = get_smeared_EM(it, partp4);
			//else cout<<"wrong "<<pid<<endl;
			}
		}

                // boost particles
                switch(ttype)
                {
                        case (TauType::A):
                        partp4 = boost_particle_lab2Breit(partp4);
                        break;
                        case (TauType::B):
                        partp4 = boost_particle_lab2Breit(partp4);
                        break;
                        case (TauType::C):
                        partp4 = boost_particle_lab2CM(partp4);
                        break;
                        default:
                        break;
                }
                // add it to the particle momentum containter 
                m_mommap.insert(make_pair((int) it,partp4));
                m_yes.insert(make_pair((int) it,false));
        }

        set_q_beam(ttype, EType::MCReco);
        set_q_jet(ttype, EType::MCReco);
        set_yes();
        set_pj(ttype); // call get_pjq()/get_pjb() right after get_tau1_true()
        tau1 = compute_tau1(EType::MCReco);
        if(ttype==TauType::A) set_jet(EType::MCReco);
        set_suppression_factor(ttype); // call get_suppression_factor() right after get_tau1_reco()

        return tau1;

	return 0.;
}

void EventShape::set_suppression_factor(TauType ttype)
{
        double Q2= m_true_Q2;
	TLorentzVector qj_lb;
	TLorentzVector qj_br;
	TLorentzVector qj_cm;
	TLorentzVector qb_lb;
	TLorentzVector qb_br;
	TLorentzVector qb_cm;
	switch(ttype)
	{
		case (TauType::A):	
		qj_lb = m_qj_a;
		qj_br = boost_particle_lab2Breit(qj_lb);
		qj_cm = boost_particle_lab2CM(qj_lb);
                qb_lb = m_qb_a;
                qb_br = boost_particle_lab2Breit(qb_lb);
                qb_cm = boost_particle_lab2CM(qb_lb);

		break; 
		case (TauType::B):
		qj_lb = m_qj_b;
		qj_br = boost_particle_lab2Breit(qj_lb);
		qj_cm = boost_particle_lab2CM(qj_lb);
                qb_lb = m_qb_b;
                qb_br = boost_particle_lab2Breit(qb_lb);
                qb_cm = boost_particle_lab2CM(qb_lb);

		break;
		case (TauType::C):
		qj_lb = m_qj_c;
		qj_br = boost_particle_lab2Breit(qj_lb);
		qj_cm = boost_particle_lab2CM(qj_lb);
                qb_lb = m_qb_c;
                qb_br = boost_particle_lab2Breit(qb_lb);
                qb_cm = boost_particle_lab2CM(qb_lb);

		break;
	}
	double tau1_lb=0.; double tau1_br=0.; double tau1_cm=0.;
	double tau1_lb_inacc=0.; double tau1_br_inacc=0.; double tau1_cm_inacc=0.;
	double tau1_q=0.; 
        for (std::map<int,TLorentzVector>::iterator iip=m_particles.begin();
                iip!=m_particles.end(); ++iip)
        {
                TLorentzVector _pp = iip->second;
		TLorentzVector _pp_br = boost_particle_lab2Breit(_pp);
		TLorentzVector _pp_cm = boost_particle_lab2CM(_pp);
                if(verbosity>7) cout<<"p4 "<<_pp.Px()<<" "<<_pp.Py()<<" "<<_pp.Pz()<<" "<<_pp.E()<<endl;
		if(verbosity>7) cout<<"qj "<<qj_lb.Pz()<<" "<<qj_lb.E()<<endl;
		if(verbosity>7) cout<<"qb "<<qb_lb.Pz()<<" "<<qb_lb.E()<<endl;
                double eta=_pp.Eta(); //double phi=_pp.Phi();
                bool passed = ((eta > -3.5) && (eta<4.0)) ? true:false;
		double itau1 = (m_yes[iip->first]) ? qj_lb.Dot(_pp) : qb_lb.Dot(_pp);
		double itau1_br = (m_yes[iip->first]) ? qj_br.Dot(_pp_br) : qb_br.Dot(_pp_br);
		double itau1_cm = (m_yes[iip->first]) ? qj_cm.Dot(_pp_cm) : qb_cm.Dot(_pp_cm);
		itau1 *= (2./Q2);
		itau1_br *= (2./Q2);
		itau1_cm *= (2./Q2);
		tau1_lb += itau1;
		tau1_br += itau1_br;
		tau1_cm += itau1_cm;
		if(verbosity>7)cout<<"eta "<<eta<<" passed "<<passed<<" tau1 "<<tau1_lb<<" type "<<ttype<<endl;
                if (passed)
                {
		tau1_lb_inacc += itau1;
		tau1_br_inacc += itau1_br;
		tau1_cm_inacc += itau1_cm;
                }
		double itau1_q = (m_yes[iip->first]) ? qj_lb.Dot(_pp): 0.;
		itau1_q *= (2./Q2);
		tau1_q += itau1_q;
	}	
	
	switch(ttype)
	{
		case (TauType::A):
		m_sf_a_lb = tau1_lb_inacc/tau1_lb;
		m_sf_a_br = tau1_br_inacc/tau1_br;
		m_sf_a_cm = tau1_cm_inacc/tau1_cm;
		m_qf_a = tau1_q/tau1_lb;
		break;
		case (TauType::B):
		m_sf_b_lb = tau1_lb_inacc/tau1_lb;
		m_sf_b_br = tau1_br_inacc/tau1_br;
		m_sf_b_cm = tau1_cm_inacc/tau1_cm;
		m_qf_b = tau1_q/tau1_lb;
		break;
		case (TauType::C):
		m_sf_c_lb = tau1_lb_inacc/tau1_lb;
		m_sf_c_br = tau1_br_inacc/tau1_br;
		m_sf_c_cm = tau1_cm_inacc/tau1_cm;
		m_qf_c = tau1_q/tau1_lb;
		break;
		default:
		break;
	}
	if (verbosity>7) cout<<"sf_a_lb "<<m_sf_a_lb<<" sf_b_br "<<m_sf_b_br<<endl;
}

void EventShape::set_q_jet(TauType ttype, EType etype)
{	
	double parton_x = (etype == EType::MCTrue) ? m_true_parton_x : m_reco_parton_x;
	// A: jet axis, B: q+xP, C: k
	switch(ttype)
	{
		case (TauType::A):
		find_qj(etype);
		break;
		case (TauType::B):
		m_qj = m_beamlepton - m_scatlepton + parton_x*m_beamhadron;
		m_qj.SetE(m_qj.P());
		m_qj_b = m_qj;
		m_qj = boost_particle_lab2Breit(m_qj);
		break;
		case (TauType::C):
		m_qj = m_beamlepton;
		m_qj.SetE(m_qj.P());
		m_qj_c = m_qj;
		m_qj = boost_particle_lab2CM(m_qj);
		//m_qj.SetPz(-m_qj.Pz());
		break;
		default:
		break;
	}
        if (verbosity) cout<<"qj (px,py,pz,e) true "<< (etype==EType::MCTrue)<<" " << m_qj.Px()<<" "<<
        m_qj.Py()<<" "<<m_qj.Pz()<<" " << m_qj.E()<<endl;

}

void EventShape::set_q_beam(TauType ttype, EType etype)
{
        double parton_x = (etype == EType::MCTrue) ? m_true_parton_x : m_reco_parton_x;
	// A: xP, B: xP, C: P
        switch(ttype)
        {
                case (TauType::A):
                m_qb = parton_x*m_beamhadron;
                m_qb.SetE(m_qb.P());
		m_qb_a = m_qb;
		m_qb = boost_particle_lab2Breit(m_qb);
                break;
                case (TauType::B):
                m_qb = parton_x*m_beamhadron;
                m_qb.SetE(m_qb.P());
		m_qb_b = m_qb;
                m_qb = boost_particle_lab2Breit(m_qb);
                break;
                case (TauType::C):
                m_qb = m_beamhadron;
                m_qb.SetE(m_qb.P());
		m_qb_c = m_qb;
                m_qb = boost_particle_lab2CM(m_qb);
		//m_qb.SetPz(-m_qb.Pz());
                break;
                default:
                break;
        }
	if (verbosity) cout<<"qb (px,py,pz,e)  "<< m_qb.Px()<<" "<<
	m_qb.Py()<<" "<<m_qb.Pz()<<" "<<m_qb.E()<<endl;
}


void EventShape::find_qj(EType etype)
{
        double parton_x = (etype == EType::MCTrue)? m_true_parton_x : m_reco_parton_x;
	// find seeds: truth, anti-kt, 
	// For truth jets, use xP+q as seed jet direction   
	TLorentzVector seedqj = parton_x*m_beamhadron + m_beamlepton - m_scatlepton;
	if(verbosity>7)cout<<"(xP+q)^2= "<< seedqj.M2()<< ", (xP+q)_0= "<< seedqj.E()<<endl;
	seedqj.SetE(seedqj.P());
	m_qj = seedqj;
        m_qj = boost_particle_lab2Breit(m_qj);
	int niter=1;
	double tau1_jet; double tau1_jet_new;
	// Once jet algorithm is implemented, use clustered jet as seeds.
	if(m_do_minimization)
	{	int n=0;
		double costheta=cos(m_qj.Theta()); double phi=m_qj.Phi();

		while(true)
		{
		//m_qj = seedqj;
		set_yes();
			if (n==0) tau1_jet=compute_tau1_jet(etype);
			else
			{
				if(n==niter) break;
		 		tau1_jet_new = compute_tau1_jet(etype);
		 		if(tau1_jet_new<tau1_jet) break;	
				tau1_jet = tau1_jet_new;
			}
		//Fit parameters: theta, phi.
		double mag=m_qj.E();
	
        	TMinuit* gMinuit = new TMinuit(2);
		gMinuit->SetFCN(CostFn);
		gMinuit->Command("SET PRINT -1");//  PRINT -1/ NOW 0
        	Double_t arglist[10];
        	Int_t ierflg = 0;
        	arglist[0] = 1; // 1 for chi2, 0.5 for MLL
        	gMinuit->mnexcm("SET ERR",arglist,1,ierflg);

        	// Set starting values and step sizes for parameters
                double pi = TMath::Pi();
                double twopi = 2*pi;
        	Double_t step=pi/100.;
        	static Double_t vstart[2]={costheta,phi};// cos(thteta) = 1. sin(phi)=0.
        	static Double_t stepsize[2] = {2./200,2*step};
        	gMinuit->mnparm(0, "a1", vstart[0],stepsize[0],-1.,1.,ierflg);
        	gMinuit->mnparm(1, "a2", vstart[1],stepsize[1],0.,twopi,ierflg);

        	arglist[0] =10000000; // max calls
        	arglist[1] = 1.; // tolerance

        	gMinuit->mnexcm("SIMplex",arglist,2,ierflg);
        	gMinuit->mnexcm("MIGRAD", arglist,2,ierflg);// symmetric errors, not valid for bound fit
		gMinuit->Command("SET PRINT -1");
        	// Print results
        	Double_t amin,edm,errdef;
        	Int_t nvpar,nparx,icstat;
        	gMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);
        	double fitpar[2]; double dummy;
        	for (int i = 0;i<2; ++i)
		{
        	gMinuit->GetParameter(i,fitpar[i],dummy);
		}
		costheta= fitpar[0]; phi=fitpar[1];
		double sintheta = sqrt(1.-costheta*costheta);
		m_qj.SetPxPyPzE(mag*sintheta*cos(phi), mag*sintheta*sin(phi),mag*costheta,mag);
		++n;
		}
	}
	m_qj_a = boost_particle_Breit2lab(m_qj);	

}

void EventShape::set_yes()
{

        for (std::map<int,TLorentzVector>::iterator iit=m_mommap.begin(); iit!=m_mommap.end(); ++iit){
        	TLorentzVector pp = iit->second;
        	if(verbosity>7) cout<<"p4 "<<pp.Px()<<" "<<pp.Py()<<" "<<pp.Pz()<<" "<<pp.E()<<endl;
        	double qbdotp =m_qb.Dot(pp); double qjdotp = m_qj.Dot(pp);
        	if (qbdotp > qjdotp)
        	{ 
		if(qjdotp<0)cout<<"negative "<< qjdotp<<" p4 "<<pp.P() <<"," <<pp.E()<<endl;
                m_yes.find(iit->first)->second = true;
                if (verbosity>7) cout<<"yes "<<m_yes[iit->first]<<endl;
        	}
        	else
       	 	{
		m_yes.find(iit->first)->second = false;
        	if (verbosity>7) cout<<"yes "<<m_yes[iit->first]<<endl;
        	}
	}
}

void EventShape::set_passed()
{        
	for (std::map<int,TLorentzVector>::iterator iip=m_particles.begin();
                iip!=m_particles.end(); ++iip)
		{
                TLorentzVector pp = iip->second;
                if(verbosity>7) cout<<"p4 "<<pp.Px()<<" "<<pp.Py()<<" "<<pp.Pz()<<" "<<pp.E()<<endl;
                double eta=pp.Eta(); //double phi=pp.Phi();
		bool passed = (fabs(eta<3.5)) ? true:false;
                if (passed)
                {
                m_passed.find(iip->first)->second = true;
                if (verbosity>7) cout<<"yes "<<m_passed[iip->first]<<endl;
                }
                else
                {
                m_passed.find(iip->first)->second = false;
                if (verbosity>7) cout<<"yes "<<m_passed[iip->first]<<endl;
                }
        }
}

double EventShape::compute_tau1(EType etype)
{
        double Q2= (etype==EType::MCTrue) ? m_true_Q2 : m_reco_Q2;
	double tau1=0;
	for (std::map<int,TLorentzVector>::iterator it=m_mommap.begin(); it!=m_mommap.end(); ++it)
	{
		TLorentzVector _p4 = it->second;
		if(verbosity>7) cout<<"p4 "<<_p4.Px()<<" "<<_p4.Py()<<" "<<_p4.Pz()<<" "<<_p4.E()<<endl;
		double term = (m_yes[it->first]) ? m_qj.Dot(_p4) : m_qb.Dot(_p4);
		if(verbosity>7) cout<<"term " <<term<<endl;
		tau1 += (2.*term/Q2);
	}	

	return tau1;
}

double EventShape::compute_tau1_jet(EType etype)
{
	double Q2= (etype==EType::MCTrue) ? m_true_Q2 : m_reco_Q2;
        double tau1=0;
        for (std::map<int,TLorentzVector>::iterator it=m_mommap.begin(); it!=m_mommap.end(); ++it)
        {
                TLorentzVector _p4 = it->second;
                if(verbosity>7) cout<<"p4 "<<_p4.Px()<<" "<<_p4.Py()<<" "<<_p4.Pz()<<" "<<_p4.E()<<endl;
                double term = (m_yes[it->first]) ? m_qj.Dot(_p4) : 0.;
                if(verbosity>7) cout<<"term " <<term<<endl;
                tau1 += (2.*term/Q2);
        }

        return tau1;
}

void EventShape::set_boost(Frame f, EType etype)
{
 	cout<<fixed<<setprecision(5);
	// parton
        double parton_x = (etype == EType::MCTrue)? m_true_parton_x : m_reco_parton_x;	   
	TLorentzVector ppart=parton_x*m_beamhadron;
	// gamma, electron, proton
	TLorentzVector pgnew; TLorentzVector penew; TLorentzVector ppnew;
	switch(f)
	{
		case (Frame::CM) :
		set_boost_lab2CM(m_beamlepton, m_beamhadron);

	if(verbosity>5)
	{
        penew = boost_particle_lab2CM(m_beamlepton);
        ppnew = boost_particle_lab2CM(m_beamhadron);
        cout<<"elec "<<m_beamlepton.E()<<" , "<<m_beamlepton.Px()<<" , "<<m_beamlepton.Py()<<" , "<<m_beamlepton.Pz()<<endl;
        cout<<"hadron "<<m_beamhadron.E()<<" , "<<m_beamhadron.Px()<<" , "
		<<m_beamhadron.Py()<<" , "<<m_beamhadron.Pz()<<endl;
	cout<<"elec boost"<<penew.E()<<" , "<<penew.Px()<<" , "<<penew.Py()<<" , "<<penew.Pz()<<endl;
        cout<<"hadron boost"<<ppnew.E()<<" , "<<ppnew.Px()<<" , "<<ppnew.Py()<<" , "<<ppnew.Pz()<<endl;
	}
		break;	
		
		case (Frame::Breit) :
                set_boost_lab2Breit(m_boson, ppart);		
        if(verbosity>5)
        {
        pgnew = boost_particle_lab2Breit(m_boson);
        ppnew = boost_particle_lab2Breit(ppart);
	TLorentzVector pginv = boost_particle_Breit2lab(pgnew);
        cout<<"gamma "<<m_boson.E()<<" , "<<m_boson.Px()<<" , "<<m_boson.Py()<<" , "<<m_boson.Pz()<<endl;
        cout<<"part "<<ppart.E()<<" , "<<ppart.Px()<<" , "<<ppart.Py()<<" , "<<ppart.Pz()<<endl;
        cout<<"gamma boost "<<pgnew.E()<<" , "<<pgnew.Px()<<" , "<<pgnew.Py()<<" , "<<pgnew.Pz()<<endl;
        cout<<"part boost "<<ppnew.E()<<" , "<<ppnew.Px()<<" , "<<ppnew.Py()<<" , "<<ppnew.Pz()
	<<" pz ratio "<<ppnew.Pz()/pgnew.Pz()<<endl;
	cout<<"gamma back "<<pginv.E()<<" , "<<pginv.Px()<<" , "<<pginv.Py()<<" , "<<pginv.Pz()<<endl;	
        }

		break;
		default:
		break;
	}
}


void EventShape::set_boost_lab2CM(const TLorentzVector& _p1, const TLorentzVector& _p2)
{
	TLorentzVector p1(_p1); TLorentzVector p2(_p2);
	TLorentzVector sum = p1+p2;
	m_b1_cm = -sum.BoostVector();
	p1.Boost(m_b1_cm);
	p2.Boost(m_b1_cm);
	TVector3 p1V= p1.Vect(); TVector3 p2V=p2.Vect();	
	// z axis: p1 going direction
	TVector3 zaxis = p2V.Unit();
	TVector3 yaxis =zaxis.Cross(p1V.Unit());
	yaxis = yaxis.Unit();
	TVector3 xaxis=yaxis.Cross(zaxis);
	xaxis = xaxis.Unit();
	
	m_r1_cm.SetZAxis(zaxis,xaxis);
	m_r1_cm = m_r1_cm.Inverse();
}

void EventShape::set_boost_lab2Breit(const TLorentzVector& _p1, const TLorentzVector& _p2)
{
	TLorentzVector p1(_p1); TLorentzVector p2(_p2);
	if(verbosity>7)
	{
	double Q=-p1.M();
		cout<<"Q "<<Q<<" M "<<p1.M()<<endl;
	// 	p1 (0,0,Q,0) p2 (0,0,-Q/2,sqrt(Q2/4+m2)) 
	TLorentzVector p3(0,0,Q,0); TLorentzVector p4(0,0,-Q/2.,sqrt(pow(Q/2,2)+p2.M2()));
	cout<<"p3-p1 "<<p3.E()-p1.E()<<", "<<p3.Px()-p1.Px()<<", "<<p3.Py()-p1.Py()<<", "<<p3.Pz()-p1.Pz()<<endl;
	cout<<"p2-p4 "<<p2.E()-p4.E()<<", "<<p2.Px()-p4.Px()<<", "<<p2.Py()-p4.Py()<<", "<<p2.Pz()-p4.Pz()<<endl;
	}

	TLorentzVector sum= p1+2.*p2;
	m_b1_br = -sum.BoostVector();

	p1.Boost(m_b1_br);
	p2.Boost(m_b1_br);
	sum.Boost(m_b1_br);

	TVector3 p1V=p1.Vect(); 
	TVector3 p2V=p2.Vect();
	TVector3 zaxis = p1V.Unit();
	TVector3 yaxis = zaxis.Cross(p2V.Unit());
	yaxis = yaxis.Unit();
	TVector3 xaxis=yaxis.Cross(zaxis);
	xaxis=xaxis.Unit();
	
	m_r1_br.SetZAxis(zaxis,xaxis);
	m_r1_br = m_r1_br.Inverse();

/* If boosted to CM first, we need a second boost, in which case boost vector 
 * can become larger than the speed of light. 
	p1V.Transform(m_r1_br);
	p2V.Transform(m_r1_br);	
	p1.SetPx(p1V.Px()); p1.SetPy(p1V.Py()); p1.SetPz(p1V.Pz());
	p2.SetPx(p2V.Px()); p2.SetPy(p2V.Py()); p2.SetPz(p2V.Pz());

	double invm2 = p1.E()+p2.E(); invm2 *= invm2;
	double mp1sq=pow(p1.E(),2)-pow(p1.Vect().Mag(),2);
	double mp2sq=p2.M();  mp2sq *= mp2sq;
	double q=sqrt((pow(mp2sq - invm2 + mp1sq,2)-4*mp1sq*mp2sq)/
			(2*(mp2sq+invm2)-mp1sq));
	double e=sqrt(invm2+q*q/4);
	
	TLorentzVector b2 = TLorentzVector(0,0,q/2.,e);
	m_b2_br =b2.BoostVector();
*/
}

TLorentzVector EventShape::boost_particle_lab2CM(const TLorentzVector& porig)
{	
	TLorentzVector ptrans = porig;
	ptrans.Boost(m_b1_cm);
	TVector3 vec = ptrans.Vect();
	vec.Transform(m_r1_cm);
	ptrans.SetPx(vec.Px()); ptrans.SetPy(vec.Py()); ptrans.SetPz(vec.Pz());
	return ptrans; 
}

TLorentzVector EventShape::boost_particle_lab2Breit(const TLorentzVector& porig)
{
	TLorentzVector ptrans =porig;
	ptrans.Boost(m_b1_br);
	TVector3 vec = ptrans.Vect();
	vec.Transform(m_r1_br);
	ptrans.SetPx(vec.Px()); ptrans.SetPy(vec.Py()); ptrans.SetPz(vec.Pz());
	//ptrans.Boost(m_b2_br);
	return ptrans;
}

TLorentzVector EventShape::boost_particle_Breit2lab(const TLorentzVector& porig)
{
	TLorentzVector ptrans = porig;
	TVector3 vec = ptrans.Vect();
	vec.Transform(m_r1_br.Inverse());
	ptrans.SetPx(vec.Px()); ptrans.SetPy(vec.Py()); ptrans.SetPz(vec.Pz());
	ptrans.Boost(-m_b1_br);
	return ptrans;
}

TLorentzVector EventShape::boost_particle_CM2lab(const TLorentzVector& porig)
{
        TLorentzVector ptrans = porig;
        TVector3 vec = ptrans.Vect();
        vec.Transform(m_r1_cm.Inverse());
        ptrans.SetPx(vec.Px()); ptrans.SetPy(vec.Py()); ptrans.SetPz(vec.Pz());
        ptrans.Boost(-m_b1_cm);
        return ptrans;
}


void EventShape::set_kinematics(EType etype)
{
	switch(etype)
	{
		case (EType::MCTrue):
			set_parton_x(etype);
                        set_Q2(etype);
			set_boson(m_tevent->ExchangeBoson()->PxPyPzE());
			set_beam_hadron(m_tevent->BeamHadron()->PxPyPzE()); 
			set_beam_lepton(m_tevent->BeamLepton()->PxPyPzE());
			set_scattered_lepton(m_tevent->ScatteredLepton()->PxPyPzE());
// GetBeamPartonTheta,GetLeptonPhi,GetHardS/T/U,GetHardQ2,
// GetTrueY,GetTrueQ2,GetTrueW2, GetTrueNu,GetTrueR
		break;
		case (EType::MCReco):
		        set_parton_x(etype);
                        set_Q2(etype);
			// photon momentum is not measured, and not used.
                        set_boson(m_tevent->ExchangeBoson()->PxPyPzE());
			// beam/lepton momenta are the same as truth
                        set_beam_hadron(m_tevent->BeamHadron()->PxPyPzE());
                        set_beam_lepton(m_tevent->BeamLepton()->PxPyPzE());
                        set_scattered_lepton(m_sevent->ScatteredLepton()->PxPyPzE());
		
		break;
		default:
		cout<<"Set Valid Event Type! "<<endl;
		break;
	}
	//TLorentzVector q= m_beamlepton - m_scatlepton;
	//m_Q2 = -q.Dot(q);
	if(verbosity>7) cout<<"Q2 "<<m_true_Q2 <<endl;
}

void EventShape::set_parton_x(EType etype)
{       if (etype == MCTrue) m_true_parton_x = m_tevent->GetTrueX();
        else if (etype == MCReco)
	{
	 switch(m_xQ2_method)
	 {
	  case 1: // Null momentum (electron)
	  m_reco_parton_x = m_sevent->GetX(); break;
	  case 2: // Jacquet Blondel (hadron)
	  m_reco_parton_x = m_sevent->GetXJacquetBlondel(); break;
	  case 3: // Double Angle 
	  m_reco_parton_x = m_sevent->GetXDoubleAngle(); break;
	  default: // Electron
	  m_reco_parton_x = m_sevent->GetX(); break;
	 }
	}
}

double EventShape::get_parton_x(EType etype)
{       
	double x;
	if(etype == MCTrue) x = m_tevent->GetTrueX();
	else if (etype == MCReco) 
	{
         switch(m_xQ2_method)
	 {
          case 1: // Null Momentum (electron)
          x = m_sevent->GetX(); break;
          case 2: // Jacquet Blondel (hadron)
          x = m_sevent->GetXJacquetBlondel(); break;
          case 3: // Double Angle
          x = m_sevent->GetXDoubleAngle(); break;
          default: 
          x = m_sevent->GetX(); break;
         }
	}
        return x;
}

void EventShape::set_Q2(EType etype)
{       if (etype==MCTrue) m_true_Q2 = m_tevent->GetTrueQ2();
        else if (etype==MCReco)
	{
         switch(m_xQ2_method)
         {
          case 1: // Null Momentum (electron)
          m_reco_Q2 = m_sevent->GetQ2(); break;
          case 2: // Jacquet Blondel (hadron)
          m_reco_Q2 = m_sevent->GetQ2JacquetBlondel(); break;
          case 3: // Double Angle
          m_reco_Q2 = m_sevent->GetQ2DoubleAngle(); break;
          default:
          m_reco_Q2 = m_sevent->GetQ2(); break;
         }
	}
}

double EventShape::get_Q2(EType etype)
{       double Q2;
	if (etype == MCTrue) Q2 = m_tevent->GetTrueQ2();
	else if (etype == MCReco)
	{
         switch(m_xQ2_method)
         {
          case 1: // Null Momentum (electron)
          Q2 = m_sevent->GetQ2(); break;
          case 2: // Jacquet Blondel (hadron)
          Q2 = m_sevent->GetQ2JacquetBlondel(); break;
          case 3: // Double Angle
          Q2 = m_sevent->GetQ2DoubleAngle(); break;
          default:
          Q2 = m_sevent->GetQ2(); break;
         }
	}
        return Q2;
}

void EventShape::set_scattered_lepton(Frame f, EType etype)
{       TLorentzVector p4_e= (etype == MCTrue) ? m_tevent->ScatteredLepton()->Get4Vector()
        : m_sevent->ScatteredLepton()->Get4Vector();
        switch(f){
                case (Frame::Breit):
                if(etype == MCTrue) m_true_e_br = boost_particle_lab2Breit(p4_e);
                else m_reco_e_br = boost_particle_lab2Breit(p4_e);
                break;
                case (Frame::CM):
                if(etype == MCTrue) m_true_e_cm = boost_particle_lab2CM(p4_e);
                else m_reco_e_cm = boost_particle_lab2CM(p4_e);
                break;
                default:
                m_true_e_lb = p4_e; 
		m_reco_e_lb = p4_e;
                break;
        }
}

TLorentzVector EventShape::get_scattered_lepton(Frame f, EType etype)
{
        TLorentzVector p4_e;
        switch(f){
                case (Frame::Breit):
                p4_e = (etype==MCTrue) ? m_true_e_br :  m_reco_e_br;
                break;
                case (Frame::CM):
                p4_e = (etype==MCTrue) ? m_true_e_cm : m_reco_e_cm;
                break;
                default:
                p4_e = (etype==MCTrue) ? m_true_e_lb : m_reco_e_lb;
                break;
        }
        return p4_e;
}

double EventShape::get_e_y(Frame f, EType etype)
{
        TLorentzVector p4_e = get_scattered_lepton(f,etype);
        return p4_e.Rapidity();
}
double EventShape::get_e_pt(Frame f, EType etype)
{
        TLorentzVector p4_e = get_scattered_lepton(f,etype);
        return p4_e.Pt();
}
double EventShape::get_e_eta(Frame f, EType etype)
{
        TLorentzVector p4_e = get_scattered_lepton(f,etype);
        return p4_e.Eta();
}

void EventShape::set_jet(EType etype)
{

        TLorentzVector pjet(0.,0.,0.,0.);
        for(std::map<int,TLorentzVector>::iterator itt=m_mommap.begin(); itt!=m_mommap.end();++itt)
        {
                TLorentzVector pp = itt->second;
		if (m_yes[itt->first]) pjet += pp;
        }
	
	// Always set in Breit frame
        if(etype == MCTrue) m_true_jet_br = pjet;
        else m_reco_jet_br = pjet;

	TLorentzVector pjet_lb = boost_particle_Breit2lab(pjet);
        if(etype == MCTrue) m_true_jet_lb = pjet_lb;
	else m_reco_jet_lb = pjet_lb;

        TLorentzVector pjet_cm = boost_particle_lab2CM(pjet_lb);
        if(etype == MCTrue) m_true_jet_cm = pjet_cm;
        else m_reco_jet_cm = pjet_cm;
	if(verbosity>3)cout<<"Jet_pt in Breit "<<pjet.Pt()<<" lab "<<pjet_lb.Pt()<<" CM "<<pjet_cm.Pt()<<endl;
        if(verbosity>3)cout<<"Jet_p in Breit "<<pjet.P()<<" lab "<<pjet_lb.P()<<" CM "<<pjet_cm.P()<<endl;
}

TLorentzVector EventShape::get_jet(Frame f, EType etype)
{
        TLorentzVector p4_jet;
        switch(f){
                case (Frame::Breit):
                p4_jet = (etype==MCTrue) ? m_true_jet_br :  m_reco_jet_br;
                break;
                case (Frame::CM):
                p4_jet = (etype==MCTrue) ? m_true_jet_cm : m_reco_jet_cm;
                break;
                default:
                p4_jet = (etype==MCTrue) ? m_true_jet_lb : m_reco_jet_lb;
                break;
        }
        return p4_jet;
}

double EventShape::get_jet_y(Frame f, EType etype)
{       TLorentzVector p4_jet = get_jet(f,etype);
        return p4_jet.Rapidity();
}

double EventShape::get_jet_pt(Frame f, EType etype)
{       TLorentzVector p4_jet = get_jet(f,etype);
        return p4_jet.Pt();
}

double EventShape::get_jet_eta(Frame f, EType etype)
{       TLorentzVector p4_jet = get_jet(f,etype);
        return p4_jet.Eta();
}

double EventShape::get_jet_p(Frame f, EType etype)
{	TLorentzVector p4_jet = get_jet(f,etype);
	return p4_jet.P();
}
double EventShape::get_jet_e(Frame f, EType etype)
{	TLorentzVector p4_jet = get_jet(f,etype);
	return p4_jet.E();
}
double EventShape::get_jet_mass(Frame f, EType etype)
{	TLorentzVector p4_jet = get_jet(f,etype);
	return p4_jet.M();
}

double EventShape::get_suppression_factor(TauType ttype, Frame f)
{
        double sf=0.;
        switch(ttype){
		case (TauType::A):
		 if (f == Frame::Lab) sf = m_sf_a_lb;
		 else if (f == Frame::Breit) sf = m_sf_a_br;
		 else if (f == Frame::CM) sf = m_sf_a_cm;
		break;
                case (TauType::B):
		 if (f == Frame::Lab) sf =  m_sf_b_lb;
                 else if (f == Frame::Breit) sf = m_sf_b_br;
		 else if (f == Frame::CM) sf = m_sf_b_cm;
                break;
                case (TauType::C):
                 if (f == Frame::Lab) sf = m_sf_c_lb;
		 else if (f == Frame::Breit) sf = m_sf_c_br;
		 else if (f == Frame::CM) sf = m_sf_c_cm;
                break;
                default:
                break;
        }
        return sf;	
}

double EventShape::get_tau1_q_fraction(TauType ttype)
{
	double qf=0.;
	switch(ttype)
	{
		case (TauType::A):
		qf = m_qf_a; break;
		case (TauType::B):
		qf = m_qf_b; break;
		case (TauType::C):
		qf = m_qf_c; break;
		default: break;
	}
	return qf;
}

TLorentzVector EventShape::get_qj(TauType ttype, Frame f)
{
	TLorentzVector p4(0.,0.,0.,0.);
	switch(ttype)
	{
		case (TauType::A):
			switch(f)
			{
			case (Frame::Lab):
			p4 = m_qj_a; break;
			case (Frame::Breit):
			p4 = boost_particle_lab2Breit(m_qj_a); break;
			case (Frame::CM): 
			p4 = boost_particle_lab2CM(m_qj_a); break;
			default: break;
			}
		break;

		case (TauType::B):
	        	switch(f)
	        	{
	                case (Frame::Lab):
	                p4 = m_qj_b; break;
	                case (Frame::Breit):
	                p4 = boost_particle_lab2Breit(m_qj_b); break;
	                case (Frame::CM):
	                p4 = boost_particle_lab2CM(m_qj_b); break;
	                default: break;
	        	}
		break;

		case (TauType::C):
		        switch(f)
	        	{
	                case (Frame::Lab):
	                p4 = m_qj_c; break;
	                case (Frame::Breit):
	                p4 = boost_particle_lab2Breit(m_qj_c); break;
	                case (Frame::CM):
	                p4 = boost_particle_lab2CM(m_qj_c); break;
	                default: break;
	        	}
		break;

		default: break;
	}
	return p4;
}

TLorentzVector EventShape::get_qb(TauType ttype, Frame f)
{
        TLorentzVector p4(0.,0.,0.,0.);
        switch(ttype)
        {
                case (TauType::A):
                        switch(f)
                        {
                        case (Frame::Lab):
                        p4 = m_qb_a; break;
                        case (Frame::Breit):
                        p4 = boost_particle_lab2Breit(m_qb_a); break;
                        case (Frame::CM):
                        p4 = boost_particle_lab2CM(m_qb_a); break;
                        default: break;
                        }
                break;

                case (TauType::B):
                        switch(f)
                        {
                        case (Frame::Lab):
                        p4 = m_qb_b; break;
                        case (Frame::Breit):
                        p4 = boost_particle_lab2Breit(m_qb_b); break;
                        case (Frame::CM):
                        p4 = boost_particle_lab2CM(m_qb_b); break;
                        default: break;
                        }
                break;
                case (TauType::C):
                        switch(f)
                        {
                        case (Frame::Lab):
                        p4 = m_qb_c; break;
                        case (Frame::Breit):
                        p4 = boost_particle_lab2Breit(m_qb_c); break;
                        case (Frame::CM):
                        p4 = boost_particle_lab2CM(m_qb_c); break;
                        default: break;
                        }
                break;

                default: break;
        }
        return p4;
}

void EventShape::set_particles(EType etype)
{
	m_passed.clear();
	m_particles.clear();
        unsigned int nentries= (etype == MCTrue) ? 
			m_tevent->GetNTracks() : m_sevent->GetNTracks();
        for(unsigned int it =0; it< nentries;++it)
        {
                if(it<3) continue;
		TLorentzVector partp4;
		if(etype==EType::MCTrue)
		{
                const Particle *part = m_tevent->GetTrack(it);
                if(part->GetStatus() !=1) continue;
                if(part->GetE() == m_scatlepton.E()) continue;
                //if(!pass_cuts(it, true)) continue;
                partp4 = part->PxPyPzE();
		}
		else if(etype == EType::MCReco)
		{
                const Smear::ParticleMCS *part = m_sevent->GetTrack(it);
		if (!part) continue;
		if(part->GetStatus() !=1) continue;
		double pE = part->GetE(); double pP = part->GetP();
                if((fabs(pE) < threshold) && (fabs(pP) < threshold)) continue;
                if(part->GetE() == m_scatlepton.E()) continue;
                if(!pass_cuts(it, true)) continue;
                partp4 = part->PxPyPzE();

		double mass =  m_tevent->GetTrack(it)->GetM();
                if(m_smear_mode==0)
                {
                partp4.SetPxPyPzE(m_tevent->GetTrack(it)->GetPx(),
                              m_tevent->GetTrack(it)->GetPy(),
                              m_tevent->GetTrack(it)->GetPz(),
                              m_tevent->GetTrack(it)->GetE());
                }
                else
                {
                        int pid = std::abs(m_tevent->GetTrack(it)->GetPdgCode());
		 	if(fabs(pE) < threshold)
		 	{// Tracking info only
			 double energy=0.;
                         if(m_smear_mode==3 || m_smear_mode==4)
                         { // 3 sig separation
                //partp4.SetPxPyPzE(m_tevent->GetTrack(it)->GetPx(),
                //              m_tevent->GetTrack(it)->GetPy(),
                //              m_tevent->GetTrack(it)->GetPz(),
                //              m_tevent->GetTrack(it)->GetE());
                          mass=0.139;
                          double p3 = m_tevent->GetTrack(it)->GetP();
                          double eta = m_tevent->GetTrack(it)->GetEta();
                          bool is_pkpi = ((pid==211) || (pid==321) || (pid ==2212));
                          bool is_3sig = is_pid_3sig(p3, eta);
			  if (is_3sig && is_pkpi)
			  {
			   mass= m_tevent->GetTrack(it)->GetM();
			   energy = sqrt(pow(p3,2) + pow(mass,2));
			  }
                         }
                 	 if (m_smear_mode ==2){ mass =0.139; energy = sqrt(pow(pP,2) + pow(mass,2));}
                         if (m_smear_mode ==1){ //energy = m_tevent->GetTrack(it)->GetE();}
                                               mass = m_tevent->GetTrack(it)->GetM();
                                               energy = sqrt(pow(pP,2)+pow(mass,2));}
		 	 partp4.SetE(energy);

                         if (m_tracking_mode>0)
                         {
                         double pt = m_tevent->GetTrack(it)->GetPt();
                         double eff = ftreff->Eval(pt);
                         if (eff < rnd->Uniform(0.,1.)) continue; // no track, no calo, skip
                         }

		 	}// end track only 
		 	else if(fabs(pP) < threshold)
		 	{// Calo info only (gamma or neutral h)

                        partp4 = get_smeared_calo_only(it, partp4);

		 	}// end calo only
			else if(fabs(pP) >= threshold && fabs(pE) >= threshold)
                        {// Both track and Calo
                	if(pid!=11)
                	partp4 = get_smeared_chargedH(it,partp4 );
                	else 
                	partp4 = get_smeared_EM(it, partp4);
			}
		}
		}// MCReco
		
		m_particles.insert(make_pair((int) it, partp4));
		m_passed.insert(make_pair((int)it, true));
	}
}

void EventShape::set_pj(TauType ttype)
{
        for (std::map<int,TLorentzVector>::iterator ipart=m_particles.begin(); 
		ipart!=m_particles.end(); ++ipart)
        {
		int ip=ipart->first;
                TLorentzVector pp = ipart->second;
                TLorentzVector pp_br = boost_particle_lab2Breit(pp);
                TLorentzVector pp_cm = boost_particle_lab2CM(pp);
		switch(ttype)
		{
		 case (TauType::A):
		  if(m_yes[ip]){
			m_pjq_a_lab.push_back(pp);
			m_pjq_a_breit.push_back(pp_br);
			m_pjq_a_cm.push_back(pp_cm);
			m_pjq_a_id.push_back(ip);
		  }
		  else
		  {
			m_pjb_a_lab.push_back(pp);
			m_pjb_a_breit.push_back(pp_br);
			m_pjb_a_cm.push_back(pp_cm);
			m_pjb_a_id.push_back(ip);
		  }
		 break;
		 case (TauType::B):
		  if(m_yes[ip])
		  {
                        m_pjq_b_lab.push_back(pp);
                        m_pjq_b_breit.push_back(pp_br);
                        m_pjq_b_cm.push_back(pp_cm);
		  }
		  else
		  {
                        m_pjb_b_lab.push_back(pp);
                        m_pjb_b_breit.push_back(pp_br);
                        m_pjb_b_cm.push_back(pp_cm);
		  }
		 break;
		 case (TauType::C):
		  if(m_yes[ip])
		  {
                        m_pjq_c_lab.push_back(pp);
                        m_pjq_c_breit.push_back(pp_br);
                        m_pjq_c_cm.push_back(pp_cm);
		  }
		  else
		  {
                        m_pjb_c_lab.push_back(pp);
                        m_pjb_c_breit.push_back(pp_br);
                        m_pjb_c_cm.push_back(pp_cm);
		  }
		 break;
		 default:
		 break;
		}
	}
}

void EventShape::clear_pj()
{
     m_pjq_a_id.clear();
     m_pjq_a_lab.clear();
     m_pjq_a_breit.clear();
     m_pjq_a_cm.clear();
     m_pjq_b_lab.clear();
     m_pjq_b_breit.clear();
     m_pjq_b_cm.clear();
     m_pjq_c_lab.clear();
     m_pjq_c_breit.clear();
     m_pjq_c_cm.clear();
     m_pjb_a_lab.clear();
     m_pjb_a_breit.clear();
     m_pjb_a_cm.clear();
     m_pjb_b_lab.clear();
     m_pjb_b_breit.clear();
     m_pjb_b_cm.clear();
     m_pjb_c_lab.clear();
     m_pjb_c_breit.clear();
     m_pjb_c_cm.clear();

}

std::vector<int> EventShape::get_pjq_id(TauType ttype)
{
	std::vector<int> pjqidvec;
	pjqidvec.clear();
        switch(ttype)
        {
                case (TauType::A):
			pjqidvec = m_pjq_a_id;
		break;
		default:
			pjqidvec = m_pjq_a_id;
		break;
	}
	return pjqidvec;
}
std::vector<int> EventShape::get_pjb_id(TauType ttype)
{
	std::vector<int> pjbidvec;
	pjbidvec.clear();
        switch(ttype)
        {
                case (TauType::A):
                        pjbidvec = m_pjb_a_id;
                break;
                default:
                        pjbidvec = m_pjb_a_id;
                break;
        }
        return pjbidvec;
}


std::vector<float> EventShape::get_pjq(TauType ttype, Frame f, Var var, EType etype)
{
        double Q2= (etype==EType::MCTrue) ? m_true_Q2 : m_reco_Q2;
	std::vector<TLorentzVector> pjqvec;
	pjqvec.clear();
	switch(ttype)
	{
		case (TauType::A):
			if(f==Frame::Lab) pjqvec = m_pjq_a_lab;
			else if(f==Frame::Breit) pjqvec = m_pjq_a_breit;
			else if(f==Frame::CM) pjqvec = m_pjq_a_cm;
		break;
		case (TauType::B):
			if(f==Frame::Lab) pjqvec = m_pjq_b_lab;
			else if(f==Frame::Breit) pjqvec = m_pjq_b_breit;
			else if(f==Frame::CM) pjqvec = m_pjq_b_cm;
		break;
		case (TauType::C):
			if(f==Lab) pjqvec = m_pjq_c_lab;
			else if(f==Breit) pjqvec = m_pjq_c_breit;
			else if(f==CM) pjqvec = m_pjq_c_cm;
		break;
		default: 
		break;
	}
	std::vector<float> varvec;
	varvec.clear();
	for (std::vector<TLorentzVector>::iterator vv = pjqvec.begin(); vv != pjqvec.end(); ++vv)
	{
	TLorentzVector vect = *vv;
        float val = -99.;
        switch(var)
        {
                case (Var::Dot):
                val = 2*m_qj.Dot(vect)/Q2; break;
                case (Var::Eta):
                val = vect.Eta(); break;
                case (Var::Pt):
                val = vect.Pt(); break;
                case (Var::P):
                val = vect.P(); break;
                default: break;
        }
        varvec.push_back(val);
        }
        return varvec;

}

std::vector<float> EventShape::get_pjb(TauType ttype, Frame f, Var var, EType etype)
{
        double Q2= (etype==EType::MCTrue) ? m_true_Q2 : m_reco_Q2;
	std::vector<TLorentzVector> pjbvec;
	pjbvec.clear();
        switch(ttype)
        {
                case (TauType::A):
                        if(f==Frame::Lab) pjbvec= m_pjb_a_lab;
                        else if(f==Frame::Breit) pjbvec= m_pjb_a_breit;
                        else if(f==Frame::CM) pjbvec= m_pjb_a_cm;
                break;
                case (TauType::B):
                        if(f==Frame::Lab) pjbvec= m_pjb_b_lab;
                        else if(f==Frame::Breit) pjbvec= m_pjb_b_breit;
                        else if(f==Frame::CM) pjbvec= m_pjb_b_cm;
                break;
                case (TauType::C):
                        if(f==Lab) pjbvec= m_pjb_c_lab;
                        else if(f==Breit) pjbvec= m_pjb_c_breit;
                        else if(f==CM) pjbvec= m_pjb_c_cm;
                break;
                default:
                break;
        }
	std::vector<float> varvec;
	varvec.clear();
 	for (std::vector<TLorentzVector>::iterator iv = pjbvec.begin() ; iv != pjbvec.end(); ++iv)
	{
	TLorentzVector vect = *iv;
	float val = -99.;
	switch(var)
	{
		case (Var::Dot):
		val = 2*m_qb.Dot(vect)/Q2; break;	
		case (Var::Eta):
		val = vect.Eta(); break;
		case (Var::Pt):
		val = vect.Pt(); break;
		case (Var::P):
		val = vect.P(); break;
		default: break;
	}
	varvec.push_back(val);
	}
	return varvec;
}

TLorentzVector EventShape::get_smeared_chargedH(unsigned int it, const TLorentzVector& p4)
{
		//int pid = std::abs(m_tevent->GetTrack(it)->GetPdgCode());
		TLorentzVector p(p4);
                int pid = std::abs(m_tevent->GetTrack(it)->GetPdgCode());
		if (m_smear_mode ==1)
		{//// time-like 4 vector ()
		 double mass = m_tevent->GetTrack(it)->GetM();
                 double energy = sqrt(pow(p4.P(),2)+pow(mass,2));
		 p.SetE(energy);
		 // true enery
                 //p.SetE(m_tevent->GetTrack(it)->GetE());
		}
		else if (m_smear_mode ==2)
		{
                double mass =0.139;
                double energy = sqrt(pow(p.P(),2) + pow(mass,2));
		p.SetE(energy);
		}
		else if (m_smear_mode ==3 || m_smear_mode ==4)
		{
                //p.SetPxPyPzE(m_tevent->GetTrack(it)->GetPx(),
                //              m_tevent->GetTrack(it)->GetPy(),
                //              m_tevent->GetTrack(it)->GetPz(),
                //              m_tevent->GetTrack(it)->GetE());

                double mass =0.139;
		double energy= 0.;
                double p3 = m_tevent->GetTrack(it)->GetP();
                double eta = m_tevent->GetTrack(it)->GetEta();
                bool is_pkpi = ((pid==211) || (pid==321) || (pid ==2212));
                bool is_3sig = is_pid_3sig(p3, eta);
		if(is_3sig && is_pkpi)
		{
			 mass= m_tevent->GetTrack(it)->GetM();
			 energy = sqrt(pow(p.P(),2) + pow(mass,2));
		}
                p.SetE(energy);
		}

                if (m_tracking_mode>0)
                {
                double pt = m_tevent->GetTrack(it)->GetPt();
                double eff = ftreff->Eval(pt);
                if (eff < rnd->Uniform(0.,1.))  // no track, p3 from calo
		 {
		  p = get_smeared_calo_only(it, p);
		 }
                }

		if(fabs(p.P())>fabs(p.E()))cout<<"chargedH pid "<<pid<<" space-like particle"<<endl;
		return p;
}

TLorentzVector EventShape::get_smeared_EM(unsigned int it, const TLorentzVector& p4)
{
		TLorentzVector p(p4);
                int pid = std::abs(m_tevent->GetTrack(it)->GetPdgCode());
		if (m_smear_mode==1)
		{
                 double mass = m_tevent->GetTrack(it)->GetM();
                 double energy = sqrt(pow(p4.P(),2)+pow(mass,2));
                 p.SetE(energy);

		//p.SetE(m_tevent->GetTrack(it)->GetE());
		}
		else if(m_smear_mode==2)
		{
                double mass =0.;
                double energy = sqrt(pow(p.P(),2) + pow(mass,2));
                p.SetE(energy);			
		}
		else if(m_smear_mode==3 || m_smear_mode==4)
		{

                p.SetPxPyPzE(m_tevent->GetTrack(it)->GetPx(),
                              m_tevent->GetTrack(it)->GetPy(),
                              m_tevent->GetTrack(it)->GetPz(),
                              m_tevent->GetTrack(it)->GetE());
                double mass =0.;
                double energy = sqrt(pow(p.P(),2) + pow(mass,2));
                p.SetE(energy);
		}

                if (m_tracking_mode>0)
                {
                double pt = m_tevent->GetTrack(it)->GetPt();
                double eff = ftreff->Eval(pt);
                if (eff < rnd->Uniform(0.,1.))  // no track, p3 from calo
                 {
                  p = get_smeared_calo_only(it, p);
                 }
                }
                if(fabs(p.P())>fabs(p.E()))cout<<"EM pid "<<pid<<" space-like particle"<<endl;
		return p;
}

TLorentzVector EventShape::get_smeared_calo_only(unsigned int it, const TLorentzVector& p4)
{
		TLorentzVector p(p4);
                 int pid =std::abs(m_tevent->GetTrack(it)->GetPdgCode());
		double mass = m_tevent->GetTrack(it)->GetM();
		double energy = m_tevent->GetTrack(it)->GetE();
                double phi = m_tevent->GetTrack(it)->GetPhi(); 
		double theta = m_tevent->GetTrack(it)->GetTheta();
                if (m_smear_mode ==1 || m_smear_mode==2 || m_smear_mode==3)
		{
		}
                if (m_smear_mode==4)
                {
		energy = m_sevent->GetTrack(it)->GetE();
		theta = m_sevent->GetTrack(it)->GetTheta();
		phi = m_sevent->GetTrack(it)->GetPhi();
		
		//if (pid==2112 || pid==130) mass = 0.4937;
                }
                double p3 = sqrt(pow(energy,2) - pow(mass,2));
                p.SetPx(p3*sin(theta)*cos(phi));
                p.SetPy(p3*sin(theta)*sin(phi));
                p.SetPz(p3*cos(theta));
		p.SetE(energy);
                if(verbosity>7 && (fabs(p.P())-fabs(p.E())>0.00001))
		cout<<"calo only pid "<<pid<<" p "<<p.P()<<" E "<<p.E()<<endl;
		if(energy<mass) p.SetPxPyPzE(0.,0.,0.,0.);
		return p;
}

void EventShape::set_tracking_mode(int tracking_mode)
{
   m_tracking_mode = tracking_mode;
   switch(tracking_mode)
   {

    case 4: // .93 at saturation, 0.5 at 0.2 GeV/c
    ftreff->SetParameters(0.93, 0.0395111, 0.197555); break;

    case 3: // .95 at saturation, 0.7 at 0.2 GeV/c
    ftreff->SetParameters(0.95,0.0331696, 0.165848); break;
//    case 3: // .95 at saturation, 0.94 at 0.2 GeV/c
//    ftreff->SetParameters(0.95, 0.0209571, 0.104786); break;

    case 2: // .97 at saturation, 0.3 at 0.2 GeV/c
    ftreff->SetParameters(0.97,0.0476587, 0.238294); break;

    case 1: // .97 at saturation, 0.5 at 0.5GeV.c 
    ftreff->SetParameters(0.97, 0.0395111, 0.197555); break;
    default:
    ftreff->SetParameters(1.,1.,-10000.); break;
   }
}
bool EventShape::pass_cuts(unsigned int it, bool is_true_part)
{
	bool pass=true;
	//list variables, kinematic cuts first
	//GetPx,GetPy,GetPz,GetPxPyPzE,GetPt,GetEta,GetE,Id
	float pt = (is_true_part) ? m_tevent->GetTrack(it)->GetPt()
				  : m_sevent->GetTrack(it)->GetPt();
	if (pt<m_pt_cut) return false;
	return pass;
}

bool EventShape::is_pid_3sig(double p, double eta)
{
	int etabin=-9; // backward to forward
	if(eta > -3.5 && eta < -1 ) etabin = 0; 
	else if (eta >= -1 && eta < 1) etabin = 1;
	else if (eta >= 1  && eta < 2) etabin = 2;
	else if (eta >= 2  && eta < 3) etabin = 3;
	else if (eta >= 3  && eta < 3.5) etabin = 4;
	switch (etabin)
	{
		case 0:
		if(p<7. && p>0.) return true; 
		case 1:
		if(p<5. && p>0.) return true;
		case 2:
		if(p<8. && p>0.) return true; 
		case 3:
		if(p<20. && p>0.) return true; 
		case 4:
		if(p<45. && p>0.) return true; 
		default:
		return false;
	}
}

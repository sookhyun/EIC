#include "MyJetAnalysis.h"

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/Fun4AllServer.h>
#include <fun4all/PHTFileServer.h>
#include <g4main/PHG4Particle.h>
#include <g4main/PHG4TruthInfoContainer.h>

#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>

#include <g4eval/JetEvalStack.h>

//#include <trackbase_historic/SvtxTrackMap.h>

#include <g4jets/JetMap.h>
#include <g4jets/Jet.h>

#include <TDatabasePDG.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TString.h>
#include <TTree.h>
#include <TVector3.h>

#include <HepMC/GenEvent.h>
#include <HepMC/GenVertex.h>
#include <phhepmc/PHHepMCGenEvent.h>
#include <phhepmc/PHHepMCGenEventMap.h>
#include <phhepmc/PHGenIntegral.h>  // for PHGenIntegral
#include <phhepmc/PHGenIntegralv1.h>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>

using namespace std;

MyJetAnalysis::MyJetAnalysis(const std::string& recojetname, const std::string& truthjetname, const std::string& outputfilename,float match_dist)
  : SubsysReco("MyJetAnalysis_" + recojetname + "_" + truthjetname)
  , m_recoJetName(recojetname)
  , m_truthJetName(truthjetname)
  , m_outputFileName(outputfilename)
  , m_etaRange(-1, 1)
  , m_ptRange(5, 100)
  , m_trackJetMatchingRadius(.7)
  , m_T(nullptr)
  , m_event(-1)
  , m_id(-1)
  , m_nComponent(-1)
  , m_eta(numeric_limits<float>::signaling_NaN())
  , m_phi(numeric_limits<float>::signaling_NaN())
  , m_e(numeric_limits<float>::signaling_NaN())
  , m_pt(numeric_limits<float>::signaling_NaN())
  , m_truthID(-1)
  , m_truthNComponent(-1)
  , m_truthEta(numeric_limits<float>::signaling_NaN())
  , m_truthPhi(numeric_limits<float>::signaling_NaN())
  , m_truthE(numeric_limits<float>::signaling_NaN())
  , m_truthPt(numeric_limits<float>::signaling_NaN())
  , m_nMatchedTrack(-1)
  , m_matchDist(match_dist)
  , verbosity(0)
{
  m_trackdR.fill(numeric_limits<float>::signaling_NaN());
  m_trackpT.fill(numeric_limits<float>::signaling_NaN());
}

MyJetAnalysis::~MyJetAnalysis()
{
}

int MyJetAnalysis::Init(PHCompositeNode* topNode)
{
  if (Verbosity() >= MyJetAnalysis::VERBOSITY_SOME)
    cout << "MyJetAnalysis::Init - Outoput to " << m_outputFileName << endl;

  PHTFileServer::get().open(m_outputFileName, "RECREATE");

  m_T = new TTree("T", "MyJetAnalysis Tree");
  // parton
  m_T->Branch("process_id",&m_process_id,"process_id/I");
  m_T->Branch("parton_id1",&m_parton_id1,"parton_id1/I");
  m_T->Branch("parton_id2",&m_parton_id2,"parton_id2/I");
  m_T->Branch("x1", &m_x1,"x1/F");
  m_T->Branch("x2", &m_x2,"x2/F");
  m_T->Branch("parton_flavor", &m_flavor,"parton_flavor/I");
  m_T->Branch("parton_eta", &m_parton_eta,"parton_eta/F");   
  m_T->Branch("parton_phi", &m_parton_phi,"parton_phi/F"); 
  m_T->Branch("parton_zt", &m_parton_zt,"parton_zt/F"); 
  m_T->Branch("photon_eta", &m_photon_eta,"photon_eta/F");
  m_T->Branch("photon_phi", &m_photon_phi,"photon_phi/F");
  m_T->Branch("photon_pt", &m_photon_pt,"photon_pt/F");

  // reconstructed jets
  m_T->Branch("m_event", &m_event, "event/I");
  m_T->Branch("id", &m_id, "id/I");
  m_T->Branch("nComponent", &m_nComponent, "nComponent/I");
  m_T->Branch("eta", &m_eta, "eta/F");
  m_T->Branch("phi", &m_phi, "phi/F");
  m_T->Branch("e", &m_e, "e/F");
  m_T->Branch("pt", &m_pt, "pt/F");
  // truth jets
  m_T->Branch("truthID", &m_truthID, "truthID/I");
  m_T->Branch("truthNComponent", &m_truthNComponent, "truthNComponent/I");
  m_T->Branch("truthEta", &m_truthEta, "truthEta/F");
  m_T->Branch("truthPhi", &m_truthPhi, "truthPhi/F");
  m_T->Branch("truthE", &m_truthE, "truthE/F");
  m_T->Branch("truthPt", &m_truthPt, "truthPt/F");
  // particle level information
  m_T->Branch("nTruthConstituents", &m_nTruthConstituents, "nTruthConstituents/I");
  m_T->Branch("truth_part_px", m_truth_part_px.data(), "truth_part_px[nTruthConstituents]/F");
  m_T->Branch("truth_part_py", m_truth_part_py.data(), "truth_part_py[nTruthConstituents]/F");
  m_T->Branch("truth_part_pz", m_truth_part_pz.data(), "truth_part_pz[nTruthConstituents]/F");
  m_T->Branch("truth_part_e", m_truth_part_e.data(), "truth_part_e[nTruthConstituents]/F");
  m_T->Branch("truth_part_phi", m_truth_part_phi.data(), "truth_part_phi[nTruthConstituents]/F");
  m_T->Branch("truth_part_pid", m_truth_part_pid.data(), "truth_part_pid[nTruthConstituents]/I");
  m_T->Branch("truth_part_qe", m_truth_part_qe.data(), "truth_part_qe[nTruthConstituents]/I");

  //      //! number of matched tracks
  //      int m_nMatchedTrack;
  m_T->Branch("nMatchedTrack", &m_nMatchedTrack, "nMatchedTrack/I");
  //      std::array<float, kMaxMatchedTrack> m_trackdR;
  m_T->Branch("trackdR", m_trackdR.data(), "trackdR[nMatchedTrack]/F");
  //      std::array<float, kMaxMatchedTrack> m_trackpT;
  m_T->Branch("trackpT", m_trackpT.data(), "trackpT[nMatchedTrack]/F");

  return Fun4AllReturnCodes::EVENT_OK;
}

int MyJetAnalysis::End(PHCompositeNode* topNode)
{
  m_intlumi = genintegral->get_Integrated_Lumi();

  cout<<"MyJetAnalysis::End \n \t integrated lumi (m_Pythia8->info.nAccepted()/(m_Pythia8->info.sigmaGen() * 1e9) )"<< m_intlumi 
			<< "\n \t cross section " <<genintegral->get_CrossSection_Processed_Event()  
			<< "\n \t n processed events "<< genintegral->get_N_Processed_Event() 
			<< "\n \t n accepted events (m_Pythia8->info.nAccepted() ) " <<genintegral->get_N_Generator_Accepted_Event()  <<endl;
  cout << "MyJetAnalysis::End - Output to " << m_outputFileName << endl;
  PHTFileServer::get().cd(m_outputFileName);
  m_T->Write();

  return Fun4AllReturnCodes::EVENT_OK;
}

int MyJetAnalysis::InitRun(PHCompositeNode* topNode)
{
//  m_jetEvalStack = shared_ptr<JetEvalStack>(new JetEvalStack(topNode, m_recoJetName, m_truthJetName));
//  m_jetEvalStack->get_stvx_eval_stack()->set_use_initial_vertex(initial_vertex);
  genintegral = findNode::getClass<PHGenIntegral>(topNode, "PHGenIntegral");
  m_intlumi = genintegral->get_Integrated_Lumi();
  cout<<"InitRun integrated lumi "<< m_intlumi <<endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

int MyJetAnalysis::process_event(PHCompositeNode* topNode)
{
  if (Verbosity() >= MyJetAnalysis::VERBOSITY_SOME)
    cout << "MyJetAnalysis::process_event() entered" << endl;

  ////////////////
  // Get nodes
  ////////////////

  //PHG4TruthInfoContainer* _truthinfo = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");
  PHHepMCGenEventMap* geneventmap = findNode::getClass<PHHepMCGenEventMap>(topNode, "PHHepMCGenEventMap");
  if (!geneventmap)
  {
    std::cout << PHWHERE << " - Fatal error - missing node PHHepMCGenEventMap" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  int _embedding_id=1;
  PHHepMCGenEvent* genevt = geneventmap->get(_embedding_id);
  if (!genevt)
  {
    std::cout << PHWHERE << " - Fatal error - node PHHepMCGenEventMap missing subevent with embedding ID " << _embedding_id;
    std::cout << ". Print PHHepMCGenEventMap:";
    geneventmap->identify();
    return Fun4AllReturnCodes::ABORTRUN;
  }

  HepMC::GenEvent* truthevent = genevt->getEvent();
  if (!truthevent)
  {
    cout << PHWHERE << "no evt pointer under phhepmvgenevent found " << endl;
    return 0;
  }

  JetEvalStack *_jetevalstack = new JetEvalStack(topNode, m_truthJetName.c_str(), m_truthJetName.c_str());
  JetTruthEval *trutheval = _jetevalstack->get_truth_eval();

  JetMap* truth_jets = findNode::getClass<JetMap>(topNode, m_truthJetName);
  if (!truth_jets)
  {
    std::cout << PHWHERE << " - Fatal error - node " << m_truthJetName << " JetMap missing." << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  const double jet_radius = truth_jets->get_par();
  if(verbosity>10)cout<<"jet radius "<<jet_radius<<endl;
  int ijet_t = 0;
//  bool pass_event = false;
 ++m_event;
  ////////////////   
  // Parton info
  ////////////////
  HepMC::PdfInfo *pdfinfo = truthevent->pdf_info();
  m_process_id = truthevent->signal_process_id();
  m_parton_id1 = pdfinfo->id1();
  m_parton_id2 = pdfinfo->id2();
  m_x1 = pdfinfo->x1();
  m_x2 = pdfinfo->x2();
  //////////////////
  // Photon tagging 
  ////////////////// 
  m_photon_phi=-99.; m_photon_eta=-99.; m_photon_pt=-99.;
  if(m_process_id == 241 || m_process_id == 242)
  {
    int code = photon_tagging(truthevent);
    if(!code) cout<<"photon not found!"<<endl;
  }  

  if(verbosity>10)cout<<"process " <<  m_process_id << " x1 "<<m_x1<<" x2 "<<m_x2<<" partid1 "<<  m_parton_id1 <<  " partid2 "<< m_parton_id2<<endl;
  if(verbosity>10)cout<<"# of jets "<<truth_jets->size()<<endl;
  int njets = truth_jets->size();
  m_nTruthConstituents =0;
  for (JetMap::Iter iter = truth_jets->begin(); iter != truth_jets->end();
       ++iter)
  {
    if(ijet_t < njets-2){     ijet_t++; continue;	}  
    if(verbosity>10)cout<<"jet # "<<ijet_t<<endl;
    Jet* this_jet = iter->second;
    if(!this_jet)continue;
    ///////////////////
    // jet level info
    ///////////////////
    int jetid= -99;
    if(ijet_t == (njets-1)) jetid = 0;//leading jet
    else if(ijet_t == (njets-2)) jetid=1;//subleading jet
    m_truthID = jetid;
    m_truthNComponent = this_jet->size_comp();
    m_truthEta = this_jet->get_eta();
    m_truthPhi = this_jet->get_phi();
    m_truthE = this_jet->get_e();
    m_truthPt = this_jet->get_pt();
    if(verbosity>10) cout<<" jet pt "<<m_truthPt<<endl;
    if(jetid==1){
    m_parton_eta=-99.; m_parton_phi=-99.; m_parton_zt=-99.;
    const int jet_flavor = parton_tagging(this_jet, truthevent, jet_radius+m_matchDist);
    m_flavor =  jet_flavor;
    if(verbosity>10 ) cout<<" jet flavor " << jet_flavor<<" zt "<<m_parton_zt<< " partid1 "<<  m_parton_id1 <<  " partid2 "<< m_parton_id2<<endl;
    }
  ////////////////////////
  // particle level info
  ////////////////////////


  m_nTruthConstituents =0;
  int method =2;
  if(method==1)
  {
    std::set<PHG4Particle *> truthjetcomp =
        trutheval->all_truth_particles(this_jet);

    for (std::set<PHG4Particle *>::iterator iter2 = truthjetcomp.begin();
         iter2 != truthjetcomp.end();
         ++iter2)
    {
      PHG4Particle *truthpart = *iter2;
      if (!truthpart)
      {
        cout << "no truth particles in the jet??" << endl;
        break;
      }
	m_truth_part_px[m_nTruthConstituents] = truthpart->get_px(); 
        m_truth_part_py[m_nTruthConstituents] = truthpart->get_py();
        m_truth_part_pz[m_nTruthConstituents] = truthpart->get_pz();
        m_truth_part_e[m_nTruthConstituents] = truthpart->get_e();
        m_truth_part_pid[m_nTruthConstituents] = truthpart->get_pid();
	float charge = TDatabasePDG::Instance()->GetParticle(m_truth_part_pid[m_nTruthConstituents])->Charge();
	m_truth_part_qe[m_nTruthConstituents] = charge;
	m_nTruthConstituents++;
	if(verbosity>10) cout<<"m_truth_part_pid "<<m_truth_part_pid[m_nTruthConstituents]<<endl;
    }// for particles

/*
    if (m_pt < 10 || fabs(m_eta) > 5)
      continue;
    if (m_pt > m_ptRange.first && m_pt < m_ptRange.second && (m_eta) > m_etaRanga.first && (m_eta) < m_etaRange.second)
    {
      pass_event = true;
    }
    else
    {
      continue;
    }
*/
    m_T->Fill();
    ijet_t++;
  }//method 1
  else if(method==2)
  {

   for (HepMC::GenEvent::particle_const_iterator p = truthevent->particles_begin();
       p != truthevent->particles_end(); ++p)
   {
    //https://rivet.hepforge.org/code/hepmc.bak/classHepMC_1_1GenParticle.html	  
    HepMC::GenParticle* gpart = (*p);
	    
    float dR = delta_r(gpart->momentum().pseudoRapidity(), m_truthEta,
                      gpart->momentum().phi(), m_truthPhi);

    //http://home.thep.lu.se/~torbjorn/pythia81html/ParticleProperties.html    
    int status = gpart->status();
    //if ((status==23 || status==33) && dR < 1.0) cout<<"outgoing parton "<< gpart->pdg_id()<<" status "<<status<<endl;
    if (status !=1) continue; // only FS particles
    if (dR > jet_radius) continue;
    m_truth_part_px[m_nTruthConstituents] = gpart->momentum().px();
        m_truth_part_py[m_nTruthConstituents] = gpart->momentum().py();
        m_truth_part_pz[m_nTruthConstituents] = gpart->momentum().pz();
        m_truth_part_e[m_nTruthConstituents] = gpart->momentum().e();
	m_truth_part_phi[m_nTruthConstituents] = gpart->momentum().phi();
	m_truth_part_pid[m_nTruthConstituents] = gpart->pdg_id();
        float charge = TDatabasePDG::Instance()->GetParticle(m_truth_part_pid[m_nTruthConstituents])->Charge();
        m_truth_part_qe[m_nTruthConstituents] = charge;
        m_nTruthConstituents++;
        if(verbosity>10) cout<<"  m_truth_part_pid "<<m_truth_part_pid[m_nTruthConstituents]<< " px "<<gpart->momentum().px()<<" py "<<gpart->momentum().py()
			<<" pz "<< gpart->momentum().pz()<<endl;    
   }  //       for (HepMC::GenEvent::particle_const_iterator p =  
    m_T->Fill();
    ijet_t++;
  }
  if(verbosity>10) cout<<"---   nTruthConstituents " <<m_nTruthConstituents<<" m_truthNComponent "<<m_truthNComponent << endl;
  }  //   for (JetMap::Iter iter = jets->begin(); iter != jets->end(); ++iter)
  
// reco tree
/*
  m_jetEvalStack->next_event(topNode);
  JetRecoEval* recoeval = m_jetEvalStack->get_reco_eval();
  ++m_event;
  // interface to jets
  JetMap* jets = findNode::getClass<JetMap>(topNode, m_recoJetName);
  if (!jets)
  {
    cout
        << "MyJetAnalysis::process_event - Error can not find DST JetMap node "
        << m_recoJetName << endl;
    exit(-1);
  }
  // interface to tracks
  SvtxTrackMap* trackmap = findNode::getClass<SvtxTrackMap>(topNode, "SvtxTrackMap");
  if (!trackmap)
  {
    trackmap = findNode::getClass<SvtxTrackMap>(topNode, "TrackMap");
    if (!trackmap)
    {
      cout
          << "MyJetAnalysis::process_event - Error can not find DST trackmap node SvtxTrackMap" << endl;
      exit(-1);
    }
  }
  for (JetMap::Iter iter = jets->begin(); iter != jets->end(); ++iter)
  {
    Jet* jet = iter->second;
    assert(jet);
    bool eta_cut = (jet->get_eta() >= m_etaRange.first) and (jet->get_eta() <= m_etaRange.second);
    bool pt_cut = (jet->get_pt() >= m_ptRange.first) and (jet->get_pt() <= m_ptRange.second);
    if ((not eta_cut) or (not pt_cut))
    {
      if (Verbosity() >= MyJetAnalysis::VERBOSITY_MORE)
      {
        cout << "MyJetAnalysis::process_event() - jet failed acceptance cut: ";
        cout << "eta cut: " << eta_cut << ", ptcut: " << pt_cut << endl;
        cout << "jet eta: " << jet->get_eta() << ", jet pt: " << jet->get_pt() << endl;
        jet->identify();
      }
      continue;
    }
    // fill trees - jet spectrum
    Jet* truthjet = recoeval->max_truth_jet_by_energy(jet);
    m_id = jet->get_id();
    m_nComponent = jet->size_comp();
    m_eta = jet->get_eta();
    m_phi = jet->get_phi();
    m_e = jet->get_e();
    m_pt = jet->get_pt();
    m_truthID = -1;
    m_truthNComponent = -1;
    m_truthEta = NAN;
    m_truthPhi = NAN;
    m_truthE = NAN;
    m_truthPt = NAN;
    if (truthjet)
    {
      m_truthID = truthjet->get_id();
      m_truthNComponent = truthjet->size_comp();
      m_truthEta = truthjet->get_eta();
      m_truthPhi = truthjet->get_phi();
      m_truthE = truthjet->get_e();
      m_truthPt = truthjet->get_pt();
    }
    // fill trees - jet track matching
    m_nMatchedTrack = 0;
    for (SvtxTrackMap::Iter iter = trackmap->begin();
         iter != trackmap->end();
         ++iter)
    {
      SvtxTrack* track = iter->second;
      TVector3 v(track->get_px(), track->get_py(), track->get_pz());
      const double dEta = v.Eta() - m_eta;
      const double dPhi = v.Phi() - m_phi;
      const double dR = sqrt(dEta * dEta + dPhi * dPhi);
      if (dR < m_trackJetMatchingRadius)
      {
        //matched track to jet
        assert(m_nMatchedTrack < kMaxMatchedTrack);
        m_trackdR[m_nMatchedTrack] = dR;
        m_trackpT[m_nMatchedTrack] = v.Perp();
        ++m_nMatchedTrack;
      }
      if (m_nMatchedTrack >= kMaxMatchedTrack)
      {
        cout << "MyJetAnalysis::process_event() - reached max track that matching a jet. Quit iterating tracks" << endl;
        break;
      }
    }  //    for (SvtxTrackMap::Iter iter = trackmap->begin();
  }  //   for (JetMap::Iter iter = jets->begin(); iter != jets->end(); ++iter)
*/

  return Fun4AllReturnCodes::EVENT_OK;
}

int MyJetAnalysis::parton_tagging(Jet* this_jet, HepMC::GenEvent* theEvent,
                                      const double match_radius)
{
  float this_pt = this_jet->get_pt();
  float this_phi = this_jet->get_phi();
  float this_eta = this_jet->get_eta();

  int jet_flavor = 0;
  double jet_parton_zt = 0;
  float dist = 999;
  //std::cout << " truth jet #" << ijet_t << ", pt / eta / phi = " << this_pt << " / " << this_eta << " / " << this_phi << ", checking flavor" << std::endl;

  for (HepMC::GenEvent::particle_const_iterator p = theEvent->particles_begin();
       p != theEvent->particles_end(); ++p)
  {

    float dR = delta_r((*p)->momentum().pseudoRapidity(), this_eta,
                      (*p)->momentum().phi(), this_phi);
    if (dR > match_radius)
      continue;
    int status = (*p)->status();
    if(verbosity>10 && (status == 23 || status ==33) ) cout<<"outgoing parton "<< (*p)->pdg_id()<<" status "<<status<< " dR "<<dR<< endl;
    int pidabs = abs((*p)->pdg_id());
    const double zt = (*p)->momentum().perp() / this_pt;

    if (pidabs == TDatabasePDG::Instance()->GetParticle("u")->PdgCode()      //
        or pidabs == TDatabasePDG::Instance()->GetParticle("d")->PdgCode()   //
        or pidabs == TDatabasePDG::Instance()->GetParticle("s")->PdgCode()
        or pidabs == TDatabasePDG::Instance()->GetParticle("c")->PdgCode()
        or pidabs == TDatabasePDG::Instance()->GetParticle("b")->PdgCode()
 	or pidabs == TDatabasePDG::Instance()->GetParticle("g")->PdgCode()
	)  // default 0
    {
      //if (pidabs > abs(jet_flavor)) // for heavy flavor
      //if(dR < dist && (status == 23 || status ==33)) // distance method
      if (zt > jet_parton_zt && (status == 23 || status ==33)) // zt method
      {
	dist = dR; 
        jet_parton_zt = zt;
        jet_flavor = (*p)->pdg_id();
	m_parton_eta = (*p)->momentum().pseudoRapidity();
	m_parton_phi = (*p)->momentum().phi();
	m_parton_zt = zt;
      }
    }
  }  //       for (HepMC::GenEvent::particle_const_iterator p =
/*
  if (abs(jet_flavor) == TDatabasePDG::Instance()->GetParticle("g")->PdgCode())
  {
    _h2_g->Fill(this_jet->get_pt(), this_eta);
  }
  else if (abs(jet_flavor) == TDatabasePDG::Instance()->GetParticle("u")->PdgCode())
  {
    _h2_u->Fill(this_jet->get_pt(), this_eta);
  }
  else if (abs(jet_flavor) == TDatabasePDG::Instance()->GetParticle("d")->PdgCode())
  {
    _h2_d->Fill(this_jet->get_pt(), this_eta);
  }
  else if (abs(jet_flavor) == TDatabasePDG::Instance()->GetParticle("s")->PdgCode())
  {
    _h2_s->Fill(this_jet->get_pt(), this_eta);
  }  
*/
  //this_jet->set_property(static_cast<Jet::PROPERTY>(prop_JetPartonFlavor),jet_flavor);
  //this_jet->set_property(static_cast<Jet::PROPERTY>(prop_JetPartonZT),jet_parton_zt);
  //          this_jet->identify();
  if(verbosity>10)cout<<"distance "<<dist<<" zt "<<jet_parton_zt<<endl;
  return jet_flavor;
}


bool MyJetAnalysis::photon_tagging(HepMC::GenEvent* theEvent)
{
  bool found = false;
  for (HepMC::GenEvent::particle_const_iterator p = theEvent->particles_begin();
       p != theEvent->particles_end(); ++p)
  {
    int status = (*p)->status();
    int pidabs = abs((*p)->pdg_id());
    if((status == 23 || status ==33) && pidabs == 23/*TDatabasePDG::Instance()->GetParticle("gamma")->PdgCode()*/) 
    {                                 
    // WeakZ0:gmZmode - irrespective of the option used, the particle produced will always be assigned code 23 for Z^0, 
    // and open decay channels is purely dictated by what is set for the Z^0.
      if(found) cout<<"more than 1 photon found "<<endl;
      m_photon_phi = (*p)->momentum().phi();
      m_photon_eta = (*p)->momentum().pseudoRapidity();
      m_photon_pt = sqrt(pow((*p)->momentum().px(),2) + pow((*p)->momentum().py(),2));
      found=true;
      if(verbosity>10)cout<<"photon "<< (*p)->pdg_id()<<" status "<<status<< 
			" phi "<< m_photon_phi<<" eta "<<m_photon_eta<<" pt "<< m_photon_pt<<endl;

    }
  }  //       for (HepMC::GenEvent::particle_const_iterator p =  

  return found;
}

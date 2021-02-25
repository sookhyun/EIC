#include "MyJetAnalysis.h"

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/Fun4AllServer.h>
#include <fun4all/PHTFileServer.h>

#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>

#include <g4eval/JetEvalStack.h>

#include <trackbase_historic/SvtxTrackMap.h>

#include <g4jets/JetMap.h>

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

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>

using namespace std;

MyJetAnalysis::MyJetAnalysis(const std::string& recojetname, const std::string& truthjetname, const std::string& outputfilename,int flavor)
  : SubsysReco("MyJetAnalysis_" + recojetname + "_" + truthjetname)
  , m_recoJetName(recojetname)
  , m_truthJetName(truthjetname)
  , m_outputFileName(outputfilename)
  , m_etaRange(-1, 1)
  , m_ptRange(5, 100)
  , m_trackJetMatchingRadius(.7)
  , m_hInclusiveE(nullptr)
  , m_hInclusiveEta(nullptr)
  , m_hInclusivePhi(nullptr)
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
  , m_flavor(flavor)
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

  // Histograms
  m_hInclusiveE = new TH1F(
      "hInclusive_E",  //
      TString(m_recoJetName) + " inclusive jet E;Total jet energy (GeV)", 100, 0, 100);

  m_hInclusiveEta =
      new TH1F(
          "hInclusive_eta",  //
          TString(m_recoJetName) + " inclusive jet #eta;#eta;Jet energy density", 50, -1, 1);
  m_hInclusivePhi =
      new TH1F(
          "hInclusive_phi",  //
          TString(m_recoJetName) + " inclusive jet #phi;#phi;Jet energy density", 50, -M_PI, M_PI);

  //Trees
  m_T = new TTree("T", "MyJetAnalysis Tree");

  //      int m_event;
  m_T->Branch("m_event", &m_event, "event/I");
  //      int m_id;
  m_T->Branch("id", &m_id, "id/I");
  //      int m_nComponent;
  m_T->Branch("nComponent", &m_nComponent, "nComponent/I");
  //      float m_eta;
  m_T->Branch("eta", &m_eta, "eta/F");
  //      float m_phi;
  m_T->Branch("phi", &m_phi, "phi/F");
  //      float m_e;
  m_T->Branch("e", &m_e, "e/F");
  //      float m_pt;
  m_T->Branch("pt", &m_pt, "pt/F");
  //
  //      int m_truthID;
  m_T->Branch("truthID", &m_truthID, "truthID/I");
  //      int m_truthNComponent;
  m_T->Branch("truthNComponent", &m_truthNComponent, "truthNComponent/I");
  //      float m_truthEta;
  m_T->Branch("truthEta", &m_truthEta, "truthEta/F");
  //      float m_truthPhi;
  m_T->Branch("truthPhi", &m_truthPhi, "truthPhi/F");
  //      float m_truthE;
  m_T->Branch("truthE", &m_truthE, "truthE/F");
  //      float m_truthPt;
  m_T->Branch("truthPt", &m_truthPt, "truthPt/F");
  //
  //      //! number of matched tracks
  //      int m_nMatchedTrack;
  m_T->Branch("nMatchedTrack", &m_nMatchedTrack, "nMatchedTrack/I");
  //      std::array<float, kMaxMatchedTrack> m_trackdR;
  m_T->Branch("id", m_trackdR.data(), "trackdR[nMatchedTrack]/F");
  //      std::array<float, kMaxMatchedTrack> m_trackpT;
  m_T->Branch("id", m_trackpT.data(), "trackpT[nMatchedTrack]/F");

  return Fun4AllReturnCodes::EVENT_OK;
}

int MyJetAnalysis::End(PHCompositeNode* topNode)
{
  cout << "MyJetAnalysis::End - Output to " << m_outputFileName << endl;
  PHTFileServer::get().cd(m_outputFileName);

  m_hInclusiveE->Write();
  m_hInclusiveEta->Write();
  m_hInclusivePhi->Write();
  m_T->Write();

  return Fun4AllReturnCodes::EVENT_OK;
}

int MyJetAnalysis::InitRun(PHCompositeNode* topNode)
{
//  m_jetEvalStack = shared_ptr<JetEvalStack>(new JetEvalStack(topNode, m_recoJetName, m_truthJetName));
//  m_jetEvalStack->get_stvx_eval_stack()->set_use_initial_vertex(initial_vertex);
  return Fun4AllReturnCodes::EVENT_OK;
}

int MyJetAnalysis::process_event(PHCompositeNode* topNode)
{
  if (Verbosity() >= MyJetAnalysis::VERBOSITY_SOME)
    cout << "MyJetAnalysis::process_event() entered" << endl;

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

  HepMC::GenEvent* theEvent = genevt->getEvent();

  JetMap* truth_jets = findNode::getClass<JetMap>(topNode, m_truthJetName);
  if (!truth_jets)
  {
    std::cout << PHWHERE << " - Fatal error - node " << m_truthJetName << " JetMap missing." << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  const double jet_radius = truth_jets->get_par();
  cout<<"jet radius "<<jet_radius<<endl;
  int ijet_t = 0;
//  bool pass_event = false;
  for (JetMap::Iter iter = truth_jets->begin(); iter != truth_jets->end();
       ++iter)
  {
    Jet* this_jet = iter->second;
/*
    float this_pt = this_jet->get_pt();
    float this_eta = this_jet->get_eta();

    if (this_pt < 10 || fabs(this_eta) > 5)
      continue;

    //_h2all->Fill(this_jet->get_pt(), this_eta);

    if (this_pt > _pt_min && this_pt < _pt_max && (this_eta) > _eta_min && (this_eta) < _eta_max)
    {
//      pass_event = true;
      //_h2->Fill(this_jet->get_pt(), this_eta);
      //if (Verbosity() >= HFJetTruthTrigger::VERBOSITY_MORE)
        this_jet->identify();
    }
    else
    {
      continue;
    }
*/
    const int jet_flavor = parton_tagging(this_jet, theEvent, jet_radius);
    cout<<"jet flavor "<<jet_flavor<<endl;
    
    if (abs(jet_flavor) == m_flavor)
    {
      //pass_event = true;
    }

    ijet_t++;
  }

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

    // fill histograms
    assert(m_hInclusiveE);
    m_hInclusiveE->Fill(jet->get_e());
    assert(m_hInclusiveEta);
    m_hInclusiveEta->Fill(jet->get_eta());
    assert(m_hInclusivePhi);
    m_hInclusivePhi->Fill(jet->get_phi());

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

    m_T->Fill();
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

  //std::cout << " truth jet #" << ijet_t << ", pt / eta / phi = " << this_pt << " / " << this_eta << " / " << this_phi << ", checking flavor" << std::endl;

  //TODO: lack taggign scheme of gluon splitting -> QQ_bar
  for (HepMC::GenEvent::particle_const_iterator p = theEvent->particles_begin();
       p != theEvent->particles_end(); ++p)
  {
    float dR = deltaR((*p)->momentum().pseudoRapidity(), this_eta,
                      (*p)->momentum().phi(), this_phi);
    if (dR > match_radius)
      continue;

    int pidabs = abs((*p)->pdg_id());
    const double zt = (*p)->momentum().perp() / this_pt;

    if (pidabs == TDatabasePDG::Instance()->GetParticle("u")->PdgCode()      //
        or pidabs == TDatabasePDG::Instance()->GetParticle("d")->PdgCode()   //
        or pidabs == TDatabasePDG::Instance()->GetParticle("s")->PdgCode())  // handle heavy quarks only. All other favor tagged as default 0
    {
      if (pidabs > abs(jet_flavor))  // heavy quark found
      {
        jet_parton_zt = zt;
        jet_flavor = (*p)->pdg_id();
      }
      else if (pidabs == abs(jet_flavor))  // same quark mass. next compare zt
      {
        if (zt > jet_parton_zt)
        {
          jet_parton_zt = zt;
          jet_flavor = (*p)->pdg_id();
        }
      }

      if (pidabs == TDatabasePDG::Instance()->GetParticle("b")->PdgCode())
      {
        if (Verbosity() >= MyJetAnalysis::VERBOSITY_MORE)
          std::cout << __PRETTY_FUNCTION__
                    << " --BOTTOM--> pt / eta / phi = "
                    << (*p)->momentum().perp() << " / "
                    << (*p)->momentum().pseudoRapidity() << " / "
                    << (*p)->momentum().phi() << std::endl;
      }
      else if (pidabs == TDatabasePDG::Instance()->GetParticle("c")->PdgCode())
      {
        if (Verbosity() >= MyJetAnalysis::VERBOSITY_MORE)
          std::cout << __PRETTY_FUNCTION__
                    << " --CHARM --> pt / eta / phi = "
                    << (*p)->momentum().perp() << " / "
                    << (*p)->momentum().pseudoRapidity() << " / "
                    << (*p)->momentum().phi() << std::endl;
      }
    }
  }  //       for (HepMC::GenEvent::particle_const_iterator p =

  if (abs(jet_flavor) == TDatabasePDG::Instance()->GetParticle("b")->PdgCode())
  {
    //_h2_b->Fill(this_jet->get_pt(), this_eta);
  }
  else if (abs(jet_flavor) == TDatabasePDG::Instance()->GetParticle("c")->PdgCode())
  {
    //_h2_c->Fill(this_jet->get_pt(), this_eta);
  }

  //this_jet->set_property(static_cast<Jet::PROPERTY>(prop_JetPartonFlavor),jet_flavor);
  //this_jet->set_property(static_cast<Jet::PROPERTY>(prop_JetPartonZT),jet_parton_zt);
  //          this_jet->identify();

  return jet_flavor;
}



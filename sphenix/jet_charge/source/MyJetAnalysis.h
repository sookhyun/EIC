#ifndef MYJETANALYSIS_H
#define MYJETANALYSIS_H

#include <fun4all/SubsysReco.h>
#include <cmath>
#include <memory>
#include <string>
#include <utility>  // std::pair, std::make_pair

#include <array>

class PHCompositeNode;
class Jet;
class JetEvalStack;
class TTree;
class TH1;
namespace HepMC
{
  class GenEvent;
}

/// \class MyJetAnalysis
class MyJetAnalysis : public SubsysReco
{
 public:
  MyJetAnalysis(
      const std::string &recojetname = "AntiKt_Tower_r04",
      const std::string &truthjetname = "AntiKt_Truth_r04",
      const std::string &outputfilename = "myjetanalysis.root",
      int flavor = 3);

  virtual ~MyJetAnalysis();

  //! set eta range
  void
  setEtaRange(double low, double high)
  {
    m_etaRange.first = low;
    m_etaRange.second = high;
  }
  //! set eta range
  void
  setPtRange(double low, double high)
  {
    m_ptRange.first = low;
    m_ptRange.second = high;
  }
  void use_initial_vertex(const bool b = true) {initial_vertex = b;}
  int Init(PHCompositeNode *topNode);
  int InitRun(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);
  int End(PHCompositeNode *topNode);


  float
  deltaR(float eta1, float eta2, float phi1, float phi2)
  {

    float deta = eta1 - eta2;
    float dphi = phi1 - phi2;
    if (dphi > +3.14159)
      dphi -= 2 * 3.14159;
    if (dphi < -3.14159)
      dphi += 2 * 3.14159;

    return sqrt(pow(deta, 2) + pow(dphi, 2));

  }

  double
  get_eta_max() const
  {
    return _eta_max;
  }

  void
  set_eta_max(double etaMax)
  {
    _eta_max = etaMax;
  }

  double
  get_eta_min() const
  {
    return _eta_min;
  }

  void
  set_eta_min(double etaMin)
  {
    _eta_min = etaMin;
  }

  double
  get_pt_max() const
  {
    return _pt_max;
  }

  void
  set_pt_max(double ptMax)
  {
    _pt_max = ptMax;
  }

  double
  get_pt_min() const
  {
    return _pt_min;
  }

  void
  set_pt_min(double ptMin)
  {
    _pt_min = ptMin;
  }

 private:
  //! tag jet flavor by parton matching, like PRL 113, 132301 (2014)
  int
  parton_tagging(Jet * jet, HepMC::GenEvent*, const double match_radius);
  //! cache the jet evaluation modules
  std::shared_ptr<JetEvalStack> m_jetEvalStack;

  std::string m_recoJetName;
  std::string m_truthJetName;
  std::string m_outputFileName;

  //! eta range
  std::pair<double, double> m_etaRange;

  //! pT range
  std::pair<double, double> m_ptRange;

  //! flag to use initial vertex in track evaluator
  bool initial_vertex = false;

  //! max track-jet matching radius
  double m_trackJetMatchingRadius;

  //! Output histograms
  TH1 *m_hInclusiveE;
  TH1 *m_hInclusiveEta;
  TH1 *m_hInclusivePhi;

  //! Output Tree variables
  TTree *m_T;

  int m_event;
  int m_id;
  int m_nComponent;
  float m_eta;
  float m_phi;
  float m_e;
  float m_pt;

  int m_process_id;
  int m_parton_id1;
  int m_parton_id2;
  float m_x1;
  float m_x2;

  int m_truthID;
  int m_truthNComponent;
  int m_nTruthConstituents;
  float m_truthEta;
  float m_truthPhi;
  float m_truthE;
  float m_truthPt;

  //! number of matched tracks
  int m_nMatchedTrack;

  int m_flavor;
  int _maxevent;

  double _pt_min;
  double _pt_max;

  double _eta_min;
  double _eta_max;


  enum
  {
    //! max number of tracks
    kMaxMatchedTrack = 1000
  };

  std::array<float, kMaxMatchedTrack> m_trackdR;
  std::array<float, kMaxMatchedTrack> m_trackpT;
  std::array<float, kMaxMatchedTrack> m_truth_part_px;  
  std::array<float, kMaxMatchedTrack> m_truth_part_py;
  std::array<float, kMaxMatchedTrack> m_truth_part_pz;
  std::array<float, kMaxMatchedTrack> m_truth_part_e;
  std::array<float, kMaxMatchedTrack> m_truth_part_pid;
};

#endif  // MYJETANALYSIS_H

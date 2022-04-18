/**
 * @file FastJetAlgo.h
 * @author Sookhyun Lee
 * @date May 2020
 */

#ifndef _FASTJETALGO_H
#define _FASTJETALGO_H

#include <fastjet/ClusterSequence.hh>
#include <fastjet/Selector.hh>
#include <fastjet/JetDefinition.hh>
#include <fastjet/PseudoJet.hh>


class FastJetAlgo
{
  public:

   enum Frame
   {
     Lab = 0;
     CM = 1;
     Breit = 2;
   };

   FastJetAlgo(fastjet::JetAlgorithm _algo, 
	       Frame _frame, 
	       double _R) {}
   virtual ~FastJetAlgo() {}

   void set_jet_pt_min(float _jet_pt_min){jet_pt_min=_jet_pt_min;}

   void set_verbosity(int _verbosity){verbosity = _verbosity;}
  private:
    int verbosity;
    fastjet::JetAlgorithm algo;
    fastjet::JetDefinition jetdef;
    Frame frame;
    float R;


};


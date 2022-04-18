/**
 * @file FastJetAlgo.cc
 * @author Sookhyun Lee
 * @date May 2020
 */

#include "FastJetAlgo.h"


#include <iostream>

using namespace std;

FastJetAlgo::FastJetAlgo(fastjet::JetAlgorithm _algo=fastjet::antikt_algorithm,
               		 Frame _frame = Frame::Lab,
               		 double _R =1.0
			)
  : verbosity(0)
  , algo(_algo)
  , frame(_frame) 
  , R(_R)
{
  fastjet::ClusterSequence clseq;
  if (verbosity)
  {
    clseq.print_banner();
  }
  else
  {
    ostringstream nullstream;
    clseq.set_fastjet_banner_stream(&nullstream);
    clseq.print_banner();
    clseq.set_fastjet_banner_stream(&cout);
  }
}  


}

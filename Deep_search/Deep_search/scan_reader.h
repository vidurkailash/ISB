#ifndef scan_reader_class_h
#define scan_reader_class_h

#include "deep_class.h"
#include "MSReader.h"

#include <numeric>

class scan_reader {
public:
  match_lists mzml(peptide_lists& my_peptide_lists, my_parameters& my_params);

private:

  //MH: Binary search function. Could and should be overloaded for PPM instead of MZ tolerance
  int findPeakMZ(MSToolkit::Spectrum& spec, double mz, double tol);

  //MH: quick helper function for sorting from lowest to highes RT
  static bool compareRTime(const my_features& a, const my_features& b) { return a.rtime < b.rtime; }
  static bool compareRTimeD(const my_features& a, const my_features& b) { return a.d_pep_seq_rt < b.d_pep_seq_rt; }
  static bool compareSeqScan(const my_markers& a, const my_markers& b) { 
    int i=a.pep_seq.compare(b.pep_seq);
    if(i==0) return (a.spec_sn<b.spec_sn);
    else return (i<0);
  }
};

#endif
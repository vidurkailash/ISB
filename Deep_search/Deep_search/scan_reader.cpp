#include "scan_reader.h"

using namespace std;
using namespace MSToolkit;


bool scan_reader::mzml(peptide_lists& my_peptide_lists, my_parameters& my_params) {

	//READ MZML FILE FOR BOTH FULLY TRYPTIC AND MISCLEAVED PEPTIDES (MZML)
	MSReader myfile;
	Spectrum mySpec;

	myfile.setFilter(MS1);
	if(!myfile.readFile(my_params.mzml.c_str(), mySpec)) return false;


	//MH: We need to get through precursor peak extraction fast.
	//To do so, we need to do it a) using a single pass through all arrays, and
	//b) using fast lookups instead of full iteration.
	//Step 1 is to sort all PSMs by retention time
	sort(my_peptide_lists.all_psm.begin(), my_peptide_lists.all_psm.end(), compareRTime);
	size_t pepIndex = 0; // we will start from the first sorted peptide

	//MH: Start reading our mzML file.
	while (mySpec.getRTime() != 0) { //this line checks that the current spectrum is valid. The moment it is not, the end of file was reached.

		size_t i = pepIndex;
		while (i < my_peptide_lists.all_psm.size() && (my_peptide_lists.all_psm[i].xml_rtime / 60) < (mySpec.getRTime() - my_params.ret_time)) i++;
		if (i == my_peptide_lists.all_psm.size()) break; //if we've checked every peptide, stop now.
		pepIndex = i; //mark our new start point for the next iteration

		while (i < my_peptide_lists.all_psm.size() && (my_peptide_lists.all_psm[i].xml_rtime / 60) < (mySpec.getRTime() + my_params.ret_time)) {
			//compute the desired m/z value to find, and a tolerance around that value
			double mz = (my_peptide_lists.all_psm[i].pre_neutral_mass + my_peptide_lists.all_psm[i].charge * 1.007276466) / my_peptide_lists.all_psm[i].charge;
			double tolerance = my_params.ppm * mz / 1000000;
			int ret = findPeakMZ(mySpec, mz, tolerance);
			dsXIC dsx;
			dsx.rTime = mySpec.getRTime();
			dsx.tot = 0;
			if (ret > -1) {
				dsx.intensity = mySpec[ret].intensity;
			} else {
				dsx.intensity = 0;
			}
			my_peptide_lists.all_psm[i].XIC.push_back(dsx);
			i++;//MH: go to the next PSM
		}
		myfile.readFile(NULL, mySpec); //read the next spectrum
	}

	return true;

}



//MH: Function to binary search for an mz value within a specific tolerance
int scan_reader::findPeakMZ(MSToolkit::Spectrum& spec, double mz, double tol){
  int sz = spec.size();
  int lower = 0;
  int mid = sz / 2;
  int upper = sz;

  //our boundaries
  double LB=mz-tol;
  double UB=mz+tol;

  //stop now if the spectrum is empty.
  if(sz<1) return -1;

  while(lower<upper){
    if(spec[mid].mz>UB){ //too high
      upper = mid;
      mid = (lower + upper) / 2;
    } else if(spec[mid].mz<LB) { //too low
      lower = mid + 1;
      mid = (lower + upper) / 2;
    } else { //we have a match, now we must step left and right to make sure it's the max
      float max=spec[mid].intensity;
      int maxI=mid;
      int i=mid;
      while(i>0 && spec[i].mz>LB){
        i--;
        if(spec[i].intensity>max){
          max=spec[i].intensity;
          maxI=i;
        }
      }
      i=mid;
      while(i<sz-1 && spec[i].mz<UB){
        i++;
        if (spec[i].intensity > max) {
          max = spec[i].intensity;
          maxI = i;
        }
      }
      return maxI; 
    }
  }
  return -1; //-1 means can't find peak;
}


metrics deep_functions::calc1(peptide_lists& my_peptide_lists, metrics& my_metrics, my_parameters& my_params) {


	//int count2 = 0;
	//int count1 = 0;
	//float count0 = 0;

	//for (size_t i = 0; i < my_peptide_lists.xic_mc_results.size(); i++) {
	//	for (size_t j = 0; j < my_peptide_lists.xic_ft_results.size(); j++) {
	//		size_t present = my_peptide_lists.xic_mc_results[i].pep_seq.find(my_peptide_lists.xic_ft_results[j].pep_seq);
	//		if (present != string::npos && (my_peptide_lists.xic_mc_results[i].pep_seq != my_peptide_lists.xic_ft_results[j].pep_seq)) {
	//			my_peptide_lists.peptide_matches.push_back(dsPair()); 
	//			my_peptide_lists.peptide_matches.back().missIndex = i; 
	//			my_peptide_lists.peptide_matches.back().trypIndex = j;
	//			my_peptide_lists.peptide_matches.back().ft_pep_seq = my_peptide_lists.xic_ft_results[j].pep_seq;
	//			my_peptide_lists.peptide_matches.back().mc_pep_seq = my_peptide_lists.xic_mc_results[i].pep_seq;
	//			my_peptide_lists.peptide_matches.back().ft_areaXIC = my_peptide_lists.xic_ft_results[j].areaXIC;
	//			my_peptide_lists.peptide_matches.back().mc_areaXIC = my_peptide_lists.xic_mc_results[i].areaXIC;
	//		}
	//	}
	//}
	//my_metrics.find_num = my_peptide_lists.peptide_matches.size();

	//sort(my_peptide_lists.peptide_matches.begin(), my_peptide_lists.peptide_matches.end(), compareTrypIndex);
	//

	///*for (int i = 0; i < my_peptide_lists.peptide_matches.size(); i++) {
	//	cout << my_peptide_lists.peptide_matches[i].trypIndex << "  " << my_peptide_lists.peptide_matches[i].missIndex << "  " << my_peptide_lists.peptide_matches[i].ft_pep_seq << "  " << my_peptide_lists.peptide_matches[i].mc_pep_seq << endl; 
	//}
	//cout << my_peptide_lists.peptide_matches.size() << endl; */

	//
	//vector<dsPair> tmp1; 
	//
	//my_peptide_lists.peptide_matches_ftm.push_back(tmp1); 
	//my_peptide_lists.peptide_matches_ftm.back().push_back(my_peptide_lists.peptide_matches[0]);
	//for (size_t i = 1; i < my_peptide_lists.peptide_matches.size(); i++) {
	//	if (my_peptide_lists.peptide_matches[i].trypIndex != my_peptide_lists.peptide_matches[i - 1].trypIndex) {
	//		my_peptide_lists.peptide_matches_ftm.push_back(tmp1);
	//	}
	//	my_peptide_lists.peptide_matches_ftm.back().push_back(my_peptide_lists.peptide_matches[i]);
	//}

	///*for (int i = 0; i < my_peptide_lists.peptide_matches_ftm.size(); i++) {
	//	for (int j = 0; j < my_peptide_lists.peptide_matches_ftm[i].size(); j++) {
	//		cout << my_peptide_lists.peptide_matches_ftm[i].size() << "  " << my_peptide_lists.peptide_matches_ftm[i][j].trypIndex << "  " << my_peptide_lists.peptide_matches_ftm[i][j].missIndex << endl;
	//	}
	//}*/

	//
	//for (size_t i = 0; i < my_peptide_lists.peptide_matches_ftm.size(); i++) {
	//	/*cout << my_peptide_lists.peptide_matches_ftm[i].size() << endl;*/ 
	//	if (my_peptide_lists.peptide_matches_ftm[i].size() > 1) {
	//		float transient = 0;
	//		for (size_t j = 0; j < my_peptide_lists.peptide_matches_ftm[i].size(); j++) transient += my_peptide_lists.xic_mc_results[my_peptide_lists.peptide_matches_ftm[i][j].missIndex].areaXIC;
	//		float hold = my_peptide_lists.xic_ft_results[my_peptide_lists.peptide_matches_ftm[i][0].trypIndex].areaXIC / transient; 
	//		my_peptide_lists.peptide_matches_ftm[i].back().ratio = hold; 
	//	}
	//	if (my_peptide_lists.peptide_matches_ftm[i].size() == 1) {
	//		float hold1 = my_peptide_lists.xic_ft_results[my_peptide_lists.peptide_matches_ftm[i][0].trypIndex].areaXIC / my_peptide_lists.xic_mc_results[my_peptide_lists.peptide_matches_ftm[i][0].missIndex].areaXIC; 
	//		my_peptide_lists.peptide_matches_ftm[i][0].ratio = hold1; 
	//		hold1 = 0; 
	//	}
	//}


	//vector<float> rat;
	//float rat1 = 0;
	//int num = 0;
	//for (size_t i = 0; i < my_peptide_lists.peptide_matches_ftm.size(); i++) {
	//	for (size_t j = 0; j < my_peptide_lists.peptide_matches_ftm[i].size(); j++) {
	//		rat.push_back(my_peptide_lists.peptide_matches_ftm[i][j].ratio);
	//		rat1 += my_peptide_lists.peptide_matches_ftm[i][j].ratio;
	//		if (my_peptide_lists.peptide_matches_ftm[i][j].ratio != 0) num++;
	//	}
	//}

	//float rat2 = rat1 / num;  
	//my_metrics.intensity_final = rat2;

	//float var = 0;
	//for (size_t i = 0; i < rat.size(); ++i) {
	//	var += pow(rat[i] - rat2, 2);
	//}
	//var = var / rat.size();
	//float stdv = sqrt(var);
	//my_metrics.stdv_final = stdv;



	//sort(my_peptide_lists.peptide_matches.begin(), my_peptide_lists.peptide_matches.end(), compareMissIndex);

	///*for (int i = 0; i < my_peptide_lists.peptide_matches.size(); i++) {
	//	cout << my_peptide_lists.peptide_matches[i].missIndex << "  " << my_peptide_lists.peptide_matches[i].trypIndex << "  " << my_peptide_lists.peptide_matches[i].mc_pep_seq << "  " << my_peptide_lists.peptide_matches[i].ft_pep_seq << endl;
	//}
	//cout << my_peptide_lists.peptide_matches.size() << endl; */

	//my_peptide_lists.peptide_matches_mcm.push_back(tmp1);
	//my_peptide_lists.peptide_matches_mcm.back().push_back(my_peptide_lists.peptide_matches[0]);
	//for (size_t i = 1; i < my_peptide_lists.peptide_matches.size(); i++) {
	//	if (my_peptide_lists.peptide_matches[i].missIndex != my_peptide_lists.peptide_matches[i - 1].missIndex) {
	//		my_peptide_lists.peptide_matches_mcm.push_back(tmp1);
	//	}
	//	my_peptide_lists.peptide_matches_mcm.back().push_back(my_peptide_lists.peptide_matches[i]);
	//}

	//// sanity checks 
	//for (size_t i = 0; i < my_peptide_lists.peptide_matches_mcm.size(); i++) {
	//	for (size_t j = 0; j < my_peptide_lists.peptide_matches_mcm[i].size(); j++) {
	//		if (my_peptide_lists.xic_mc_results[my_peptide_lists.peptide_matches_mcm[i][j].missIndex].miss_cleaves == 2) {
	//			count2++;
	//		}
	//		if (my_peptide_lists.xic_mc_results[my_peptide_lists.peptide_matches_mcm[i][j].missIndex].miss_cleaves == 1) {
	//			count1++;
	//		}
	//	}
	//}

	//my_metrics.twice_mc = count2;
	//my_metrics.once_mc = count1;



	//
	//

	//return my_metrics;

}

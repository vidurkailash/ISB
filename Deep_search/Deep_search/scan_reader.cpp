#include "scan_reader.h"

using namespace std;
using namespace MSToolkit;

//MH: This function should be rewritten to first, make a list of all target peptides (whether miscleaved or not),
//then extract the precursor information into an array for each peptide.
match_lists scan_reader::mzml(peptide_lists& my_peptide_lists, my_parameters& my_params) {

	//READ MZML FILE FOR BOTH FULLY TRYPTIC AND MISCLEAVED PEPTIDES (MZML)
	MSReader myfile;
	Spectrum mySpec;

	myfile.setFilter(MS1);
	myfile.readFile(my_params.mzml.c_str(), mySpec);


  //MH: We need to get through precursor peak extraction fast.
  //To do so, we need to do it a) using a single pass through all arrays, and
  //b) using fast lookups instead of full iteration.
  //Step 1 is to sort my_peptide_lists.dlist
  sort(my_peptide_lists.d_list.begin(),my_peptide_lists.d_list.end(),compareRTime);
  size_t pepIndex=0; // we will start from the first sorted peptide

	match_lists my_match_lists;

	int b = 0;
	int c = 0;
	while (mySpec.getRTime() != 0) {

    //MH: Instead of iterating over all peptides, just iterate over the ones related to this spectrum
    //First we adjust our start point to the next peptide that falls within our RT window
    size_t i=pepIndex;
    while(i<my_peptide_lists.d_list.size() && (my_peptide_lists.d_list[i].rtime/60) < (mySpec.getRTime()-2) ) i++;
    if(i==my_peptide_lists.d_list.size()) break; //if we've checked every peptide, stop now.
    pepIndex=i; //mark our new start point for the next iteration

    //MH: From here, iterate over all peptides until we reach the end of our RT window
    while(i < my_peptide_lists.d_list.size() && (my_peptide_lists.d_list[i].rtime / 60) < (mySpec.getRTime() + 2) ){

      //MH: Now we have to check for the peak itself in this spectrum.
      //Instead of iterating over the whole spectrum, just do a binary search to
      //rapidly find the position in the array.
      int ret = findPeakMZ(mySpec, my_peptide_lists.d_list[i].mz,0.01);
      if(ret>-1) {
        my_match_lists.results.push_back(my_markers()); //input 
        my_match_lists.results.back().spec_rt = mySpec.getRTime(); //input
        my_match_lists.results.back().spec_sn = mySpec.getScanNumber(); //input
        my_match_lists.results.back().spec_size = mySpec.size(); //input
        my_match_lists.results.back().spec_mz = mySpec[ret].mz; //input
        my_match_lists.results.back().spec_intensity = mySpec[ret].intensity; //input
        my_match_lists.results.back().mzml_mz = my_peptide_lists.d_list[i].mz; //input   
        my_match_lists.results.back().mzml_rt = my_peptide_lists.d_list[i].rtime / 60; // input
        my_match_lists.results.back().pep_seq = my_peptide_lists.d_list[i].pep_seq; //input 
        my_match_lists.results.back().miss_cleaves = my_peptide_lists.d_list[i].miss_cleaves;
      }

      //MH: go to the next peptide
      i++;
    }

		myfile.readFile(NULL, mySpec);
	}

	cout << "pass 1 done" << "\n" << endl;
  
  //MH: When deleting duplicates, it helps to sort so that duplicates are next to each other
  //Then you need to make only a single pass. Furthermore, don't use erase, simply copy
  //the information you want to keep to a temporary vector, then copy it back to your
  //desired vector
  sort(my_match_lists.results.begin(), my_match_lists.results.end(),compareSeqScan);
  vector<my_markers> tmp;
  tmp.push_back(my_match_lists.results[0]); //First entry is always novel
  for(size_t i=1;i<my_match_lists.results.size();i++){
    //skip duplicates
    if( (my_match_lists.results[i].pep_seq == my_match_lists.results[i-1].pep_seq) && (my_match_lists.results[i].spec_sn == my_match_lists.results[i-1].spec_sn)) continue;
    //keep novel entries
    tmp.push_back(my_match_lists.results[i]);
  }
  my_match_lists.results=tmp; //copy the vector back


  //MH: Now resort the list for d_pep
  sort(my_peptide_lists.d_list.begin(), my_peptide_lists.d_list.end(), compareRTimeD);
  pepIndex = 0; // reset our index

  //MH: Restart reading our mzML file. Note, this could all be done in one pass, instead of two.
  myfile.readFile(my_params.mzml.c_str(), mySpec);
  while (mySpec.getRTime() != 0) {

    size_t i = pepIndex;
    while (i < my_peptide_lists.d_list.size() && (my_peptide_lists.d_list[i].d_pep_seq_rt / 60) < (mySpec.getRTime() - 2)) i++;
    if (i == my_peptide_lists.d_list.size()) break; //if we've checked every peptide, stop now.
    pepIndex = i; //mark our new start point for the next iteration

    while (i < my_peptide_lists.d_list.size() && (my_peptide_lists.d_list[i].d_pep_seq_rt / 60) < (mySpec.getRTime() + 2)) {
      int ret = findPeakMZ(mySpec, my_peptide_lists.d_list[i].d_pep_seq_mz, 0.01);
      if (ret > -1) {
        my_match_lists.results1.push_back(my_markers()); //input 
        my_match_lists.results1.back().spec_rt = mySpec.getRTime(); //input
        my_match_lists.results1.back().spec_sn = mySpec.getScanNumber(); //input
        my_match_lists.results1.back().spec_size = mySpec.size(); //input
        my_match_lists.results1.back().spec_mz = mySpec[ret].mz; //input
        my_match_lists.results1.back().spec_intensity = mySpec[ret].intensity; //input
        my_match_lists.results1.back().mzml_mz = my_peptide_lists.d_list[i].d_pep_seq_mz; //input   
        my_match_lists.results1.back().mzml_rt = my_peptide_lists.d_list[i].d_pep_seq_rt / 60; // input
        my_match_lists.results1.back().pep_seq = my_peptide_lists.d_list[i].d_pep_seq; //input 
        my_match_lists.results1.back().miss_cleaves = my_peptide_lists.d_list[i].d_miss_cleaves;
      }
      i++;//MH: go to the next peptide
    }
    myfile.readFile(NULL, mySpec);
  }

  cout << "pass 2 done" << "\n" << endl;

  //MH: Delete duplicates again.
  sort(my_match_lists.results1.begin(), my_match_lists.results1.end(), compareSeqScan);
  tmp.clear(); //clear our old tmp vector
  tmp.push_back(my_match_lists.results1[0]);
  for (size_t i = 1; i < my_match_lists.results1.size(); i++) {
    if ((my_match_lists.results1[i].pep_seq == my_match_lists.results1[i - 1].pep_seq) && (my_match_lists.results1[i].spec_sn == my_match_lists.results1[i - 1].spec_sn)) continue;
    tmp.push_back(my_match_lists.results1[i]);
  }
  my_match_lists.results1 = tmp;

	return my_match_lists; 

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
      //if(mid==0) break;
      upper = mid - 1;
      mid = (lower + upper) / 2;
    } else if(spec[mid].mz<LB) { //too low
      lower = mid + 1;
      mid = (lower + upper) / 2;
      //if(mid==sz) break;
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



	vector<vector<my_compare>> bigone;
	vector<vector<my_compare>> bigone1;
	vector<my_compare> end;
	vector<my_compare> end1; 
	vector<my_compare> test; 
	int testing = 0; 
	int iteration = 0;
	int iteration1 = 0;
	int count2 = 0;
	int count1 = 0;
	float count0 = 0;

	for (int i = 0; i < my_peptide_lists.final.size(); i++) {
		for (int j = 0; j < my_peptide_lists.final1.size(); j++) {
			size_t present = my_peptide_lists.final[i].seq.find(my_peptide_lists.final1[j].seq);
			if (present != string::npos && (my_peptide_lists.final[i].seq != my_peptide_lists.final1[j].seq)) {
				end.push_back(my_compare());
				end[iteration].miss_cleave_seq = my_peptide_lists.final[i].seq;
				end[iteration].tryp_seq = my_peptide_lists.final1[j].seq;
				end[iteration].mc_tot = my_peptide_lists.final[i].tot;
				end[iteration].tp_tot = my_peptide_lists.final1[j].tot;
				end[iteration].mc_mc = my_peptide_lists.final[i].mc;
				end[iteration].tp_mc = my_peptide_lists.final1[j].mc;
				iteration++;
			}
		}
	}

	end1 = end;


	for (int i = 0; i < end.size() - 1; i++) {
		bigone.push_back(test);
		bigone[testing].push_back(end[i]);
		for (int j = i + 1; j < end.size(); j++) {
			if (end[i].miss_cleave_seq == end[j].miss_cleave_seq) {
				bigone[testing].push_back(end[j]);
				end.erase(end.begin() + j);
				j = (i + 1) - 1;
			}
		}
		testing++;
	}

	bigone.push_back(test);
	bigone[testing].push_back(end[end.size() - 1]);


	// sanity checks 
	for (int i = 0; i < bigone.size(); i++) {
		for (int j = 0; j < bigone[i].size(); j++) {
			if (bigone[i][j].mc_mc == 2) {
				count2++;
			}
			if (bigone[i][j].mc_mc == 1) {
				count1++;
			}
			if (bigone[i][j].tp_mc == 0) {
				count0++;
			}
		}
		bigone[i][bigone[i].size() - 1].zero_frac = count0 / float(bigone[i].size());
		bigone[i][bigone[i].size() - 1].matches = (int)bigone[i].size();
		count0 = 0;
	}
	
  //vector<double> zero;
  double zero=0;
	for (int i = 0; i < bigone.size(); i++) {
		for (int j = 0; j < bigone[i].size(); j++) {
			zero+=bigone[i][j].zero_frac;
		}
	}
	double b = zero / bigone.size();
	my_metrics.zero = b; 
	my_metrics.twice_mc = count2; 
	my_metrics.once_mc = count1; 


	//useful metrics 
	vector<my_compare> temp;
	int temp_it = 0;
	for (int i = 0; i < end1.size() - 1; i++) {
		temp.push_back(my_compare());
		temp[temp_it].miss_cleave_seq = end1[i].miss_cleave_seq;
		temp[temp_it].tryp_seq = end1[i].tryp_seq;
		temp[temp_it].mc_tot = end1[i].mc_tot;
		temp[temp_it].tp_tot = end1[i].tp_tot;
		temp[temp_it].mc_mc = end1[i].mc_mc;
		temp[temp_it].tp_mc = end1[i].tp_mc;
		temp_it++;
		for (int j = i + 1; j < end1.size(); j++) {
			if (end1[i].tryp_seq == end1[j].tryp_seq) {
				temp.push_back(my_compare());
				temp[temp_it].miss_cleave_seq = end1[j].miss_cleave_seq;
				temp[temp_it].tryp_seq = end1[j].tryp_seq;
				temp[temp_it].mc_tot = end1[j].mc_tot;
				temp[temp_it].tp_tot = end1[j].tp_tot;
				temp[temp_it].mc_mc = end1[j].mc_mc;
				temp[temp_it].tp_mc = end1[j].tp_mc;
				temp_it++;
				end1.erase(end1.begin() + j); 
				j = (i + 1) - 1;
			}
		}
	}

	temp.push_back(my_compare());
	temp[temp_it].miss_cleave_seq = end[end.size() - 1].miss_cleave_seq;
	temp[temp_it].tryp_seq = end[end.size() - 1].tryp_seq;
	temp[temp_it].mc_tot = end[end.size() - 1].mc_tot;
	temp[temp_it].tp_tot = end[end.size() - 1].tp_tot;
	temp[temp_it].mc_mc = end[end.size() - 1].mc_mc;
	temp[temp_it].tp_mc = end[end.size() - 1].tp_mc;

	for (int i = 0; i < temp.size() - 1; i++) {
		bigone1.push_back(test);
		bigone1[iteration1].push_back(temp[i]);
		for (int j = i + 1; j < temp.size(); j++) {
			if (temp[i].tryp_seq == temp[j].tryp_seq) {
				bigone1[iteration1].push_back(temp[j]);
				temp.erase(temp.begin() + j);
				j = (i + 1) - 1; 
			}
		}
		iteration1++;
	}

	bigone1.push_back(test); 
	bigone1[iteration1].push_back(temp[temp.size() - 1]);

	//vector<float> transient; 
	for (int i = 0; i < bigone1.size(); i++) {
		if (bigone1[i].size() > 1 && bigone1[i][0].tp_mc == 0) {
			float transient=0;
			for (int j = 0; j < bigone1[i].size(); j++) transient+=bigone1[i][j].mc_tot;
			
			float hold1 = bigone1[i][0].tp_tot / transient;
			bigone1[i][bigone1[i].size() - 1].ratio = hold1;
		}
		if (bigone1[i].size() == 1 && bigone1[i][0].tp_mc == 0) {
			float hold2 = bigone1[i][0].tp_tot / bigone1[i][0].mc_tot;
			bigone1[i][0].ratio = hold2; 
			hold2 = 0; 
		}
	}

	vector<double> rat; 
	double rat1=0;
	int asdf = 0;
	for (int i = 0; i < bigone1.size(); i++) {
		for (int j = 0; j < bigone1[i].size(); j++) {
			rat.push_back(bigone1[i][j].ratio);
			rat1+= bigone1[i][j].ratio;
			if(bigone1[i][j].ratio!=0) asdf++;
		}
	}
	
	double rat2 = rat1 / asdf;  //MH: I'm not sure this is being calculated correctly
	my_metrics.intensity_final = rat2;
 
	double var = 0;  
	for (int i = 0; i < rat.size(); ++i) { 
		var += pow(rat[i] - rat2, 2); 
	}
	var = var / rat.size(); 
	double stdv = sqrt(var); 
	my_metrics.stdv_final = stdv; 

	vector<my_compare> last_one;
	vector<my_compare> last_one1;
	for (size_t i = 0; i < bigone1.size(); i++) {
		for (size_t j = 0; j < bigone1[i].size(); j++) {
			last_one.push_back(bigone1[i][j]);
		}
	}
	sort(last_one.begin(), last_one.end(), compareInten); 
	last_one1.push_back(last_one[0]); 
	for (size_t i = 1; i < last_one.size(); i++) {
		if (last_one[i].tryp_seq == last_one[i - 1].tryp_seq) continue;
		last_one1.push_back(last_one[i]); 
	}
	last_one = last_one1; 
	for (size_t i = 0; i < 20; i++) {
		cout << last_one[i].miss_cleave_seq << "  " << last_one[i].mc_tot << "  " << last_one[i].tryp_seq << "  " << last_one[i].tp_tot << "\t" << last_one[i].tp_tot/last_one[i].mc_tot << "\tXX" << endl; 
	}



	//PRINT FUNCTIONS 
	/*int meme = 0; 
	for (int i = 0; i < bigone.size(); i++) {
		for (int j = 0; j < bigone[i].size(); j++) {
			cout << bigone[i][j].miss_cleave_seq << "   " << bigone[i][j].mc_tot << "  " << bigone[i][j].mc_mc << "   " << bigone[i][j].tryp_seq << "   " << bigone[i][j].tp_tot << "   " << bigone[i][j].tp_mc << "    " << bigone[i][j].matches << "    " <<  bigone[i][j].zero_frac << endl; 
			meme++; 
		}
	}
	

	int meme1 = 0;
	for (int i = 0; i < bigone1.size(); i++) {
		for (int j = 0; j < bigone1[i].size(); j++) {
			cout << bigone1[i][j].miss_cleave_seq << "   " << bigone1[i][j].mc_tot << "  " << bigone1[i][j].mc_mc << "   " << bigone1[i][j].tryp_seq << "   " << bigone1[i][j].tp_tot << "   " << bigone1[i][j].tp_mc << "    " << bigone1[i][j].ratio << endl;
			meme1++;
		}
	}*/
	/*for (int i = 0; i < bigone1.size(); i++) {
		cout << bigone1[i].size() << endl; 
	}*/
	


	
	/*for (int i = 0; i < temp.size(); i++) {
		cout << temp[i].miss_cleave_seq << "   " << temp[i].mc_tot << "   " << temp[i].mc_mc << "   " << temp[i].tryp_seq << "   " << temp[i].tp_tot << "   " << temp[i].tp_mc << endl;
	}
	cout << temp.size() << endl; 

	for (int i = 0; i < temp.size(); i++) {
		if (temp[i].tryp_seq == "QFEMEELILELAAQVLEDK") {
			cout << "gotem ;) " << endl; 
		}
	}*/
	/*for (int i = 0; i < end.size(); i++) {
		cout << end[i].miss_cleave_seq << "   " << end[i].mc_tot << "   " << end[i].mc_mc << "   " << end[i].tryp_seq << "   " << end[i].tp_tot << "   " << end[i].tp_mc << endl;
	}
	cout << end.size() << endl;*/

	//for (int i = 0; i < my_peptide_lists.master[649].size(); i++) {
	//	/*for (int j = 0; j < my_peptide_lists.master[i].size(); j++) {
	//		cout << my_peptide_lists.master[i][j].seq << "    " << my_peptide_lists.master[i][j].x << "    " << my_peptide_lists.master[i][j].y << "   " << my_peptide_lists.master[i][j].tot << endl;
	//	}*/
	//	cout << my_peptide_lists.master[649][i].seq << "    " << my_peptide_lists.master[649][i].x << "    " << my_peptide_lists.master[649][i].y << "    " << my_peptide_lists.master[649][i].tot << endl;
	//}

	//cout << my_peptide_lists.master.size() << endl;


	//for (int i = 0; i < my_peptide_lists.master1[588].size(); i++) {
	//	/*for (int j = 0; j < my_peptide_lists.master1[i].size(); j++) {
	//		cout << my_peptide_lists.master1[i][j].seq << "    " << my_peptide_lists.master1[i][j].x << "    " << my_peptide_lists.master1[i][j].y << "   " << my_peptide_lists.master1[i][j].tot << endl;
	//	}*/
	//	cout << my_peptide_lists.master1[588][i].seq << "    " << my_peptide_lists.master1[588][i].x << "    " << my_peptide_lists.master1[588][i].y << "     " << my_peptide_lists.master1[588][i].tot << endl;
	//}

	//cout << my_peptide_lists.master1.size() << endl;





	/*cout << my_peptide_lists.final[649].seq << "   " << my_peptide_lists.final[649].tot << "   " << my_peptide_lists.final[649].mc << endl;

	cout << my_peptide_lists.final1[588].seq << "   " << my_peptide_lists.final1[588].tot << "    " << my_peptide_lists.final1[588].mc << endl;*/



	return my_metrics; 
}
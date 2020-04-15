#include "rapidxml.hpp"
#include "deep_class.h"

using namespace rapidxml;
using namespace std;

peptide_lists deep_functions::xml_parse(my_parameters& my_params) {
	xml_document<> doc;

	// Read the xml file into a vector
	ifstream theFile(my_params.filename);
	vector<char> buffer((istreambuf_iterator<char>(theFile)), istreambuf_iterator<char>());
	buffer.push_back('\0');
	// Parse the buffer using the xml file parsing library into doc 
	doc.parse<0>(&buffer[0]);
	// Find root node
	xml_node<>* root_node = doc.first_node("msms_pipeline_analysis");
	xml_node<>* spectrum_node = root_node->first_node("msms_run_summary");

	peptide_lists my_peptide_lists;

	int c = 0;
	int d = 0;

	for (spectrum_node; spectrum_node; spectrum_node = spectrum_node->next_sibling())
	{
		if (string(spectrum_node->value()) == "msms_run_summary") {
			continue;
		}

		switch (my_params.ident) {
		//IP
		case 'a':
			for (xml_node<>* sample_node = spectrum_node->first_node("spectrum_query"); sample_node; sample_node = sample_node->next_sibling())
			{
				for (xml_node<>* secondary_node = sample_node->first_node("search_result")->first_node("search_hit")->first_node("analysis_result"); secondary_node; secondary_node = secondary_node->next_sibling())
				{
					if (/*bool(secondary_node->first_attribute("analysis")) == 1 &&*/ string(secondary_node->first_attribute("analysis")->value()) == "interprophet") {
						if (atof(secondary_node->first_node("interprophet_result")->first_attribute("probability")->value()) >= 0.00) {
							d++;
						}
						if (atof(secondary_node->first_node("interprophet_result")->first_attribute("probability")->value()) >= my_params.ipro_prob) {
							my_peptide_lists.all.push_back(my_features());
							my_peptide_lists.all[c].pep_seq = string(sample_node->first_node("search_result")->first_node("search_hit")->first_attribute("peptide")->value());
							my_peptide_lists.all[c].charge = atoi(sample_node->first_attribute("assumed_charge")->value());
							my_peptide_lists.all[c].mass = atof(sample_node->first_attribute("precursor_neutral_mass")->value());
							my_peptide_lists.all[c].rtime = (float)atof(sample_node->first_attribute("retention_time_sec")->value());
							my_peptide_lists.all[c].prev_aa = string(sample_node->first_node("search_result")->first_node("search_hit")->first_attribute("peptide_prev_aa")->value());
							my_peptide_lists.all[c].next_aa = string(sample_node->first_node("search_result")->first_node("search_hit")->first_attribute("peptide_next_aa")->value());
							c++;
						}
					}
				}
			}
			break; 
		//PP
		case 'b':
			for (xml_node<>* sample_node = spectrum_node->first_node("spectrum_query"); sample_node; sample_node = sample_node->next_sibling())
			{
				if (atof(sample_node->first_node("search_result")->first_node("search_hit")->first_node("analysis_result")->first_node("peptideprophet_result")->first_attribute("probability")->value()) >= 0.00) {
					d++;
				}
				if (atof(sample_node->first_node("search_result")->first_node("search_hit")->first_node("analysis_result")->first_node("peptideprophet_result")->first_attribute("probability")->value()) >= my_params.pep_prob) {
					my_peptide_lists.all.push_back(my_features());
					my_peptide_lists.all[c].pep_seq = string(sample_node->first_node("search_result")->first_node("search_hit")->first_attribute("peptide")->value());
					my_peptide_lists.all[c].charge = atoi(sample_node->first_attribute("assumed_charge")->value());
					my_peptide_lists.all[c].mass = atof(sample_node->first_attribute("precursor_neutral_mass")->value());
					my_peptide_lists.all[c].rtime = (float)atof(sample_node->first_attribute("retention_time_sec")->value());
					my_peptide_lists.all[c].prev_aa = string(sample_node->first_node("search_result")->first_node("search_hit")->first_attribute("peptide_prev_aa")->value());
					my_peptide_lists.all[c].next_aa = string(sample_node->first_node("search_result")->first_node("search_hit")->first_attribute("peptide_next_aa")->value());
					c++;
				}
			}
			break;
		default:
			cout << "error" << endl; 
			exit(1);
			break; 
		}

	}
	my_peptide_lists.total = d;
                                  
	return my_peptide_lists;

}

bool deep_functions::tryptic_calc(peptide_lists& my_peptide_lists) {

	int c = 0;
	int d = 0;

	for (int i = 0; i < my_peptide_lists.all.size(); i++) {

		if ((my_peptide_lists.all[i].prev_aa == "R" || my_peptide_lists.all[i].prev_aa == "K" || my_peptide_lists.all[i].prev_aa == "-") &&
			(my_peptide_lists.all[i].pep_seq.back() == 'R' || my_peptide_lists.all[i].pep_seq.back() == 'K' || my_peptide_lists.all[i].next_aa == "-")) {

			my_peptide_lists.tryptic.push_back(my_features());
			my_peptide_lists.tryptic[c].pep_seq = my_peptide_lists.all[i].pep_seq;
			my_peptide_lists.tryptic[c].charge = my_peptide_lists.all[i].charge;
			my_peptide_lists.tryptic[c].mass = my_peptide_lists.all[i].mass;
			my_peptide_lists.tryptic[c].prev_aa = my_peptide_lists.all[i].prev_aa;
			my_peptide_lists.tryptic[c].rtime = my_peptide_lists.all[i].rtime;
			my_peptide_lists.tryptic[c].next_aa = my_peptide_lists.all[i].next_aa;

			c++;
		}
		else {

			my_peptide_lists.non_tryptic.push_back(my_features());
			my_peptide_lists.non_tryptic[d].pep_seq = my_peptide_lists.all[i].pep_seq;
			my_peptide_lists.non_tryptic[d].charge = my_peptide_lists.all[i].charge;
			my_peptide_lists.non_tryptic[d].mass = my_peptide_lists.all[i].mass;
			my_peptide_lists.non_tryptic[d].prev_aa = my_peptide_lists.all[i].prev_aa;
			my_peptide_lists.non_tryptic[d].rtime = my_peptide_lists.all[i].rtime;
			my_peptide_lists.non_tryptic[d].next_aa = my_peptide_lists.all[i].next_aa;

			d++;
		}
		
	}

	return true;

}

bool deep_functions::miss_cleave(peptide_lists& my_peptide_lists) {

	int c = 0;
	for (int i = 0; i < my_peptide_lists.tryptic.size(); i++) {
		for (int j = 0; j < my_peptide_lists.tryptic[i].pep_seq.size() - 1; j++) {
			int k = j + 1;
			if ((my_peptide_lists.tryptic[i].pep_seq[j] == 'K' || my_peptide_lists.tryptic[i].pep_seq[j] == 'R') && my_peptide_lists.tryptic[i].pep_seq[k] != 'P') {
				c++;
			}
		}
		my_peptide_lists.tryptic[i].miss_cleaves = c;
		c = 0;
	}

	return true;

}

bool deep_functions::delete_dup(peptide_lists& my_peptide_lists) {

  //MH: I think this function is supposed to find the unique set of peptide sequences
  //from the complete list of tryptic PSMs. delete_dup is a misnomer. Nothing is deleted.
  //Instead new arrays of unique sequences are created.
  size_t i;

  //MH: Sort the tryptic array by sequence and then charge state
  sort(my_peptide_lists.tryptic.begin(),my_peptide_lists.tryptic.end(),compareSeqZ);

  //MH: Iterate over all tryptic PSMs, keeping each unique instance of sequence and charge
  my_peptide_lists.tryp_unique_z.push_back(my_peptide_lists.tryptic[0]);
  for(i=1;i<my_peptide_lists.tryptic.size();i++){
    if(my_peptide_lists.tryp_unique_z.back().pep_seq.compare(my_peptide_lists.tryptic[i].pep_seq)==0 &&
      my_peptide_lists.tryp_unique_z.back().charge == my_peptide_lists.tryptic[i].charge) continue;
    my_peptide_lists.tryp_unique_z.push_back(my_peptide_lists.tryptic[i]);
  }

  //MH: Repeat the process, keeping each unique instance of sequence;
  my_peptide_lists.tryp_unique.push_back(my_peptide_lists.tryptic[0]);
  for (i = 1; i < my_peptide_lists.tryptic.size(); i++) {
    if (my_peptide_lists.tryp_unique.back().pep_seq.compare(my_peptide_lists.tryptic[i].pep_seq) == 0) continue;
    my_peptide_lists.tryp_unique.push_back(my_peptide_lists.tryptic[i]);
  }

  //MH: Now of the unique peptide sequences, make a subset of miscleaved peptides
  i=0;
  while(my_peptide_lists.tryp_unique[i].miss_cleaves==0) i++; //no boundary checks here
  my_peptide_lists.miss_unique.push_back(my_peptide_lists.tryp_unique[i]);
  for (i = i+1; i < my_peptide_lists.tryp_unique.size(); i++) {
    if (my_peptide_lists.tryp_unique[i].miss_cleaves==0) continue;
    my_peptide_lists.miss_unique.push_back(my_peptide_lists.tryp_unique[i]);
  }

	return true;

}

bool deep_functions::lcd(peptide_lists& my_peptide_lists) {
	
  
  //MH: This new method finds all miscleavages for each fully cleaved peptide, and stores them in a
  //new array.
  size_t i,j;
  for(i=0;i<my_peptide_lists.tryp_unique.size();i++){

    //MH: skip any peptides that have miscleavages
    if(my_peptide_lists.tryp_unique[i].miss_cleaves>0) continue;

    for(j=0;j<my_peptide_lists.miss_unique.size();j++){

      //MH: we have a match, make a new entry
      size_t found = my_peptide_lists.miss_unique[j].pep_seq.find(my_peptide_lists.tryp_unique[i].pep_seq);
      if(found != string::npos && my_peptide_lists.miss_unique[j].pep_seq.compare(my_peptide_lists.tryp_unique[i].pep_seq)!=0) {
        my_features f= my_peptide_lists.miss_unique[j];
        f.d_pep_seq = my_peptide_lists.tryp_unique[i].pep_seq;
        f.d_pep_seq_rt = my_peptide_lists.tryp_unique[i].rtime;
        f.d_miss_cleaves = my_peptide_lists.tryp_unique[i].miss_cleaves;
        f.d_pep_seq_mass = my_peptide_lists.tryp_unique[i].mass;
        f.d_pep_seq_charge = my_peptide_lists.tryp_unique[i].charge;
        if (found == 0) {
          f.cleave_loc = 'R';
          f.cleave_pos = 0;
        } else {
          f.cleave_loc = 'L';
          f.cleave_pos = (int)found;
        }
        f.mz = (f.mass + (f.charge * 1.00727)) / f.charge;
        f.d_pep_seq_mz = (f.d_pep_seq_mass + (f.d_pep_seq_charge * 1.00727)) / f.d_pep_seq_charge;
        my_peptide_lists.d_list.push_back(f);
      }

    }

  }

  //MH: I think this approach is flawed. MISCLVDRSEQENCEK is a miscleavage of both MISCLVDR and SEQENCEK
  //Therefore it isn't possible to have only a single tryptic representation for each miscleavage.
	///*for (int i = 0; i < my_peptide_lists.miss_unique.size(); i++) {
	//	for (int j = 0; j < my_peptide_lists.tryp_unique.size(); j++) {
	//		size_t found = my_peptide_lists.miss_unique[i].pep_seq.find(my_peptide_lists.tryp_unique[j].pep_seq);
	//		if (found != string::npos && (my_peptide_lists.miss_unique[i].pep_seq != my_peptide_lists.tryp_unique[j].pep_seq)) {
	//			my_peptide_lists.miss_unique[i].d_pep_seq = my_peptide_lists.tryp_unique[j].pep_seq;
	//			my_peptide_lists.miss_unique[i].d_pep_seq_rt = my_peptide_lists.tryp_unique[j].rtime;
	//			my_peptide_lists.miss_unique[i].d_miss_cleaves = my_peptide_lists.tryp_unique[j].miss_cleaves;
	//			if (found == 0) {
	//				my_peptide_lists.miss_unique[i].cleave_loc = 'R';
	//				my_peptide_lists.miss_unique[i].cleave_pos = (int)found;

	//			}
	//			else {
	//				my_peptide_lists.miss_unique[i].cleave_loc = 'L';
	//				my_peptide_lists.miss_unique[i].cleave_pos = (int)found;

	//			}
	//			my_peptide_lists.miss_unique[i].d_pep_seq_mass = my_peptide_lists.tryp_unique[j].mass;
	//			my_peptide_lists.miss_unique[i].d_pep_seq_charge = my_peptide_lists.tryp_unique[j].charge;
	//		}

	//	}

	//	my_peptide_lists.miss_unique[i].mz = (my_peptide_lists.miss_unique[i].mass + ((my_peptide_lists.miss_unique[i].charge) * 1.00727)) / my_peptide_lists.miss_unique[i].charge;
	//	my_peptide_lists.miss_unique[i].d_pep_seq_mz = (my_peptide_lists.miss_unique[i].d_pep_seq_mass + ((my_peptide_lists.miss_unique[i].d_pep_seq_charge) * 1.00727)) / my_peptide_lists.miss_unique[i].d_pep_seq_charge;
	//}*/




	/*for (int i = 0; i < my_peptide_lists.tryp_unique.size(); i++) {
		cout << my_peptide_lists.tryp_unique[i].pep_seq << "  " << my_peptide_lists.tryp_unique[i].charge << "   " << my_peptide_lists.tryp_unique[i].mass << "  " << my_peptide_lists.tryp_unique[i].mz << "  " << my_peptide_lists.tryp_unique[i].cleave_loc << "   " << my_peptide_lists.tryp_unique[i].cleave_pos  << "   " << my_peptide_lists.tryp_unique[i].d_pep_seq << "  " << my_peptide_lists.tryp_unique[i].d_pep_seq_charge << "  " << my_peptide_lists.tryp_unique[i].d_pep_seq_mass << "  " << my_peptide_lists.tryp_unique[i].d_pep_seq_mz << endl;
	}*/

	/*cout << my_peptide_lists.tryp_unique.size() << endl;*/



	return true; 

}

peptide_lists deep_functions::new_list(peptide_lists& my_peptide_lists) {

	my_peptide_lists.d_list = my_peptide_lists.miss_unique; 
	
	for (int i = 0; i < my_peptide_lists.d_list.size(); i++) {
		if (my_peptide_lists.d_list[i].d_pep_seq.size() == 0 ){
			my_peptide_lists.d_list.erase(my_peptide_lists.d_list.begin() + i);
			i = i - 1;
		}

	}



	/*for (int i = 0; i < my_peptide_lists.d_list.size(); i++) {
		cout << my_peptide_lists.d_list[i].pep_seq << "  " << my_peptide_lists.d_list[i].charge << "   " << my_peptide_lists.d_list[i].mass << "  " << my_peptide_lists.d_list[i].mz << "  " << my_peptide_lists.d_list[i].cleave_loc << "   " << my_peptide_lists.d_list[i].cleave_pos << " ||||| " << my_peptide_lists.d_list[i].d_pep_seq << "  " << my_peptide_lists.d_list[i].d_pep_seq_charge << "  " << my_peptide_lists.d_list[i].d_pep_seq_mass << "  " << my_peptide_lists.d_list[i].d_pep_seq_mz << endl;
	}*/


	cout << my_peptide_lists.d_list.size() << endl;

	return my_peptide_lists;

}



metrics deep_functions::calc(peptide_lists& my_peptide_lists, metrics& my_metrics) {


	vector<int> miss_psm;
	int p = 0;
	for (int i = 0; i < my_peptide_lists.tryptic.size(); i++) {
		miss_psm.push_back(my_peptide_lists.tryptic[i].miss_cleaves);
		if (my_peptide_lists.tryptic[i].miss_cleaves > 0) {
			p++;
		}
	}
	double r = accumulate(miss_psm.begin(), miss_psm.end(), 0);


	vector<int> number;
	vector<int> miss_pep;
	double h = 0;
	for (size_t i = 0; i < my_peptide_lists.tryp_unique.size(); i++) {
		number.push_back((int)my_peptide_lists.tryp_unique[i].pep_seq.size());
		miss_pep.push_back(my_peptide_lists.tryp_unique[i].miss_cleaves);
		if (my_peptide_lists.tryp_unique[i].miss_cleaves > 0) {
			h++;
		}
	}
	double f = accumulate(number.begin(), number.end(), 0);
	double g = accumulate(miss_pep.begin(), miss_pep.end(), 0);



	my_metrics.total_psm = my_peptide_lists.total;
	my_metrics.psm_num = my_peptide_lists.all.size();
	my_metrics.tryptic_num = my_peptide_lists.tryptic.size();
	my_metrics.nontryptic_num = my_peptide_lists.non_tryptic.size();
	my_metrics.unique_pep_charge = my_peptide_lists.tryp_unique_z.size();
	my_metrics.unique_peptides = my_peptide_lists.tryp_unique.size();
	my_metrics.avg_pep_length = f / my_metrics.unique_peptides;
	my_metrics.tryp_frac = (double)my_metrics.tryptic_num / my_metrics.psm_num;
	my_metrics.nontryp_frac = (double)my_metrics.nontryptic_num / my_metrics.psm_num;
	my_metrics.pep_frac = (double)my_peptide_lists.tryp_unique.size() / my_metrics.psm_num;
	my_metrics.miss_cleave_rate_psm = r / my_metrics.tryptic_num;
	my_metrics.miss_cleave_rate_pep = g / my_metrics.unique_peptides;
	my_metrics.num_miss_cleave_pep = h; 
	my_metrics.num_miss_cleave_total = g; 
	my_metrics.miss_cleave_avg = g / h;
	my_metrics.num_miss_cleave_psm = p; 
	my_metrics.num_tot_miss_cleave_psm = r; 
	my_metrics.golden_stat_psm = my_metrics.num_miss_cleave_psm / my_metrics.tryptic_num;
	my_metrics.golden_stat_unique_pep = my_metrics.num_miss_cleave_pep / my_metrics.unique_peptides;
	my_metrics.find_num = my_peptide_lists.d_list.size(); 

	return my_metrics;

}


//MH: Instead of all of this complication, just find the boundaries of the array that you want to keep
//And you don't need to copy the data to a new array whenever you iterate over it. Just iterate over
//the exisiting array.
bool deep_functions::cleanNoise(std::vector<my_intensities>& v) {
  size_t i;
  size_t a,b; //MH: these will be your array boundaries

  //MH: apex information
  float max=0;
  size_t maxIndex=0;

  //MH: First step, find the maximum intensity (Assumed to be the peak apex)
  for(i=0;i<v.size();i++){
    if(v[i].y>max){
      max=v[i].y;
      maxIndex=i;
    }
  }

  //MH: set out threshold
  max/=10;

  //MH: Now step to the left of the apex to find your lower bound
  a=maxIndex;
  while(a>0){
    if(v[a-1].y>max) a--;
    else break;
  }

  //MH: Now step to the right of the apex to find your upper bound
  b=maxIndex;
  while(b<v.size()-1){
    if(v[b+1].y>max) b++;
    else break;
  }

  //MH: With bounds set, create a subset array
  vector<my_intensities> tmp;
  for(i=a;i<=b;i++){
    tmp.push_back(v[i]);
  }

  //MH: overwrite our original bloated array with the trimmed subset array
  v=tmp;

  //MH: report success
  return true;


	//vector<float> querry;
	//vector<float> tmp;
	// 
	//for (size_t i = 0; i < v.size(); i++) {
	//	querry.push_back(v[i].y);
	//}
	//int index = max_element(querry.begin(), querry.end()) - querry.begin();
	//float val = *max_element(querry.begin(), querry.end());
	//int i = index; 
	//
	//while (i > 0) {
	//	if (querry[i] > 0.1 * val) {
	//		tmp.push_back(querry[i]);
	//		i--;
	//	}
	//	else { break; }
	//}
	//i = index;
	//while (i < querry.size()) {
	//	if (querry[i] > 0.1 * val) {
	//		tmp.push_back(querry[i]);
	//		i++;
	//	}
	//	else { break; }
	//}

	//querry = tmp; 
	//tmp.clear(); 


	//return querry;
}

bool deep_functions::reader(peptide_lists& my_peptide_lists, metrics& my_metrics, match_lists& my_match_lists) {

  vector<my_intensities> qc;

  //MH: Again, sorting is your friend here. Although these vectors are already sorted, so we won't do that.
  //But because they are sorted, you can get through in a single pass. And you won't need to erase
  //which will save tons of time.
  for (size_t i = 0; i < my_match_lists.results.size(); i++) {
    my_peptide_lists.spectra.push_back(my_intensities()); //add first one
    my_peptide_lists.spectra.back().x = my_match_lists.results[i].spec_rt;
    my_peptide_lists.spectra.back().y = my_match_lists.results[i].spec_intensity;
    my_peptide_lists.spectra.back().seq = my_match_lists.results[i].pep_seq;
    my_peptide_lists.spectra.back().mc = my_match_lists.results[i].miss_cleaves;
  }
  my_peptide_lists.master.push_back(qc);
  my_peptide_lists.master.back().push_back(my_peptide_lists.spectra[0]);
  for (size_t i = 1; i < my_peptide_lists.spectra.size(); i++) {
    if (my_peptide_lists.spectra[i].seq != my_peptide_lists.spectra[i - 1].seq) {
      my_peptide_lists.master.push_back(qc);
    }
    my_peptide_lists.master.back().push_back(my_peptide_lists.spectra[i]);
  }

  //MH: Repeat for the other list
  for (size_t i = 0; i < my_match_lists.results1.size(); i++) {
    my_peptide_lists.spectra1.push_back(my_intensities()); //add first one
    my_peptide_lists.spectra1.back().x = my_match_lists.results1[i].spec_rt;
    my_peptide_lists.spectra1.back().y = my_match_lists.results1[i].spec_intensity;
    my_peptide_lists.spectra1.back().seq = my_match_lists.results1[i].pep_seq;
    my_peptide_lists.spectra1.back().mc = my_match_lists.results1[i].miss_cleaves;
  }
  my_peptide_lists.master1.push_back(qc);
  my_peptide_lists.master1.back().push_back(my_peptide_lists.spectra1[0]);
  for (size_t i = 1; i < my_peptide_lists.spectra1.size(); i++) {
    if (my_peptide_lists.spectra1[i].seq != my_peptide_lists.spectra1[i - 1].seq) {
      my_peptide_lists.master1.push_back(qc);
    }
    my_peptide_lists.master1.back().push_back(my_peptide_lists.spectra1[i]);
  }


  //DELETE NOISE (INTENSITIES LESS THAN 10% OF THE MAX) (NOISE)
  //MH: The faster, simpler approach
  for(size_t i=0;i<my_peptide_lists.master.size();i++){
    if(!cleanNoise(my_peptide_lists.master[i])){
      //handle error
    }
  }
  for (size_t i = 0; i < my_peptide_lists.master1.size(); i++) {
    if (!cleanNoise(my_peptide_lists.master1[i])) {
      //handle error
    }
  }

  ///*for (size_t i = 0; i < my_peptide_lists.master.size(); i++) {
  //  vector<float> rn;
  //  rn = cleanNoise(my_peptide_lists.master[i]);

  //  for (size_t n = 0; n < my_peptide_lists.master[i].size() - 1; n++) {
  //    for (size_t p = n + 1; p < my_peptide_lists.master[i].size(); p++) {
  //      vector<float>::iterator it;
  //      it = find(rn.begin(), rn.end(), my_peptide_lists.master[i][p].y);
  //      if (it == rn.end()) {
  //        my_peptide_lists.master[i].erase(my_peptide_lists.master[i].begin() + p);
  //        p = (n + 1) - 1;
  //      }
  //    }
  //  }
  //  for (size_t n = 0; n < my_peptide_lists.master[i].size(); n++) {
  //    vector<float>::iterator it;
  //    it = find(rn.begin(), rn.end(), my_peptide_lists.master[i][n].y);
  //    if (it == rn.end()) {
  //      my_peptide_lists.master[i].erase(my_peptide_lists.master[i].begin() + n);
  //    }
  //  }
  //  rn.clear();
  //}

  //for (size_t i = 0; i < my_peptide_lists.master1.size(); i++) {
  //  vector<float> rn1;
  //  rn1 = cleanNoise(my_peptide_lists.master1[i]);

  //  for (size_t n = 0; n < my_peptide_lists.master1[i].size() - 1; n++) {
  //    for (size_t p = n + 1; p < my_peptide_lists.master1[i].size(); p++) {
  //      vector<float>::iterator it;
  //      it = find(rn1.begin(), rn1.end(), my_peptide_lists.master1[i][p].y);
  //      if (it == rn1.end()) {
  //        my_peptide_lists.master1[i].erase(my_peptide_lists.master1[i].begin() + p);
  //        p = (n + 1) - 1;
  //      }
  //    }
  //  }
  //  for (size_t n = 0; n < my_peptide_lists.master1[i].size(); n++) {
  //    vector<float>::iterator it;
  //    it = find(rn1.begin(), rn1.end(), my_peptide_lists.master1[i][n].y);
  //    if (it == rn1.end()) {
  //      my_peptide_lists.master1[i].erase(my_peptide_lists.master1[i].begin() + n);
  //    }
  //  }
  //  rn1.clear();
  //}*/

  // END NOISE FUNCTION 


  //delete peptides with only one data point
  vector<vector<my_intensities>> filter;
  for (size_t i = 0; i < my_peptide_lists.master.size(); i++) {
    if (my_peptide_lists.master[i].size() > 1) {
      filter.push_back(my_peptide_lists.master[i]);
    }
  }
  my_peptide_lists.master = filter;
  filter.clear();
  for (size_t i = 0; i < my_peptide_lists.master1.size(); i++) {
    if (my_peptide_lists.master1[i].size() > 1) {
      filter.push_back(my_peptide_lists.master1[i]);
    }
  }
  my_peptide_lists.master1 = filter;
  filter.clear();



  cout << "noise deleted" << "\n" << endl;


  //metric calculations
  vector<float> peak;
  float avg=0;
  //MH: this array isn't necessary
  //vector<float> peak2; 
  for (int i = 0; i < my_peptide_lists.master.size(); i++) {
    for (int j = 0; j < my_peptide_lists.master[i].size(); j++) {
      peak.push_back(my_peptide_lists.master[i][j].y);
    }
    avg+= *max_element(peak.begin(), peak.end());
    //peak2.push_back(*max_element(peak.begin(), peak.end()));
    peak.clear();
  }


  //float sum = 0; for (int i = 0; i < peak2.size(); ++i) { sum += peak2[i]; }
  //float sumi = sum / peak2.size();
  my_metrics.miss_avg_high = avg/ my_peptide_lists.master.size();

  peak.clear();
  avg=0;
  //vector<float> peak1;
  //vector<float> peak3;
  for (int i = 0; i < my_peptide_lists.master1.size(); i++) {
    for (int j = 0; j < my_peptide_lists.master1[i].size(); j++) {
      peak.push_back(my_peptide_lists.master1[i][j].y);
    }
    avg += *max_element(peak.begin(), peak.end());
    //peak3.push_back(*max_element(peak1.begin(), peak1.end()));
    peak.clear();
  }


  //float sum1 = 0; for (int i = 0; i < peak3.size(); ++i) { sum1 += peak3[i]; }
  //float sumi1 = sum1 / peak3.size();
  my_metrics.tryp_avg_high = avg/ my_peptide_lists.master1.size();

  //MH: Calculating area is something performed on multiple arrays. So create a function
  // that does it on any array you pass to it. Also, have it return the total rather than
  // iterate the array again to count the total.
  for (size_t i = 0; i < my_peptide_lists.master.size(); i++) {
    my_peptide_lists.master[i].back().tot = calcPeakArea(my_peptide_lists.master[i]);
  }
  for (size_t i = 0; i < my_peptide_lists.master1.size(); i++) {
    my_peptide_lists.master1[i].back().tot = calcPeakArea(my_peptide_lists.master1[i]);
  }


  vector<my_intensities> final;
  int cycle = 0;
  for (int i = 0; i < my_peptide_lists.master.size(); i++) {
    final.push_back(my_intensities());
    final[cycle].seq = my_peptide_lists.master[i][my_peptide_lists.master[i].size() - 1].seq;
    final[cycle].tot = my_peptide_lists.master[i][my_peptide_lists.master[i].size() - 1].tot;
    final[cycle].mc = my_peptide_lists.master[i][my_peptide_lists.master[i].size() - 1].mc;
    cycle++;
  }

  vector<my_intensities> final1;
  int cycle1 = 0;
  for (int i = 0; i < my_peptide_lists.master1.size(); i++) {
    final1.push_back(my_intensities());
    final1[cycle1].seq = my_peptide_lists.master1[i][my_peptide_lists.master1[i].size() - 1].seq;
    final1[cycle1].tot = my_peptide_lists.master1[i][my_peptide_lists.master1[i].size() - 1].tot;
    final1[cycle1].mc = my_peptide_lists.master1[i][my_peptide_lists.master1[i].size() - 1].mc;
    cycle1++;
  }

  my_peptide_lists.final = final;
  my_peptide_lists.final1 = final1;

  // END AREA FUNCTION


  cout << "intensity calculated" << "\n" << endl;






  return true;

}


//MH: This function replicates computing (and storing) the area between every two points.
// It then returns the sum of the area of the peak. Note that the boundaries are assumed
// to be the first and last points in the array.
float deep_functions::calcPeakArea(std::vector<my_intensities>& v){
  float w;
  float hRect;
  float hTri;
  float total=0;
  for (size_t i = 0; i < v.size()-1; i++) {
    w = v[i + 1].x - v[i].x;
    hTri = abs(v[i + 1].y - v[i].y);
    if (v[i + 1].y < v[i].y) hRect = v[i + 1].y;
    else hRect = v[i].y;
    v[i].tot = w * hRect + w * hTri / 2;
    total+=v[i].tot;
  }
  return total;
}




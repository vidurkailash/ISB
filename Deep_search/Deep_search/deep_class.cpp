#include <string.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <vector>
#include "rapidxml.hpp"
#include <map> 
#include <algorithm>
#include <list>
#include <iterator>
#include <cstdio>
#include<numeric>
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
							my_peptide_lists.all[c].rtime = atof(sample_node->first_attribute("retention_time_sec")->value());
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
					my_peptide_lists.all[c].rtime = atof(sample_node->first_attribute("retention_time_sec")->value());
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

peptide_lists deep_functions::tryptic_calc(peptide_lists& my_peptide_lists) {

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

	return my_peptide_lists;

}

peptide_lists deep_functions::miss_cleave(peptide_lists& my_peptide_lists) {

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

	return my_peptide_lists;

}

peptide_lists deep_functions::delete_dup(peptide_lists& my_peptide_lists) {

	my_peptide_lists.tryp_unique_z = my_peptide_lists.tryptic;
	for (int i = 0; i < my_peptide_lists.tryp_unique_z.size() - 1; i++) {
		for (int j = i + 1; j < my_peptide_lists.tryp_unique_z.size(); j++) {
			if ((my_peptide_lists.tryp_unique_z[i].pep_seq == my_peptide_lists.tryp_unique_z[j].pep_seq) && (my_peptide_lists.tryp_unique_z[i].charge == my_peptide_lists.tryp_unique_z[j].charge)) {
				my_peptide_lists.tryp_unique_z.erase(my_peptide_lists.tryp_unique_z.begin() + j);
				j = (i + 1) - 1;
			}
			else { continue; }
		}
	}


	my_peptide_lists.tryp_unique = my_peptide_lists.tryptic;
	for (int i = 0; i < my_peptide_lists.tryp_unique.size() - 1; i++) {
		for (int j = i + 1; j < my_peptide_lists.tryp_unique.size(); j++) {
			if ((my_peptide_lists.tryp_unique[i].pep_seq == my_peptide_lists.tryp_unique[j].pep_seq)) {
				my_peptide_lists.tryp_unique.erase(my_peptide_lists.tryp_unique.begin() + j);
				j = (i + 1) - 1;
			}
			else { continue; }
		}
	}


	my_peptide_lists.miss_unique = my_peptide_lists.tryp_unique; 
	for (int i = 0; i < my_peptide_lists.miss_unique.size(); i++) {
		if (my_peptide_lists.miss_unique[i].miss_cleaves == 0) {
			my_peptide_lists.miss_unique.erase(my_peptide_lists.miss_unique.begin() + i); 
			i = i - 1; 
		}
	}



	return my_peptide_lists;

}



peptide_lists deep_functions::lcd(peptide_lists& my_peptide_lists) {
	

	/*for (int i = 0; i < my_peptide_lists.tryp_unique.size(); i++) {
		for (int j = 0; j < my_peptide_lists.miss_unique.size(); j++) {
			size_t found = my_peptide_lists.tryp_unique[i].pep_seq.find(my_peptide_lists.miss_unique[j].pep_seq);
			if (found != string::npos && (my_peptide_lists.miss_unique[j].pep_seq != my_peptide_lists.tryp_unique[i].pep_seq)) {
				my_peptide_lists.tryp_unique[i].d_pep_seq = my_peptide_lists.miss_unique[j].pep_seq;
				my_peptide_lists.tryp_unique[i].d_pep_seq_rt = my_peptide_lists.miss_unique[j].rtime;
				if (found == 0) {
					my_peptide_lists.tryp_unique[i].cleave_loc = 'R';
					my_peptide_lists.tryp_unique[i].cleave_pos = found;

				}
				else {
					my_peptide_lists.tryp_unique[i].cleave_loc = 'L';
					my_peptide_lists.tryp_unique[i].cleave_pos = found;

				}
				my_peptide_lists.tryp_unique[i].d_pep_seq_mass = my_peptide_lists.miss_unique[j].mass;
				my_peptide_lists.tryp_unique[i].d_pep_seq_charge = my_peptide_lists.miss_unique[j].charge;
			}
			
			
		}

		my_peptide_lists.tryp_unique[i].mz = (my_peptide_lists.tryp_unique[i].mass + ((my_peptide_lists.tryp_unique[i].charge) * 1.00727)) / my_peptide_lists.tryp_unique[i].charge;
		my_peptide_lists.tryp_unique[i].d_pep_seq_mz = (my_peptide_lists.tryp_unique[i].d_pep_seq_mass + ((my_peptide_lists.tryp_unique[i].d_pep_seq_charge) * 1.00727)) / my_peptide_lists.tryp_unique[i].d_pep_seq_charge;
	}*/



	for (int i = 0; i < my_peptide_lists.miss_unique.size(); i++) {
		for (int j = 0; j < my_peptide_lists.tryp_unique.size(); j++) {
			size_t found = my_peptide_lists.miss_unique[i].pep_seq.find(my_peptide_lists.tryp_unique[j].pep_seq);
			if (found != string::npos && (my_peptide_lists.miss_unique[i].pep_seq != my_peptide_lists.tryp_unique[j].pep_seq)) {
				my_peptide_lists.miss_unique[i].d_pep_seq = my_peptide_lists.tryp_unique[j].pep_seq;
				my_peptide_lists.miss_unique[i].d_pep_seq_rt = my_peptide_lists.tryp_unique[j].rtime;
				my_peptide_lists.miss_unique[i].d_miss_cleaves = my_peptide_lists.tryp_unique[j].miss_cleaves;
				if (found == 0) {
					my_peptide_lists.miss_unique[i].cleave_loc = 'R';
					my_peptide_lists.miss_unique[i].cleave_pos = found;

				}
				else {
					my_peptide_lists.miss_unique[i].cleave_loc = 'L';
					my_peptide_lists.miss_unique[i].cleave_pos = found;

				}
				my_peptide_lists.miss_unique[i].d_pep_seq_mass = my_peptide_lists.tryp_unique[j].mass;
				my_peptide_lists.miss_unique[i].d_pep_seq_charge = my_peptide_lists.tryp_unique[j].charge;
			}

		}

		my_peptide_lists.miss_unique[i].mz = (my_peptide_lists.miss_unique[i].mass + ((my_peptide_lists.miss_unique[i].charge) * 1.00727)) / my_peptide_lists.miss_unique[i].charge;
		my_peptide_lists.miss_unique[i].d_pep_seq_mz = (my_peptide_lists.miss_unique[i].d_pep_seq_mass + ((my_peptide_lists.miss_unique[i].d_pep_seq_charge) * 1.00727)) / my_peptide_lists.miss_unique[i].d_pep_seq_charge;
	}




	/*for (int i = 0; i < my_peptide_lists.tryp_unique.size(); i++) {
		cout << my_peptide_lists.tryp_unique[i].pep_seq << "  " << my_peptide_lists.tryp_unique[i].charge << "   " << my_peptide_lists.tryp_unique[i].mass << "  " << my_peptide_lists.tryp_unique[i].mz << "  " << my_peptide_lists.tryp_unique[i].cleave_loc << "   " << my_peptide_lists.tryp_unique[i].cleave_pos  << "   " << my_peptide_lists.tryp_unique[i].d_pep_seq << "  " << my_peptide_lists.tryp_unique[i].d_pep_seq_charge << "  " << my_peptide_lists.tryp_unique[i].d_pep_seq_mass << "  " << my_peptide_lists.tryp_unique[i].d_pep_seq_mz << endl;
	}*/

	/*cout << my_peptide_lists.tryp_unique.size() << endl;*/



	return my_peptide_lists; 

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


	/*cout << my_peptide_lists.d_list.size() << endl;*/

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
	for (int i = 0; i < my_peptide_lists.tryp_unique.size(); i++) {
		number.push_back(my_peptide_lists.tryp_unique[i].pep_seq.size());
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
	my_metrics.tryp_frac = my_metrics.tryptic_num / my_metrics.psm_num;
	my_metrics.nontryp_frac = my_metrics.nontryptic_num / my_metrics.psm_num;
	my_metrics.pep_frac = my_peptide_lists.tryp_unique.size() / my_metrics.psm_num;
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
#include "protein_reader.h"

using namespace std; 


bool protein_reader::prot_stats(peptide_lists& my_peptide_lists) {



	vector<dsPeptide> temp;
	my_peptide_lists.xic_total_results = my_peptide_lists.xic_ft_results; 
	
	int c = 0; 

	for (size_t i = 0; i < my_peptide_lists.xic_mc_results.size(); i++) {
		my_peptide_lists.xic_total_results.push_back(my_peptide_lists.xic_mc_results[i]);
	}
	
	for (size_t i = 0; i < my_peptide_lists.xic_total_results.size(); i++) {
		if (my_peptide_lists.xic_total_results[i].proteotypic == 0) continue;
		temp.push_back(my_peptide_lists.xic_total_results[i]);
	}
	my_peptide_lists.xic_total_results = temp;
	
	
	
	sort(my_peptide_lists.xic_total_results.begin(), my_peptide_lists.xic_total_results.end(), compareProt);
	
	
	temp.clear(); 
	my_peptide_lists.prot_f_real.push_back(temp);
	my_peptide_lists.prot_f_real.back().push_back(my_peptide_lists.xic_total_results[0]);
	for (size_t i = 1; i < my_peptide_lists.xic_total_results.size(); i++) {
		if (my_peptide_lists.xic_total_results[i].prot_seq != my_peptide_lists.xic_total_results[i - 1].prot_seq) {
			my_peptide_lists.prot_f_real.push_back(temp);
		}
		my_peptide_lists.prot_f_real.back().push_back(my_peptide_lists.xic_total_results[i]);
	}


	
	char f = 'f';
	char m = 'm';

	for (size_t i = 0; i < my_peptide_lists.prot_f_real.size(); i++) {
		my_peptide_lists.prot_f.push_back(dsProtein());
		my_peptide_lists.prot_f.back().sumTryp = prot_calc(my_peptide_lists.prot_f_real[i], f);
		my_peptide_lists.prot_f.back().sumMiss = prot_calc(my_peptide_lists.prot_f_real[i], m);
		my_peptide_lists.prot_f.back().prot_seq = my_peptide_lists.prot_f_real[i][0].prot_seq; 
		my_peptide_lists.prot_f.back().trypPeptides = index_counter(my_peptide_lists.prot_f_real[i], f); 
		my_peptide_lists.prot_f.back().missPeptides = index_counter(my_peptide_lists.prot_f_real[i], m);
		my_peptide_lists.prot_f.back().sumTryp = prot_calc(my_peptide_lists.prot_f_real[i], f);
		my_peptide_lists.prot_f.back().percentMiss = my_peptide_lists.prot_f.back().sumMiss / (my_peptide_lists.prot_f.back().sumMiss + my_peptide_lists.prot_f.back().sumTryp); 
	}

	


	return true; 
}

float protein_reader::prot_calc(std::vector<dsPeptide>& v, char cleave) {

	float val = 0; 

	switch (cleave) {
	case 'f':
		for (size_t i = 0; i < v.size(); i++) {
			if (v[i].miss_cleaves == 0) {
				val += v[i].areaXIC;
			}
		}
		break; 
	case 'm':
		for (size_t i = 0; i < v.size(); i++) {
			if (v[i].miss_cleaves != 0) {
				val += v[i].areaXIC;
			}
		}
		break; 
	default:
		break; 

	}

	return val; 
}


vector<string> protein_reader::index_counter(std::vector<dsPeptide>& v, char cleave) {

	vector<string> tmp; 
	tmp.clear(); 

	switch (cleave) {
	case 'f':
		for (size_t i = 0; i < v.size(); i++) {
			if (v[i].miss_cleaves == 0) {
				tmp.push_back(v[i].pep_seq);
			}
		}
		break;
	case 'm':
		for (size_t i = 0; i < v.size(); i++) {
			if (v[i].miss_cleaves != 0) {
				tmp.push_back(v[i].pep_seq);
			}
		}
		break;
	default:
		break;

	}

	return tmp; 
}

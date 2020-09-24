#include "protein_reader.h"

using namespace std; 


bool protein_reader::prot_stats(peptide_lists& my_peptide_lists) {



	/*vector<dsPeptide> temp;
	my_peptide_lists.xic_total_results = my_peptide_lists.xic_ft_results; 
	
	int c = 0; 

	for (size_t i = 0; i < my_peptide_lists.xic_mc_results.size(); i++) {
		my_peptide_lists.xic_total_results.push_back(my_peptide_lists.xic_mc_results[i]);
	}
	
	for (size_t i = 0; i < my_peptide_lists.xic_total_results.size(); i++) {
		if (my_peptide_lists.xic_total_results[i].proteotypic == 0) continue;
		temp.push_back(my_peptide_lists.xic_total_results[i]);
	}
	my_peptide_lists.xic_total_results = temp;*/
	
	
	
	sort(my_peptide_lists.enzymatic_unique.begin(), my_peptide_lists.enzymatic_unique.end(), compareProt);

	vector<dsPeptide> tmp;
	my_peptide_lists.prot_f_real.push_back(tmp);
		my_peptide_lists.prot_f_real.back().push_back(my_peptide_lists.enzymatic_unique[0]);
		for (size_t i = 1; i < my_peptide_lists.enzymatic_unique.size(); i++) {
			if (my_peptide_lists.enzymatic_unique[i].prot_seq != my_peptide_lists.enzymatic_unique[i - 1].prot_seq) {
				my_peptide_lists.prot_f_real.push_back(tmp);
			}
			my_peptide_lists.prot_f_real.back().push_back(my_peptide_lists.enzymatic_unique[i]);
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
			my_peptide_lists.prot_f.back().percentMiss = my_peptide_lists.prot_f.back().sumMiss / (my_peptide_lists.prot_f.back().sumMiss + my_peptide_lists.prot_f.back().sumTryp);
			my_peptide_lists.prot_f.back().total = my_peptide_lists.prot_f.back().sumMiss + my_peptide_lists.prot_f.back().sumTryp;
		}
	
	

	/*sort(my_peptide_lists.prot_f.begin(), my_peptide_lists.prot_f.end(), compareTotal);
	for (int i = 0; i < 500; i++) {
		cout << my_peptide_lists.prot_f[i].total << endl;
	}*/
	
	vector<string> greater_than_40;
	vector<string> greater_than_15; 
	vector<string> total; 

	/*for (int i = 0; i < my_peptide_lists.prot_f.size(); i++) {
		if (my_peptide_lists.prot_f[i].percentMiss > 0.40) {
			greater_than_40.push_back(my_peptide_lists.prot_f[i].prot_seq);
		}
	}
	cout << "greater than 40: " << greater_than_40.size() << endl; 
	for (int i = 0; i < greater_than_40.size(); i++) {
		cout << greater_than_40[i] << endl; 
	}*/


	/*for (int i = 0; i < my_peptide_lists.prot_f.size(); i++) {
		if (my_peptide_lists.prot_f[i].percentMiss > 0.15) {
			greater_than_15.push_back(my_peptide_lists.prot_f[i].prot_seq);
		}
	}
	cout << "greater than 15: " << greater_than_15.size() << endl;
	for (int i = 0; i < greater_than_15.size(); i++) {
		cout << greater_than_15[i] << endl;
	}*/
	

	/*for (int i = 0; i < my_peptide_lists.prot_f.size(); i++) {
		if (my_peptide_lists.prot_f[i].percentMiss >= 0) {
			total.push_back(my_peptide_lists.prot_f[i].prot_seq);
		}
	}
	cout << "total " << total.size() << endl;
	for (int i = 0; i < total.size(); i++) {
		cout << total[i] << endl;
	}*/
	
	return true; 
}

//bool protein_reader::prot_stats(peptide_lists& my_peptide_lists) {
//
//
//
//	vector<dsPeptide> temp;
//	my_peptide_lists.xic_total_results = my_peptide_lists.xic_ft_results;
//
//	int c = 0;
//
//	for (size_t i = 0; i < my_peptide_lists.xic_mc_results.size(); i++) {
//		my_peptide_lists.xic_total_results.push_back(my_peptide_lists.xic_mc_results[i]);
//	}
//
//	for (size_t i = 0; i < my_peptide_lists.xic_total_results.size(); i++) {
//		if (my_peptide_lists.xic_total_results[i].proteotypic == 0) continue;
//		temp.push_back(my_peptide_lists.xic_total_results[i]);
//	}
//	my_peptide_lists.xic_total_results = temp;
//
//
//
//	sort(my_peptide_lists.xic_total_results.begin(), my_peptide_lists.xic_total_results.end(), compareProt);
//
//
//	temp.clear();
//	my_peptide_lists.prot_f_real.push_back(temp);
//	my_peptide_lists.prot_f_real.back().push_back(my_peptide_lists.xic_total_results[0]);
//	for (size_t i = 1; i < my_peptide_lists.xic_total_results.size(); i++) {
//		if (my_peptide_lists.xic_total_results[i].prot_seq != my_peptide_lists.xic_total_results[i - 1].prot_seq) {
//			my_peptide_lists.prot_f_real.push_back(temp);
//		}
//		my_peptide_lists.prot_f_real.back().push_back(my_peptide_lists.xic_total_results[i]);
//	}
//
//
//
//	char f = 'f';
//	char m = 'm';
//
//	for (size_t i = 0; i < my_peptide_lists.prot_f_real.size(); i++) {
//		my_peptide_lists.prot_f.push_back(dsProtein());
//		my_peptide_lists.prot_f.back().sumTryp = prot_calc(my_peptide_lists.prot_f_real[i], f);
//		my_peptide_lists.prot_f.back().sumMiss = prot_calc(my_peptide_lists.prot_f_real[i], m);
//		my_peptide_lists.prot_f.back().prot_seq = my_peptide_lists.prot_f_real[i][0].prot_seq;
//		my_peptide_lists.prot_f.back().trypPeptides = index_counter(my_peptide_lists.prot_f_real[i], f);
//		my_peptide_lists.prot_f.back().missPeptides = index_counter(my_peptide_lists.prot_f_real[i], m);
//		my_peptide_lists.prot_f.back().percentMiss = my_peptide_lists.prot_f.back().sumMiss / (my_peptide_lists.prot_f.back().sumMiss + my_peptide_lists.prot_f.back().sumTryp);
//		my_peptide_lists.prot_f.back().total = my_peptide_lists.prot_f.back().sumMiss + my_peptide_lists.prot_f.back().sumTryp;
//	}
//
//	/*sort(my_peptide_lists.prot_f.begin(), my_peptide_lists.prot_f.end(), compareTotal);
//	for (int i = 0; i < 500; i++) {
//		cout << my_peptide_lists.prot_f[i].total << endl;
//	}*/
//
//	vector<string> greater_than_40;
//	vector<string> greater_than_15;
//	vector<string> total;
//
//	/*for (int i = 0; i < my_peptide_lists.prot_f.size(); i++) {
//		if (my_peptide_lists.prot_f[i].percentMiss > 0.40) {
//			greater_than_40.push_back(my_peptide_lists.prot_f[i].prot_seq);
//		}
//	}
//	cout << "greater than 40: " << greater_than_40.size() << endl;
//	for (int i = 0; i < greater_than_40.size(); i++) {
//		cout << greater_than_40[i] << endl;
//	}*/
//
//
//	/*for (int i = 0; i < my_peptide_lists.prot_f.size(); i++) {
//		if (my_peptide_lists.prot_f[i].percentMiss > 0.15) {
//			greater_than_15.push_back(my_peptide_lists.prot_f[i].prot_seq);
//		}
//	}
//	cout << "greater than 15: " << greater_than_15.size() << endl;
//	for (int i = 0; i < greater_than_15.size(); i++) {
//		cout << greater_than_15[i] << endl;
//	}*/
//
//
//	/*for (int i = 0; i < my_peptide_lists.prot_f.size(); i++) {
//		if (my_peptide_lists.prot_f[i].percentMiss >= 0) {
//			total.push_back(my_peptide_lists.prot_f[i].prot_seq);
//		}
//	}
//	cout << "total " << total.size() << endl;
//	for (int i = 0; i < total.size(); i++) {
//		cout << total[i] << endl;
//	}*/
//
//	return true;
//}

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


vector<size_t> protein_reader::index_counter(std::vector<dsPeptide>& v, char cleave) {

	vector<size_t> tmp; 
	tmp.clear(); 

	switch (cleave) {
	case 'f':
		for (size_t i = 0; i < v.size(); i++) {
			if (v[i].miss_cleaves == 0) {
				tmp.push_back(i);
			}
		}
		break;
	case 'm':
		for (size_t i = 0; i < v.size(); i++) {
			if (v[i].miss_cleaves != 0) {
				tmp.push_back(i);
			}
		}
		break;
	default:
		break;

	}

	return tmp; 
}

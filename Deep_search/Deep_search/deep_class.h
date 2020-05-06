#ifndef deep_class_h
#define deep_class_h

#include "deep_structs.h"
#include <fstream>
#include <iostream>
#include<numeric>
#include <vector> 
#include <algorithm>

class peptide_lists {
public:
	int total = 0;
	//std::vector<my_features> all;
	//std::vector<my_features> tryptic;
	//std::vector<my_features> non_tryptic;
	//std::vector<my_features> tryp_unique_z;
	//std::vector<my_features> tryp_unique;
	//std::vector<my_features> fully_tryp_unique;
	//std::vector<my_features> miss_unique; 

	///*std::vector<my_features> d_list;*/ 

	//std::vector<my_intensities> spectra;
	//std::vector<my_intensities> spectra1;
	//std::vector<std::vector<my_intensities>> master; 
	//std::vector<std::vector<my_intensities>> master1;
	//std::string file; 
	//std::vector<my_intensities> final; 
	//std::vector<my_intensities> final1;
	// 
	//

	// REAL VECTORS

	std::vector<dsPeptide> all_real; 
	std::vector<dsPeptide> tryptic_real;
	std::vector<dsPeptide> non_tryptic_real;
	std::vector<dsPeptide> tryp_unique_z_real;
	std::vector<dsPeptide> tryp_unique_real;
	std::vector<dsPeptide> miss_unique_real;
	std::vector<dsPeptide> fully_tryp_unique_real;

	std::vector<dsPeptide> xic_ft_results;
	std::vector<dsPeptide> xic_mc_results; 
	
	std::vector<dsPair> peptide_matches; 
	std::vector<dsPair> peptide_matches1;
	std::vector<std::vector<dsPair>> peptide_matches_ftm;
	std::vector<std::vector<dsPair>> peptide_matches_mcm;

	std::vector<dsPeptide> xic_total_results;
	std::vector<std::vector<dsPeptide>> prot_f_real;
	std::vector<dsProtein> prot_f;
	std::vector<dsProtein> prot_f1;

};

//class match_lists {
//public:
//	std::vector<my_markers> results; 
//	std::vector<my_markers> results1;
//	
//};

class deep_functions {
public:
	peptide_lists xml_parse(my_parameters& my_params);
	bool tryptic_calc(peptide_lists& my_peptide_lists);
	bool miss_cleave(peptide_lists& my_peptide_lists);
	bool delete_dup(peptide_lists& my_peptide_lists);
	bool reader(peptide_lists& my_peptide_lists, metrics& my_metrics);
	metrics calc(peptide_lists& my_peptide_lists, metrics& my_metrics);
	bool lcd(peptide_lists& my_peptide_lists, my_parameters& my_params); 
	/*peptide_lists new_list(peptide_lists& my_peptide_lists); */
	metrics calc1(peptide_lists& my_peptide_lists, metrics& my_metircs, my_parameters& my_params);

	bool cleanNoise(std::vector<dsXIC>& v);
	float calcPeakArea(std::vector<dsXIC>& v);
	void print(peptide_lists& my_peptide_lists, metrics& my_metrics);

private: 
	static bool compareInten(const dsPair& a, const dsPair& b) { return a.ft_areaXIC > b.ft_areaXIC; };
	static bool compareSeqZ(const dsPeptide& a, const dsPeptide& b) {
		int i = a.pep_seq.compare(b.pep_seq);
		if (i == 0) return (a.charge < b.charge);
		else return (i < 0);
	}
	static bool compareTrypIndex(const dsPair& a, const dsPair& b) { return a.trypIndex < b.trypIndex; }
	static bool compareMissIndex(const dsPair& a, const dsPair& b) { return a.missIndex < b.missIndex; }
	static bool comparePercentMiss(const dsProtein& a, const dsProtein& b) { return a.percentMiss < b.percentMiss;  }
}; 


#endif
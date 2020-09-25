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
	int total_psms_in_xml;
	std::vector<dsPSM> all_psm; 
	std::vector<dsPeptide> all_peptides; 
	//std::vector<std::vector<dsPeptide>> prot_f_real;
	std::vector<dsProtein> all_proteins;


	std::vector<dsPeptide> tryptic_real;
	std::vector<dsPeptide> non_tryptic_real;
	std::vector<dsPeptide> semi_tryptic_real; 
	std::vector<dsPeptide> tryp_unique_z_real;
	std::vector<dsPeptide> semi_tryptic_unique_z_real;
	std::vector<dsPeptide> tryp_unique_real;
	std::vector<dsPeptide> miss_unique_real;
	std::vector<dsPeptide> semi_tryptic_unique_real;
	std::vector<dsPeptide> fully_tryp_unique_real;

	std::vector<dsPeptide> xic_ft_results;
	std::vector<dsPeptide> xic_mc_results; 
	std::vector<dsPeptide> xic_st_results; 
	
	std::vector<dsPair> peptide_matches; 
	std::vector<dsPair> peptide_matches1;
	std::vector<std::vector<dsPair>> peptide_matches_ftm;
	std::vector<std::vector<dsPair>> peptide_matches_mcm;

	std::vector<dsPeptide> xic_total_results;
	
	
	std::vector<dsProtein> prot_f1;


	//Functions
	bool delete_dup();
	bool enzymatic_calc(my_parameters& my_params);
	void json(std::string fn);
	bool miss_cleave(my_parameters& my_params);
	bool prot_stats();
	bool reader();
	void xml_parse(my_parameters& my_params);

private:
	bool cleanNoise(std::vector<dsXIC>& v);
	float calcPeakArea(std::vector<dsXIC>& v);
	dsPeptide convert_PSM_to_peptide(dsPSM& psm);

	static bool compareInten(const dsPair& a, const dsPair& b) { return a.ft_areaXIC > b.ft_areaXIC; };
	static bool compareSeqZ(const dsPSM& a, const dsPSM& b) {
		int i = a.pep_seq.compare(b.pep_seq);
		if (i == 0) return (a.charge < b.charge);
		else return (i < 0);
	}
	static bool compareTrypIndex(const dsPair& a, const dsPair& b) { return a.trypIndex < b.trypIndex; }
	static bool compareMissIndex(const dsPair& a, const dsPair& b) { return a.missIndex < b.missIndex; }
	static bool compareTotal(const dsProtein& a, const dsProtein& b) { return a.total > b.total; }
	static bool compareProt(const dsPeptide& a, const dsPeptide& b) { return a.prot_seq < b.prot_seq; }

};

//class match_lists {
//public:
//	std::vector<my_markers> results; 
//	std::vector<my_markers> results1;
//	
//};

class deep_functions {
public:
	//peptide_lists xml_parse(my_parameters& my_params);
	//bool enzymatic_calc(peptide_lists& my_peptide_lists, my_parameters& my_params);
	bool semi_enzymatic_calc(peptide_lists& my_peptide_lists, my_parameters& my_params);
	//bool miss_cleave(peptide_lists& my_peptide_lists, my_parameters& my_params);
	//bool delete_dup(peptide_lists& my_peptide_lists);
	//bool reader(peptide_lists& my_peptide_lists, metrics& my_metrics);
	metrics calc(peptide_lists& my_peptide_lists, metrics& my_metrics);
	//bool lcd(peptide_lists& my_peptide_lists, my_parameters& my_params); 
	/*peptide_lists new_list(peptide_lists& my_peptide_lists); */
	metrics calc1(peptide_lists& my_peptide_lists, metrics& my_metircs, my_parameters& my_params);

	//bool cleanNoise(std::vector<dsXIC>& v);
	//float calcPeakArea(std::vector<dsXIC>& v);
	void print(peptide_lists& my_peptide_lists, metrics& my_metrics);
	//void json(peptide_lists& my_peptide_lists, std::string fn); 

//private: 
//	static bool compareInten(const dsPair& a, const dsPair& b) { return a.ft_areaXIC > b.ft_areaXIC; };
//	static bool compareSeqZ(const dsPSM& a, const dsPSM& b) {
//		int i = a.pep_seq.compare(b.pep_seq);
//		if (i == 0) return (a.charge < b.charge);
//		else return (i < 0);
//	}
//	static bool compareTrypIndex(const dsPair& a, const dsPair& b) { return a.trypIndex < b.trypIndex; }
//	static bool compareMissIndex(const dsPair& a, const dsPair& b) { return a.missIndex < b.missIndex; }
//	static bool compareTotal(const dsProtein& a, const dsProtein& b) { return a.total > b.total;  }
}; 


#endif
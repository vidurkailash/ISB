#ifndef deep_class_h
#define deep_class_h

#include "deep_structs.h"
#include <vector> 

class peptide_lists {
public:
	int total = 0;
	std::vector<my_features> all;
	std::vector<my_features> tryptic;
	std::vector<my_features> non_tryptic;
	std::vector<my_features> tryp_unique_z;
	std::vector<my_features> tryp_unique;
	std::vector<my_features> miss_unique; 
	std::vector<my_features> d_list; 
	std::vector<my_intensities> spectra;
	std::vector<my_intensities> spectra1;
	std::vector<std::vector<my_intensities>> master; 
	std::vector<std::vector<my_intensities>> master1;
	std::string file; 
	std::vector<my_intensities> final; 
	std::vector<my_intensities> final1;
	

};

class match_lists {
public:
	std::vector<my_markers> results; 
	std::vector<my_markers> results1;
	
};

class deep_functions {
public:
	peptide_lists xml_parse(my_parameters& my_params);
	peptide_lists tryptic_calc(peptide_lists& my_peptide_lists);
	peptide_lists miss_cleave(peptide_lists& my_peptide_lists);
	peptide_lists delete_dup(peptide_lists& my_peptide_lists);
	//match_lists mzml(peptide_lists& my_peptide_lists, my_parameters& my_params);
	peptide_lists reader(peptide_lists& my_peptide_lists, metrics& my_metrics, match_lists& my_match_lists);
	metrics calc(peptide_lists& my_peptide_lists, metrics& my_metrics);
	peptide_lists lcd(peptide_lists& my_peptide_lists); 
	peptide_lists new_list(peptide_lists& my_peptide_lists); 
	metrics calc1(peptide_lists& my_peptide_lists, metrics& my_metircs, my_parameters& my_params);

}; 


#endif
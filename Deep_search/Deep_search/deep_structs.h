#ifndef deep_structs
#define deep_structs

#include <string>

typedef struct parameters {

	double pep_prob; //peptide prophet probability parameter 
	double ipro_prob; //iprophet probability parameter 
	std::string filename; //filename parameter
	char ident; //identifieer for ipro or peptide prophet 
	std::string mzml; //spectra file


} my_parameters;

//MH: Good rules to know about data types. 
//1. When you need discrete values, use int (i.e. charge states)
//2. When you store array indexes, use size_t (non-negative, scales with architecture)
//3. When storing m/z values, use double to have maximum precision with your decimal places
//   This is very important when working with parts-per-million and parts-per-billion precisions.
//4. When decimal place precision isn't so important (i.e. retention time, peak intensity) then use float.
//   Or if you are not sure, don't use float and just use double.
typedef struct features {

	std::string pep_seq;
	int charge;
	float rtime;
	double mass; 
	double mz; 
	int miss_cleaves;
	std::string prev_aa;
	std::string next_aa;

	std::string d_pep_seq;
	int d_pep_seq_charge;
	float d_pep_seq_rt;
	double d_pep_seq_mass;
	double d_pep_seq_mz; 
	int d_miss_cleaves;
	char cleave_loc;
	int cleave_pos; 

	
} my_features;



typedef struct intensities {
	float x;
	float y; 
	std::string seq;
	float tot;
	int mc;

} my_intensities;


typedef struct compare {
	std::string miss_cleave_seq;
	std::string tryp_seq;
	float mc_tot;
	float tp_tot;
	int mc_mc;
	int tp_mc;
	int matches; 
	float zero_frac; 
	float ratio; 

} my_compare;

//MH: commenting this out as it seems to be kruft.
//typedef struct results {
//	float miss_tot;
//	float tryp_tot; 
//
//} my_results;

typedef struct markers {

	float spec_rt;
	int spec_sn;
	int spec_size; 
	double spec_mz; 
	float mzml_rt;
	double mzml_mz;
	float spec_intensity;
	std::string pep_seq; 
	int miss_cleaves; 

} my_markers; 


typedef struct metrics {

	size_t total_psm; //check
	size_t psm_num; //check
	size_t tryptic_num; //check
	size_t nontryptic_num; //check
	size_t unique_pep_charge; //check
	size_t unique_peptides; //check
	double avg_pep_length; //check
	double tryp_frac; //check
	double nontryp_frac; //check
	double pep_frac; //check 
	double miss_cleave_rate_psm; //check
	double miss_cleave_rate_pep; //check
	double miss_cleave_avg;
	double num_miss_cleave_pep; 
	double num_miss_cleave_total; 
	double num_miss_cleave_psm; 
	double num_tot_miss_cleave_psm; 
	double golden_stat_psm;
	double golden_stat_unique_pep; 
	size_t find_num;
	double miss_avg_high; 
	double tryp_avg_high;
	int twice_mc; 
	int once_mc; 
	double zero; 
	double intensity_final; 
	double stdv_final; 

} my_metrics;


#endif
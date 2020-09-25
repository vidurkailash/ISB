#ifndef deep_structs
#define deep_structs

#include <string>
#include <vector>

typedef struct dsXIC {
  float rTime;
  float intensity;
  float tot; 
} dsXIC;


typedef struct dsPSM {
  std::string pep_seq;
  std::string prot_seq;
	bool decoy;
  bool proteotypic = 0;
  int charge;
  float xml_rtime;
  double pre_neutral_mass;
  double calc_neutral_mass;
  double xml_mz;
  double probability;
  int miss_cleaves; 
  char prev_aa;     
  char next_aa;   
	bool non_enzymatic = 0; 
	bool semi_enzymatic = 0; 
	int psm_count;
	std::vector<dsXIC> XIC;  //stands for eXtracted Ion Chromatogram
	double areaXIC;
} dsPSM;


typedef struct dsPeptide {
	//basic peptide identifiers
  std::string pep_seq;
  std::string prot_seq;
	bool decoy;
	double calc_neutral_mass;

	//describes peptide type
  bool proteotypic = 0;
  bool non_enzymatic = 0;
	int miss_cleaves; //this variable has great dual functionality. If it is 0, then fully tryptic peptide. If >0, then it must be miscleaved.

	//peptide quantifiers
  double areaXIC;
  int psm_count;

} dsPeptide;


typedef struct dsPair {
  size_t trypIndex;
  size_t missIndex;
  std::string ft_pep_seq; 
  std::string mc_pep_seq; 
  double ratio;
  double ft_areaXIC;
  double mc_areaXIC;
} dsPair;


typedef struct dsProtein {
	std::string prot_seq;
	std::vector<size_t> peptides; 
	double sumEnz;
	double sumMiss;
	double sumNonSp;
	double total; //total signal of all peptides
	int sumPSM;
} dsProtein;



typedef struct parameters {

	double probability;   //peptide prophet probability parameter 
	std::string filename; //filename parameter
	bool iprophet;        //identifieer for ipro or peptide prophet 
	std::string mzml;     //spectra file
	float ret_time; 
	double ppm; 
	std::string cleave_loc;
	std::string hyphen; 
	std::string anti_cleave_loc;
	std::string decoy;

} my_parameters;



//MH: Good rules to know about data types. 
//1. When you need discrete values, use int (i.e. charge states)
//2. When you store array indexes, use size_t (non-negative, scales with architecture)
//3. When storing m/z values, use double to have maximum precision with your decimal places
//   This is very important when working with parts-per-million and parts-per-billion precisions.
//4. When decimal place precision isn't so important (i.e. retention time, peak intensity) then use float.
//   Or if you are not sure, don't use float and just use double.



typedef struct metrics {

	//PSM statistics
	int psm_total; //number of target PSMs above probability threshold read in from the results
	int psm_enzymatic;
	int psm_miscleave;
	int psm_nonspecific;
	int psm_unique;

	//peptide statistics
	int pep_count; //number of peptides
	double pep_total; //total signal
	double pep_enzymatic;
	double pep_miscleave;
	double pep_nonspecific;

	//protein statistics
	int prot_count;
	double prot_avg_enzymatic;
	double prot_avg_miscleave;
	double prot_avg_nonspecific;

	//size_t psm_num; //check
	//size_t tryptic_num; //check
	//size_t nontryptic_num; //check
	//size_t unique_pep_charge; //check
	//size_t unique_peptides; //check
	//double avg_pep_length; //check
	//double tryp_frac; //check
	//double nontryp_frac; //check
	//double pep_frac; //check 
	//double miss_cleave_rate_psm; //check
	//double miss_cleave_rate_pep; //check
	//double miss_cleave_avg;
	//double num_miss_cleave_pep; 
	//double num_miss_cleave_total; 
	//double num_miss_cleave_psm; 
	//double num_tot_miss_cleave_psm; 
	//double golden_stat_psm;
	//double golden_stat_unique_pep; 
	//size_t find_num;
	//double miss_avg_high; 
	//double tryp_avg_high;
	//int twice_mc; 
	//int once_mc; 
	//float intensity_final; 
	//float stdv_final; 
	//float protein_final; 
	//float protein_stdv; 
	//float total_intensity; 

} my_metrics;


#endif
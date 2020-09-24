#ifndef protein_reader_class_h
#define protein_reader_class_h

#include "deep_class.h" 


#include <numeric> 
#include <vector> 

class protein_reader {

public: 
	bool prot_stats(peptide_lists& my_peptide_lists);
	float prot_calc(std::vector<dsPeptide>& v, char cleave);
	std::vector<size_t> index_counter(std::vector<dsPeptide>& v, char cleave); 

private: 
	static bool compareProt(const dsPeptide& a, const dsPeptide& b) { return a.prot_seq < b.prot_seq; }
	static bool compareTotal(const dsProtein& a, const dsProtein& b) { return a.total > b.total; }
};


#endif
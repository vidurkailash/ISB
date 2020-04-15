// deep_search_run.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <string.h>
#include <stdio.h>
#include <fstream>
#include <vector>
#include "deep_class.h"
#include "scan_reader.h"
#include <algorithm>
#include <list>
#include <iterator>
#include <cstdio>
#include <numeric>


using namespace std;

my_parameters cmd_input(int argc, char* argv[]);
void warnings(my_parameters& my_params);

// COMMAND LINE PARSE FUNCTION 
my_parameters cmd_input(int argc, char* argv[]) {

	double input_prob = 10;
	string input_filename = "q";
	string input_mzml = "spec";
  my_parameters my_params;

  //MH: set default values
  my_params.ident='b';

	for (int i = 1; i < argc; i++) {
		if (string(argv[i]) == "--threshold" || string(argv[i]) == "-t") {
			if (i + 1 < argc) {
				input_prob = atof(argv[++i]);
			}
			else {
				cout << "probability entered is invalid" << endl;
				exit(1);
				// call warning function 
			}
		}
    else if (string(argv[i]) == "--iprophet" || string(argv[i]) == "-i") {
      my_params.ident = 'a';
    }
		else if (string(argv[i]) == "--filename" || string(argv[i]) == "-f") {
			if (i + 1 < argc) {
				input_filename = string(argv[++i]);
			}
			else {
				cout << "xml filename entered is invalid" << endl;
				exit(1);
				// call warning function 
			}
		}
		else if (string(argv[i]) == "--mzml" || string(argv[i]) == "-m") {
			if (i + 1 < argc) {
				input_mzml = string(argv[++i]);
			}
			else {
				cout << "mzml filename entered is invalid" << endl;
				exit(1);
				// call warning function 
			}
		}
		else {
			cout << "one or more parameters entered are invalid, please consult documentation and re-enter valid parameters" << endl;
			exit(1);
		}
	}

	//MAKES SURE THESE PARAMETERS ARE ENTERED AND ARE VALID
	if (input_prob == 10) {
		cout << "probability entered is invalid" << endl; exit(1);
	}
	if (input_filename == "q") {
		cout << "xml filename entered is invalid" << endl; exit(1);
	}
	if (input_mzml == "spec") {
		cout << "mzml filename entered is invalid" << endl; exit(1);
	}


	//INPUT PARAMETERS IN DATA STRUCTURE
  //MH: no longer need for user intervention.
	if(my_params.ident == 'a') my_params.ipro_prob = input_prob;
	else my_params.pep_prob = input_prob;	
	my_params.filename = input_filename;
	my_params.mzml = input_mzml; 

	return my_params;

}


void warnings(my_parameters& my_params) {

	if (my_params.filename.end()[-1] != 'l') {
		cout << "filename is incorrect or not found" << endl;
		exit(1);
	}

	switch (my_params.ident) {
	case 'a':
		if (my_params.ipro_prob < 0 || my_params.ipro_prob > 1) {
			cout << "threshold iprophet probability value invalid" << endl;
			exit(1);
		}
		if (my_params.ipro_prob < 0.9) {
			char ans;
			cout << "Are you sure you want this probability?: (Y/N)" << endl;
			cin >> ans;
			if (ans == 'N') {
				cout << "Please re enter valid probability" << endl;
				exit(1);
			}
		}
		break;
	case 'b':
		if (my_params.pep_prob < 0 || my_params.pep_prob > 1) {
			cout << "threshold peptide probability value invalid" << endl;
			exit(1);
		}
		if (my_params.pep_prob < 0.9) {
			char ans;
			cout << "Are you sure you want this probability?: (Y/N)" << endl;
			cin >> ans;
			if (ans == 'N') {
				cout << "Please re enter valid probability" << endl;
				exit(1);
			}
		}
		break;
	default:
		cout << "please enter a probability" << endl;
		break;
	}

  cout << "Parameters good." << endl;

}




int main(int argc, char* argv[])
{

  //TODO: Add usage statement
  if(argc==1){
  }


	deep_functions my_deep_functions;
  scan_reader sr;

	my_parameters my_params;
	my_params = cmd_input(argc, argv);
	warnings(my_params);

	peptide_lists my_peptide_lists;
	match_lists my_match_lists;
  metrics my_metrics;
	
  //MH: Parse the pepXML file and return a list of PSMs
  my_peptide_lists = my_deep_functions.xml_parse(my_params);
	cout << "XML parsed, " << my_peptide_lists.all.size() << " of " << my_peptide_lists.total << " PSMs above probability threshold." << endl; 
	
  //MH: Count the tryptic and non-tryptic PSMs
	if(my_deep_functions.tryptic_calc(my_peptide_lists)){
    cout << my_peptide_lists.tryptic.size() << " tryptic PSMs, and " << my_peptide_lists.non_tryptic.size() << " non-tryptic PSMs" << endl;
  } else {
    //MH: Right now, the function can never return false
  }
	
  //MH: Count the miscleaved PSMs
  if(my_deep_functions.miss_cleave(my_peptide_lists)){
    int mc=0;
    for(size_t i=0;i<my_peptide_lists.tryptic.size();i++){
      if(my_peptide_lists.tryptic[i].miss_cleaves>0) mc++;
    }
    cout << mc << " of " << my_peptide_lists.tryptic.size() << " tryptic PSMs are miscleaved." << endl;
  } else {
    //MH: Right now, the function can never return false
  }

  //MH: Trim the list from PSMs to peptides
	if(my_deep_functions.delete_dup(my_peptide_lists)){
    cout << my_peptide_lists.tryp_unique_z.size() << " unique tryptic PSMs." << endl;
    cout << "  " << my_peptide_lists.tryp_unique.size() << " are unique tryptic peptides," << endl;
    cout << "  and " << my_peptide_lists.miss_unique.size() << " of those have miscleavages." << endl;
  } else {
    //MH: Right now, the function can never return false
  }

  //MH: This function iterates over all peptides to determine pairs of tryptic and miscleaved to target.
	if(my_deep_functions.lcd(my_peptide_lists)){
    cout << my_peptide_lists.d_list.size() << " unique tryptic-miscleaved peptide pairs." << endl;
    //for(size_t i=0;i<my_peptide_lists.d_list.size();i++){
    //  cout << my_peptide_lists.d_list[i].d_pep_seq << "\t" << my_peptide_lists.d_list[i].pep_seq << endl;
    //}
  } else {
    //MH: Right now, the function can never return false
  }
  //MH: Obsolete
	//my_peptide_lists = my_deep_functions.new_list(my_peptide_lists); 

  //MH: check if there are peptide pairs to analyze
	if (my_peptide_lists.d_list.size() == 0) {
    cout << "no matches between micleaved and trypic sequences were found" << endl; 
		return 1;
	}


  cout << my_peptide_lists.d_list.size() << " matches found... next step (it's a big one)" << endl;
	my_match_lists = sr.mzml(my_peptide_lists, my_params);
	cout << "mzml file cross checked... next step " << "\n" << endl;

	if(my_deep_functions.reader(my_peptide_lists, my_metrics, my_match_lists)){
  } else {
      //MH: Right now, the function can never return false
  }


		cout << "almost there... calculating metrics... results will be printed shortly" << "\n" << endl;

		my_metrics = my_deep_functions.calc(my_peptide_lists, my_metrics);
		my_metrics = my_deep_functions.calc1(my_peptide_lists, my_metrics, my_params);
		
	

	cout << "\n" << "----------- XML DATA -----------" << "\n" << endl;
	cout << "-- SANTIY CHECKS --" << "\n" << endl; 
	cout << "Total PSM in file:\t\t\t\t\t\t\t\t" << my_metrics.total_psm << "\n" << endl;
	cout << "# of PSM above entered threshold value:\t\t\t\t\t\t" << my_metrics.psm_num << "\n" << endl;
	cout << "# of tryptic PSM above entered threshold value:\t\t\t\t\t" << my_metrics.tryptic_num << "\n" << endl;
	cout << "# of non tryptic PSM above entered threshold value:\t\t\t\t" << my_metrics.nontryptic_num << "\n" << endl;
	cout << "Tryptic peptides with unique charges:\t\t\t\t\t\t" << my_metrics.unique_pep_charge << "\n" << endl;
	cout << "Tryptic peptides with unique sequences:\t\t\t\t\t\t" << my_metrics.unique_peptides << "\n" << endl;
	cout << "Avg seuqnce length of unique tryptic peptide:\t\t\t\t\t" << my_metrics.avg_pep_length << "\n" << endl;
	cout << "Tryptic psm / psm above threshold:\t\t\t\t\t\t" << my_metrics.tryp_frac << "\n" << endl;
	cout << "Non tryptic psm / psm above threshold:\t\t\t\t\t\t" << my_metrics.nontryp_frac << "\n" << endl;
	cout << "Unique peptides / psm above threshold:\t\t\t\t\t\t" << my_metrics.pep_frac << "\n" << endl;
	cout << "Avg # of misscleaves per unique peptide that is misclevaed:\t\t\t" << my_metrics.miss_cleave_avg << "\n" << endl;

	cout << "-- USEFUL STATS --" << "\n" << endl; 
	cout << "# of unique pepides that are miss cleaved:\t\t\t\t\t" << my_metrics.num_miss_cleave_pep << "\n" << endl;
	cout << "# of total misscleaves amognst unique peptides:\t\t\t\t\t" << my_metrics.num_miss_cleave_total << "\n" << endl;
	cout << "# of total misscleaves amongst unique peptides / total # of unique peptides:\t" << my_metrics.miss_cleave_rate_pep << "\n" << endl;
	cout << "# of misscleaved tryptic psm above threshold: \t\t\t\t\t" << my_metrics.num_miss_cleave_psm << "\n" << endl;
	cout << "# of total miss cleaves amongst tryptic psm:\t\t\t\t\t" << my_metrics.num_tot_miss_cleave_psm << "\n" << endl;
	cout << "# of total misscleaves amongst tryptic psm / total # tryptic psm aove threshold:\t" << my_metrics.miss_cleave_rate_psm << "\n" << endl;
	cout << "# of misscleaved unique peptides / total # of unique peptides:\t\t\t" << my_metrics.golden_stat_unique_pep << " ** " << "\n" << endl;
	cout << "# of miss cleaved tryptic psm / total tryptic psm above threshold:\t\t" << my_metrics.golden_stat_psm << " ** " << "\n" << endl; 

	cout << "\n" << "---------- MZML DATA ---------" << "\n" << endl;
	cout << "-- SANTIY CHECKS --" << "\n" << endl; 
	cout << "# of fully tryptic peptides found that has a match to a miscleaved peptide:\t" << my_metrics.find_num << "\n" << endl; 
	cout << "avg highest peak intensity for miscleaved peptide: \t\t\t\t" << my_metrics.miss_avg_high << "\n" << endl; 
	cout << "avg highest peak intensity for trypitc peptide: \t\t\t\t\t" << my_metrics.tryp_avg_high << "\n" << endl; 
	cout << "# of miscleaved peptides that have fully tryptic matches that have 2 miscleaves:\t " << my_metrics.twice_mc << "\n" << endl; 
	cout << "# of miscleaved peptides that have fully tryptic matches that have 1 miscleave:\t " << my_metrics.once_mc << "\n" << endl;

	cout << "-- USEFUL STATS --" << "\n" << endl; 
	cout << "fraction of tryptic matches that have zero miscleaves:\t\t\t " << my_metrics.zero << "\n" << endl; 
	cout << "avg ratio of intensities: 0 miscleave full tryptic peptides / 1 or 2 miscleave peptides:\t" << my_metrics.intensity_final << " ** " << "\n" << endl; 
	cout << "stdv of ratio of intensities: \t\t\t\t\t\t\t\t" << my_metrics.stdv_final << " ** " << "\n" << endl; 

	return 0;

}

// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu

// Tips for Getting Started: 
//   1. Use the Solution Explorer window to add/manage files
//   2. Use the Team Explorer window to connect to source control
//   3. Use the Output window to see build output and other messages
//   4. Use the Error List window to view errors
//   5. Go to Project > Add New Item to create new code files, or Project > Add Existing Item to add existing code files to the project
//   6. In the future, to open this project again, go to File > Open > Project and select the .sln file

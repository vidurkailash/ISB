// deep_search_run.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <string.h>
#include <stdio.h>
#include <fstream>
#include <vector>
#include "deep_class.h"
#include "scan_reader.h"
#include "protein_reader.h"
#include <algorithm>
#include <list>
#include <iterator>
#include <cstdio>
#include <numeric>


using namespace std;

my_parameters cmd_input(int argc, char* argv[]);
void warnings(my_parameters& my_params);
void info_description(); 

// COMMAND LINE PARSE FUNCTION 
my_parameters cmd_input(int argc, char* argv[]) {

	double input_prob = 10;
	string input_filename = "q";
	string input_mzml = "spec";
	float input_rtime = 2; 
	double input_ppm = 10; 
	string input_loc = "KR"; 
	string input_anti_loc = "P"; 
	
 
    my_parameters my_params;

  //MH: set default values
  my_params.ident='b';

	for (int i = 1; i < argc; i++) {
		if (string(argv[i]) == "--threshold" || string(argv[i]) == "-t") {
			if (i + 1 < argc) {
				input_prob = atof(argv[++i]);
			}
			else {
				/*cout << "probability entered is invalid" << endl;
				exit(1);*/
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
				/*cout << "xml filename entered is invalid" << endl;
				exit(1);*/
				// call warning function 
			}
		}
		else if (string(argv[i]) == "--mzml" || string(argv[i]) == "-m") {
			if (i + 1 < argc) {
				input_mzml = string(argv[++i]);
			}
			else {
				/*cout << "mzml filename entered is invalid" << endl;
				exit(1);*/
				// call warning function 
			}
		}
		else if (string(argv[i]) == "--rtime" || string(argv[i]) == "-r") {
			if (i + 1 < argc) {
				input_rtime = (float)atof(argv[++i]);
			}
		}
		else if (string(argv[i]) == "--ppm" || string(argv[i]) == "-p") {
			if (i + 1 < argc) {
				input_ppm = (float)atof(argv[++i]);
		    }
		}
		else if (string(argv[i]) == "--loc" || string(argv[i]) == "-c") {
			if (i + 1 < argc) {
				input_loc = string(argv[++i]);
			}
		}
		else if (string(argv[i]) == "--anti" || string(argv[i]) == "-a") {
			if (i + 1 < argc) {
				input_anti_loc = string(argv[++i]);
			}
		}

		/*else {
			cout << "one or more parameters entered are invalid, please consult documentation and re-enter valid parameters" << endl;
			exit(1);
		}*/
	}

	//MAKES SURE THESE PARAMETERS ARE ENTERED AND ARE VALID
	if (input_prob == 10) {
		/*cout << "probability entered is invalid" << endl; exit(1);*/
		info_description(); 
	}
	if (input_filename == "q") {
		/*cout << "xml filename entered is invalid" << endl; exit(1);*/
		info_description();
	}
	if (input_mzml == "spec") {
		/*cout << "mzml filename entered is invalid" << endl; exit(1);*/
		info_description();
	}


	//INPUT PARAMETERS IN DATA STRUCTURE
  //MH: no longer need for user intervention.
	if(my_params.ident == 'a') my_params.ipro_prob = input_prob;
	else my_params.pep_prob = input_prob;	
	my_params.filename = input_filename;
	my_params.mzml = input_mzml; 
	my_params.run_time = input_rtime;
	my_params.ppm = input_ppm; 
	input_loc.push_back('-'); 
	my_params.cleave_loc = input_loc; 
	my_params.anti_cleave_loc = input_anti_loc;

	return my_params;

}


void info_description() {
	
	cout << "QUICK TUTORIAL" << "\n" << endl; 
	cout << "Example:" << endl; 
	cout << ">deep_search_run.exe -f <file.pep.xml> -t <0-1> -m <file.mzml> [-i -r <float> -p <float> -c <string> -a <string> ]" << "\n" << endl;
	cout << "The three mandatory parameters are the xml file(-f), the threshold value(-t), and the corresponding mzml file(-m)." << endl; 
	cout << "The optional parameters include the iprophet tag(-i), rtime search window(-r), ppm threshold(-p), cleave sites (-c), and their exceptions (-a). If not specified, the defaults are peptide prophet search, 2 min, 10 ppm, K/R, and P." << "\n" << endl;  
	cout << "Some common errors: parameter tags are not correct, xml file and mzml file do not corresponding with each other, threshold value, rtime value, and/or ppm value out of bounds." << "\n" << endl; 
	cout << "Consult documentation for further clarification." << endl; 

	exit(1); 

}


void warnings(my_parameters& my_params) {

	//if (my_params.filename.end()[-1] != 'l') {
	//	cout << "Error: xml and/or mzml filename is incorrect or not found. Consult documentation." << endl;
	//	info_description();
	//	/*exit(1);*/
	//}
	//if (my_params.mzml.end()[-1] != 'l') {
	//	cout << "Error: xml and/or mzml filename is incorrect or not found. Consult documentation." << endl;
	//	info_description();
	//	/*exit(1);*/
	//}

	switch (my_params.ident) {
	case 'a':
		if (my_params.ipro_prob < 0 || my_params.ipro_prob > 1) {
			cout << "Error: threshold iprophet probability value invalid. Consult documentation." << endl;
			info_description();
			/*exit(1);*/
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
			cout << "Error: threshold peptide probability value invalid. Consult documentation." << endl;
			info_description();
			/*exit(1);*/
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

	deep_functions my_deep_functions;
	scan_reader sr;
	protein_reader pr; 

	my_parameters my_params;
	my_params = cmd_input(argc, argv);
	warnings(my_params);

	peptide_lists my_peptide_lists;
	
	metrics my_metrics;
	
  //Parse the pepXML file and return a list of PSMs
  my_peptide_lists = my_deep_functions.xml_parse(my_params);
	cout << "XML parsed, " << my_peptide_lists.all_real.size() << " of " << my_peptide_lists.total << " PSMs above probability threshold." << endl; 



	
  //Count the tryptic and non-tryptic PSMs
	if(my_deep_functions.tryptic_calc(my_peptide_lists, my_params)){
    cout << my_peptide_lists.tryptic_real.size() << " tryptic PSMs, and " << my_peptide_lists.non_tryptic_real.size() << " non-tryptic PSMs" << endl;
  } else {
    //Right now, the function can never return false
  }
	
	
	

  //Count the miscleaved PSMs
  if(my_deep_functions.miss_cleave(my_peptide_lists, my_params)){
    int mc=0;
    for(size_t i=0;i<my_peptide_lists.tryptic_real.size();i++){
      if(my_peptide_lists.tryptic_real[i].miss_cleaves>0) mc++;
    }
    cout << mc << " of " << my_peptide_lists.tryptic_real.size() << " tryptic PSMs are miscleaved." << endl;
  } else {
    //Right now, the function can never return false
  }

 


	if(my_deep_functions.delete_dup(my_peptide_lists)){
    cout << my_peptide_lists.tryp_unique_z_real.size() << " unique tryptic PSMs." << endl;
    cout << my_peptide_lists.tryp_unique_real.size() << " are unique tryptic peptides," << endl;
    cout << "and " << my_peptide_lists.miss_unique_real.size() << " of those have miscleavages." << endl;
	} 
	else {
    //Right now, the function can never return false
	}



 //MH: This function iterates over all peptides to determine pairs of tryptic and miscleaved to target.
	if(my_deep_functions.lcd(my_peptide_lists, my_params)){
		cout << "\n" << "m/z calculated for both lists" << "\n" << endl; 
	} else {
    //MH: Right now, the function can never return false
	}

  

	my_peptide_lists = sr.mzml(my_peptide_lists, my_params);
	cout << "mzml file cross checked... next step " << "\n" << endl;

	
	



	if(my_deep_functions.reader(my_peptide_lists, my_metrics)){} 
	else {
      //MH: Right now, the function can never return false
	}


	
	cout << "computing protein level stats" << endl; 
	if(pr.prot_stats(my_peptide_lists)){}
	else{}


	cout << "almost there... calculating metrics... results will be printed shortly" << "\n" << endl;

	my_metrics = my_deep_functions.calc(my_peptide_lists, my_metrics);
	my_metrics = my_deep_functions.calc1(my_peptide_lists, my_metrics, my_params);
		
	my_deep_functions.print(my_peptide_lists, my_metrics); 

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

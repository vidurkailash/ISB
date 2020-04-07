#include <string.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <map> 
#include <algorithm>
#include <list>
#include <iterator>
#include <cstdio>
#include<numeric>
#include "deep_class.h"
#include "MSReader.h"

using namespace std;
using namespace MSToolkit;



match_lists deep_functions::mzml(peptide_lists& my_peptide_lists, my_parameters& my_params) {

	//READ MZML FILE FOR BOTH FULLY TRYPTIC AND MISCLEAVED PEPTIDES (MZML)
	cout << "starting pass 1" << "\n" << endl; 

	MSReader myfile;
	Spectrum mySpec;

	myfile.setFilter(MS1);
	myfile.readFile(my_params.mzml.c_str(), mySpec);


	match_lists my_match_lists;


	int b = 0;
	int c = 0;
	while (mySpec.getRTime() != 0) {

		/*cout << "-------------" << mySpec.getRTime() << "-------------" << endl; */

		for (int z = 0; z < mySpec.size(); z++) {
			/*cout << "-------------" << mySpec.getScanNumber() << "----------" << endl; */
			/*cout << mySpec[z].mz << "\t" << mySpec[z].intensity << endl;*/
			for (int i = 0; i < my_peptide_lists.d_list.size(); i++) {
				/*cout << my_peptide_lists.d_list[i].rtime / 60 << "      " << my_peptide_lists.d_list[i].d_pep_seq_rt / 60 << endl;*/

				if (((my_peptide_lists.d_list[i].rtime / 60) <= mySpec.getRTime() + 2) && ((my_peptide_lists.d_list[i].rtime / 60) >= mySpec.getRTime() - 2) && (my_peptide_lists.d_list[i].mz <= mySpec[z].mz + 0.01) && (my_peptide_lists.d_list[i].mz >= mySpec[z].mz - 0.01)) {

					my_match_lists.results.push_back(my_markers()); //input 
					/*cout << "-------------" << mySpec.getRTime() << "-------------" << endl;*/
					my_match_lists.results[c].spec_rt = mySpec.getRTime(); //input
					/*cout << "-------------" << mySpec.getScanNumber() << "----------" << mySpec.size() << "----------" << endl;*/
					my_match_lists.results[c].spec_sn = mySpec.getScanNumber(); //input
					my_match_lists.results[c].spec_size = mySpec.size(); //input
					/*cout << mySpec[z].mz << "\t" << mySpec[z].intensity << endl;
					cout << my_peptide_lists.d_list[i].rtime / 60 << "   " << my_peptide_lists.d_list[i].mz << endl;*/
					my_match_lists.results[c].spec_mz = mySpec[z].mz; //input
					my_match_lists.results[c].spec_intensity = mySpec[z].intensity; //input
					my_match_lists.results[c].mzml_mz = my_peptide_lists.d_list[i].mz; //input   
					my_match_lists.results[c].mzml_rt = my_peptide_lists.d_list[i].rtime / 60; // input
					my_match_lists.results[c].pep_seq = my_peptide_lists.d_list[i].pep_seq; //input 
					my_match_lists.results[c].miss_cleaves = my_peptide_lists.d_list[i].miss_cleaves;
					b++;
					c++;
					/*cout << "NEXT MATCH" << "     " << b << "    " << my_match_lists.results.size() << endl;*/
					continue;
				}
			}
			/*if (z < (0.25 * mySpec.size()) + 0.5 && z >((0.25 * mySpec.size()) - 0.5)) {
				cout << "25% done with pass 1" << "\n" << endl;
			}
			if (z < ((0.5 * mySpec.size()) + 0.5) && z >((0.5 * mySpec.size()) - 0.5)) {
				
			}
			if (z < ((0.75 * mySpec.size()) + 0.5) && z >((0.75 * mySpec.size()) - 0.5)) {
				cout << "75% done with pass 1" << "\n" << endl;
			}*/

		}

		/*if (mySpec.getRTime() < 75.01 && mySpec.getRTime() > 74.99) {
			cout << "50% done with pass 1" << "\n" << endl;
		}*/

		myfile.readFile(NULL, mySpec);
	}


	cout << "pass 1 done" << "\n" << endl;


	MSReader myfile1;
	Spectrum mySpec1;

	myfile1.setFilter(MS1);
	myfile1.readFile(my_params.mzml.c_str(), mySpec1);


	int k = 0;
	int l = 0;
	while (mySpec1.getRTime() != 0) {

		/*cout << "-------------" << mySpec.getRTime() << "-------------" << endl; */

		for (int z = 0; z < mySpec1.size(); z++) {
			/*cout << "-------------" << mySpec.getScanNumber() << "----------" << endl; */
			/*cout << mySpec[z].mz << "\t" << mySpec[z].intensity << endl;*/
			for (int i = 0; i < my_peptide_lists.d_list.size(); i++) {
				/*cout << my_peptide_lists.d_list[i].rtime / 60 << "      " << my_peptide_lists.d_list[i].d_pep_seq_rt / 60 << endl;*/


				if (((my_peptide_lists.d_list[i].d_pep_seq_rt / 60) <= mySpec1.getRTime() + 2) && ((my_peptide_lists.d_list[i].d_pep_seq_rt / 60) >= mySpec1.getRTime() - 2) && (my_peptide_lists.d_list[i].d_pep_seq_mz <= mySpec1[z].mz + 0.01) && (my_peptide_lists.d_list[i].d_pep_seq_mz >= mySpec1[z].mz - 0.01)) {

					my_match_lists.results1.push_back(my_markers()); //input 
					/*cout << "-------------" << mySpec.getRTime() << "-------------" << endl;*/
					my_match_lists.results1[l].spec_rt = mySpec1.getRTime(); //input
					/*cout << "-------------" << mySpec.getScanNumber() << "----------" << mySpec.size() << "----------" << endl;*/
					my_match_lists.results1[l].spec_sn = mySpec1.getScanNumber(); //input
					my_match_lists.results1[l].spec_size = mySpec1.size(); //input
					/*cout << mySpec[z].mz << "\t" << mySpec[z].intensity << endl;
					cout << my_peptide_lists.d_list[i].rtime / 60 << "   " << my_peptide_lists.d_list[i].mz << endl;*/
					my_match_lists.results1[l].spec_mz = mySpec1[z].mz; //input
					my_match_lists.results1[l].spec_intensity = mySpec1[z].intensity; //input
					my_match_lists.results1[l].mzml_mz = my_peptide_lists.d_list[i].d_pep_seq_mz; //input
					my_match_lists.results1[l].mzml_rt = my_peptide_lists.d_list[i].d_pep_seq_rt / 60; //input
					my_match_lists.results1[l].pep_seq = my_peptide_lists.d_list[i].d_pep_seq; //input 
					my_match_lists.results1[l].miss_cleaves = my_peptide_lists.d_list[i].d_miss_cleaves;
					k++;
					l++;
					/*cout << "NEXT MATCH" << "     " << b << "    " << my_match_lists.results.size() << endl;*/
					continue;
				}

			}
			/*if (z < (0.25 * mySpec1.size()) + 0.5 && z > ((0.25 * mySpec1.size()) - 0.5)) {
				cout << "25% done with pass 1" << "\n" << endl;
			}
			if (z < ((0.5 * mySpec1.size()) + 0.5) && z > ((0.5 * mySpec1.size()) - 0.5)) {
				cout << "50% done with pass 1" << "\n" << endl;
			}
			if (z < ((0.75 * mySpec1.size()) + 0.5) && z > ((0.75 * mySpec1.size()) - 0.5)) {
				cout << "75% done with pass 1" << "\n" << endl;
			}*/

		}

		/*if (mySpec1.getRTime() < 75.01 && mySpec1.getRTime() > 74.99) {
			cout << "50% done with pass 2" << "\n" << endl;
		}*/
		myfile1.readFile(NULL, mySpec1);
	}


	// END MZML FUNCTION

	cout << "pass 2 done" << "\n" << endl; 
	cout << my_match_lists.results.size() << endl;

	// DELETE DUPLICATES (DD)
	for (int i = 0; i < my_match_lists.results.size() - 1; i++) {
		for (int j = i + 1; j < my_match_lists.results.size(); j++) {
			if ((my_match_lists.results[i].pep_seq == my_match_lists.results[j].pep_seq) && (my_match_lists.results[i].spec_sn == my_match_lists.results[j].spec_sn)) {
				my_match_lists.results.erase(my_match_lists.results.begin() + j);
				j = (i + 1) - 1;
			}
		}
	}
	cout << my_match_lists.results.size() << endl; 
	cout << my_match_lists.results1.size() << endl;

	for (int i = 0; i < my_match_lists.results1.size() - 1; i++) {
		for (int j = i + 1; j < my_match_lists.results1.size(); j++) {
			if ((my_match_lists.results1[i].pep_seq == my_match_lists.results1[j].pep_seq) && (my_match_lists.results1[i].spec_sn == my_match_lists.results1[j].spec_sn)) {
				my_match_lists.results1.erase(my_match_lists.results1.begin() + j);
				j = (i + 1) - 1;
			}
		}
	}
	cout << my_match_lists.results1.size() << endl;

	// END DD FUNCTION


	return my_match_lists; 

}




peptide_lists deep_functions::reader(peptide_lists& my_peptide_lists, metrics& my_metrics, match_lists& my_match_lists) {


	//cout << "sort started" << "\n" << endl;
	//cout << my_match_lists.results.size() << endl; 
	///*for (int i = 0; i < my_match_lists.results.size(); i++) {
	//	cout << my_match_lists.results[i].pep_seq << endl;
	//}*/
	//cout << "print done" << endl; 

	//ORDER PEPTIDES AND PUT INTO VECOTR OF VECTORS (VV)
	int q = 0;
	for (int i = 0; i < my_match_lists.results.size() - 1; i++) {
		my_peptide_lists.spectra.push_back(my_intensities());
		my_peptide_lists.spectra[q].x = my_match_lists.results[i].spec_rt;
		my_peptide_lists.spectra[q].y = my_match_lists.results[i].spec_intensity;
		my_peptide_lists.spectra[q].seq = my_match_lists.results[i].pep_seq;
		my_peptide_lists.spectra[q].mc = my_match_lists.results[i].miss_cleaves;
		q++;
		for (int j = i + 1; j < my_match_lists.results.size(); j++) {
			if (my_match_lists.results[i].pep_seq == my_match_lists.results[j].pep_seq) {
				my_peptide_lists.spectra.push_back(my_intensities());
				my_peptide_lists.spectra[q].x = my_match_lists.results[j].spec_rt;
				my_peptide_lists.spectra[q].y = my_match_lists.results[j].spec_intensity;
				my_peptide_lists.spectra[q].seq = my_match_lists.results[i].pep_seq;
				my_peptide_lists.spectra[q].mc = my_match_lists.results[i].miss_cleaves;
				q++;
				my_match_lists.results.erase(my_match_lists.results.begin() + j);
				j = (i + 1) - 1;
			}
		}
		/*q++;*/
	}



	int num = 0;
	for (int i = 0; i < my_match_lists.results1.size() - 1; i++) {
		my_peptide_lists.spectra1.push_back(my_intensities());
		my_peptide_lists.spectra1[num].x = my_match_lists.results1[i].spec_rt;
		my_peptide_lists.spectra1[num].y = my_match_lists.results1[i].spec_intensity;
		my_peptide_lists.spectra1[num].seq = my_match_lists.results1[i].pep_seq;
		my_peptide_lists.spectra1[num].mc = my_match_lists.results1[i].miss_cleaves;
		num++;
		for (int j = i + 1; j < my_match_lists.results1.size(); j++) {
			if (my_match_lists.results1[i].pep_seq == my_match_lists.results1[j].pep_seq) {
				my_peptide_lists.spectra1.push_back(my_intensities());
				my_peptide_lists.spectra1[num].x = my_match_lists.results1[j].spec_rt;
				my_peptide_lists.spectra1[num].y = my_match_lists.results1[j].spec_intensity;
				my_peptide_lists.spectra1[num].seq = my_match_lists.results1[i].pep_seq;
				my_peptide_lists.spectra1[num].mc = my_match_lists.results1[i].miss_cleaves;
				num++;
				my_match_lists.results1.erase(my_match_lists.results1.begin() + j);
				j = (i + 1) - 1;
			}
		}
		/*q++;*/
	}

	

	vector<my_intensities> qc;

	int inc = 0;
	for (int i = 0; i < my_peptide_lists.spectra.size() - 1; i++) {
		my_peptide_lists.master.push_back(qc);
		my_peptide_lists.master[inc].push_back(my_peptide_lists.spectra[i]);
		for (int j = i + 1; j < my_peptide_lists.spectra.size(); j++) {
			if (my_peptide_lists.spectra[i].seq == my_peptide_lists.spectra[j].seq) {
				my_peptide_lists.master[inc].push_back(my_peptide_lists.spectra[j]);
				my_peptide_lists.spectra.erase(my_peptide_lists.spectra.begin() + j);
				j = (i + 1) - 1;
			}
		}
		inc++;
	}

	 

	vector<my_intensities> qc1;

	int inc1 = 0;
	for (int i = 0; i < my_peptide_lists.spectra1.size() - 1; i++) {
		my_peptide_lists.master1.push_back(qc1);
		my_peptide_lists.master1[inc1].push_back(my_peptide_lists.spectra1[i]);
		for (int j = i + 1; j < my_peptide_lists.spectra1.size(); j++) {
			if (my_peptide_lists.spectra1[i].seq == my_peptide_lists.spectra1[j].seq) {
				my_peptide_lists.master1[inc1].push_back(my_peptide_lists.spectra1[j]);
				my_peptide_lists.spectra1.erase(my_peptide_lists.spectra1.begin() + j);
				j = (i + 1) - 1;
			}
		}
		inc1++;
	}


	cout << "vector querry done" << "\n" << endl; 

	// END VV FUNCTION


	//DELETE NOISE (INTENSITIES LESS THAN 10% OF THE MAX) (NOISE)
	vector<float> querry;
	for (int i = 0; i < my_peptide_lists.master.size(); i++) {
		for (int j = 0; j < my_peptide_lists.master[i].size(); j++) {
			querry.push_back(my_peptide_lists.master[i][j].y);
		}
		for (int k = 0; k < querry.size() - 1; k++) {
			for (int m = k + 1; m < querry.size(); m++) {
				if (querry[m] < 0.1 * (*max_element(querry.begin(), querry.end()))) {
					querry.erase(querry.begin() + m);
					m = (k + 1) - 1;
				}
			}
		}
		for (int k = 0; k < querry.size(); k++) {
			if (querry[k] < 0.1 * (*max_element(querry.begin(), querry.end()))) {
				querry.erase(querry.begin() + k);
			}
		}
		for (int n = 0; n < my_peptide_lists.master[i].size() - 1; n++) {
			for (int p = n + 1; p < my_peptide_lists.master[i].size(); p++) {
				vector<float>::iterator it;
				it = find(querry.begin(), querry.end(), my_peptide_lists.master[i][p].y);
				if (it == querry.end()) {
					my_peptide_lists.master[i].erase(my_peptide_lists.master[i].begin() + p);
					p = (n + 1) - 1;
				}
			}
		}
		for (int n = 0; n < my_peptide_lists.master[i].size(); n++) {
			vector<float>::iterator it;
			it = find(querry.begin(), querry.end(), my_peptide_lists.master[i][n].y);
			if (it == querry.end()) {
				my_peptide_lists.master[i].erase(my_peptide_lists.master[i].begin() + n);
			}
		}

		querry.clear();
	}


	vector<float> querry1;
	for (int i = 0; i < my_peptide_lists.master1.size(); i++) {
		for (int j = 0; j < my_peptide_lists.master1[i].size(); j++) {
			querry1.push_back(my_peptide_lists.master1[i][j].y);
		}
		for (int k = 0; k < querry1.size() - 1; k++) {
			for (int m = k + 1; m < querry1.size(); m++) {
				if (querry1[m] < 0.1 * (*max_element(querry1.begin(), querry1.end()))) {
					querry1.erase(querry1.begin() + m);
					m = (k + 1) - 1;
				}
			}
		}
		for (int k = 0; k < querry1.size(); k++) {
			if (querry1[k] < 0.1 * (*max_element(querry1.begin(), querry1.end()))) {
				querry1.erase(querry1.begin() + k);
			}
		}
		for (int n = 0; n < my_peptide_lists.master1[i].size() - 1; n++) {
			for (int p = n + 1; p < my_peptide_lists.master1[i].size(); p++) {
				vector<float>::iterator it1;
				it1 = find(querry1.begin(), querry1.end(), my_peptide_lists.master1[i][p].y);
				if (it1 == querry1.end()) {
					my_peptide_lists.master1[i].erase(my_peptide_lists.master1[i].begin() + p);
					p = (n + 1) - 1;
				}
			}
		}
		for (int n = 0; n < my_peptide_lists.master1[i].size(); n++) {
			vector<float>::iterator it1;
			it1 = find(querry1.begin(), querry1.end(), my_peptide_lists.master1[i][n].y);
			if (it1 == querry1.end()) {
				my_peptide_lists.master1[i].erase(my_peptide_lists.master1[i].begin() + n);
			}
		}

		querry1.clear();
	}

	// END NOISE FUNCTION 


	cout << "noise deleted" << "\n" << endl; 

	//metric calculations
	vector<float> peak; 
	vector<float> peak2;
	for (int i = 0; i < my_peptide_lists.master.size(); i++) {
		for (int j = 0; j < my_peptide_lists.master[i].size(); j++) {
			peak.push_back(my_peptide_lists.master[i][j].y);
		}
		peak2.push_back(*max_element(peak.begin(), peak.end()));
		peak.clear(); 
	}


	float sum = 0; for (int i = 0; i < peak2.size(); ++i) { sum += peak2[i]; }
	float sumi = sum / peak2.size(); 
	my_metrics.miss_avg_high = sumi;

	vector<float> peak1; 
	vector<float> peak3;
	for (int i = 0; i < my_peptide_lists.master1.size(); i++) {
		for (int j = 0; j < my_peptide_lists.master1[i].size(); j++) {
			peak1.push_back(my_peptide_lists.master1[i][j].y);
		}
		peak3.push_back(*max_element(peak1.begin(), peak1.end()));
		peak1.clear(); 
	}
	

	float sum1 = 0; for (int i = 0; i < peak3.size(); ++i) { sum1 += peak3[i]; }
	float sumi1 = sum1 / peak3.size();
	my_metrics.tryp_avg_high = sumi1; 

	
	



	//INTEGRATE INTENSITY FUNCTION (AREA)
	float u = 0;
	float v = 0;
	float h = 0;

	for (int i = 0; i < my_peptide_lists.master.size(); i++) {
		for (int j = 0; j < my_peptide_lists.master[i].size(); j++) {
			if (my_peptide_lists.master[i][j + 1].y > my_peptide_lists.master[i][j].y) {
				u = ((my_peptide_lists.master[i][j + 1].x - my_peptide_lists.master[i][j].x) * (my_peptide_lists.master[i][j].y)) + (((my_peptide_lists.master[i][j + 1].x - my_peptide_lists.master[i][j].x) * (my_peptide_lists.master[i][j + 1].y - my_peptide_lists.master[i][j].y)) / 2) + u;
			}
			else {
				v = ((my_peptide_lists.master[i][j + 1].x - my_peptide_lists.master[i][j].x) * (my_peptide_lists.master[i][j + 1].y)) + (((my_peptide_lists.master[i][j + 1].x - my_peptide_lists.master[i][j].x) * (my_peptide_lists.master[i][j].y - my_peptide_lists.master[i][j + 1].y)) / 2) + v;
			}
			h = u + v;
			my_peptide_lists.master[i][j].tot = h; 
			u = 0;
			v = 0;
		}
		
	}

	float u1 = 0;
	float v1 = 0;
	float h1 = 0;

	for (int i = 0; i < my_peptide_lists.master1.size(); i++) {
		for (int j = 0; j < my_peptide_lists.master1[i].size(); j++) {
			if (my_peptide_lists.master1[i][j + 1].y > my_peptide_lists.master1[i][j].y) {
				u1 = ((my_peptide_lists.master1[i][j + 1].x - my_peptide_lists.master1[i][j].x) * (my_peptide_lists.master1[i][j].y)) + (((my_peptide_lists.master1[i][j + 1].x - my_peptide_lists.master1[i][j].x) * (my_peptide_lists.master1[i][j + 1].y - my_peptide_lists.master1[i][j].y)) / 2) + u1;
			}
			else {
				v1 = ((my_peptide_lists.master1[i][j + 1].x - my_peptide_lists.master1[i][j].x) * (my_peptide_lists.master1[i][j + 1].y)) + (((my_peptide_lists.master1[i][j + 1].x - my_peptide_lists.master1[i][j].x) * (my_peptide_lists.master1[i][j].y - my_peptide_lists.master1[i][j + 1].y)) / 2) + v1;
			}
			h1 = u1 + v1;
			my_peptide_lists.master1[i][j].tot = h1;
			u1 = 0;
			v1 = 0;
		}

	}

	for (int i = 0; i < my_peptide_lists.master.size(); i++) {
		my_peptide_lists.master[i][my_peptide_lists.master[i].size() - 1].tot = 0; 
	}

	vector<float> total; 
	for (int i = 0; i < my_peptide_lists.master.size(); i++) {
		for (int j = 0; j < my_peptide_lists.master[i].size(); j++) {
			total.push_back(my_peptide_lists.master[i][j].tot);
		}
		float a = accumulate(total.begin(), total.end(), 0);
		my_peptide_lists.master[i][my_peptide_lists.master[i].size() - 1].tot = a; 
		a = 0; 
		total.clear(); 
	}

	for (int i = 0; i < my_peptide_lists.master.size(); i++) {
		for (int j = 0; j < my_peptide_lists.master[i].size(); j++) {
			if (my_peptide_lists.master[i][j].tot == 0) {
				my_peptide_lists.master[i][j].tot = my_peptide_lists.master[i][j].y;
			}
		}
	}



	for (int i = 0; i < my_peptide_lists.master1.size(); i++) {
		my_peptide_lists.master1[i][my_peptide_lists.master1[i].size() - 1].tot = 0;
	}

	vector<float> total1;
	for (int i = 0; i < my_peptide_lists.master1.size(); i++) {
		for (int j = 0; j < my_peptide_lists.master1[i].size(); j++) {
			total1.push_back(my_peptide_lists.master1[i][j].tot);
		}
		float a1 = accumulate(total1.begin(), total1.end(), 0);
		my_peptide_lists.master1[i][my_peptide_lists.master1[i].size() - 1].tot = a1;
		a1 = 0;
		total1.clear();
	}

	for (int i = 0; i < my_peptide_lists.master1.size(); i++) {
		for (int j = 0; j < my_peptide_lists.master1[i].size(); j++) {
			if (my_peptide_lists.master1[i][j].tot == 0) {
				my_peptide_lists.master1[i][j].tot = my_peptide_lists.master1[i][j].y;
			}
		}
	}

	vector<my_intensities> final;
	int cycle = 0;
	for (int i = 0; i < my_peptide_lists.master.size(); i++) {
		final.push_back(my_intensities());
		final[cycle].seq = my_peptide_lists.master[i][my_peptide_lists.master[i].size() - 1].seq;
		final[cycle].tot = my_peptide_lists.master[i][my_peptide_lists.master[i].size() - 1].tot;
		final[cycle].mc = my_peptide_lists.master[i][my_peptide_lists.master[i].size() - 1].mc;
		cycle++;
	}

	vector<my_intensities> final1;
	int cycle1 = 0;
	for (int i = 0; i < my_peptide_lists.master1.size(); i++) {
		final1.push_back(my_intensities());
		final1[cycle1].seq = my_peptide_lists.master1[i][my_peptide_lists.master1[i].size() - 1].seq;
		final1[cycle1].tot = my_peptide_lists.master1[i][my_peptide_lists.master1[i].size() - 1].tot;
		final1[cycle1].mc = my_peptide_lists.master1[i][my_peptide_lists.master1[i].size() - 1].mc;
		cycle1++;
	}

	my_peptide_lists.final = final;
	my_peptide_lists.final1 = final1; 

	// END AREA FUNCTION

	
	cout << "intensity calculated" << "\n" << endl; 

	
	
		


	return my_peptide_lists;

}



metrics deep_functions::calc1(peptide_lists& my_peptide_lists, metrics& my_metrics, my_parameters& my_params) {



	vector<vector<my_compare>> bigone;
	vector<vector<my_compare>> bigone1;
	vector<my_compare> end;
	vector<my_compare> end1; 
	vector<my_compare> test; 
	int testing = 0; 
	int iteration = 0;
	int iteration1 = 0;
	int count2 = 0;
	int count1 = 0;
	float count0 = 0;

	for (int i = 0; i < my_peptide_lists.final.size(); i++) {
		for (int j = 0; j < my_peptide_lists.final1.size(); j++) {
			size_t present = my_peptide_lists.final[i].seq.find(my_peptide_lists.final1[j].seq);
			if (present != string::npos && (my_peptide_lists.final[i].seq != my_peptide_lists.final1[j].seq)) {
				end.push_back(my_compare());
				end[iteration].miss_cleave_seq = my_peptide_lists.final[i].seq;
				end[iteration].tryp_seq = my_peptide_lists.final1[j].seq;
				end[iteration].mc_tot = my_peptide_lists.final[i].tot;
				end[iteration].tp_tot = my_peptide_lists.final1[j].tot;
				end[iteration].mc_mc = my_peptide_lists.final[i].mc;
				end[iteration].tp_mc = my_peptide_lists.final1[j].mc;
				iteration++;
			}
		}
	}

	end1 = end;


	for (int i = 0; i < end.size() - 1; i++) {
		bigone.push_back(test);
		bigone[testing].push_back(end[i]);
		for (int j = i + 1; j < end.size(); j++) {
			if (end[i].miss_cleave_seq == end[j].miss_cleave_seq) {
				bigone[testing].push_back(end[j]);
				end.erase(end.begin() + j);
				j = (i + 1) - 1;
			}
		}
		testing++;
	}

	bigone.push_back(test);
	bigone[testing].push_back(end[end.size() - 1]);


	// sanity checks 
	for (int i = 0; i < bigone.size(); i++) {
		for (int j = 0; j < bigone[i].size(); j++) {
			if (bigone[i][j].mc_mc == 2) {
				count2++;
			}
			if (bigone[i][j].mc_mc == 1) {
				count1++;
			}
			if (bigone[i][j].tp_mc == 0) {
				count0++;
			}
		}
		bigone[i][bigone[i].size() - 1].zero_frac = count0 / float(bigone[i].size());
		bigone[i][bigone[i].size() - 1].matches = bigone[i].size();
		count0 = 0;
	}
	vector<double> zero;
	for (int i = 0; i < bigone.size(); i++) {
		for (int j = 0; j < bigone[i].size(); j++) {
			zero.push_back(bigone[i][j].zero_frac);
		}
	}
	double a = accumulate(zero.begin(), zero.end(), 0);
	double b = a / bigone.size();
	my_metrics.zero = b; 
	my_metrics.twice_mc = count2; 
	my_metrics.once_mc = count1; 


	//useful metrics 
	vector<my_compare> temp;
	int temp_it = 0;
	for (int i = 0; i < end1.size() - 1; i++) {
		temp.push_back(my_compare());
		temp[temp_it].miss_cleave_seq = end1[i].miss_cleave_seq;
		temp[temp_it].tryp_seq = end1[i].tryp_seq;
		temp[temp_it].mc_tot = end1[i].mc_tot;
		temp[temp_it].tp_tot = end1[i].tp_tot;
		temp[temp_it].mc_mc = end1[i].mc_mc;
		temp[temp_it].tp_mc = end1[i].tp_mc;
		temp_it++;
		for (int j = i + 1; j < end1.size(); j++) {
			if (end1[i].tryp_seq == end1[j].tryp_seq) {
				temp.push_back(my_compare());
				temp[temp_it].miss_cleave_seq = end1[j].miss_cleave_seq;
				temp[temp_it].tryp_seq = end1[j].tryp_seq;
				temp[temp_it].mc_tot = end1[j].mc_tot;
				temp[temp_it].tp_tot = end1[j].tp_tot;
				temp[temp_it].mc_mc = end1[j].mc_mc;
				temp[temp_it].tp_mc = end1[j].tp_mc;
				temp_it++;
				end1.erase(end1.begin() + j); 
				j = (i + 1) - 1;
			}
		}
	}

	temp.push_back(my_compare());
	temp[temp_it].miss_cleave_seq = end[end.size() - 1].miss_cleave_seq;
	temp[temp_it].tryp_seq = end[end.size() - 1].tryp_seq;
	temp[temp_it].mc_tot = end[end.size() - 1].mc_tot;
	temp[temp_it].tp_tot = end[end.size() - 1].tp_tot;
	temp[temp_it].mc_mc = end[end.size() - 1].mc_mc;
	temp[temp_it].tp_mc = end[end.size() - 1].tp_mc;

	for (int i = 0; i < temp.size() - 1; i++) {
		bigone1.push_back(test);
		bigone1[iteration1].push_back(temp[i]);
		for (int j = i + 1; j < temp.size(); j++) {
			if (temp[i].tryp_seq == temp[j].tryp_seq) {
				bigone1[iteration1].push_back(temp[j]);
				temp.erase(temp.begin() + j);
				j = (i + 1) - 1; 
			}
		}
		iteration1++;
	}

	bigone1.push_back(test); 
	bigone1[iteration1].push_back(temp[temp.size() - 1]);

	vector<float> transient; 
	for (int i = 0; i < bigone1.size(); i++) {
		if (bigone1[i].size() > 1 && bigone1[i][0].tp_mc == 0) {
			for (int j = 0; j < bigone1[i].size(); j++) {
				transient.push_back(bigone1[i][j].mc_tot);
			}
			float hold = accumulate(transient.begin(), transient.end(), 0);
			
			float hold1 = bigone1[i][0].tp_tot / hold;
			bigone1[i][bigone1[i].size() - 1].ratio = hold1;
			transient.clear(); hold = 0; hold1 = 0; 
		}
		if (bigone1[i].size() == 1 && bigone1[i][0].tp_mc == 0) {
			float hold2 = bigone1[i][0].tp_tot / bigone1[i][0].mc_tot;
			bigone1[i][0].ratio = hold2; 
			hold2 = 0; 
		}
	}

	vector<double> rat; 
	for (int i = 0; i < bigone1.size(); i++) {
		for (int j = 0; j < bigone1[i].size(); j++) {
			rat.push_back(bigone1[i][j].ratio);
		}
	}
	
	double rat1 = accumulate(rat.begin(), rat.end(), 0);
	int asdf = 0;
	for (int i = 0; i < rat.size(); i++) {
		if (rat[i] != 0) {
			asdf++;
		}
	}
	double rat2 = rat1 / asdf;
	my_metrics.intensity_final = rat2;
 
	double var = 0;  for (int i = 0; i < rat.size(); ++i) { var += pow(rat[i] - rat2, 2); }
	var = var / rat.size(); 
	double stdv = 0; stdv = sqrt(var); 
	my_metrics.stdv_final = stdv; 




	//PRINT FUNCTIONS 
	/*int meme = 0; 
	for (int i = 0; i < bigone.size(); i++) {
		for (int j = 0; j < bigone[i].size(); j++) {
			cout << bigone[i][j].miss_cleave_seq << "   " << bigone[i][j].mc_tot << "  " << bigone[i][j].mc_mc << "   " << bigone[i][j].tryp_seq << "   " << bigone[i][j].tp_tot << "   " << bigone[i][j].tp_mc << "    " << bigone[i][j].matches << "    " <<  bigone[i][j].zero_frac << endl; 
			meme++; 
		}
	}
	

	int meme1 = 0;
	for (int i = 0; i < bigone1.size(); i++) {
		for (int j = 0; j < bigone1[i].size(); j++) {
			cout << bigone1[i][j].miss_cleave_seq << "   " << bigone1[i][j].mc_tot << "  " << bigone1[i][j].mc_mc << "   " << bigone1[i][j].tryp_seq << "   " << bigone1[i][j].tp_tot << "   " << bigone1[i][j].tp_mc << "    " << bigone1[i][j].ratio << endl;
			meme1++;
		}
	}*/
	/*for (int i = 0; i < bigone1.size(); i++) {
		cout << bigone1[i].size() << endl; 
	}*/
	


	
	/*for (int i = 0; i < temp.size(); i++) {
		cout << temp[i].miss_cleave_seq << "   " << temp[i].mc_tot << "   " << temp[i].mc_mc << "   " << temp[i].tryp_seq << "   " << temp[i].tp_tot << "   " << temp[i].tp_mc << endl;
	}
	cout << temp.size() << endl; 

	for (int i = 0; i < temp.size(); i++) {
		if (temp[i].tryp_seq == "QFEMEELILELAAQVLEDK") {
			cout << "gotem ;) " << endl; 
		}
	}*/
	/*for (int i = 0; i < end.size(); i++) {
		cout << end[i].miss_cleave_seq << "   " << end[i].mc_tot << "   " << end[i].mc_mc << "   " << end[i].tryp_seq << "   " << end[i].tp_tot << "   " << end[i].tp_mc << endl;
	}
	cout << end.size() << endl;*/

	//for (int i = 0; i < my_peptide_lists.master[649].size(); i++) {
	//	/*for (int j = 0; j < my_peptide_lists.master[i].size(); j++) {
	//		cout << my_peptide_lists.master[i][j].seq << "    " << my_peptide_lists.master[i][j].x << "    " << my_peptide_lists.master[i][j].y << "   " << my_peptide_lists.master[i][j].tot << endl;
	//	}*/
	//	cout << my_peptide_lists.master[649][i].seq << "    " << my_peptide_lists.master[649][i].x << "    " << my_peptide_lists.master[649][i].y << "    " << my_peptide_lists.master[649][i].tot << endl;
	//}

	//cout << my_peptide_lists.master.size() << endl;


	//for (int i = 0; i < my_peptide_lists.master1[588].size(); i++) {
	//	/*for (int j = 0; j < my_peptide_lists.master1[i].size(); j++) {
	//		cout << my_peptide_lists.master1[i][j].seq << "    " << my_peptide_lists.master1[i][j].x << "    " << my_peptide_lists.master1[i][j].y << "   " << my_peptide_lists.master1[i][j].tot << endl;
	//	}*/
	//	cout << my_peptide_lists.master1[588][i].seq << "    " << my_peptide_lists.master1[588][i].x << "    " << my_peptide_lists.master1[588][i].y << "     " << my_peptide_lists.master1[588][i].tot << endl;
	//}

	//cout << my_peptide_lists.master1.size() << endl;





	/*cout << my_peptide_lists.final[649].seq << "   " << my_peptide_lists.final[649].tot << "   " << my_peptide_lists.final[649].mc << endl;

	cout << my_peptide_lists.final1[588].seq << "   " << my_peptide_lists.final1[588].tot << "    " << my_peptide_lists.final1[588].mc << endl;*/



	return my_metrics; 
}
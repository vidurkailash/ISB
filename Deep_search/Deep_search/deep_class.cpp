#include "rapidxml.hpp"
#include "deep_class.h"
#include "rapidjson/document.h"
#include "rapidjson/prettywriter.h"
#include "rapidjson/stringbuffer.h"


using namespace rapidxml;
using namespace std;
using namespace rapidjson;



void peptide_lists::xml_parse(my_parameters& my_params) {
  xml_document<> doc;

  // Read the xml file into a vector
  ifstream theFile(my_params.filename);
  vector<char> buffer((istreambuf_iterator<char>(theFile)), istreambuf_iterator<char>());
  buffer.push_back('\0');
  // Parse the buffer using the xml file parsing library into doc 
  doc.parse<0>(&buffer[0]);
  // Find root node
  xml_node<>* root_node = doc.first_node("msms_pipeline_analysis");
  xml_node<>* spectrum_node = root_node->first_node("msms_run_summary");

  total_psms_in_xml = 0;

  for (spectrum_node; spectrum_node; spectrum_node = spectrum_node->next_sibling()){
      
    //skip past the run summary?
    if (string(spectrum_node->value()) == "msms_run_summary") continue;

    //grab the spectrum query
    for (xml_node<>* sample_node = spectrum_node->first_node("spectrum_query"); sample_node; sample_node = sample_node->next_sibling()) {
      total_psms_in_xml++;
      
      //grab on the first hit of the first search result.
      xml_node<>* search_hit = sample_node->first_node("search_result")->first_node("search_hit");
      
      //go to next PSM if below probability threshold.
      xml_node<>* analysis_node = search_hit->first_node("analysis_result");
      string aType=analysis_node->first_attribute("analysis")->value();
      double probability;
      if(my_params.iprophet){
        while(aType.compare("interprophet")!=0){
          analysis_node = analysis_node->next_sibling("analysis_result");
          if(analysis_node ==NULL) {
            cout << " iProphet analysis not found in pepXML. Exiting." << endl;
            exit(12);
          }
          aType = analysis_node->first_attribute("analysis")->value();
        }
        probability = atof(analysis_node->first_node("interprophet_result")->first_attribute("probability")->value());
      } else {
        while (aType.compare("peptideprophet") != 0) {
          analysis_node = analysis_node->next_sibling("analysis_result");
          if (analysis_node == NULL) {
            cout << " PeptideProphet analysis not found in pepXML. Exiting." << endl;
            exit(12);
          }
          aType = analysis_node->first_attribute("analysis")->value();
        }
        probability = atof(analysis_node->first_node("peptideprophet_result")->first_attribute("probability")->value());
      }
      if(probability<my_params.probability) continue;

      //Keeping this PSM, so record necessary information
      dsPSM psm;
      psm.pep_seq = string(search_hit->first_attribute("peptide")->value());
      psm.charge = atoi(sample_node->first_attribute("assumed_charge")->value());
      psm.pre_neutral_mass = atof(sample_node->first_attribute("precursor_neutral_mass")->value());
      psm.xml_rtime = (float)atof(sample_node->first_attribute("retention_time_sec")->value());
      psm.prev_aa = search_hit->first_attribute("peptide_prev_aa")->value()[0];
      psm.next_aa = search_hit->first_attribute("peptide_next_aa")->value()[0];
      psm.prot_seq = string(search_hit->first_attribute("protein")->value());
      psm.calc_neutral_mass = atof(search_hit->first_attribute("calc_neutral_pep_mass")->value());
                            
      //check if peptide is proteotypic, decoy, etc.
      if (my_params.decoy.empty())  psm.decoy = false;
      else {
        psm.decoy = true;
        if (psm.prot_seq.find(my_params.decoy) == string::npos)  psm.decoy = false;
      }
      if (atoi(search_hit->first_attribute("num_tot_proteins")->value()) == 1) {
        psm.proteotypic = true;
      } else {
        psm.proteotypic = false;
        xml_node<>* ap_node = search_hit->first_node("alternative_protein");
        while (ap_node != NULL) {
          string tmp = ap_node->first_attribute("protein")->value();
          if (!my_params.decoy.empty() && tmp.find(my_params.decoy) == string::npos) {
            psm.decoy = false;
            break;
          }
          ap_node = ap_node->next_sibling("alternative_protein");
        }
      }

      //update peptide sequence if it is modified
      xml_node<>* mod_node = search_hit->first_node("modification_info");
      if (mod_node != NULL)  psm.pep_seq = mod_node->first_attribute("modified_peptide")->value();

      all_psm.push_back(psm);
    }
  }

}

bool peptide_lists::enzymatic_calc(my_parameters& my_params) {


    for (size_t i = 0; i < all_psm.size(); i++) {
        size_t found = my_params.cleave_loc.find(all_psm[i].prev_aa);
        size_t found1 = my_params.cleave_loc.find(all_psm[i].pep_seq.back());
        size_t found2 = my_params.hyphen.find(all_psm[i].next_aa);
        if (found!= string::npos && (found1!=string::npos || found2!=string::npos)) {
            all_psm[i].non_enzymatic = false; 
        } else {
          //check to see if the non-specific cut is the leading Met on the protein. If so, this is probably true enzymatic.
          //ALSO: DO THIS SMARTER - LOOK IN THE DB
          if (all_psm[i].prev_aa=='M') all_psm[i].non_enzymatic = false;
          else all_psm[i].non_enzymatic = true;
        }
        
    }
    return true;

}

bool deep_functions::semi_enzymatic_calc(peptide_lists& my_peptide_lists, my_parameters& my_params) {

    for (size_t i = 0; i < my_peptide_lists.all_psm.size(); i++) {
        size_t found = my_params.cleave_loc.find(my_peptide_lists.all_psm[i].prev_aa);
        size_t found1 = my_params.cleave_loc.find(my_peptide_lists.all_psm[i].pep_seq.back());
        size_t found2 = my_params.hyphen.find(my_peptide_lists.all_psm[i].next_aa);
        if (found != string::npos || (found1 != string::npos || found2 != string::npos)) {
            my_peptide_lists.all_psm[i].semi_enzymatic = 1;  
        }
        if (found != string::npos && (found1 != string::npos || found2 != string::npos) && my_peptide_lists.all_psm[i].semi_enzymatic == 1) {
            my_peptide_lists.all_psm[i].semi_enzymatic = 0;
        }
    }
  

    return true; 
}


bool peptide_lists::miss_cleave(my_parameters& my_params) {

  for (int i = 0; i < all_psm.size(); i++) {
    all_psm[i].miss_cleaves=0;
    for (int j = 0; j < all_psm[i].pep_seq.size() - 1; j++) {
      size_t found = my_params.cleave_loc.find(all_psm[i].pep_seq[j]);
      size_t found1 = my_params.anti_cleave_loc.find(all_psm [i].pep_seq[j + 1]);
      if (found!=string::npos && found1==string::npos) {
        all_psm[i].miss_cleaves++;
      }
    }
  }

  return true;

}

dsPeptide peptide_lists::convert_PSM_to_peptide(dsPSM& psm){
  dsPeptide pep;
  pep.pep_seq=psm.pep_seq;
  pep.prot_seq=psm.prot_seq;
  pep.decoy=psm.decoy;
  pep.miss_cleaves=psm.miss_cleaves;
  pep.non_enzymatic=psm.non_enzymatic;
  pep.proteotypic=psm.proteotypic;
  pep.calc_neutral_mass=psm.calc_neutral_mass;
  pep.areaXIC=psm.areaXIC;
  pep.psm_count=psm.psm_count;
  return pep;
}

//This function combines redundant PSMs into a single PSM entry with the lowest retention time.
bool peptide_lists::delete_dup() {

    size_t i;

    //All PSMs are sorted so that they are in order of sequence and charge
    sort(all_psm.begin(), all_psm.end(), compareSeqZ);
    vector<dsPSM> v;
    v.push_back(all_psm[0]);
    v.back().psm_count = 1;

    for (i = 1; i < all_psm.size(); i++) {

      //if we have a redundant PSM, increase the count and adjust RT.
      if (v.back().pep_seq.compare(all_psm[i].pep_seq) == 0 &&
        v.back().charge == all_psm[i].charge) {
        v.back().psm_count++;
        //always store lowest RT
        if (all_psm[i].xml_rtime < v.back().xml_rtime) v.back().xml_rtime = all_psm[i].xml_rtime;

      //otherwise, add this new PSM to the vector
      } else {
        v.push_back(all_psm[i]);
        v.back().psm_count = 1;
      }
    }

    all_psm=v; //copy the vector over the existing PSMs
   
    return true;

}


bool peptide_lists::reader() {

  //iterate over all extracted precursor signals to find the peak and compute area
  for (size_t i = 0; i < all_psm.size(); i++) {
    cleanNoise(all_psm[i].XIC);
    all_psm[i].areaXIC = calcPeakArea(all_psm[i].XIC);
  }

  //MH: Create a peptide list from the PSM list
  for (size_t i = 0; i < all_psm.size(); i++) {
    //check if there is already a peptide for this PSM
    size_t j;
    for (j = 0; j < all_peptides.size(); j++) {
      if (all_peptides[j].pep_seq.compare(all_psm[i].pep_seq) == 0) break;
    }

    if (j == all_peptides.size()) { //add new peptide
      dsPeptide pep=convert_PSM_to_peptide(all_psm[i]);
      all_peptides.push_back(pep);
    } else { //otherwise, sum areas (and psms)
      all_peptides[j].areaXIC += all_psm[i].areaXIC;
      all_peptides[j].psm_count += all_psm[i].psm_count;
    }
  }

  return true;

}


metrics deep_functions::calc(peptide_lists& my_peptide_lists, metrics& my_metrics) {


	vector<int> miss_psm;
	int p = 0;
	for (int i = 0; i < my_peptide_lists.tryptic_real.size(); i++) {
		miss_psm.push_back(my_peptide_lists.tryptic_real[i].miss_cleaves);
		if (my_peptide_lists.tryptic_real[i].miss_cleaves > 0) {
			p++;
		}
	}
	double r = accumulate(miss_psm.begin(), miss_psm.end(), 0);


	vector<int> number;
	vector<int> miss_pep;
	double h = 0;
	for (size_t i = 0; i < my_peptide_lists.tryp_unique_real.size(); i++) {
		number.push_back((int)my_peptide_lists.tryp_unique_real[i].pep_seq.size());
		miss_pep.push_back(my_peptide_lists.tryp_unique_real[i].miss_cleaves);
		if (my_peptide_lists.tryp_unique_real[i].miss_cleaves > 0) {
			h++;
		}
	}
	double f = accumulate(number.begin(), number.end(), 0);
	double g = accumulate(miss_pep.begin(), miss_pep.end(), 0);



	//my_metrics.total_psm = my_peptide_lists.total;
	//my_metrics.psm_num = my_peptide_lists.all_psm.size();
	//my_metrics.tryptic_num = my_peptide_lists.tryptic_real.size();
	//my_metrics.nontryptic_num = my_peptide_lists.non_tryptic_real.size();
	//my_metrics.unique_pep_charge = my_peptide_lists.tryp_unique_z_real.size();
	//my_metrics.unique_peptides = my_peptide_lists.tryp_unique_real.size();
	//my_metrics.avg_pep_length = f / my_metrics.unique_peptides;
	//my_metrics.tryp_frac = (double)my_metrics.tryptic_num / my_metrics.psm_num;
	//my_metrics.nontryp_frac = (double)my_metrics.nontryptic_num / my_metrics.psm_num;
	//my_metrics.pep_frac = (double)my_peptide_lists.tryp_unique_real.size() / my_metrics.psm_num;
	//my_metrics.miss_cleave_rate_psm = r / my_metrics.tryptic_num;
	//my_metrics.miss_cleave_rate_pep = g / my_metrics.unique_peptides;
	//my_metrics.num_miss_cleave_pep = h; 
	//my_metrics.num_miss_cleave_total = g; 
	//my_metrics.miss_cleave_avg = g / h;
	//my_metrics.num_miss_cleave_psm = p; 
	//my_metrics.num_tot_miss_cleave_psm = r; 
	//my_metrics.golden_stat_psm = my_metrics.num_miss_cleave_psm / my_metrics.tryptic_num;
	//my_metrics.golden_stat_unique_pep = my_metrics.num_miss_cleave_pep / my_metrics.unique_peptides;
	

	return my_metrics;

}



bool peptide_lists::cleanNoise(std::vector<dsXIC>& v) {
  size_t i;
  size_t a,b; //MH: these will be your array boundaries

  //MH: apex information
  float max=0;
  size_t maxIndex=0;

  //MH: First step, find the maximum intensity (Assumed to be the peak apex)
  for(i=0;i<v.size();i++){
    if(v[i].intensity>max){
      max=v[i].intensity;
      maxIndex=i;
    }
  }

  //MH: set out threshold
  max/=10;

  //MH: Now step to the left of the apex to find your lower bound
  a=maxIndex;
  while(a>0){
    a--;
    if (v[a].intensity < max) break;
  }

  //MH: Now step to the right of the apex to find your upper bound
  b=maxIndex;
  while(b<v.size()-1){
    b++;
    if (v[b].intensity < max) break;
  }

  //MH: With bounds set, create a subset array
  vector<dsXIC> tmp;
  for(i=a;i<=b;i++){
    tmp.push_back(v[i]);
  }

  //MH: overwrite our original bloated array with the trimmed subset array
  v=tmp;

  //MH: report success
  return true;
	
}




//MH: This function replicates computing (and storing) the area between every two points.
// It then returns the sum of the area of the peak. Note that the boundaries are assumed
// to be the first and last points in the array.
float peptide_lists::calcPeakArea(std::vector<dsXIC>& v){
  float w;
  float hRect;
  float hTri;
  float total=0;
  for (size_t i = 0; i < v.size()-1; i++) {
    w = v[i + 1].rTime - v[i].rTime;
    hTri = abs(v[i + 1].intensity - v[i].intensity);
    if (v[i + 1].intensity < v[i].intensity) hRect = v[i + 1].intensity;
    else hRect = v[i].intensity;
    v[i].tot = w * hRect + w * hTri / 2;
    total+=v[i].tot;
  }
  return total;
}


void deep_functions::print(peptide_lists& my_peptide_lists, metrics& my_metrics) {


  //  cout << "\n" << "----------- XML STATS -----------" << "\n" << endl;
  //  cout << "Total PSM in file:\t\t\t\t\t\t\t\t\t\t" << my_metrics.total_psm << "\n" << endl;
  //  cout << "# of PSM above entered threshold value:\t\t\t\t\t\t\t\t" << my_metrics.psm_num << "\n" << endl;
  //  
  //  cout << "\n" << "----------- PEPTIDE STATS -----------" << "\n" << endl;
  //  cout << "-- Sanity Checks --" << "\n" << endl;
  //  cout << "# of tryptic PSM above entered threshold value:\t\t\t\t\t\t\t" << my_metrics.tryptic_num << "\n" << endl;
  //  cout << "# of non tryptic PSM above entered threshold value:\t\t\t\t\t\t" << my_metrics.nontryptic_num << "\n" << endl;
  //  cout << "Tryptic peptides with unique charges:\t\t\t\t\t\t\t\t" << my_metrics.unique_pep_charge << "\n" << endl;
  //  cout << "Tryptic peptides with unique sequences:\t\t\t\t\t\t\t\t" << my_metrics.unique_peptides << "\n" << endl;
  //  cout << "Avg seuqnce length of unique tryptic peptide:\t\t\t\t\t\t\t" << my_metrics.avg_pep_length << "\n" << endl;
  //  cout << "Tryptic psm / psm above threshold:\t\t\t\t\t\t\t\t" << my_metrics.tryp_frac << "\n" << endl;
  //  cout << "Non tryptic psm / psm above threshold:\t\t\t\t\t\t\t\t" << my_metrics.nontryp_frac << "\n" << endl;
  //  cout << "Unique peptides / psm above threshold:\t\t\t\t\t\t\t\t" << my_metrics.pep_frac << "\n" << endl;
  //  cout << "Avg # of misscleaves per unique peptide that is misclevaed:\t\t\t\t\t" << my_metrics.miss_cleave_avg << "\n" << endl;
  //  cout << "# of fully tryptic peptides found that has a match to a miscleaved peptide:\t\t\t" << my_metrics.find_num << "\n" << endl;
  //  cout << "avg highest peak intensity for miscleaved peptide: \t\t\t\t\t\t" << my_metrics.miss_avg_high << "\n" << endl;
  //  cout << "avg highest peak intensity for trypitc peptide: \t\t\t\t\t\t" << my_metrics.tryp_avg_high << "\n" << endl;
  //  cout << "# of miscleaved peptides that have fully tryptic matches that have 2 miscleaves:\t\t " << my_metrics.twice_mc << "\n" << endl;
  //  cout << "# of miscleaved peptides that have fully tryptic matches that have 1 miscleave:\t\t\t " << my_metrics.once_mc << "\n" << endl;

  //  cout << "-- Useful Stats --" << "\n" << endl;
  //  cout << "# of unique pepides that are miss cleaved:\t\t\t\t\t\t\t" << my_metrics.num_miss_cleave_pep << "\n" << endl;
  //  cout << "# of total misscleaves amognst unique peptides:\t\t\t\t\t\t\t" << my_metrics.num_miss_cleave_total << "\n" << endl;
  //  cout << "# of total misscleaves amongst unique peptides / total # of unique peptides:\t\t\t" << my_metrics.miss_cleave_rate_pep << "\n" << endl;
  //  cout << "# of misscleaved tryptic psm above threshold: \t\t\t\t\t\t\t" << my_metrics.num_miss_cleave_psm << "\n" << endl;
  //  cout << "# of total miss cleaves amongst tryptic psm:\t\t\t\t\t\t\t" << my_metrics.num_tot_miss_cleave_psm << "\n" << endl;
  //  cout << "# of total misscleaves amongst tryptic psm / total # tryptic psm above threshold:\t\t" << my_metrics.miss_cleave_rate_psm << "\n" << endl;
  //  cout << "# of misscleaved unique peptides / total # of unique peptides:\t\t\t\t\t" << my_metrics.golden_stat_unique_pep << " ** " << "\n" << endl;
  //  cout << "# of miss cleaved tryptic psm / total tryptic psm above threshold:\t\t\t\t" << my_metrics.golden_stat_psm << " ** " << "\n" << endl;
  //  cout << "avg ratio of intensities: 0 miscleave full tryptic peptides / 1 or 2 miscleave peptides:\t" << my_metrics.intensity_final << " ** " << "\n" << endl;
  //  cout << "stdv of ratio of intensities: \t\t\t\t\t\t\t\t\t" << my_metrics.stdv_final << " ** " << "\n" << endl;
  //  cout << "Global Stat: total miscleaved peptide intensity / sum of all peptide intensities: \t\t" << my_metrics.total_intensity << " ** " << "\n" << endl;
  // 
  //  vector<float> rat;
  //  float rat1 = 0;
  //  for (size_t i = 0; i < my_peptide_lists.prot_f.size(); i++) {
  //      rat.push_back(my_peptide_lists.prot_f[i].percentMiss);
  //      rat1 += my_peptide_lists.prot_f[i].percentMiss;
  //  }

  //  float rat2 = rat1 / my_peptide_lists.prot_f.size();
  //  my_metrics.protein_final = rat2;

  //  float var = 0;
  //  for (size_t i = 0; i < rat.size(); ++i) {
  //      var += pow(rat[i] - rat2, 2);
  //  }
  //  var = var / rat.size();
  //  float stdv = sqrt(var);
  //  my_metrics.protein_stdv = stdv;


  //  cout << "\n" << "----------- PROTEIN STATS -----------" << "\n" << endl;
  //  cout << "avg percent miss within protein: misscleaved peptides sum intensity / total sum intensity:\t" << my_metrics.protein_final << " ** " << "\n" << endl;
  //  cout << "stdv of percent miss: \t\t\t\t\t\t\t\t\t\t" << my_metrics.protein_stdv << " ** " << "\n" << endl;

  // 
  //  




  //  sort(my_peptide_lists.peptide_matches.begin(), my_peptide_lists.peptide_matches.end(), compareInten);


  //  my_peptide_lists.peptide_matches1.push_back(my_peptide_lists.peptide_matches[0]);
  //  for (size_t i = 1; i < my_peptide_lists.peptide_matches.size(); i++) {
  //      if (my_peptide_lists.peptide_matches[i].ft_pep_seq == my_peptide_lists.peptide_matches[i - 1].ft_pep_seq) continue;
  //      my_peptide_lists.peptide_matches1.push_back(my_peptide_lists.peptide_matches[i]);
  //  }

  //  cout << "\n" << "-- 30 TWENTY MOST ABUNDANT PEPTIDES --" << "\n" << endl;
  //  my_peptide_lists.peptide_matches = my_peptide_lists.peptide_matches1;
  //  if (my_peptide_lists.peptide_matches.size() < 30) {
  //      for (size_t i = 0; i < my_peptide_lists.peptide_matches.size(); i++) {
  //          cout << i + 1 << ": " << " T sequence: " << my_peptide_lists.peptide_matches[i].ft_pep_seq << " || T intensity: " << my_peptide_lists.peptide_matches[i].ft_areaXIC << " || MC sequence: " << my_peptide_lists.peptide_matches[i].mc_pep_seq << " || MC intensity: " << my_peptide_lists.peptide_matches[i].mc_areaXIC << " || intensity ratio: " << my_peptide_lists.peptide_matches[i].ft_areaXIC / my_peptide_lists.peptide_matches[i].mc_areaXIC << endl;
  //      }
  //  }
  //  else {
  //      for (size_t i = 0; i < 30; i++) {
  //          cout << i + 1 << ": " << " T sequence: " << my_peptide_lists.peptide_matches[i].ft_pep_seq << " || T intensity: " << my_peptide_lists.peptide_matches[i].ft_areaXIC << " || MC sequence: " << my_peptide_lists.peptide_matches[i].mc_pep_seq << " || MC intensity: " << my_peptide_lists.peptide_matches[i].mc_areaXIC << " || intensity ratio: " << my_peptide_lists.peptide_matches[i].ft_areaXIC / my_peptide_lists.peptide_matches[i].mc_areaXIC << endl;
  //      }
  //  }
  //  

  //

  // 

  // 


  //  sort(my_peptide_lists.prot_f.begin(), my_peptide_lists.prot_f.end(), compareTotal); 

 


  //  cout << "\n" << "-- 50 TWENTY MOST ABUNDANT PROTEINS--" << "\n" << endl;
  //  
  //  if (my_peptide_lists.prot_f.size() < 50) {
  //      for (size_t i = 0; i < my_peptide_lists.prot_f.size(); i++) {
  //          cout << i + 1 << ":  " << "Prot sequence: " <<  my_peptide_lists.prot_f[i].prot_seq << " || Tot intensity: " << my_peptide_lists.prot_f[i].total << " || amount MC peptides: " << my_peptide_lists.prot_f[i].percentMiss <<  endl;
  //      }
  //  }
  //  else {
  //      for (size_t i = 0; i < 50; i++) {
  //          cout << i + 1 << ":  " << "Prot sequence: " << my_peptide_lists.prot_f[i].prot_seq << " || Tot intensity: " << my_peptide_lists.prot_f[i].total << " || amount MC peptides: " << my_peptide_lists.prot_f[i].percentMiss << endl;
  //      }
  //  }

  //  int count =0; 
  //  for (size_t i = 0; i < my_peptide_lists.prot_f.size(); i++) {
  //      if (my_peptide_lists.prot_f[i].percentMiss == 0 /*&& my_peptide_lists.prot_f[i].trypPeptides.size() >= 5*/) {
  //          count++; 
  //      }
  //  }

  //  cout << "\n" << count << "  proteins out of " << my_peptide_lists.prot_f.size() << " are entirely made up of tryptic peptides" << "\n" << endl; 

  // /* for (int i = 0; i < my_peptide_lists.prot_f.size(); i++) {
  //      cout << my_peptide_lists.prot_f[i].prot_seq << "   " << my_peptide_lists.prot_f[i].trypPeptides.size() << "   " << my_peptide_lists.prot_f[i].missPeptides.size() << "   " << my_peptide_lists.prot_f[i].sumTryp << "   " << my_peptide_lists.prot_f[i].sumMiss << "   " << my_peptide_lists.prot_f[i].percentMiss << endl;
  //  }*/







  //  cout << "------------------------------------------------------------------" << endl; 

  //  cout << "\n" << "*** MOST IMPORTANT DIAGNOSTIC METRICS ***" << "\n" << endl; 

  //  cout << "-- PEPTIDE STATS --" << "\n" << endl; 
  //  cout << "# of misscleaved unique peptides / total # of unique peptides:\t\t\t\t\t" << my_metrics.golden_stat_unique_pep << " ** " << "\n" << endl;
  //  cout << "# of miss cleaved tryptic psm / total tryptic psm above threshold:\t\t\t\t" << my_metrics.golden_stat_psm << " ** " << "\n" << endl;
  //  cout << "avg ratio of intensities: 0 miscleave full tryptic peptides / 1 or 2 miscleave peptides:\t" << my_metrics.intensity_final << " ** " << "\n" << endl;
  //  cout << "stdv of ratio of intensities: \t\t\t\t\t\t\t\t\t" << my_metrics.stdv_final << " ** " << "\n" << endl;
  //  cout << "Global Stat: total miscleaved peptide intensity / sum of all peptide intensities: \t\t" << my_metrics.total_intensity << " ** " << "\n" << endl;

  //  cout << "\n" << "-- PROTEIN STATS --" << "\n" << endl;
  //  cout << "avg percent miss within protein: misscleaved peptides sum intensity / total sum intensity:\t" << my_metrics.protein_final << " ** " << "\n" << endl;
  //  cout << "stdv of percent miss: \t\t\t\t\t\t\t\t\t\t" << my_metrics.protein_stdv << " ** " << "\n" << endl;





  // /* for (size_t i = 0; i < my_peptide_lists.prot_f.size(); i++) {
  //      cout << i + 1 << "-  " << "protein sequence- " << my_peptide_lists.prot_f[i].prot_seq << "    -      total intensity- " << my_peptide_lists.prot_f[i].total << "    -    amount made up of miss cleaved peptides- " << my_peptide_lists.prot_f[i].percentMiss << endl;
  //  }*/




  ///*  for (int i = 0; i < 20; i++) {
  //      cout << my_peptide_lists.prot_f[i].missPeptides.size();
  //      cout << "  " << my_peptide_lists.prot_f[i].trypPeptides.size() << endl;
  //     
  //  }*/



  //
}

void peptide_lists::json(string fn) {
    
  StringBuffer s;
  PrettyWriter<StringBuffer> writer(s);
     
  writer.StartObject();
  writer.Key("Proteins");
  writer.StartArray();
  for (int i = 0; i <all_proteins.size(); i++) {
    writer.StartObject();
    writer.Key("Protein Name");
            
    writer.String(all_proteins[i].prot_seq.c_str());
    writer.Key("Total Intensity");
    writer.Double(all_proteins[i].total);
    writer.Key("Total Enzymatic Peptide Intensity");
    writer.Double(all_proteins[i].sumEnz);
    writer.Key("Total Mis-cleaved Peptide Intensity");
    writer.Double(all_proteins[i].sumMiss);
    writer.Key("Total Nonspecific Peptide Intensity");
    writer.Double(all_proteins[i].sumNonSp);
    writer.Key("Percent Enzymatic");
    writer.Double(all_proteins[i].sumEnz / all_proteins[i].total * 100);
    writer.Key("Percent Mis-Cleaved");
    writer.Double(all_proteins[i].sumMiss / all_proteins[i].total * 100);
    writer.Key("Percent Nonspecific");
    writer.Double(all_proteins[i].sumNonSp / all_proteins[i].total*100);
    writer.Key("Peptides");
    writer.StartArray();
    for (int j = 0; j < all_proteins[i].peptides.size(); j++) {
      size_t index = all_proteins[i].peptides[j];
      writer.StartObject();
      writer.Key("Sequence");
      writer.String(all_peptides[index].pep_seq.c_str());
      writer.Key("Abundance");
      writer.Double(all_peptides[index].areaXIC);
      writer.Key("PSM Count");
      writer.Int(all_peptides[index].psm_count);
      writer.Key("Mis-cleaved");
      if (all_peptides[index].miss_cleaves > 0) writer.Bool(true);
      else writer.Bool(false);
      writer.Key("Nonspecific");
      if (all_peptides[index].non_enzymatic) writer.Bool(true);
      else writer.Bool(false);
      writer.EndObject();
    }
    writer.EndArray();
    writer.EndObject();
  }
  writer.EndArray();
  writer.EndObject();

  FILE* f=fopen(fn.c_str(),"wt");
  fprintf(f,"%s",s.GetString());
  fclose(f);   
 
}

bool peptide_lists::prot_stats() {

  //Sort peptides by protein name - note that peptides cannot be sorted ever again to maintain indexes.
  sort(all_peptides.begin(), all_peptides.end(), compareProt);

  //Make the first proteotypic protein from the first peptide
  size_t i=0;
  while(!all_peptides[i].proteotypic) i++;
  dsProtein prot;
  prot.prot_seq=all_peptides[i].prot_seq;
  all_proteins.push_back(prot);
  all_proteins.back().peptides.push_back(i);

  //iterate over all remaining peptides, combining them into proteins to store in the protein array
  for (i = i+1; i < all_peptides.size(); i++) {

    //skip non-proteotypic peptides
    if (!all_peptides[i].proteotypic) continue;

    //if peptide belongs to current protein, add it to that protein's peptide list
    if(all_peptides[i].prot_seq.compare(all_proteins.back().prot_seq)==0){
      all_proteins.back().peptides.push_back(i);
    
    //otherwise, start a new protein
    } else {
      prot.prot_seq=all_peptides[i].prot_seq;
      all_proteins.push_back(prot);
      all_proteins.back().peptides.push_back(i);
    }

  }


  //iterate over protein list to sum up all peptide signals
  //this could be done above, but is done here for ease of reading the code
  for(i=0;i<all_proteins.size();i++){
    all_proteins[i].sumMiss=0;
    all_proteins[i].sumNonSp=0;
    all_proteins[i].sumEnz=0;
    all_proteins[i].sumPSM=0;
    all_proteins[i].total=0;
    for(size_t j=0;j<all_proteins[i].peptides.size();j++){
      size_t index=all_proteins[i].peptides[j];
      if(all_peptides[index].miss_cleaves>0) all_proteins[i].sumMiss += all_peptides[index].areaXIC;
      if (all_peptides[index].non_enzymatic) all_proteins[i].sumNonSp += all_peptides[index].areaXIC;
      if(all_peptides[index].miss_cleaves==0 && !all_peptides[index].non_enzymatic) all_proteins[i].sumEnz+= all_peptides[index].areaXIC;
      all_proteins[i].total += all_peptides[index].areaXIC;
      all_proteins[i].sumPSM += all_peptides[index].psm_count;
    }
  }

  return true;
}

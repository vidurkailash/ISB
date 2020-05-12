#include "rapidxml.hpp"
#include "deep_class.h"

using namespace rapidxml;
using namespace std;



peptide_lists deep_functions::xml_parse(my_parameters& my_params) {
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

    peptide_lists my_peptide_lists;

    int c = 0;
    int d = 0;
    int e = 0;

    for (spectrum_node; spectrum_node; spectrum_node = spectrum_node->next_sibling())
    {
        if (string(spectrum_node->value()) == "msms_run_summary") {
            continue;
        }

        switch (my_params.ident) {
            //IP
        case 'a':
            for (xml_node<>* sample_node = spectrum_node->first_node("spectrum_query"); sample_node; sample_node = sample_node->next_sibling())
            {
                for (xml_node<>* secondary_node = sample_node->first_node("search_result")->first_node("search_hit")->first_node("analysis_result"); secondary_node; secondary_node = secondary_node->next_sibling())
                {
                    if (/*bool(secondary_node->first_attribute("analysis")) == 1 &&*/ string(secondary_node->first_attribute("analysis")->value()) == "interprophet") {
                        if (atof(secondary_node->first_node("interprophet_result")->first_attribute("probability")->value()) >= 0.00) {
                            d++;
                        }
                        if (atof(secondary_node->first_node("interprophet_result")->first_attribute("probability")->value()) >= my_params.ipro_prob) {
                            my_peptide_lists.all_real.push_back(dsPeptide());
                            my_peptide_lists.all_real[c].pep_seq = string(sample_node->first_node("search_result")->first_node("search_hit")->first_attribute("peptide")->value());
                            my_peptide_lists.all_real[c].charge = atoi(sample_node->first_attribute("assumed_charge")->value());
                            my_peptide_lists.all_real[c].pre_neutral_mass = atof(sample_node->first_attribute("precursor_neutral_mass")->value());
                            my_peptide_lists.all_real[c].xml_rtime = (float)atof(sample_node->first_attribute("retention_time_sec")->value());
                            my_peptide_lists.all_real[c].prev_aa = (string)sample_node->first_node("search_result")->first_node("search_hit")->first_attribute("peptide_prev_aa")->value();
                            my_peptide_lists.all_real[c].next_aa = string(sample_node->first_node("search_result")->first_node("search_hit")->first_attribute("peptide_next_aa")->value());
                            my_peptide_lists.all_real[c].prot_seq = string(sample_node->first_node("search_result")->first_node("search_hit")->first_attribute("protein")->value());
                            my_peptide_lists.all_real[c].calc_neutral_mass = atof(sample_node->first_node("search_result")->first_node("search_hit")->first_attribute("calc_neutral_pep_mass")->value());
                            if (atoi(sample_node->first_node("search_result")->first_node("search_hit")->first_attribute("num_tot_proteins")->value()) == 1) {
                                my_peptide_lists.all_real[c].proteotypic = 1;
                            }
                            c++;
                        }
                    }
                }
            }
            break;
            //PP
        case 'b':
            for (xml_node<>* sample_node = spectrum_node->first_node("spectrum_query"); sample_node; sample_node = sample_node->next_sibling())
            {
                if (atof(sample_node->first_node("search_result")->first_node("search_hit")->first_node("analysis_result")->first_node("peptideprophet_result")->first_attribute("probability")->value()) >= 0.00) {
                    d++;
                }
                if (atof(sample_node->first_node("search_result")->first_node("search_hit")->first_node("analysis_result")->first_node("peptideprophet_result")->first_attribute("probability")->value()) >= my_params.pep_prob) {
                    my_peptide_lists.all_real.push_back(dsPeptide());
                    my_peptide_lists.all_real[c].pep_seq = string(sample_node->first_node("search_result")->first_node("search_hit")->first_attribute("peptide")->value());
                    my_peptide_lists.all_real[c].charge = atoi(sample_node->first_attribute("assumed_charge")->value());
                    my_peptide_lists.all_real[c].pre_neutral_mass = atof(sample_node->first_attribute("precursor_neutral_mass")->value());
                    my_peptide_lists.all_real[c].xml_rtime = (float)atof(sample_node->first_attribute("retention_time_sec")->value());
                    my_peptide_lists.all_real[c].prev_aa = string(sample_node->first_node("search_result")->first_node("search_hit")->first_attribute("peptide_prev_aa")->value());
                    my_peptide_lists.all_real[c].next_aa = string(sample_node->first_node("search_result")->first_node("search_hit")->first_attribute("peptide_next_aa")->value());
                    my_peptide_lists.all_real[c].prot_seq = string(sample_node->first_node("search_result")->first_node("search_hit")->first_attribute("protein")->value());
                    my_peptide_lists.all_real[c].calc_neutral_mass = atof(sample_node->first_node("search_result")->first_node("search_hit")->first_attribute("calc_neutral_pep_mass")->value());
                    if (atoi(sample_node->first_node("search_result")->first_node("search_hit")->first_attribute("num_tot_proteins")->value()) == 1) {
                        my_peptide_lists.all_real[c].proteotypic = 1;
                    }
                    c++;
                }
            }
            break;
        default:
            cout << "error" << endl;
            exit(1);
            break;
        }

    }
    my_peptide_lists.total = d;

    return my_peptide_lists;

}

bool deep_functions::tryptic_calc(peptide_lists& my_peptide_lists, my_parameters& my_params) {

   
    //NEW METHOD
    int c = 0;
    int d = 0;

    for (size_t i = 0; i < my_peptide_lists.all_real.size(); i++) {
        size_t found = my_params.cleave_loc.find(my_peptide_lists.all_real[i].prev_aa);
        size_t found1 = my_params.cleave_loc.find(my_peptide_lists.all_real[i].pep_seq.back());
        if (found!= string::npos && found1!=string::npos) {

            my_peptide_lists.tryptic_real.push_back(dsPeptide());
            my_peptide_lists.tryptic_real[c].pep_seq = my_peptide_lists.all_real[i].pep_seq;
            my_peptide_lists.tryptic_real[c].charge = my_peptide_lists.all_real[i].charge;
            my_peptide_lists.tryptic_real[c].pre_neutral_mass = my_peptide_lists.all_real[i].pre_neutral_mass;
            my_peptide_lists.tryptic_real[c].prev_aa = my_peptide_lists.all_real[i].prev_aa;
            my_peptide_lists.tryptic_real[c].xml_rtime = my_peptide_lists.all_real[i].xml_rtime;
            my_peptide_lists.tryptic_real[c].next_aa = my_peptide_lists.all_real[i].next_aa;
            my_peptide_lists.tryptic_real[c].prot_seq = my_peptide_lists.all_real[i].prot_seq;
            my_peptide_lists.tryptic_real[c].proteotypic = my_peptide_lists.all_real[i].proteotypic;
            my_peptide_lists.tryptic_real[c].calc_neutral_mass = my_peptide_lists.all_real[i].calc_neutral_mass;

            c++;
        }
        else {

            my_peptide_lists.non_tryptic_real.push_back(dsPeptide());
            my_peptide_lists.non_tryptic_real[d].pep_seq = my_peptide_lists.all_real[i].pep_seq;
            my_peptide_lists.non_tryptic_real[d].charge = my_peptide_lists.all_real[i].charge;
            my_peptide_lists.non_tryptic_real[d].pre_neutral_mass = my_peptide_lists.all_real[i].pre_neutral_mass;
            my_peptide_lists.non_tryptic_real[d].prev_aa = my_peptide_lists.all_real[i].prev_aa;
            my_peptide_lists.non_tryptic_real[d].xml_rtime = my_peptide_lists.all_real[i].xml_rtime;
            my_peptide_lists.non_tryptic_real[d].next_aa = my_peptide_lists.all_real[i].next_aa;
            my_peptide_lists.non_tryptic_real[d].prot_seq = my_peptide_lists.all_real[i].prot_seq;
            my_peptide_lists.non_tryptic_real[d].proteotypic = my_peptide_lists.all_real[i].proteotypic;
            my_peptide_lists.non_tryptic_real[d].calc_neutral_mass = my_peptide_lists.all_real[i].calc_neutral_mass;

            d++;
        }
    }



   
    //OLD METHOD
    /*for (size_t i = 0; i < my_peptide_lists.all_real.size(); i++) {
        if ((my_peptide_lists.all_real[i].prev_aa == "R" || my_peptide_lists.all_real[i].prev_aa == "K" || my_peptide_lists.all_real[i].prev_aa == "-") &&
            (my_peptide_lists.all_real[i].pep_seq.back() == 'R' || my_peptide_lists.all_real[i].pep_seq.back() == 'K' || my_peptide_lists.all_real[i].next_aa == "-")) {

            my_peptide_lists.tryptic_real.push_back(dsPeptide());
            my_peptide_lists.tryptic_real[c].pep_seq = my_peptide_lists.all_real[i].pep_seq;
            my_peptide_lists.tryptic_real[c].charge = my_peptide_lists.all_real[i].charge;
            my_peptide_lists.tryptic_real[c].pre_neutral_mass = my_peptide_lists.all_real[i].pre_neutral_mass;
            my_peptide_lists.tryptic_real[c].prev_aa = my_peptide_lists.all_real[i].prev_aa;
            my_peptide_lists.tryptic_real[c].xml_rtime = my_peptide_lists.all_real[i].xml_rtime;
            my_peptide_lists.tryptic_real[c].next_aa = my_peptide_lists.all_real[i].next_aa;
            my_peptide_lists.tryptic_real[c].prot_seq = my_peptide_lists.all_real[i].prot_seq;
            my_peptide_lists.tryptic_real[c].proteotypic = my_peptide_lists.all_real[i].proteotypic;
            my_peptide_lists.tryptic_real[c].calc_neutral_mass = my_peptide_lists.all_real[i].calc_neutral_mass;

            c++;
        }
        else {

            my_peptide_lists.non_tryptic_real.push_back(dsPeptide());
            my_peptide_lists.non_tryptic_real[d].pep_seq = my_peptide_lists.all_real[i].pep_seq;
            my_peptide_lists.non_tryptic_real[d].charge = my_peptide_lists.all_real[i].charge;
            my_peptide_lists.non_tryptic_real[d].pre_neutral_mass = my_peptide_lists.all_real[i].pre_neutral_mass;
            my_peptide_lists.non_tryptic_real[d].prev_aa = my_peptide_lists.all_real[i].prev_aa;
            my_peptide_lists.non_tryptic_real[d].xml_rtime = my_peptide_lists.all_real[i].xml_rtime;
            my_peptide_lists.non_tryptic_real[d].next_aa = my_peptide_lists.all_real[i].next_aa;
            my_peptide_lists.non_tryptic_real[d].prot_seq = my_peptide_lists.all_real[i].prot_seq;
            my_peptide_lists.non_tryptic_real[d].proteotypic = my_peptide_lists.all_real[i].proteotypic;
            my_peptide_lists.non_tryptic_real[d].calc_neutral_mass = my_peptide_lists.all_real[i].calc_neutral_mass;

            d++;
        }
    }*/




    return true;

}

bool deep_functions::miss_cleave(peptide_lists& my_peptide_lists, my_parameters& my_params) {

   //NEW METHOD
    int c = 0;
    for (int i = 0; i < my_peptide_lists.tryptic_real.size(); i++) {
        for (int j = 0; j < my_peptide_lists.tryptic_real[i].pep_seq.size() - 1; j++) {
            size_t found = my_params.cleave_loc.find(my_peptide_lists.tryptic_real[i].pep_seq[j]);
            size_t found1 = my_params.anti_cleave_loc.find(my_peptide_lists.tryptic_real[i].pep_seq[j + 1]);
            if (found!=string::npos && found1==string::npos) {
                c++;
            }
        }
        my_peptide_lists.tryptic_real[i].miss_cleaves = c;
        c = 0;
    }
    
    
    //OLD METHOD 
    /*int c = 0;
    for (int i = 0; i < my_peptide_lists.tryptic_real.size(); i++) {
        for (int j = 0; j < my_peptide_lists.tryptic_real[i].pep_seq.size() - 1; j++) {
            int k = j + 1;
            if ((my_peptide_lists.tryptic_real[i].pep_seq[j] == 'K' || my_peptide_lists.tryptic_real[i].pep_seq[j] == 'R') && my_peptide_lists.tryptic_real[i].pep_seq[k] != 'P') {
                c++;
            }
        }
        my_peptide_lists.tryptic_real[i].miss_cleaves = c;
        c = 0;
    }*/

    return true;

}

bool deep_functions::delete_dup(peptide_lists& my_peptide_lists) {


    size_t i;

    //Sort the tryptic array by sequence and then charge state
    sort(my_peptide_lists.tryptic_real.begin(), my_peptide_lists.tryptic_real.end(), compareSeqZ);

    //Iterate over all tryptic PSMs, keeping each unique instance of sequence and charge
    my_peptide_lists.tryp_unique_z_real.push_back(my_peptide_lists.tryptic_real[0]);
    for (i = 1; i < my_peptide_lists.tryptic_real.size(); i++) {
        if (my_peptide_lists.tryp_unique_z_real.back().pep_seq.compare(my_peptide_lists.tryptic_real[i].pep_seq) == 0 &&
            my_peptide_lists.tryp_unique_z_real.back().charge == my_peptide_lists.tryptic_real[i].charge) continue;
        my_peptide_lists.tryp_unique_z_real.push_back(my_peptide_lists.tryptic_real[i]);
    }

    //Repeat the process, keeping each unique instance of sequence;
    my_peptide_lists.tryp_unique_real.push_back(my_peptide_lists.tryptic_real[0]);
    for (i = 1; i < my_peptide_lists.tryptic_real.size(); i++) {
        if (my_peptide_lists.tryp_unique_real.back().pep_seq.compare(my_peptide_lists.tryptic_real[i].pep_seq) == 0) continue;
        my_peptide_lists.tryp_unique_real.push_back(my_peptide_lists.tryptic_real[i]);
    }

    //Now of the unique peptide sequences, make a subset of miscleaved peptides
    i = 0;
    while (my_peptide_lists.tryp_unique_real[i].miss_cleaves == 0) i++; //no boundary checks here
    my_peptide_lists.miss_unique_real.push_back(my_peptide_lists.tryp_unique_real[i]);
    for (i = i + 1; i < my_peptide_lists.tryp_unique_real.size(); i++) {
        if (my_peptide_lists.tryp_unique_real[i].miss_cleaves == 0) continue;
        my_peptide_lists.miss_unique_real.push_back(my_peptide_lists.tryp_unique_real[i]);
    }

    for (i = 0; i < my_peptide_lists.tryp_unique_real.size(); i++) {
        if (my_peptide_lists.tryp_unique_real[i].miss_cleaves != 0) continue;
        my_peptide_lists.fully_tryp_unique_real.push_back(my_peptide_lists.tryp_unique_real[i]);
    }



    return true;

}

bool deep_functions::lcd(peptide_lists& my_peptide_lists, my_parameters& my_params) {


    //MH: This new method finds all miscleavages for each fully cleaved peptide, and stores them in a
    //new array.


    for (size_t i = 0; i < my_peptide_lists.fully_tryp_unique_real.size(); i++) {
        my_peptide_lists.fully_tryp_unique_real[i].xml_mz = (my_peptide_lists.fully_tryp_unique_real[i].pre_neutral_mass + (my_peptide_lists.fully_tryp_unique_real[i].charge * 1.00727)) / my_peptide_lists.fully_tryp_unique_real[i].charge;
        my_peptide_lists.fully_tryp_unique_real[i].tolerance = (my_params.ppm / (1000000)) * my_peptide_lists.fully_tryp_unique_real[i].calc_neutral_mass;
    }

    for (size_t i = 0; i < my_peptide_lists.miss_unique_real.size(); i++) {
        my_peptide_lists.miss_unique_real[i].xml_mz = (my_peptide_lists.miss_unique_real[i].pre_neutral_mass + (my_peptide_lists.miss_unique_real[i].charge * 1.00727)) / my_peptide_lists.miss_unique_real[i].charge;
        my_peptide_lists.miss_unique_real[i].tolerance = (my_params.ppm / (1000000)) * my_peptide_lists.miss_unique_real[i].calc_neutral_mass;
    }





    //size_t i,j;
    //for(i=0;i<my_peptide_lists.tryp_unique.size();i++){

    //  //MH: skip any peptides that have miscleavages
    //  if(my_peptide_lists.tryp_unique[i].miss_cleaves>0) continue;

    //  for(j=0;j<my_peptide_lists.miss_unique.size();j++){

    //    //MH: we have a match, make a new entry
    //    size_t found = my_peptide_lists.miss_unique[j].pep_seq.find(my_peptide_lists.tryp_unique[i].pep_seq);
    //    if(found != string::npos && my_peptide_lists.miss_unique[j].pep_seq.compare(my_peptide_lists.tryp_unique[i].pep_seq)!=0) {
    //      my_features f= my_peptide_lists.miss_unique[j];
    //      f.d_pep_seq = my_peptide_lists.tryp_unique[i].pep_seq;
    //      f.d_pep_seq_rt = my_peptide_lists.tryp_unique[i].rtime;
    //      f.d_miss_cleaves = my_peptide_lists.tryp_unique[i].miss_cleaves;
    //      f.d_pep_seq_mass = my_peptide_lists.tryp_unique[i].mass;
    //      f.d_pep_seq_charge = my_peptide_lists.tryp_unique[i].charge;
    //      if (found == 0) {
    //        f.cleave_loc = 'R';
    //        f.cleave_pos = 0;
    //      } else {
    //        f.cleave_loc = 'L';
    //        f.cleave_pos = (int)found;
    //      }
    //      f.mz = (f.mass + (f.charge * 1.00727)) / f.charge;
    //      f.d_pep_seq_mz = (f.d_pep_seq_mass + (f.d_pep_seq_charge * 1.00727)) / f.d_pep_seq_charge;
    //      my_peptide_lists.d_list.push_back(f);
    //    }

    //  }

    //}


    return true; 

}

bool deep_functions::reader(peptide_lists& my_peptide_lists, metrics& my_metrics) {

    dsXIC x;
    int a = 0; 
    int b = 0; 
    for (size_t i = 0; i < my_peptide_lists.xic_ft_results.size(); i++) {
        if (my_peptide_lists.xic_ft_results[i].pep_seq == my_peptide_lists.xic_ft_results[i+1].pep_seq) {
            x.rTime = my_peptide_lists.xic_ft_results[i].spec_rt;
            x.intensity = my_peptide_lists.xic_ft_results[i].spec_intensity;
            my_peptide_lists.xic_ft_results[a].XIC.push_back(x);
        }
        else if (my_peptide_lists.xic_ft_results[i].pep_seq != my_peptide_lists.xic_ft_results[i + 1].pep_seq) {
            x.rTime = my_peptide_lists.xic_ft_results[i].spec_rt;
            x.intensity = my_peptide_lists.xic_ft_results[i].spec_intensity;
            my_peptide_lists.xic_ft_results[a].XIC.push_back(x);
            a = i + 1; 
        }
      
    }
    vector<dsPeptide> tmp;
    tmp.push_back(my_peptide_lists.xic_ft_results[0]);
    for (size_t i = 1; i < my_peptide_lists.xic_ft_results.size(); i++) {
        if (my_peptide_lists.xic_ft_results[i].pep_seq == my_peptide_lists.xic_ft_results[i - 1].pep_seq) continue;
        tmp.push_back(my_peptide_lists.xic_ft_results[i]);
    }
    my_peptide_lists.xic_ft_results = tmp; 
    tmp.clear();



    for (size_t i = 0; i < my_peptide_lists.xic_mc_results.size(); i++) {
        if (my_peptide_lists.xic_mc_results[i].pep_seq == my_peptide_lists.xic_mc_results[i + 1].pep_seq) {
            x.rTime = my_peptide_lists.xic_mc_results[i].spec_rt;
            x.intensity = my_peptide_lists.xic_mc_results[i].spec_intensity;
            my_peptide_lists.xic_mc_results[b].XIC.push_back(x);
        }
        else if (my_peptide_lists.xic_mc_results[i].pep_seq != my_peptide_lists.xic_mc_results[i + 1].pep_seq) {
            x.rTime = my_peptide_lists.xic_mc_results[i].spec_rt;
            x.intensity = my_peptide_lists.xic_mc_results[i].spec_intensity;
            my_peptide_lists.xic_mc_results[b].XIC.push_back(x);
            b = i + 1;
        }

    }
    tmp.push_back(my_peptide_lists.xic_mc_results[0]);
    for (size_t i = 1; i < my_peptide_lists.xic_mc_results.size(); i++) {
        if (my_peptide_lists.xic_mc_results[i].pep_seq == my_peptide_lists.xic_mc_results[i - 1].pep_seq) continue;
        tmp.push_back(my_peptide_lists.xic_mc_results[i]);
    }
    my_peptide_lists.xic_mc_results = tmp;





    //DELETE NOISE (INTENSITIES LESS THAN 10% OF THE MAX) (NOISE)

    for (size_t i = 0; i < my_peptide_lists.xic_ft_results.size(); i++) {
        if (!cleanNoise(my_peptide_lists.xic_ft_results[i].XIC)) {
            //handle error
        }
    }
    for (size_t i = 0; i < my_peptide_lists.xic_mc_results.size(); i++) {
        if (!cleanNoise(my_peptide_lists.xic_mc_results[i].XIC)) {
            //handle error
        }
    }

   


    //delete peptides with only one data point
    vector<dsPeptide> filter;
    for (size_t i = 0; i < my_peptide_lists.xic_ft_results.size(); i++) {
        if (my_peptide_lists.xic_ft_results[i].XIC.size() > 1) {
            filter.push_back(my_peptide_lists.xic_ft_results[i]);
        }
    }
    my_peptide_lists.xic_ft_results = filter;
    filter.clear();
    for (size_t i = 0; i < my_peptide_lists.xic_mc_results.size(); i++) {
        if (my_peptide_lists.xic_mc_results[i].XIC.size() > 1) {
            filter.push_back(my_peptide_lists.xic_mc_results[i]);
        }
    }
    my_peptide_lists.xic_mc_results = filter;
    filter.clear();



    cout << "noise deleted" << "\n" << endl;


    //metric calculations
    vector<float> peak;
    float avg = 0;
    for (int i = 0; i < my_peptide_lists.xic_ft_results.size(); i++) {
        for (int j = 0; j < my_peptide_lists.xic_ft_results[i].XIC.size(); j++) {
            peak.push_back(my_peptide_lists.xic_ft_results[i].XIC[j].intensity);
        }
        avg += *max_element(peak.begin(), peak.end());
        peak.clear();
    }


    my_metrics.tryp_avg_high = avg / my_peptide_lists.xic_ft_results.size();

    peak.clear();
    avg = 0;
    for (int i = 0; i < my_peptide_lists.xic_mc_results.size(); i++) {
        for (int j = 0; j < my_peptide_lists.xic_mc_results[i].XIC.size(); j++) {
            peak.push_back(my_peptide_lists.xic_mc_results[i].XIC[j].intensity);
        }
        avg += *max_element(peak.begin(), peak.end());
        peak.clear();
    }


    my_metrics.miss_avg_high = avg / my_peptide_lists.xic_mc_results.size();

   
    float big = 0; 
    float big1 = 0; 
    //CALC AREA 
    for (size_t i = 0; i < my_peptide_lists.xic_ft_results.size(); i++) {
        my_peptide_lists.xic_ft_results[i].areaXIC = calcPeakArea(my_peptide_lists.xic_ft_results[i].XIC);
        big += my_peptide_lists.xic_ft_results[i].areaXIC; 
    }
    for (size_t i = 0; i < my_peptide_lists.xic_mc_results.size(); i++) {
        my_peptide_lists.xic_mc_results[i].areaXIC = calcPeakArea(my_peptide_lists.xic_mc_results[i].XIC);
        big1 += my_peptide_lists.xic_mc_results[i].areaXIC; 
    }

    my_metrics.total_intensity = big1 / (big + big1); 
   
    // END AREA FUNCTION


    cout << "intensity calculated" << "\n" << endl;



    for (size_t i = 0; i < my_peptide_lists.xic_ft_results.size(); i++) {

    }



    return true;

}

//XML

//XML

//peptide_lists deep_functions::xml_parse(my_parameters& my_params) {
//	xml_document<> doc;
//
//	// Read the xml file into a vector
//	ifstream theFile(my_params.filename);
//	vector<char> buffer((istreambuf_iterator<char>(theFile)), istreambuf_iterator<char>());
//	buffer.push_back('\0');
//	// Parse the buffer using the xml file parsing library into doc 
//	doc.parse<0>(&buffer[0]);
//	// Find root node
//	xml_node<>* root_node = doc.first_node("msms_pipeline_analysis");
//	xml_node<>* spectrum_node = root_node->first_node("msms_run_summary");
//
//	peptide_lists my_peptide_lists;
//
//	int c = 0;
//	int d = 0;
//    int e = 0; 
//
//	for (spectrum_node; spectrum_node; spectrum_node = spectrum_node->next_sibling())
//	{
//		if (string(spectrum_node->value()) == "msms_run_summary") {
//			continue;
//		}
//
//		switch (my_params.ident) {
//		//IP
//		case 'a':
//			for (xml_node<>* sample_node = spectrum_node->first_node("spectrum_query"); sample_node; sample_node = sample_node->next_sibling())
//			{
//				for (xml_node<>* secondary_node = sample_node->first_node("search_result")->first_node("search_hit")->first_node("analysis_result"); secondary_node; secondary_node = secondary_node->next_sibling())
//				{
//					if (/*bool(secondary_node->first_attribute("analysis")) == 1 &&*/ string(secondary_node->first_attribute("analysis")->value()) == "interprophet") {
//						if (atof(secondary_node->first_node("interprophet_result")->first_attribute("probability")->value()) >= 0.00) {
//							d++;
//						}
//						if (atof(secondary_node->first_node("interprophet_result")->first_attribute("probability")->value()) >= my_params.ipro_prob) {
//							my_peptide_lists.all.push_back(my_features());
//							my_peptide_lists.all[c].pep_seq = string(sample_node->first_node("search_result")->first_node("search_hit")->first_attribute("peptide")->value());
//							my_peptide_lists.all[c].charge = atoi(sample_node->first_attribute("assumed_charge")->value());
//							my_peptide_lists.all[c].mass = atof(sample_node->first_attribute("precursor_neutral_mass")->value());
//							my_peptide_lists.all[c].rtime = (float)atof(sample_node->first_attribute("retention_time_sec")->value());
//							my_peptide_lists.all[c].prev_aa = string(sample_node->first_node("search_result")->first_node("search_hit")->first_attribute("peptide_prev_aa")->value());
//							my_peptide_lists.all[c].next_aa = string(sample_node->first_node("search_result")->first_node("search_hit")->first_attribute("peptide_next_aa")->value());
//                            my_peptide_lists.all[c].prot_seq = string(sample_node->first_node("search_result")->first_node("search_hit")->first_attribute("protein")->value());
//                            my_peptide_lists.all[c].calc_mass = atof(sample_node->first_node("search_result")->first_node("search_hit")->first_attribute("calc_neutral_pep_mass")->value());
//                            if (atoi(sample_node->first_node("search_result")->first_node("search_hit")->first_attribute("num_tot_proteins")->value()) == 1) {
//                                my_peptide_lists.all[c].proteotypic = 1;
//                            }
//                            c++; 
//						}
//					}
//				}
//			}
//			break; 
//		//PP
//		case 'b':
//			for (xml_node<>* sample_node = spectrum_node->first_node("spectrum_query"); sample_node; sample_node = sample_node->next_sibling())
//			{
//				if (atof(sample_node->first_node("search_result")->first_node("search_hit")->first_node("analysis_result")->first_node("peptideprophet_result")->first_attribute("probability")->value()) >= 0.00) {
//					d++;
//				}
//				if (atof(sample_node->first_node("search_result")->first_node("search_hit")->first_node("analysis_result")->first_node("peptideprophet_result")->first_attribute("probability")->value()) >= my_params.pep_prob) {
//					my_peptide_lists.all.push_back(my_features());
//					my_peptide_lists.all[c].pep_seq = string(sample_node->first_node("search_result")->first_node("search_hit")->first_attribute("peptide")->value());
//					my_peptide_lists.all[c].charge = atoi(sample_node->first_attribute("assumed_charge")->value());
//					my_peptide_lists.all[c].mass = atof(sample_node->first_attribute("precursor_neutral_mass")->value());
//					my_peptide_lists.all[c].rtime = (float)atof(sample_node->first_attribute("retention_time_sec")->value());
//					my_peptide_lists.all[c].prev_aa = string(sample_node->first_node("search_result")->first_node("search_hit")->first_attribute("peptide_prev_aa")->value());
//					my_peptide_lists.all[c].next_aa = string(sample_node->first_node("search_result")->first_node("search_hit")->first_attribute("peptide_next_aa")->value());
//                    my_peptide_lists.all[c].prot_seq = string(sample_node->first_node("search_result")->first_node("search_hit")->first_attribute("protein")->value());
//                    my_peptide_lists.all[c].calc_mass = atof(sample_node->first_node("search_result")->first_node("search_hit")->first_attribute("calc_neutral_pep_mass")->value());
//                    if (atoi(sample_node->first_node("search_result")->first_node("search_hit")->first_attribute("num_tot_proteins")->value()) == 1) {
//                        my_peptide_lists.all[c].proteotypic = 1;
//                    }
//					c++;
//				}
//			}
//			break;
//		default:
//			cout << "error" << endl; 
//			exit(1);
//			break; 
//		}
//
//	}
//	my_peptide_lists.total = d;
//                                  
//	return my_peptide_lists;
//
//}

//TRYPTIC_CALC

//bool deep_functions::tryptic_calc(peptide_lists& my_peptide_lists) {
//
//	int c = 0;
//	int d = 0;
//
//	for (int i = 0; i < my_peptide_lists.all.size(); i++) {
//
//		if ((my_peptide_lists.all[i].prev_aa == "R" || my_peptide_lists.all[i].prev_aa == "K" || my_peptide_lists.all[i].prev_aa == "-") &&
//			(my_peptide_lists.all[i].pep_seq.back() == 'R' || my_peptide_lists.all[i].pep_seq.back() == 'K' || my_peptide_lists.all[i].next_aa == "-")) {
//
//			my_peptide_lists.tryptic.push_back(my_features());
//			my_peptide_lists.tryptic[c].pep_seq = my_peptide_lists.all[i].pep_seq;
//			my_peptide_lists.tryptic[c].charge = my_peptide_lists.all[i].charge;
//			my_peptide_lists.tryptic[c].mass = my_peptide_lists.all[i].mass;
//			my_peptide_lists.tryptic[c].prev_aa = my_peptide_lists.all[i].prev_aa;
//			my_peptide_lists.tryptic[c].rtime = my_peptide_lists.all[i].rtime;
//			my_peptide_lists.tryptic[c].next_aa = my_peptide_lists.all[i].next_aa;
//            my_peptide_lists.tryptic[c].prot_seq = my_peptide_lists.all[i].prot_seq;
//            my_peptide_lists.tryptic[c].proteotypic = my_peptide_lists.all[i].proteotypic;
//            my_peptide_lists.tryptic[c].calc_mass = my_peptide_lists.all[i].calc_mass;
//
//			c++;
//		}
//		else {
//
//			my_peptide_lists.non_tryptic.push_back(my_features());
//			my_peptide_lists.non_tryptic[d].pep_seq = my_peptide_lists.all[i].pep_seq;
//			my_peptide_lists.non_tryptic[d].charge = my_peptide_lists.all[i].charge;
//			my_peptide_lists.non_tryptic[d].mass = my_peptide_lists.all[i].mass;
//			my_peptide_lists.non_tryptic[d].prev_aa = my_peptide_lists.all[i].prev_aa;
//			my_peptide_lists.non_tryptic[d].rtime = my_peptide_lists.all[i].rtime;
//			my_peptide_lists.non_tryptic[d].next_aa = my_peptide_lists.all[i].next_aa;
//            my_peptide_lists.non_tryptic[d].prot_seq = my_peptide_lists.all[i].prot_seq;
//            my_peptide_lists.non_tryptic[d].proteotypic = my_peptide_lists.all[i].proteotypic;
//            my_peptide_lists.non_tryptic[d].calc_mass = my_peptide_lists.all[i].calc_mass;
//
//			d++;
//		}
//		
//	}
//
//	return true;
//
//}

//MISS_CLEAVE

//bool deep_functions::miss_cleave(peptide_lists& my_peptide_lists) {
//
//	int c = 0;
//	for (int i = 0; i < my_peptide_lists.tryptic.size(); i++) {
//		for (int j = 0; j < my_peptide_lists.tryptic[i].pep_seq.size() - 1; j++) {
//			int k = j + 1;
//			if ((my_peptide_lists.tryptic[i].pep_seq[j] == 'K' || my_peptide_lists.tryptic[i].pep_seq[j] == 'R') && my_peptide_lists.tryptic[i].pep_seq[k] != 'P') {
//				c++;
//			}
//		}
//		my_peptide_lists.tryptic[i].miss_cleaves = c;
//		c = 0;
//	}
//
//	return true;
//
//}

//DELETE_DUP

//bool deep_functions::delete_dup(peptide_lists& my_peptide_lists) {
//
//  //MH: I think this function is supposed to find the unique set of peptide sequences
//  //from the complete list of tryptic PSMs. delete_dup is a misnomer. Nothing is deleted.
//  //Instead new arrays of unique sequences are created.
//  size_t i;
//
//  //MH: Sort the tryptic array by sequence and then charge state
//  sort(my_peptide_lists.tryptic.begin(),my_peptide_lists.tryptic.end(),compareSeqZ);
//
//  //MH: Iterate over all tryptic PSMs, keeping each unique instance of sequence and charge
//  my_peptide_lists.tryp_unique_z.push_back(my_peptide_lists.tryptic[0]);
//  for(i=1;i<my_peptide_lists.tryptic.size();i++){
//    if(my_peptide_lists.tryp_unique_z.back().pep_seq.compare(my_peptide_lists.tryptic[i].pep_seq)==0 &&
//      my_peptide_lists.tryp_unique_z.back().charge == my_peptide_lists.tryptic[i].charge) continue;
//    my_peptide_lists.tryp_unique_z.push_back(my_peptide_lists.tryptic[i]);
//  }
//
//  //MH: Repeat the process, keeping each unique instance of sequence;
//  my_peptide_lists.tryp_unique.push_back(my_peptide_lists.tryptic[0]);
//  for (i = 1; i < my_peptide_lists.tryptic.size(); i++) {
//    if (my_peptide_lists.tryp_unique.back().pep_seq.compare(my_peptide_lists.tryptic[i].pep_seq) == 0) continue;
//    my_peptide_lists.tryp_unique.push_back(my_peptide_lists.tryptic[i]);
//  }
//
//  //MH: Now of the unique peptide sequences, make a subset of miscleaved peptides
//  i=0;
//  while(my_peptide_lists.tryp_unique[i].miss_cleaves==0) i++; //no boundary checks here
//  my_peptide_lists.miss_unique.push_back(my_peptide_lists.tryp_unique[i]);
//  for (i = i+1; i < my_peptide_lists.tryp_unique.size(); i++) {
//    if (my_peptide_lists.tryp_unique[i].miss_cleaves==0) continue;
//    my_peptide_lists.miss_unique.push_back(my_peptide_lists.tryp_unique[i]);
//  }
//
//  for (i = 0; i < my_peptide_lists.tryp_unique.size(); i++) {
//      if (my_peptide_lists.tryp_unique[i].miss_cleaves != 0) continue; 
//      my_peptide_lists.fully_tryp_unique.push_back(my_peptide_lists.tryp_unique[i]); 
//  }
//
//
//  
//
//
//
//
//	return true;
//
//}

//LCD

//bool deep_functions::lcd(peptide_lists& my_peptide_lists) {
//	
//  
//  //MH: This new method finds all miscleavages for each fully cleaved peptide, and stores them in a
//  //new array.
//
//
//    for (size_t i = 0; i < my_peptide_lists.fully_tryp_unique.size(); i++) {
//        my_peptide_lists.fully_tryp_unique[i].mz = (my_peptide_lists.fully_tryp_unique[i].mass + (my_peptide_lists.fully_tryp_unique[i].charge * 1.00727)) / my_peptide_lists.fully_tryp_unique[i].charge;
//    }
//
//    for (size_t i = 0; i < my_peptide_lists.miss_unique.size(); i++) {
//        my_peptide_lists.miss_unique[i].mz = (my_peptide_lists.miss_unique[i].mass + (my_peptide_lists.miss_unique[i].charge * 1.00727)) / my_peptide_lists.miss_unique[i].charge;
//    }
//
//
//
//  //size_t i,j;
//  //for(i=0;i<my_peptide_lists.tryp_unique.size();i++){
//
//  //  //MH: skip any peptides that have miscleavages
//  //  if(my_peptide_lists.tryp_unique[i].miss_cleaves>0) continue;
//
//  //  for(j=0;j<my_peptide_lists.miss_unique.size();j++){
//
//  //    //MH: we have a match, make a new entry
//  //    size_t found = my_peptide_lists.miss_unique[j].pep_seq.find(my_peptide_lists.tryp_unique[i].pep_seq);
//  //    if(found != string::npos && my_peptide_lists.miss_unique[j].pep_seq.compare(my_peptide_lists.tryp_unique[i].pep_seq)!=0) {
//  //      my_features f= my_peptide_lists.miss_unique[j];
//  //      f.d_pep_seq = my_peptide_lists.tryp_unique[i].pep_seq;
//  //      f.d_pep_seq_rt = my_peptide_lists.tryp_unique[i].rtime;
//  //      f.d_miss_cleaves = my_peptide_lists.tryp_unique[i].miss_cleaves;
//  //      f.d_pep_seq_mass = my_peptide_lists.tryp_unique[i].mass;
//  //      f.d_pep_seq_charge = my_peptide_lists.tryp_unique[i].charge;
//  //      if (found == 0) {
//  //        f.cleave_loc = 'R';
//  //        f.cleave_pos = 0;
//  //      } else {
//  //        f.cleave_loc = 'L';
//  //        f.cleave_pos = (int)found;
//  //      }
//  //      f.mz = (f.mass + (f.charge * 1.00727)) / f.charge;
//  //      f.d_pep_seq_mz = (f.d_pep_seq_mass + (f.d_pep_seq_charge * 1.00727)) / f.d_pep_seq_charge;
//  //      my_peptide_lists.d_list.push_back(f);
//  //    }
//
//  //  }
//
//  //}
//
//
//  
//
//
//
//  
//  //MH: I think this approach is flawed. MISCLVDRSEQENCEK is a miscleavage of both MISCLVDR and SEQENCEK
//  //Therefore it isn't possible to have only a single tryptic representation for each miscleavage.
//	///*for (int i = 0; i < my_peptide_lists.miss_unique.size(); i++) {
//	//	for (int j = 0; j < my_peptide_lists.tryp_unique.size(); j++) {
//	//		size_t found = my_peptide_lists.miss_unique[i].pep_seq.find(my_peptide_lists.tryp_unique[j].pep_seq);
//	//		if (found != string::npos && (my_peptide_lists.miss_unique[i].pep_seq != my_peptide_lists.tryp_unique[j].pep_seq)) {
//	//			my_peptide_lists.miss_unique[i].d_pep_seq = my_peptide_lists.tryp_unique[j].pep_seq;
//	//			my_peptide_lists.miss_unique[i].d_pep_seq_rt = my_peptide_lists.tryp_unique[j].rtime;
//	//			my_peptide_lists.miss_unique[i].d_miss_cleaves = my_peptide_lists.tryp_unique[j].miss_cleaves;
//	//			if (found == 0) {
//	//				my_peptide_lists.miss_unique[i].cleave_loc = 'R';
//	//				my_peptide_lists.miss_unique[i].cleave_pos = (int)found;
//
//	//			}
//	//			else {
//	//				my_peptide_lists.miss_unique[i].cleave_loc = 'L';
//	//				my_peptide_lists.miss_unique[i].cleave_pos = (int)found;
//
//	//			}
//	//			my_peptide_lists.miss_unique[i].d_pep_seq_mass = my_peptide_lists.tryp_unique[j].mass;
//	//			my_peptide_lists.miss_unique[i].d_pep_seq_charge = my_peptide_lists.tryp_unique[j].charge;
//	//		}
//
//	//	}
//
//	//	my_peptide_lists.miss_unique[i].mz = (my_peptide_lists.miss_unique[i].mass + ((my_peptide_lists.miss_unique[i].charge) * 1.00727)) / my_peptide_lists.miss_unique[i].charge;
//	//	my_peptide_lists.miss_unique[i].d_pep_seq_mz = (my_peptide_lists.miss_unique[i].d_pep_seq_mass + ((my_peptide_lists.miss_unique[i].d_pep_seq_charge) * 1.00727)) / my_peptide_lists.miss_unique[i].d_pep_seq_charge;
//	//}*/
//
//
//
//
//	/*for (int i = 0; i < my_peptide_lists.tryp_unique.size(); i++) {
//		cout << my_peptide_lists.tryp_unique[i].pep_seq << "  " << my_peptide_lists.tryp_unique[i].charge << "   " << my_peptide_lists.tryp_unique[i].mass << "  " << my_peptide_lists.tryp_unique[i].mz << "  " << my_peptide_lists.tryp_unique[i].cleave_loc << "   " << my_peptide_lists.tryp_unique[i].cleave_pos  << "   " << my_peptide_lists.tryp_unique[i].d_pep_seq << "  " << my_peptide_lists.tryp_unique[i].d_pep_seq_charge << "  " << my_peptide_lists.tryp_unique[i].d_pep_seq_mass << "  " << my_peptide_lists.tryp_unique[i].d_pep_seq_mz << endl;
//	}*/
//
//	/*cout << my_peptide_lists.tryp_unique.size() << endl;*/
//
//
//
//	return true; 
//
//}

//NEW_LIST -> OBSOLETE

//peptide_lists deep_functions::new_list(peptide_lists& my_peptide_lists) {
//
//	my_peptide_lists.d_list = my_peptide_lists.miss_unique; 
//	
//	for (int i = 0; i < my_peptide_lists.d_list.size(); i++) {
//		if (my_peptide_lists.d_list[i].d_pep_seq.size() == 0 ){
//			my_peptide_lists.d_list.erase(my_peptide_lists.d_list.begin() + i);
//			i = i - 1;
//		}
//
//	}
//
//
//
//	/*for (int i = 0; i < my_peptide_lists.d_list.size(); i++) {
//		cout << my_peptide_lists.d_list[i].pep_seq << "  " << my_peptide_lists.d_list[i].charge << "   " << my_peptide_lists.d_list[i].mass << "  " << my_peptide_lists.d_list[i].mz << "  " << my_peptide_lists.d_list[i].cleave_loc << "   " << my_peptide_lists.d_list[i].cleave_pos << " ||||| " << my_peptide_lists.d_list[i].d_pep_seq << "  " << my_peptide_lists.d_list[i].d_pep_seq_charge << "  " << my_peptide_lists.d_list[i].d_pep_seq_mass << "  " << my_peptide_lists.d_list[i].d_pep_seq_mz << endl;
//	}*/
//
//
//	cout << my_peptide_lists.d_list.size() << endl;
//
//	return my_peptide_lists;
//
//}

//READER

//bool deep_functions::reader(peptide_lists& my_peptide_lists, metrics& my_metrics, match_lists& my_match_lists) {
//
//  vector<my_intensities> qc;
//
//  //MH: Again, sorting is your friend here. Although these vectors are already sorted, so we won't do that.
//  //But because they are sorted, you can get through in a single pass. And you won't need to erase
//  //which will save tons of time.
//  for (size_t i = 0; i < my_match_lists.results.size(); i++) {
//    my_peptide_lists.spectra.push_back(my_intensities()); //add first one
//    my_peptide_lists.spectra.back().x = my_match_lists.results[i].spec_rt;
//    my_peptide_lists.spectra.back().y = my_match_lists.results[i].spec_intensity;
//    my_peptide_lists.spectra.back().seq = my_match_lists.results[i].pep_seq;
//    my_peptide_lists.spectra.back().mc = my_match_lists.results[i].miss_cleaves;
//    my_peptide_lists.spectra.back().prot_seq = my_match_lists.results[i].prot_seq;
//    my_peptide_lists.spectra.back().proteotypic = my_match_lists.results[i].proteotypic;
//
//  }
//
//
//  /*for (int i = 0; i < 10; i++) {
//      cout << my_peptide_lists.spectra[i].seq << "   " << my_peptide_lists.spectra[i].prot_seq << "  " << my_peptide_lists.spectra[i].proteotypic << endl; 
//  }
//
//  cout << "======" << endl; */
//
//  my_peptide_lists.master.push_back(qc);
//  my_peptide_lists.master.back().push_back(my_peptide_lists.spectra[0]);
//  for (size_t i = 1; i < my_peptide_lists.spectra.size(); i++) {
//    if (my_peptide_lists.spectra[i].seq != my_peptide_lists.spectra[i - 1].seq) {
//      my_peptide_lists.master.push_back(qc);
//    }
//    my_peptide_lists.master.back().push_back(my_peptide_lists.spectra[i]);
//  }
//
//  //MH: Repeat for the other list
//  for (size_t i = 0; i < my_match_lists.results1.size(); i++) {
//    my_peptide_lists.spectra1.push_back(my_intensities()); //add first one
//    my_peptide_lists.spectra1.back().x = my_match_lists.results1[i].spec_rt;
//    my_peptide_lists.spectra1.back().y = my_match_lists.results1[i].spec_intensity;
//    my_peptide_lists.spectra1.back().seq = my_match_lists.results1[i].pep_seq;
//    my_peptide_lists.spectra1.back().mc = my_match_lists.results1[i].miss_cleaves;
//    my_peptide_lists.spectra1.back().prot_seq = my_match_lists.results1[i].prot_seq;
//    my_peptide_lists.spectra1.back().proteotypic = my_match_lists.results1[i].proteotypic;
//  }
//
// /* for (int i = 0; i < 10; i++) {
//      cout << my_peptide_lists.spectra1[i].seq << "   " << my_peptide_lists.spectra1[i].prot_seq << "  " << my_peptide_lists.spectra1[i].proteotypic << endl;
//  }*/
//
//
//  my_peptide_lists.master1.push_back(qc);
//  my_peptide_lists.master1.back().push_back(my_peptide_lists.spectra1[0]);
//  for (size_t i = 1; i < my_peptide_lists.spectra1.size(); i++) {
//    if (my_peptide_lists.spectra1[i].seq != my_peptide_lists.spectra1[i - 1].seq) {
//      my_peptide_lists.master1.push_back(qc);
//    }
//    my_peptide_lists.master1.back().push_back(my_peptide_lists.spectra1[i]);
//  }
//
//
//  //DELETE NOISE (INTENSITIES LESS THAN 10% OF THE MAX) (NOISE)
//  //MH: The faster, simpler approach
//  for(size_t i=0;i<my_peptide_lists.master.size();i++){
//    if(!cleanNoise(my_peptide_lists.master[i])){
//      //handle error
//    }
//  }
//  for (size_t i = 0; i < my_peptide_lists.master1.size(); i++) {
//    if (!cleanNoise(my_peptide_lists.master1[i])) {
//      //handle error
//    }
//  }
//
//  ///*for (size_t i = 0; i < my_peptide_lists.master.size(); i++) {
//  //  vector<float> rn;
//  //  rn = cleanNoise(my_peptide_lists.master[i]);
//
//  //  for (size_t n = 0; n < my_peptide_lists.master[i].size() - 1; n++) {
//  //    for (size_t p = n + 1; p < my_peptide_lists.master[i].size(); p++) {
//  //      vector<float>::iterator it;
//  //      it = find(rn.begin(), rn.end(), my_peptide_lists.master[i][p].y);
//  //      if (it == rn.end()) {
//  //        my_peptide_lists.master[i].erase(my_peptide_lists.master[i].begin() + p);
//  //        p = (n + 1) - 1;
//  //      }
//  //    }
//  //  }
//  //  for (size_t n = 0; n < my_peptide_lists.master[i].size(); n++) {
//  //    vector<float>::iterator it;
//  //    it = find(rn.begin(), rn.end(), my_peptide_lists.master[i][n].y);
//  //    if (it == rn.end()) {
//  //      my_peptide_lists.master[i].erase(my_peptide_lists.master[i].begin() + n);
//  //    }
//  //  }
//  //  rn.clear();
//  //}
//
//  //for (size_t i = 0; i < my_peptide_lists.master1.size(); i++) {
//  //  vector<float> rn1;
//  //  rn1 = cleanNoise(my_peptide_lists.master1[i]);
//
//  //  for (size_t n = 0; n < my_peptide_lists.master1[i].size() - 1; n++) {
//  //    for (size_t p = n + 1; p < my_peptide_lists.master1[i].size(); p++) {
//  //      vector<float>::iterator it;
//  //      it = find(rn1.begin(), rn1.end(), my_peptide_lists.master1[i][p].y);
//  //      if (it == rn1.end()) {
//  //        my_peptide_lists.master1[i].erase(my_peptide_lists.master1[i].begin() + p);
//  //        p = (n + 1) - 1;
//  //      }
//  //    }
//  //  }
//  //  for (size_t n = 0; n < my_peptide_lists.master1[i].size(); n++) {
//  //    vector<float>::iterator it;
//  //    it = find(rn1.begin(), rn1.end(), my_peptide_lists.master1[i][n].y);
//  //    if (it == rn1.end()) {
//  //      my_peptide_lists.master1[i].erase(my_peptide_lists.master1[i].begin() + n);
//  //    }
//  //  }
//  //  rn1.clear();
//  //}*/
//
//  // END NOISE FUNCTION 
//
//
//  //delete peptides with only one data point
//  vector<vector<my_intensities>> filter;
//  for (size_t i = 0; i < my_peptide_lists.master.size(); i++) {
//    if (my_peptide_lists.master[i].size() > 1) {
//      filter.push_back(my_peptide_lists.master[i]);
//    }
//  }
//  my_peptide_lists.master = filter;
//  filter.clear();
//  for (size_t i = 0; i < my_peptide_lists.master1.size(); i++) {
//    if (my_peptide_lists.master1[i].size() > 1) {
//      filter.push_back(my_peptide_lists.master1[i]);
//    }
//  }
//  my_peptide_lists.master1 = filter;
//  filter.clear();
//
//
//
//
//  
//
//
//
//  cout << "noise deleted" << "\n" << endl;
//
//
//  //metric calculations
//  vector<float> peak;
//  float avg=0;
//  //MH: this array isn't necessary
//  //vector<float> peak2; 
//  for (int i = 0; i < my_peptide_lists.master.size(); i++) {
//    for (int j = 0; j < my_peptide_lists.master[i].size(); j++) {
//      peak.push_back(my_peptide_lists.master[i][j].y);
//    }
//    avg+= *max_element(peak.begin(), peak.end());
//    //peak2.push_back(*max_element(peak.begin(), peak.end()));
//    peak.clear();
//  }
//
//
//  //float sum = 0; for (int i = 0; i < peak2.size(); ++i) { sum += peak2[i]; }
//  //float sumi = sum / peak2.size();
//  my_metrics.miss_avg_high = avg/ my_peptide_lists.master.size();
//
//  peak.clear();
//  avg=0;
//  //vector<float> peak1;
//  //vector<float> peak3;
//  for (int i = 0; i < my_peptide_lists.master1.size(); i++) {
//    for (int j = 0; j < my_peptide_lists.master1[i].size(); j++) {
//      peak.push_back(my_peptide_lists.master1[i][j].y);
//    }
//    avg += *max_element(peak.begin(), peak.end());
//    //peak3.push_back(*max_element(peak1.begin(), peak1.end()));
//    peak.clear();
//  }
//
//
//  //float sum1 = 0; for (int i = 0; i < peak3.size(); ++i) { sum1 += peak3[i]; }
//  //float sumi1 = sum1 / peak3.size();
//  my_metrics.tryp_avg_high = avg/ my_peptide_lists.master1.size();
//
//  //MH: Calculating area is something performed on multiple arrays. So create a function
//  // that does it on any array you pass to it. Also, have it return the total rather than
//  // iterate the array again to count the total.
//  for (size_t i = 0; i < my_peptide_lists.master.size(); i++) {
//    my_peptide_lists.master[i].back().tot = calcPeakArea(my_peptide_lists.master[i]);
//  }
//  for (size_t i = 0; i < my_peptide_lists.master1.size(); i++) {
//    my_peptide_lists.master1[i].back().tot = calcPeakArea(my_peptide_lists.master1[i]);
//  }
//
//
//  vector<my_intensities> final;
//  int cycle = 0;
//  for (int i = 0; i < my_peptide_lists.master.size(); i++) {
//    final.push_back(my_intensities());
//    final[cycle].seq = my_peptide_lists.master[i][my_peptide_lists.master[i].size() - 1].seq;
//    final[cycle].tot = my_peptide_lists.master[i][my_peptide_lists.master[i].size() - 1].tot;
//    final[cycle].mc = my_peptide_lists.master[i][my_peptide_lists.master[i].size() - 1].mc;
//    final[cycle].prot_seq = my_peptide_lists.master[i][my_peptide_lists.master[i].size() - 1].prot_seq;
//    final[cycle].proteotypic = my_peptide_lists.master[i][my_peptide_lists.master[i].size() - 1].proteotypic;
//    cycle++;
//  }
//
//  vector<my_intensities> final1;
//  int cycle1 = 0;
//  for (int i = 0; i < my_peptide_lists.master1.size(); i++) {
//    final1.push_back(my_intensities());
//    final1[cycle1].seq = my_peptide_lists.master1[i][my_peptide_lists.master1[i].size() - 1].seq;
//    final1[cycle1].tot = my_peptide_lists.master1[i][my_peptide_lists.master1[i].size() - 1].tot;
//    final1[cycle1].mc = my_peptide_lists.master1[i][my_peptide_lists.master1[i].size() - 1].mc;
//    final1[cycle1].prot_seq = my_peptide_lists.master1[i][my_peptide_lists.master1[i].size() - 1].prot_seq;
//    final1[cycle1].proteotypic = my_peptide_lists.master1[i][my_peptide_lists.master1[i].size() - 1].proteotypic;
//    cycle1++;
//  }
//
//  my_peptide_lists.final = final;
//  my_peptide_lists.final1 = final1;
//
//  // END AREA FUNCTION
//
//
//  cout << "intensity calculated" << "\n" << endl;
//
//
//
//
//  return true;
//
//}

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



	my_metrics.total_psm = my_peptide_lists.total;
	my_metrics.psm_num = my_peptide_lists.all_real.size();
	my_metrics.tryptic_num = my_peptide_lists.tryptic_real.size();
	my_metrics.nontryptic_num = my_peptide_lists.non_tryptic_real.size();
	my_metrics.unique_pep_charge = my_peptide_lists.tryp_unique_z_real.size();
	my_metrics.unique_peptides = my_peptide_lists.tryp_unique_real.size();
	my_metrics.avg_pep_length = f / my_metrics.unique_peptides;
	my_metrics.tryp_frac = (double)my_metrics.tryptic_num / my_metrics.psm_num;
	my_metrics.nontryp_frac = (double)my_metrics.nontryptic_num / my_metrics.psm_num;
	my_metrics.pep_frac = (double)my_peptide_lists.tryp_unique_real.size() / my_metrics.psm_num;
	my_metrics.miss_cleave_rate_psm = r / my_metrics.tryptic_num;
	my_metrics.miss_cleave_rate_pep = g / my_metrics.unique_peptides;
	my_metrics.num_miss_cleave_pep = h; 
	my_metrics.num_miss_cleave_total = g; 
	my_metrics.miss_cleave_avg = g / h;
	my_metrics.num_miss_cleave_psm = p; 
	my_metrics.num_tot_miss_cleave_psm = r; 
	my_metrics.golden_stat_psm = my_metrics.num_miss_cleave_psm / my_metrics.tryptic_num;
	my_metrics.golden_stat_unique_pep = my_metrics.num_miss_cleave_pep / my_metrics.unique_peptides;
	

	return my_metrics;

}



bool deep_functions::cleanNoise(std::vector<dsXIC>& v) {
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
    if(v[a-1].intensity>max) a--;
    else break;
  }

  //MH: Now step to the right of the apex to find your upper bound
  b=maxIndex;
  while(b<v.size()-1){
    if(v[b+1].intensity>max) b++;
    else break;
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
float deep_functions::calcPeakArea(std::vector<dsXIC>& v){
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



//float deep_functions::calcPeakArea(std::vector<my_intensities>& v) {
//    float w;
//    float hRect;
//    float hTri;
//    float total = 0;
//    for (size_t i = 0; i < v.size() - 1; i++) {
//        w = v[i + 1].x - v[i].x;
//        hTri = abs(v[i + 1].y - v[i].y);
//        if (v[i + 1].y < v[i].y) hRect = v[i + 1].y;
//        else hRect = v[i].y;
//        v[i].tot = w * hRect + w * hTri / 2;
//        total += v[i].tot;
//    }
//    return total;
//}


void deep_functions::print(peptide_lists& my_peptide_lists, metrics& my_metrics) {


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
    cout << "avg ratio of intensities: 0 miscleave full tryptic peptides / 1 or 2 miscleave peptides:\t" << my_metrics.intensity_final << " ** " << "\n" << endl;
    cout << "stdv of ratio of intensities: \t\t\t\t\t\t\t\t" << my_metrics.stdv_final << " ** " << "\n" << endl;
    cout << "total miscleaved peptide intensity / sum of all peptide intensities: \t\t" << my_metrics.total_intensity << " ** " << "\n" << endl; 




    sort(my_peptide_lists.peptide_matches.begin(), my_peptide_lists.peptide_matches.end(), compareInten);


    my_peptide_lists.peptide_matches1.push_back(my_peptide_lists.peptide_matches[0]);
    for (size_t i = 1; i < my_peptide_lists.peptide_matches.size(); i++) {
        if (my_peptide_lists.peptide_matches[i].ft_pep_seq == my_peptide_lists.peptide_matches[i - 1].ft_pep_seq) continue;
        my_peptide_lists.peptide_matches1.push_back(my_peptide_lists.peptide_matches[i]);
    }

    cout << "\n" << "-- 20 TWENTY MOST ABUNDANT PEPTIDES --" << "\n" << endl;
    my_peptide_lists.peptide_matches = my_peptide_lists.peptide_matches1;
    if (my_peptide_lists.peptide_matches.size() < 20) {
        for (size_t i = 0; i < my_peptide_lists.peptide_matches.size(); i++) {
            cout << i+1 << ": " << my_peptide_lists.peptide_matches[i].mc_pep_seq << "  " << my_peptide_lists.peptide_matches[i].mc_areaXIC << "  " << my_peptide_lists.peptide_matches[i].ft_pep_seq << "  " << my_peptide_lists.peptide_matches[i].ft_areaXIC << "  " << my_peptide_lists.peptide_matches[i].ft_areaXIC / my_peptide_lists.peptide_matches[i].mc_areaXIC << endl;
        }
    }
    else {
        for (size_t i = 0; i < 20; i++) {
            cout << i+1 << ": " << my_peptide_lists.peptide_matches[i].mc_pep_seq << "  " << my_peptide_lists.peptide_matches[i].mc_areaXIC << "  " << my_peptide_lists.peptide_matches[i].ft_pep_seq << "  " << my_peptide_lists.peptide_matches[i].ft_areaXIC << "  " << my_peptide_lists.peptide_matches[i].ft_areaXIC / my_peptide_lists.peptide_matches[i].mc_areaXIC << endl;
        }
    }
    

   

    vector<float> rat;
    float rat1 = 0;
    for (size_t i = 0; i < my_peptide_lists.prot_f.size(); i++) {
        rat.push_back(my_peptide_lists.prot_f[i].percentMiss);
        rat1 += my_peptide_lists.prot_f[i].percentMiss;
    }

    float rat2 = rat1 / my_peptide_lists.prot_f.size();
    my_metrics.protein_final = rat2;

    float var = 0;
    for (size_t i = 0; i < rat.size(); ++i) {
        var += pow(rat[i] - rat2, 2);
    }
    var = var / rat.size();
    float stdv = sqrt(var);
    my_metrics.protein_stdv = stdv;


    cout << "\n" << "-- PROTEIN STATS --" << "\n" << endl;
    cout << "avg percent miss within protein: misscleaved peptides sum intensity / total sum intensity:\t" << my_metrics.protein_final << " ** " << "\n" << endl;
    cout << "stdv of percent miss: \t\t\t\t\t\t\t\t" << my_metrics.protein_stdv << " ** " << "\n" << endl;


    sort(my_peptide_lists.prot_f.begin(), my_peptide_lists.prot_f.end(), compareTotal); 

 


    cout << "\n" << "-- 20 TWENTY MOST ABUNDANT PROTEINS--" << "\n" << endl;
    
    if (my_peptide_lists.prot_f.size() < 20) {
        for (size_t i = 0; i < my_peptide_lists.prot_f.size(); i++) {
            cout << i + 1 << ": " << my_peptide_lists.prot_f[i].prot_seq << "  " << my_peptide_lists.prot_f[i].total << "  " << my_peptide_lists.prot_f[i].percentMiss <<  endl;
        }
    }
    else {
        for (size_t i = 0; i < 20; i++) {
            cout << i + 1 << ": " << my_peptide_lists.prot_f[i].prot_seq << "  " << my_peptide_lists.prot_f[i].total << "  " << my_peptide_lists.prot_f[i].percentMiss << endl;
        }
    }

    int count =0; 
    for (size_t i = 0; i < my_peptide_lists.prot_f.size(); i++) {
        if (my_peptide_lists.prot_f[i].percentMiss == 0 /*&& my_peptide_lists.prot_f[i].trypPeptides.size() >= 5*/) {
            count++; 
        }
    }

    cout << "\n" << count << "  proteins out of " << my_peptide_lists.prot_f.size() << " are entirely made up of tryptic peptides" << endl; 

   /* for (int i = 0; i < my_peptide_lists.prot_f.size(); i++) {
        cout << my_peptide_lists.prot_f[i].prot_seq << "   " << my_peptide_lists.prot_f[i].trypPeptides.size() << "   " << my_peptide_lists.prot_f[i].missPeptides.size() << "   " << my_peptide_lists.prot_f[i].sumTryp << "   " << my_peptide_lists.prot_f[i].sumMiss << "   " << my_peptide_lists.prot_f[i].percentMiss << endl;
    }*/

  
}
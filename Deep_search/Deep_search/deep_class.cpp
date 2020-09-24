#include "rapidxml.hpp"
#include "deep_class.h"
#include "rapidjson/document.h"
#include "rapidjson/prettywriter.h"
#include "rapidjson/stringbuffer.h"


using namespace rapidxml;
using namespace std;
using namespace rapidjson;



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
                            my_peptide_lists.all_psm.push_back(dsPSM());
                            my_peptide_lists.all_psm[c].pep_seq = string(sample_node->first_node("search_result")->first_node("search_hit")->first_attribute("peptide")->value());
                            my_peptide_lists.all_psm[c].charge = atoi(sample_node->first_attribute("assumed_charge")->value());
                            my_peptide_lists.all_psm[c].pre_neutral_mass = atof(sample_node->first_attribute("precursor_neutral_mass")->value());
                            my_peptide_lists.all_psm[c].xml_rtime = (float)atof(sample_node->first_attribute("retention_time_sec")->value());
                            my_peptide_lists.all_psm[c].prev_aa = (string)sample_node->first_node("search_result")->first_node("search_hit")->first_attribute("peptide_prev_aa")->value();
                            my_peptide_lists.all_psm[c].next_aa = string(sample_node->first_node("search_result")->first_node("search_hit")->first_attribute("peptide_next_aa")->value());
                            my_peptide_lists.all_psm[c].prot_seq = string(sample_node->first_node("search_result")->first_node("search_hit")->first_attribute("protein")->value());
                            my_peptide_lists.all_psm[c].calc_neutral_mass = atof(sample_node->first_node("search_result")->first_node("search_hit")->first_attribute("calc_neutral_pep_mass")->value());
                            if (atoi(sample_node->first_node("search_result")->first_node("search_hit")->first_attribute("num_tot_proteins")->value()) == 1) {
                                my_peptide_lists.all_psm[c].proteotypic = 1;
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
                    my_peptide_lists.all_psm.push_back(dsPSM());
                    my_peptide_lists.all_psm[c].pep_seq = string(sample_node->first_node("search_result")->first_node("search_hit")->first_attribute("peptide")->value());
                    my_peptide_lists.all_psm[c].charge = atoi(sample_node->first_attribute("assumed_charge")->value());
                    my_peptide_lists.all_psm[c].pre_neutral_mass = atof(sample_node->first_attribute("precursor_neutral_mass")->value());
                    my_peptide_lists.all_psm[c].xml_rtime = (float)atof(sample_node->first_attribute("retention_time_sec")->value());
                    my_peptide_lists.all_psm[c].prev_aa = string(sample_node->first_node("search_result")->first_node("search_hit")->first_attribute("peptide_prev_aa")->value());
                    my_peptide_lists.all_psm[c].next_aa = string(sample_node->first_node("search_result")->first_node("search_hit")->first_attribute("peptide_next_aa")->value());
                    my_peptide_lists.all_psm[c].prot_seq = string(sample_node->first_node("search_result")->first_node("search_hit")->first_attribute("protein")->value());
                    my_peptide_lists.all_psm[c].calc_neutral_mass = atof(sample_node->first_node("search_result")->first_node("search_hit")->first_attribute("calc_neutral_pep_mass")->value());
                    if (atoi(sample_node->first_node("search_result")->first_node("search_hit")->first_attribute("num_tot_proteins")->value()) == 1) {
                        my_peptide_lists.all_psm[c].proteotypic = 1;
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


    /*for (int i = 0; i < my_peptide_lists.all_real.size(); i++) {
        cout << my_peptide_lists.all_real[i].pep_seq << "  " << my_peptide_lists.all_real[i].xml_rtime << endl; 
    }*/



    return my_peptide_lists;

}

bool deep_functions::enzymatic_calc(peptide_lists& my_peptide_lists, my_parameters& my_params) {


    for (size_t i = 0; i < my_peptide_lists.all_psm.size(); i++) {
        size_t found = my_params.cleave_loc.find(my_peptide_lists.all_psm[i].prev_aa);
        size_t found1 = my_params.cleave_loc.find(my_peptide_lists.all_psm[i].pep_seq.back());
        size_t found2 = my_params.hyphen.find(my_peptide_lists.all_psm[i].next_aa);
        if (found!= string::npos && (found1!=string::npos || found2!=string::npos)) {
            my_peptide_lists.all_psm[i].enzymatic = 1; 
        }
        if (found == string::npos && found1 == string::npos) {
            my_peptide_lists.all_psm[i].non_enzymatic = 1;
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


bool deep_functions::miss_cleave(peptide_lists& my_peptide_lists, my_parameters& my_params) {

   
    int c = 0;
    for (int i = 0; i < my_peptide_lists.all_psm.size(); i++) {
        for (int j = 0; j < my_peptide_lists.all_psm[i].pep_seq.size() - 1; j++) {
            size_t found = my_params.cleave_loc.find(my_peptide_lists.all_psm[i].pep_seq[j]);
            size_t found1 = my_params.anti_cleave_loc.find(my_peptide_lists.all_psm [i].pep_seq[j + 1]);
            if (found!=string::npos && found1==string::npos) {
                c++;
            }
        }
        my_peptide_lists.all_psm[i].miss_cleaves = c;
        c = 0;
    }

    return true;

}

bool deep_functions::delete_dup(peptide_lists& my_peptide_lists) {

    size_t i;
    vector<dsPSM> tmp;
    vector<dsPeptide> tmp1;
    for (i = 0; i < my_peptide_lists.all_psm.size(); i++) {
        if (my_peptide_lists.all_psm[i].enzymatic == 1) {
            tmp.push_back(my_peptide_lists.all_psm[i]);
        }
    }

    sort(tmp.begin(), tmp.end(), compareSeqZ);
    
    for (i = 0; i < tmp.size(); i++) {
        tmp1.push_back(dsPeptide()); 
        tmp1[i].pep_seq = tmp[i].pep_seq;
        tmp1[i].prot_seq = tmp[i].prot_seq;
        tmp1[i].proteotypic = tmp[i].proteotypic;
        tmp1[i].enzymatic = tmp[i].enzymatic;
        tmp1[i].non_enzymatic = tmp[i].non_enzymatic;
        tmp1[i].semi_enzymatic = tmp[i].semi_enzymatic;
        tmp1[i].charge = tmp[i].charge;
        tmp1[i].xml_rtime = tmp[i].xml_rtime;
        tmp1[i].pre_neutral_mass = tmp[i].pre_neutral_mass;
        tmp1[i].calc_neutral_mass = tmp[i].calc_neutral_mass;
        tmp1[i].xml_mz = tmp[i].xml_mz;
        tmp1[i].miss_cleaves = tmp[i].miss_cleaves;
        tmp1[i].prev_aa = tmp[i].prev_aa;
        tmp1[i].next_aa = tmp[i].next_aa;
        
    }
   
    my_peptide_lists.enzymatic_unique.push_back(tmp1[0]);
    for (i = 1; i < tmp1.size(); i++) {
        if (my_peptide_lists.enzymatic_unique.back().pep_seq.compare(tmp1[i].pep_seq) == 0) continue;
        my_peptide_lists.enzymatic_unique.push_back(tmp1[i]);
    }
    
    vector<string> temp; 
    for (int i = 0; i < tmp1.size(); i++) {
        temp.push_back(tmp1[i].pep_seq); 
    }
    
    for (int i = 0; i < my_peptide_lists.enzymatic_unique.size(); i++) {
        int c = count(temp.begin(), temp.end(), my_peptide_lists.enzymatic_unique[i].pep_seq);
        my_peptide_lists.enzymatic_unique[i].psm_count = c; 
        c = 0; 
    }

   
    return true;

}


//bool deep_functions::delete_dup(peptide_lists& my_peptide_lists) {
//
//    size_t i;
//    vector<dsPSM> tmp, tmp1;
//    for (i = 0; i < my_peptide_lists.all_psm.size(); i++) {
//        if (my_peptide_lists.all_psm[i].enzymatic == 1) {
//            tmp.push_back(my_peptide_lists.all_psm[i]);
//        }
//    }
//
//    //Sort the tryptic array by sequence and then charge state
//    sort(tmp.begin(), tmp.end(), compareSeqZ);
//
//
//    //Iterate over all tryptic PSMs, keeping each unique instance of sequence and charge
//    my_peptide_lists.tryp_unique_z_real.push_back(tmp[0]);
//    for (i = 1; i < my_peptide_lists.tryptic_real.size(); i++) {
//        if (my_peptide_lists.tryp_unique_z_real.back().pep_seq.compare(my_peptide_lists.tryptic_real[i].pep_seq) == 0 &&
//            my_peptide_lists.tryp_unique_z_real.back().charge == my_peptide_lists.tryptic_real[i].charge) continue;
//        my_peptide_lists.tryp_unique_z_real.push_back(my_peptide_lists.tryptic_real[i]);
//    }
//
//
//    //Repeat the process, keeping each unique instance of sequence;
//    my_peptide_lists.tryp_unique_real.push_back(my_peptide_lists.tryptic_real[0]);
//    for (i = 1; i < my_peptide_lists.tryptic_real.size(); i++) {
//        if (my_peptide_lists.tryp_unique_real.back().pep_seq.compare(my_peptide_lists.tryptic_real[i].pep_seq) == 0) continue;
//        my_peptide_lists.tryp_unique_real.push_back(my_peptide_lists.tryptic_real[i]);
//    }
//
//
//
//    //Now of the unique peptide sequences, make a subset of miscleaved peptides
//    i = 0;
//    while (my_peptide_lists.tryp_unique_real[i].miss_cleaves == 0) i++; //no boundary checks here
//    my_peptide_lists.miss_unique_real.push_back(my_peptide_lists.tryp_unique_real[i]);
//    for (i = i + 1; i < my_peptide_lists.tryp_unique_real.size(); i++) {
//        if (my_peptide_lists.tryp_unique_real[i].miss_cleaves == 0) continue;
//        my_peptide_lists.miss_unique_real.push_back(my_peptide_lists.tryp_unique_real[i]);
//    }
//
//    for (i = 0; i < my_peptide_lists.tryp_unique_real.size(); i++) {
//        if (my_peptide_lists.tryp_unique_real[i].miss_cleaves != 0) continue;
//        my_peptide_lists.fully_tryp_unique_real.push_back(my_peptide_lists.tryp_unique_real[i]);
//    }
//
//
//    if (my_peptide_lists.semi_tryptic_real.size() > 0) {
//        sort(my_peptide_lists.semi_tryptic_real.begin(), my_peptide_lists.semi_tryptic_real.end(), compareSeqZ);
//
//        //semi tryptic 
//        my_peptide_lists.semi_tryptic_unique_z_real.push_back(my_peptide_lists.semi_tryptic_real[0]);
//        for (i = 1; i < my_peptide_lists.semi_tryptic_real.size(); i++) {
//            if (my_peptide_lists.semi_tryptic_unique_z_real.back().pep_seq.compare(my_peptide_lists.semi_tryptic_real[i].pep_seq) == 0 &&
//                my_peptide_lists.semi_tryptic_unique_z_real.back().charge == my_peptide_lists.semi_tryptic_real[i].charge) continue;
//            my_peptide_lists.semi_tryptic_unique_z_real.push_back(my_peptide_lists.semi_tryptic_real[i]);
//        }
//
//        // semi tryptic 
//        my_peptide_lists.semi_tryptic_unique_real.push_back(my_peptide_lists.semi_tryptic_real[0]);
//        for (i = 1; i < my_peptide_lists.semi_tryptic_real.size(); i++) {
//            if (my_peptide_lists.semi_tryptic_unique_real.back().pep_seq.compare(my_peptide_lists.semi_tryptic_real[i].pep_seq) == 0) continue;
//            my_peptide_lists.semi_tryptic_unique_real.push_back(my_peptide_lists.semi_tryptic_real[i]);
//        }
//
//    }
//    return true;
//
//}

bool deep_functions::lcd(peptide_lists& my_peptide_lists, my_parameters& my_params) {


    for (size_t i = 0; i < my_peptide_lists.enzymatic_unique.size(); i++) {
        my_peptide_lists.enzymatic_unique[i].xml_mz = (my_peptide_lists.enzymatic_unique[i].pre_neutral_mass + (my_peptide_lists.enzymatic_unique[i].charge * 1.00727)) / my_peptide_lists.enzymatic_unique[i].charge;
        my_peptide_lists.enzymatic_unique[i].tolerance = (my_params.ppm / (1000000)) * my_peptide_lists.enzymatic_unique[i].calc_neutral_mass;
    }

    return true; 

}

//bool deep_functions::lcd(peptide_lists& my_peptide_lists, my_parameters& my_params) {
//
//
//    for (size_t i = 0; i < my_peptide_lists.fully_tryp_unique_real.size(); i++) {
//        my_peptide_lists.fully_tryp_unique_real[i].xml_mz = (my_peptide_lists.fully_tryp_unique_real[i].pre_neutral_mass + (my_peptide_lists.fully_tryp_unique_real[i].charge * 1.00727)) / my_peptide_lists.fully_tryp_unique_real[i].charge;
//        my_peptide_lists.fully_tryp_unique_real[i].tolerance = (my_params.ppm / (1000000)) * my_peptide_lists.fully_tryp_unique_real[i].calc_neutral_mass;
//    }
//
//    for (size_t i = 0; i < my_peptide_lists.miss_unique_real.size(); i++) {
//        my_peptide_lists.miss_unique_real[i].xml_mz = (my_peptide_lists.miss_unique_real[i].pre_neutral_mass + (my_peptide_lists.miss_unique_real[i].charge * 1.00727)) / my_peptide_lists.miss_unique_real[i].charge;
//        my_peptide_lists.miss_unique_real[i].tolerance = (my_params.ppm / (1000000)) * my_peptide_lists.miss_unique_real[i].calc_neutral_mass;
//    }
//
//    for (size_t i = 0; i < my_peptide_lists.semi_tryptic_unique_real.size(); i++) {
//        my_peptide_lists.semi_tryptic_unique_real[i].xml_mz = (my_peptide_lists.semi_tryptic_unique_real[i].pre_neutral_mass + (my_peptide_lists.semi_tryptic_unique_real[i].charge * 1.00727)) / my_peptide_lists.semi_tryptic_unique_real[i].charge;
//        my_peptide_lists.semi_tryptic_unique_real[i].tolerance = (my_params.ppm / (1000000)) * my_peptide_lists.semi_tryptic_unique_real[i].calc_neutral_mass;
//    }
//
//
//    return true;
//
//}

bool deep_functions::reader(peptide_lists& my_peptide_lists, metrics& my_metrics) {

    dsXIC x;
    int a = 0; 
    /*int b = 0; */
    for (size_t i = 0; i < my_peptide_lists.enzymatic_unique.size(); i++) {
        if (my_peptide_lists.enzymatic_unique[i].pep_seq == my_peptide_lists.enzymatic_unique[i+1].pep_seq) {
            x.rTime = my_peptide_lists.enzymatic_unique[i].spec_rt;
            x.intensity = my_peptide_lists.enzymatic_unique[i].spec_intensity;
            my_peptide_lists.enzymatic_unique[a].XIC.push_back(x);
        }
        else if (my_peptide_lists.enzymatic_unique[i].pep_seq != my_peptide_lists.enzymatic_unique[i + 1].pep_seq) {
            x.rTime = my_peptide_lists.enzymatic_unique[i].spec_rt;
            x.intensity = my_peptide_lists.enzymatic_unique[i].spec_intensity;
            my_peptide_lists.enzymatic_unique[a].XIC.push_back(x);
            a = i + 1; 
        }
      
    }
    vector<dsPeptide> tmp;
    tmp.push_back(my_peptide_lists.enzymatic_unique[0]);
    for (size_t i = 1; i < my_peptide_lists.enzymatic_unique.size(); i++) {
        if (my_peptide_lists.enzymatic_unique[i].pep_seq == my_peptide_lists.enzymatic_unique[i - 1].pep_seq) continue;
        tmp.push_back(my_peptide_lists.enzymatic_unique[i]);
    }
    my_peptide_lists.enzymatic_unique = tmp; 
    tmp.clear();



    /*for (size_t i = 0; i < my_peptide_lists.xic_mc_results.size(); i++) {
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
    my_peptide_lists.xic_mc_results = tmp;*/





    //DELETE NOISE (INTENSITIES LESS THAN 10% OF THE MAX) (NOISE)

    for (size_t i = 0; i < my_peptide_lists.enzymatic_unique.size(); i++) {
        if (!cleanNoise(my_peptide_lists.enzymatic_unique[i].XIC)) {
            //handle error
        }
    }
    //for (size_t i = 0; i < my_peptide_lists.xic_mc_results.size(); i++) {
    //    if (!cleanNoise(my_peptide_lists.xic_mc_results[i].XIC)) {
    //        //handle error
    //    }
    //}

   


    //delete peptides with only one data point
    vector<dsPeptide> filter;
    for (size_t i = 0; i < my_peptide_lists.enzymatic_unique.size(); i++) {
        if (my_peptide_lists.enzymatic_unique[i].XIC.size() > 1) {
            filter.push_back(my_peptide_lists.enzymatic_unique[i]);
        }
    }
    my_peptide_lists.enzymatic_unique = filter;
    /*filter.clear();
    for (size_t i = 0; i < my_peptide_lists.xic_mc_results.size(); i++) {
        if (my_peptide_lists.xic_mc_results[i].XIC.size() > 1) {
            filter.push_back(my_peptide_lists.xic_mc_results[i]);
        }
    }
    my_peptide_lists.xic_mc_results = filter;
    filter.clear();*/



    cout << "noise deleted" << "\n" << endl;


    ////metric calculations
    //vector<float> peak;
    //float avg = 0;
    //for (int i = 0; i < my_peptide_lists.xic_ft_results.size(); i++) {
    //    for (int j = 0; j < my_peptide_lists.xic_ft_results[i].XIC.size(); j++) {
    //        peak.push_back(my_peptide_lists.xic_ft_results[i].XIC[j].intensity);
    //    }
    //    avg += *max_element(peak.begin(), peak.end());
    //    peak.clear();
    //}


    //my_metrics.tryp_avg_high = avg / my_peptide_lists.xic_ft_results.size();

    //peak.clear();
    //avg = 0;
    //for (int i = 0; i < my_peptide_lists.xic_mc_results.size(); i++) {
    //    for (int j = 0; j < my_peptide_lists.xic_mc_results[i].XIC.size(); j++) {
    //        peak.push_back(my_peptide_lists.xic_mc_results[i].XIC[j].intensity);
    //    }
    //    avg += *max_element(peak.begin(), peak.end());
    //    peak.clear();
    //}


    //my_metrics.miss_avg_high = avg / my_peptide_lists.xic_mc_results.size();

   
   /* float big = 0; 
    float big1 = 0; */
    //CALC AREA 
    for (size_t i = 0; i < my_peptide_lists.enzymatic_unique.size(); i++) {
        my_peptide_lists.enzymatic_unique[i].areaXIC = calcPeakArea(my_peptide_lists.enzymatic_unique[i].XIC);
       /* big += my_peptide_lists.xic_ft_results[i].areaXIC; */
    }
    //for (size_t i = 0; i < my_peptide_lists.xic_mc_results.size(); i++) {
    //    my_peptide_lists.xic_mc_results[i].areaXIC = calcPeakArea(my_peptide_lists.xic_mc_results[i].XIC);
    //  /*  big1 += my_peptide_lists.xic_mc_results[i].areaXIC; */
    //}

  /*  my_metrics.total_intensity = big1 / (big + big1); */
   
    // END AREA FUNCTION


    cout << "intensity calculated" << "\n" << endl;




    return true;

}
//bool deep_functions::reader(peptide_lists& my_peptide_lists, metrics& my_metrics) {
//
//    dsXIC x;
//    int a = 0;
//    int b = 0;
//    for (size_t i = 0; i < my_peptide_lists.xic_ft_results.size(); i++) {
//        if (my_peptide_lists.xic_ft_results[i].pep_seq == my_peptide_lists.xic_ft_results[i + 1].pep_seq) {
//            x.rTime = my_peptide_lists.xic_ft_results[i].spec_rt;
//            x.intensity = my_peptide_lists.xic_ft_results[i].spec_intensity;
//            my_peptide_lists.xic_ft_results[a].XIC.push_back(x);
//        }
//        else if (my_peptide_lists.xic_ft_results[i].pep_seq != my_peptide_lists.xic_ft_results[i + 1].pep_seq) {
//            x.rTime = my_peptide_lists.xic_ft_results[i].spec_rt;
//            x.intensity = my_peptide_lists.xic_ft_results[i].spec_intensity;
//            my_peptide_lists.xic_ft_results[a].XIC.push_back(x);
//            a = i + 1;
//        }
//
//    }
//    vector<dsPeptide> tmp;
//    tmp.push_back(my_peptide_lists.xic_ft_results[0]);
//    for (size_t i = 1; i < my_peptide_lists.xic_ft_results.size(); i++) {
//        if (my_peptide_lists.xic_ft_results[i].pep_seq == my_peptide_lists.xic_ft_results[i - 1].pep_seq) continue;
//        tmp.push_back(my_peptide_lists.xic_ft_results[i]);
//    }
//    my_peptide_lists.xic_ft_results = tmp;
//    tmp.clear();
//
//
//
//    for (size_t i = 0; i < my_peptide_lists.xic_mc_results.size(); i++) {
//        if (my_peptide_lists.xic_mc_results[i].pep_seq == my_peptide_lists.xic_mc_results[i + 1].pep_seq) {
//            x.rTime = my_peptide_lists.xic_mc_results[i].spec_rt;
//            x.intensity = my_peptide_lists.xic_mc_results[i].spec_intensity;
//            my_peptide_lists.xic_mc_results[b].XIC.push_back(x);
//        }
//        else if (my_peptide_lists.xic_mc_results[i].pep_seq != my_peptide_lists.xic_mc_results[i + 1].pep_seq) {
//            x.rTime = my_peptide_lists.xic_mc_results[i].spec_rt;
//            x.intensity = my_peptide_lists.xic_mc_results[i].spec_intensity;
//            my_peptide_lists.xic_mc_results[b].XIC.push_back(x);
//            b = i + 1;
//        }
//
//    }
//    tmp.push_back(my_peptide_lists.xic_mc_results[0]);
//    for (size_t i = 1; i < my_peptide_lists.xic_mc_results.size(); i++) {
//        if (my_peptide_lists.xic_mc_results[i].pep_seq == my_peptide_lists.xic_mc_results[i - 1].pep_seq) continue;
//        tmp.push_back(my_peptide_lists.xic_mc_results[i]);
//    }
//    my_peptide_lists.xic_mc_results = tmp;
//
//
//
//
//
//    //DELETE NOISE (INTENSITIES LESS THAN 10% OF THE MAX) (NOISE)
//
//    for (size_t i = 0; i < my_peptide_lists.xic_ft_results.size(); i++) {
//        if (!cleanNoise(my_peptide_lists.xic_ft_results[i].XIC)) {
//            //handle error
//        }
//    }
//    for (size_t i = 0; i < my_peptide_lists.xic_mc_results.size(); i++) {
//        if (!cleanNoise(my_peptide_lists.xic_mc_results[i].XIC)) {
//            //handle error
//        }
//    }
//
//
//
//
//    //delete peptides with only one data point
//    vector<dsPeptide> filter;
//    for (size_t i = 0; i < my_peptide_lists.xic_ft_results.size(); i++) {
//        if (my_peptide_lists.xic_ft_results[i].XIC.size() > 1) {
//            filter.push_back(my_peptide_lists.xic_ft_results[i]);
//        }
//    }
//    my_peptide_lists.xic_ft_results = filter;
//    filter.clear();
//    for (size_t i = 0; i < my_peptide_lists.xic_mc_results.size(); i++) {
//        if (my_peptide_lists.xic_mc_results[i].XIC.size() > 1) {
//            filter.push_back(my_peptide_lists.xic_mc_results[i]);
//        }
//    }
//    my_peptide_lists.xic_mc_results = filter;
//    filter.clear();
//
//
//
//    cout << "noise deleted" << "\n" << endl;
//
//
//    //metric calculations
//    vector<float> peak;
//    float avg = 0;
//    for (int i = 0; i < my_peptide_lists.xic_ft_results.size(); i++) {
//        for (int j = 0; j < my_peptide_lists.xic_ft_results[i].XIC.size(); j++) {
//            peak.push_back(my_peptide_lists.xic_ft_results[i].XIC[j].intensity);
//        }
//        avg += *max_element(peak.begin(), peak.end());
//        peak.clear();
//    }
//
//
//    my_metrics.tryp_avg_high = avg / my_peptide_lists.xic_ft_results.size();
//
//    peak.clear();
//    avg = 0;
//    for (int i = 0; i < my_peptide_lists.xic_mc_results.size(); i++) {
//        for (int j = 0; j < my_peptide_lists.xic_mc_results[i].XIC.size(); j++) {
//            peak.push_back(my_peptide_lists.xic_mc_results[i].XIC[j].intensity);
//        }
//        avg += *max_element(peak.begin(), peak.end());
//        peak.clear();
//    }
//
//
//    my_metrics.miss_avg_high = avg / my_peptide_lists.xic_mc_results.size();
//
//
//    float big = 0;
//    float big1 = 0;
//    //CALC AREA 
//    for (size_t i = 0; i < my_peptide_lists.xic_ft_results.size(); i++) {
//        my_peptide_lists.xic_ft_results[i].areaXIC = calcPeakArea(my_peptide_lists.xic_ft_results[i].XIC);
//        big += my_peptide_lists.xic_ft_results[i].areaXIC;
//    }
//    for (size_t i = 0; i < my_peptide_lists.xic_mc_results.size(); i++) {
//        my_peptide_lists.xic_mc_results[i].areaXIC = calcPeakArea(my_peptide_lists.xic_mc_results[i].XIC);
//        big1 += my_peptide_lists.xic_mc_results[i].areaXIC;
//    }
//
//    my_metrics.total_intensity = big1 / (big + big1);
//
//    // END AREA FUNCTION
//
//
//    cout << "intensity calculated" << "\n" << endl;
//
//
//
//
//    return true;
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
	my_metrics.psm_num = my_peptide_lists.all_psm.size();
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


    cout << "\n" << "----------- XML STATS -----------" << "\n" << endl;
    cout << "Total PSM in file:\t\t\t\t\t\t\t\t\t\t" << my_metrics.total_psm << "\n" << endl;
    cout << "# of PSM above entered threshold value:\t\t\t\t\t\t\t\t" << my_metrics.psm_num << "\n" << endl;
    
    cout << "\n" << "----------- PEPTIDE STATS -----------" << "\n" << endl;
    cout << "-- Sanity Checks --" << "\n" << endl;
    cout << "# of tryptic PSM above entered threshold value:\t\t\t\t\t\t\t" << my_metrics.tryptic_num << "\n" << endl;
    cout << "# of non tryptic PSM above entered threshold value:\t\t\t\t\t\t" << my_metrics.nontryptic_num << "\n" << endl;
    cout << "Tryptic peptides with unique charges:\t\t\t\t\t\t\t\t" << my_metrics.unique_pep_charge << "\n" << endl;
    cout << "Tryptic peptides with unique sequences:\t\t\t\t\t\t\t\t" << my_metrics.unique_peptides << "\n" << endl;
    cout << "Avg seuqnce length of unique tryptic peptide:\t\t\t\t\t\t\t" << my_metrics.avg_pep_length << "\n" << endl;
    cout << "Tryptic psm / psm above threshold:\t\t\t\t\t\t\t\t" << my_metrics.tryp_frac << "\n" << endl;
    cout << "Non tryptic psm / psm above threshold:\t\t\t\t\t\t\t\t" << my_metrics.nontryp_frac << "\n" << endl;
    cout << "Unique peptides / psm above threshold:\t\t\t\t\t\t\t\t" << my_metrics.pep_frac << "\n" << endl;
    cout << "Avg # of misscleaves per unique peptide that is misclevaed:\t\t\t\t\t" << my_metrics.miss_cleave_avg << "\n" << endl;
    cout << "# of fully tryptic peptides found that has a match to a miscleaved peptide:\t\t\t" << my_metrics.find_num << "\n" << endl;
    cout << "avg highest peak intensity for miscleaved peptide: \t\t\t\t\t\t" << my_metrics.miss_avg_high << "\n" << endl;
    cout << "avg highest peak intensity for trypitc peptide: \t\t\t\t\t\t" << my_metrics.tryp_avg_high << "\n" << endl;
    cout << "# of miscleaved peptides that have fully tryptic matches that have 2 miscleaves:\t\t " << my_metrics.twice_mc << "\n" << endl;
    cout << "# of miscleaved peptides that have fully tryptic matches that have 1 miscleave:\t\t\t " << my_metrics.once_mc << "\n" << endl;

    cout << "-- Useful Stats --" << "\n" << endl;
    cout << "# of unique pepides that are miss cleaved:\t\t\t\t\t\t\t" << my_metrics.num_miss_cleave_pep << "\n" << endl;
    cout << "# of total misscleaves amognst unique peptides:\t\t\t\t\t\t\t" << my_metrics.num_miss_cleave_total << "\n" << endl;
    cout << "# of total misscleaves amongst unique peptides / total # of unique peptides:\t\t\t" << my_metrics.miss_cleave_rate_pep << "\n" << endl;
    cout << "# of misscleaved tryptic psm above threshold: \t\t\t\t\t\t\t" << my_metrics.num_miss_cleave_psm << "\n" << endl;
    cout << "# of total miss cleaves amongst tryptic psm:\t\t\t\t\t\t\t" << my_metrics.num_tot_miss_cleave_psm << "\n" << endl;
    cout << "# of total misscleaves amongst tryptic psm / total # tryptic psm above threshold:\t\t" << my_metrics.miss_cleave_rate_psm << "\n" << endl;
    cout << "# of misscleaved unique peptides / total # of unique peptides:\t\t\t\t\t" << my_metrics.golden_stat_unique_pep << " ** " << "\n" << endl;
    cout << "# of miss cleaved tryptic psm / total tryptic psm above threshold:\t\t\t\t" << my_metrics.golden_stat_psm << " ** " << "\n" << endl;
    cout << "avg ratio of intensities: 0 miscleave full tryptic peptides / 1 or 2 miscleave peptides:\t" << my_metrics.intensity_final << " ** " << "\n" << endl;
    cout << "stdv of ratio of intensities: \t\t\t\t\t\t\t\t\t" << my_metrics.stdv_final << " ** " << "\n" << endl;
    cout << "Global Stat: total miscleaved peptide intensity / sum of all peptide intensities: \t\t" << my_metrics.total_intensity << " ** " << "\n" << endl;
   
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


    cout << "\n" << "----------- PROTEIN STATS -----------" << "\n" << endl;
    cout << "avg percent miss within protein: misscleaved peptides sum intensity / total sum intensity:\t" << my_metrics.protein_final << " ** " << "\n" << endl;
    cout << "stdv of percent miss: \t\t\t\t\t\t\t\t\t\t" << my_metrics.protein_stdv << " ** " << "\n" << endl;

   
    




    sort(my_peptide_lists.peptide_matches.begin(), my_peptide_lists.peptide_matches.end(), compareInten);


    my_peptide_lists.peptide_matches1.push_back(my_peptide_lists.peptide_matches[0]);
    for (size_t i = 1; i < my_peptide_lists.peptide_matches.size(); i++) {
        if (my_peptide_lists.peptide_matches[i].ft_pep_seq == my_peptide_lists.peptide_matches[i - 1].ft_pep_seq) continue;
        my_peptide_lists.peptide_matches1.push_back(my_peptide_lists.peptide_matches[i]);
    }

    cout << "\n" << "-- 30 TWENTY MOST ABUNDANT PEPTIDES --" << "\n" << endl;
    my_peptide_lists.peptide_matches = my_peptide_lists.peptide_matches1;
    if (my_peptide_lists.peptide_matches.size() < 30) {
        for (size_t i = 0; i < my_peptide_lists.peptide_matches.size(); i++) {
            cout << i + 1 << ": " << " T sequence: " << my_peptide_lists.peptide_matches[i].ft_pep_seq << " || T intensity: " << my_peptide_lists.peptide_matches[i].ft_areaXIC << " || MC sequence: " << my_peptide_lists.peptide_matches[i].mc_pep_seq << " || MC intensity: " << my_peptide_lists.peptide_matches[i].mc_areaXIC << " || intensity ratio: " << my_peptide_lists.peptide_matches[i].ft_areaXIC / my_peptide_lists.peptide_matches[i].mc_areaXIC << endl;
        }
    }
    else {
        for (size_t i = 0; i < 30; i++) {
            cout << i + 1 << ": " << " T sequence: " << my_peptide_lists.peptide_matches[i].ft_pep_seq << " || T intensity: " << my_peptide_lists.peptide_matches[i].ft_areaXIC << " || MC sequence: " << my_peptide_lists.peptide_matches[i].mc_pep_seq << " || MC intensity: " << my_peptide_lists.peptide_matches[i].mc_areaXIC << " || intensity ratio: " << my_peptide_lists.peptide_matches[i].ft_areaXIC / my_peptide_lists.peptide_matches[i].mc_areaXIC << endl;
        }
    }
    

  

   

   


    sort(my_peptide_lists.prot_f.begin(), my_peptide_lists.prot_f.end(), compareTotal); 

 


    cout << "\n" << "-- 50 TWENTY MOST ABUNDANT PROTEINS--" << "\n" << endl;
    
    if (my_peptide_lists.prot_f.size() < 50) {
        for (size_t i = 0; i < my_peptide_lists.prot_f.size(); i++) {
            cout << i + 1 << ":  " << "Prot sequence: " <<  my_peptide_lists.prot_f[i].prot_seq << " || Tot intensity: " << my_peptide_lists.prot_f[i].total << " || amount MC peptides: " << my_peptide_lists.prot_f[i].percentMiss <<  endl;
        }
    }
    else {
        for (size_t i = 0; i < 50; i++) {
            cout << i + 1 << ":  " << "Prot sequence: " << my_peptide_lists.prot_f[i].prot_seq << " || Tot intensity: " << my_peptide_lists.prot_f[i].total << " || amount MC peptides: " << my_peptide_lists.prot_f[i].percentMiss << endl;
        }
    }

    int count =0; 
    for (size_t i = 0; i < my_peptide_lists.prot_f.size(); i++) {
        if (my_peptide_lists.prot_f[i].percentMiss == 0 /*&& my_peptide_lists.prot_f[i].trypPeptides.size() >= 5*/) {
            count++; 
        }
    }

    cout << "\n" << count << "  proteins out of " << my_peptide_lists.prot_f.size() << " are entirely made up of tryptic peptides" << "\n" << endl; 

   /* for (int i = 0; i < my_peptide_lists.prot_f.size(); i++) {
        cout << my_peptide_lists.prot_f[i].prot_seq << "   " << my_peptide_lists.prot_f[i].trypPeptides.size() << "   " << my_peptide_lists.prot_f[i].missPeptides.size() << "   " << my_peptide_lists.prot_f[i].sumTryp << "   " << my_peptide_lists.prot_f[i].sumMiss << "   " << my_peptide_lists.prot_f[i].percentMiss << endl;
    }*/







    cout << "------------------------------------------------------------------" << endl; 

    cout << "\n" << "*** MOST IMPORTANT DIAGNOSTIC METRICS ***" << "\n" << endl; 

    cout << "-- PEPTIDE STATS --" << "\n" << endl; 
    cout << "# of misscleaved unique peptides / total # of unique peptides:\t\t\t\t\t" << my_metrics.golden_stat_unique_pep << " ** " << "\n" << endl;
    cout << "# of miss cleaved tryptic psm / total tryptic psm above threshold:\t\t\t\t" << my_metrics.golden_stat_psm << " ** " << "\n" << endl;
    cout << "avg ratio of intensities: 0 miscleave full tryptic peptides / 1 or 2 miscleave peptides:\t" << my_metrics.intensity_final << " ** " << "\n" << endl;
    cout << "stdv of ratio of intensities: \t\t\t\t\t\t\t\t\t" << my_metrics.stdv_final << " ** " << "\n" << endl;
    cout << "Global Stat: total miscleaved peptide intensity / sum of all peptide intensities: \t\t" << my_metrics.total_intensity << " ** " << "\n" << endl;

    cout << "\n" << "-- PROTEIN STATS --" << "\n" << endl;
    cout << "avg percent miss within protein: misscleaved peptides sum intensity / total sum intensity:\t" << my_metrics.protein_final << " ** " << "\n" << endl;
    cout << "stdv of percent miss: \t\t\t\t\t\t\t\t\t\t" << my_metrics.protein_stdv << " ** " << "\n" << endl;





   /* for (size_t i = 0; i < my_peptide_lists.prot_f.size(); i++) {
        cout << i + 1 << "-  " << "protein sequence- " << my_peptide_lists.prot_f[i].prot_seq << "    -      total intensity- " << my_peptide_lists.prot_f[i].total << "    -    amount made up of miss cleaved peptides- " << my_peptide_lists.prot_f[i].percentMiss << endl;
    }*/




  /*  for (int i = 0; i < 20; i++) {
        cout << my_peptide_lists.prot_f[i].missPeptides.size();
        cout << "  " << my_peptide_lists.prot_f[i].trypPeptides.size() << endl;
       
    }*/



  
}

void deep_functions::json(peptide_lists& my_peptide_lists, string fn) {
    
   /* sort(my_peptide_lists.prot_f.begin(), my_peptide_lists.prot_f.end(), compareTotal);
    for (int i = 0; i < 500; i++) {
        cout << my_peptide_lists.prot_f[i].total << endl; 
  }
   cout<< my_peptide_lists.prot_f.size() << endl;*/
  StringBuffer s;
  PrettyWriter<StringBuffer> writer(s);
     
  writer.StartObject();
  writer.Key("Proteins");
  writer.StartArray();
  for (int i = 0; i <my_peptide_lists.prot_f.size(); i++) {
    writer.StartObject();
    writer.Key("Protein Name");
            
    writer.String(my_peptide_lists.prot_f[i].prot_seq.c_str());
    writer.Key("Total Intensity");
    writer.Double(my_peptide_lists.prot_f[i].total);
    writer.Key("Percent Miss-Cleaved");
    writer.Double(my_peptide_lists.prot_f[i].percentMiss);
    writer.Key("Tryptic Peptides");
    writer.StartArray();
    for (int j = 0; j < my_peptide_lists.prot_f[i].trypPeptides.size(); j++) {
      writer.StartObject();
      writer.Key("Sequence");
      writer.String(my_peptide_lists.prot_f_real[i][my_peptide_lists.prot_f[i].trypPeptides[j]].pep_seq.c_str());
      writer.Key("Abundance");
      writer.Double(my_peptide_lists.prot_f_real[i][my_peptide_lists.prot_f[i].trypPeptides[j]].areaXIC);
      writer.Key("PSM Counts");
      writer.Int(my_peptide_lists.prot_f_real[i][my_peptide_lists.prot_f[i].trypPeptides[j]].psm_count);
      writer.EndObject();
    }
    writer.EndArray();
    writer.Key("Miss Cleaved Peptides");
    writer.StartArray();
    for (int j = 0; j < my_peptide_lists.prot_f[i].missPeptides.size(); j++) {
      writer.StartObject();
      writer.Key("Sequence");
      writer.String(my_peptide_lists.prot_f_real[i][my_peptide_lists.prot_f[i].missPeptides[j]].pep_seq.c_str());
      writer.Key("Abundance");
      writer.Double(my_peptide_lists.prot_f_real[i][my_peptide_lists.prot_f[i].missPeptides[j]].areaXIC);
      writer.Key("PSM Counts");
      writer.Int(my_peptide_lists.prot_f_real[i][my_peptide_lists.prot_f[i].missPeptides[j]].psm_count); 
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

  //writer.Key("a");
  //writer.StartArray();                // Between StartArray()/EndArray(),
  //for (unsigned i = 0; i < 4; i++)
  //    writer.Uint(i);                 // all values are elements of the array.
  //writer.EndArray();
       

  // {"hello":"world","t":true,"f":false,"n":null,"i":123,"pi":3.1416,"a":[0,1,2,3]}
      

  //FILE* fin;
  //FILE* fout;
  //char str[256];


  //double d;
  //int i;
  //string s;


  //string fname = "C:\\path\\data.txt";  //input file
  //string fname2 = "result.txt";  //output file


  //fin = fopen(fname.c_str(), "rt"); //Try to open input file for (r)eading, and (t)ext
  //if (fin == NULL) {
  //    cout << "Cannot find: " << fname << "\nAttempting to open from CWD." << endl;
  //    string fn = fname.substr(fname.find_last_of('\\') + 1, fname.size()); //if file cannot be opened, try CWD
  //    fin = fopen(fn.c_str(), "rt");
  //    if (fin == NULL) {
  //        cout << "Cannot find: " << fn << "\nYou're SOL. Exiting." << endl;
  //        return -1;
  //    }
  //}

  //int count = 0;
  //while (!feof(fin)) {  //continue reading our file until we reach the end
  //    if (fgets(str, 256, fin) == NULL) continue; //grab one line at a time. If a line could not be read, try again.
  //    if (strlen(str) < 1) continue;
  //    if (count == 0) d = atof(str);  //convert the first value to double-precision
  //    else if (count == 1) i = atoi(str); //convert the second value to integer
  //    else s = str; //last value is a string
  //    count++;
  //}
  //fclose(fin);

  //fout = fopen(fname2.c_str(), "wt"); //open this file for (w)rite, and (t)ext
  //fprintf(fout, "Double-precision, 2 decimal places: %.2lf\n", d);
  //fprintf(fout, "Integer: %d\n", i);
  //fprintf(fout, "String: %s\n", s.c_str());
  //fclose(fout);

    
 
}
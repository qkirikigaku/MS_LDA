#pragma once

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <numeric>
#include <random>
#include <boost/math/special_functions/digamma.hpp>
#include <boost/math/special_functions/gamma.hpp>

using namespace std;

class LDA {
    private:
        int num_topic; int num_vocab; int num_doc;
        int experiment;
        vector<vector<int> > train_doc;
        vector<vector<int> > test_doc;
        int data_type; int threshold;
        string cancer_type;
        vector<vector<vector<double> > > log_qz;
        vector<vector<double> > xi_theta; vector<vector<double> > xi_phi;
        vector<double> alpha; vector<double> beta;
        vector<vector<double> > log_sum_qz_i;
        vector<vector<double> > log_sum_qz_di;
        double temp_vlb;
        double old_vlb;
    public:
        LDA(int x, int y, int z, string w, int v);
        void run_VB();
        void initialize();
        void Update_log_qz();
        void Update_parameter();
        void Update_hyperparameter();
        void calc_vlb();
        void calc_n_dk();
        void calc_n_kv();
        void load_data();
        void write_data();
        void show_vlb();
        void Normalize(vector<double> &vec, int &length);
        double log_sum_exp(vector<double> &vec, int &length);
};

LDA::LDA(int x, int y, int z, string w, int v){
    num_topic = x; data_type = y;
    threshold = z; cancer_type = w;
    experiment = v;
}

void LDA::run_VB(){
    initialize();
    old_vlb = -10e15;
    Update_parameter();
    for (int i=0; i < 300; i++){
        cout << "iter : " << i << endl;
        Update_log_qz();
        Update_parameter();
        Update_hyperparameter();
        calc_vlb();
        show_vlb();
        if(fabs(temp_vlb - old_vlb) < 1){
            break;
        }
        old_vlb = temp_vlb;
        cout << endl;
    }
}

void LDA::initialize(){
    int d,i,k,v;
    random_device rnd;
    log_qz.resize(num_doc);
    for (d=0; d < num_doc; d++){
        log_qz[d].resize(train_doc[d].size());
        for (i=0; i < train_doc[d].size(); i++){
            log_qz[d][i].resize(num_topic, 0);
            for (k=0; k < num_topic; k++){
                log_qz[d][i][k] = (double) -rnd();
            }
            Normalize(log_qz[d][i], num_topic);
        }
    }

    mt19937 mt(rnd());
    uniform_real_distribution<double> Uniform(0.0, 1.0);
    alpha.resize(num_topic);
    for (k=0; k < num_topic; k++){
        alpha[k] = Uniform(mt);
    }

    beta.resize(num_vocab);
    for (v=0; v < num_vocab; v++){
        beta[v] = Uniform(mt);
    }

    xi_theta.resize(num_doc);
    for (d=0; d < num_doc; d++){
        xi_theta[d].resize(num_topic, 0);
    }

    xi_phi.resize(num_topic);
    for (k=0; k < num_topic; k++){
        xi_phi[k].resize(num_vocab, 0);
    }

    log_sum_qz_i.resize(num_doc);
    for (d=0; d < num_doc; d++){
        log_sum_qz_i[d].resize(num_topic, 0);
    }

    log_sum_qz_di.resize(num_topic);
    for (k=0; k < num_topic; k++){
        log_sum_qz_di[k].resize(num_vocab, -730);
    }
}

void LDA::Update_log_qz(){
    int d,v,i,k;
    
    vector<double> sum_xi_theta;
    sum_xi_theta.resize(num_doc, 0);
    for (d=0; d < num_doc; d++){
        sum_xi_theta[d] = 
            accumulate(xi_theta[d].begin(), xi_theta[d].end(), 0.0);
    }
    
    vector<double> sum_xi_phi;
    sum_xi_phi.resize(num_topic, 0);
    for (k=0; k < num_topic; k++){
        sum_xi_phi[k] = accumulate(xi_phi[k].begin(), xi_phi[k].end(), 0.0);
    }
    
    for (d=0; d < num_doc; d++){
        for (i=0; i < train_doc[d].size(); i++){
            for (k=0; k < num_topic; k++){
                log_qz[d][i][k]
                    = boost::math::digamma(xi_phi[k][train_doc[d][i]])
                    + boost::math::digamma(xi_theta[d][k])
                    - boost::math::digamma(sum_xi_phi[k])
                    - boost::math::digamma(sum_xi_theta[d]);
            }
        }
    }
    for (d=0; d < num_doc; d++){
        for (i=0; i < train_doc[d].size(); i++){
            Normalize(log_qz[d][i], num_topic);
        }
    }
}

void LDA::Update_parameter(){
    int d,v,i,k;
    calc_n_dk();
    for (d=0; d < num_doc; d++){
        for (k=0; k < num_topic; k++){
            xi_theta[d][k] = exp(log_sum_qz_i[d][k]) + alpha[k];
        }
    }
    calc_n_kv();
    for (k=0; k < num_topic; k++){
        for (v=0; v < num_vocab; v++){
            xi_phi[k][v] = exp(log_sum_qz_di[k][v]) + beta[v];
        }
    }
}

void LDA::Update_hyperparameter(){
    int d,v,i,k;
    
    calc_n_dk();
    vector<double> new_numerator_alpha;
    new_numerator_alpha.resize(num_topic, 0);
    for (k=0; k < num_topic; k++){
        for (d=0; d < num_doc; d++){
            new_numerator_alpha[k]
                += (boost::math::digamma(exp(log_sum_qz_i[d][k]) + alpha[k])
                    - boost::math::digamma(alpha[k])) * alpha[k];
        }
    }
    double new_denominator_alpha = 0;
    double sum_alpha = accumulate(alpha.begin(), alpha.end(), 0.0);
    for (d=0; d < num_doc; d++){
        new_denominator_alpha
            += boost::math::digamma(train_doc[d].size() + sum_alpha)
            - boost::math::digamma(sum_alpha);
    }
    int count = 0;
    for (k=0; k < num_topic; k++){
        if (new_numerator_alpha[k] == 0){
            count ++;
        }
    }
    if (count == 0){
        for (k=0; k < num_topic; k++){
            alpha[k] = new_numerator_alpha[k] / new_denominator_alpha;
        }
    }

    calc_n_kv();
    vector<double> new_numerator_beta;
    vector<double> new_denominator_beta;
    new_numerator_beta.resize(num_vocab, 0);
    new_denominator_beta.resize(num_vocab, 0);

    for (v=0; v < num_vocab; v++){
        for (k=0; k < num_topic; k++){
            new_numerator_beta[v]
                += (boost::math::digamma(exp(log_sum_qz_di[k][v]) + beta[v])
                    - boost::math::digamma(beta[v])) * beta[v];
        }
    }
    vector<double> sum_qz_div;
    double sum_beta = accumulate(beta.begin(), beta.end(), 0.0);
    sum_qz_div.resize(num_topic, 0);
    for (k=0; k < num_topic; k++){
        for (v=0; v < num_vocab; v++){
            sum_qz_div[k] += exp(log_sum_qz_di[k][v]) + beta[v];
        }
    }
    for (v=0; v < num_vocab; v++){
        for (k=0; k < num_topic; k++){
            new_denominator_beta[v]
                += boost::math::digamma(sum_qz_div[k])
                - boost::math::digamma(sum_beta);
        }
    }
    for (v=0; v < num_vocab; v++){
        beta[v] = new_numerator_beta[v] /new_denominator_beta[v];
    }
}

void LDA::calc_vlb(){
    int d,v,i,k;
    double first_comp = 0, second_comp = 0, third_comp = 0, 
           fourth_comp = 0, fifth_comp = 0;

    double Vbeta = 0;
    double sum_lgamma_beta = 0;
    for (v=0; v < num_vocab; v++){
        Vbeta += beta[v];
        sum_lgamma_beta += boost::math::lgamma(beta[v]);
    }
    for (k=0; k < num_topic; k++){
        double sum_xi_phi_v = 
            accumulate(xi_phi[k].begin(), xi_phi[k].end(), 0.0);
        double log_pi_gamma_xi_phi_v = 0;
        for (v=0; v < num_vocab; v++){
            log_pi_gamma_xi_phi_v += boost::math::lgamma(xi_phi[k][v]);
        }
        first_comp += boost::math::lgamma(Vbeta) - sum_lgamma_beta
                    - boost::math::lgamma(sum_xi_phi_v) + log_pi_gamma_xi_phi_v;
    }

    calc_n_kv();
    for (k=0; k < num_topic; k++){
        double sum_xi_phi_v = 
            accumulate(xi_phi[k].begin(), xi_phi[k].end(), 0.0);
        for (v=0; v < num_vocab; v++){
            second_comp += (exp(log_sum_qz_di[k][v]) + beta[v] - xi_phi[k][v])
                        * (boost::math::digamma(xi_phi[k][v]) 
                        - boost::math::digamma(sum_xi_phi_v));
        }
    }

    double sum_alpha = accumulate(alpha.begin(), alpha.end(), 0.0);
    double log_sum_gamma_alpha = 0;
    for (k=0; k < num_topic; k++){
        log_sum_gamma_alpha += boost::math::lgamma(alpha[k]);
    }
    for (d=0; d < num_doc; d++){
        double sum_xi_theta_k 
            = accumulate(xi_theta[d].begin(), xi_theta[d].end(), 0.0);
        double log_sum_gamma_xi_theta_k = 0;
        for (k=0; k < num_topic; k++){
            log_sum_gamma_xi_theta_k += boost::math::lgamma(xi_theta[d][k]);
        }
        third_comp += boost::math::lgamma(sum_alpha) - log_sum_gamma_alpha
                   - boost::math::lgamma(sum_xi_theta_k) 
                   + log_sum_gamma_xi_theta_k;
    }

    calc_n_dk();
    for (d=0; d<num_doc; d++){
        double sum_xi_theta_k 
            = accumulate(xi_theta[d].begin(), xi_theta[d].end(), 0.0);
        for (k=0; k < num_topic; k++){
            fourth_comp += (exp(log_sum_qz_i[d][k]) + alpha[k] - xi_theta[d][k])
                        * (boost::math::digamma(xi_theta[d][k])
                        - boost::math::digamma(sum_xi_theta_k));
        }
    }

    for (d=0; d < num_doc; d++){
        for (i=0; i < train_doc[d].size(); i++){
            for (k=0; k < num_topic; k++){
                fifth_comp += exp(log_qz[d][i][k]) * log_qz[d][i][k];
            }
        }
    }
    temp_vlb = first_comp + second_comp + third_comp
             + fourth_comp + fifth_comp;
}

void LDA::calc_n_dk(){
    int d,v,i,k;
    vector<vector<vector<double> > > temp_log_qz;
    temp_log_qz.resize(num_doc);
    for (d=0; d < num_doc; d++){
        temp_log_qz[d].resize(num_topic);
        for (k=0; k < num_topic; k++){
            temp_log_qz[d][k].resize(train_doc[d].size(), 0);
        }
        for (i=0; i < train_doc[d].size(); i++){
            for (k=0; k < num_topic; k++){
                temp_log_qz[d][k][i] = log_qz[d][i][k];
            }
        }
    }
    for (d=0; d < num_doc; d++){
        int doc_length = train_doc[d].size();
        for (k=0; k < num_topic; k++){
            log_sum_qz_i[d][k] = log_sum_exp(temp_log_qz[d][k], doc_length);
        }
    }
}

void LDA::calc_n_kv(){
    int d,v,i,k;
    vector<vector<vector<double> > > temp_log_qz;
    temp_log_qz.resize(num_topic);
    for (k=0; k < num_topic; k++){
        temp_log_qz[k].resize(num_vocab);
        for (d=0; d < num_doc; d++){
            for (i=0; i < train_doc[d].size(); i++){
                temp_log_qz[k][train_doc[d][i]].push_back(log_qz[d][i][k]);
            }
        }
    }
    for (k=0; k < num_topic; k++){
        for (v=0; v < num_vocab; v++){
            int vec_size = temp_log_qz[k][v].size();
            log_sum_qz_di[k][v] = log_sum_exp(temp_log_qz[k][v], vec_size);
        }
    }
}

void LDA::load_data(){
    ifstream ifs;
    string input_file_name = "data/data" + to_string(data_type) + "_o"
        + to_string(threshold) + "_" + cancer_type + ".txt";
    ifs.open(input_file_name.c_str(), ios::in);
    if(!ifs){
        cout << "Cannot open " + input_file_name << endl;
        exit(1);
    }
    char buf[1000000];
    char *temp;
    vector<vector<int> > raw_document;
    vector<int> words_number;
    ifs.getline(buf, 1000000);
    temp = strtok(buf, " ");
    num_doc = atoi(temp);
    raw_document.resize(num_doc);
    words_number.resize(num_doc, 0);
    train_doc.resize(num_doc); test_doc.resize(num_doc);
    temp = strtok(NULL, " "); num_vocab = atoi(temp);
    int temp_word_number;
    for (int d=0; d < num_doc; d++){
        ifs.getline(buf, 1000000);
        for (int v=0; v < num_vocab; v++){
            if(v == 0) temp_word_number = atoi(strtok(buf, " "));
            else temp_word_number = atoi(strtok(NULL, " "));
            for (int i=0; i < temp_word_number; i++){
                raw_document[d].push_back(v);
                words_number[d]++;
            }
        }
    }
    for (int d=0; d < num_doc; d++){
        int count = 0;
        train_doc[d].resize(words_number[d]);
        for (int i=0; i < words_number[d]; i++){
            train_doc[d][i] = raw_document[d][i];
            count ++;
        }
    }
    ifs.close();
}

void LDA::write_data(){
    ofstream ofs;
    string output_file_name = "result/data" + to_string(data_type) + "_o" +
        to_string(threshold) + "_" + cancer_type + "_" + to_string(experiment) + 
        "/result_k";
    if(num_topic < 10){
        output_file_name += "0" + to_string(num_topic) + ".txt";
    }
    else{
        output_file_name += to_string(num_topic) + ".txt";
    }
    ofs.open(output_file_name, ios::out);
    ofs << to_string(temp_vlb) << "\n";
    ofs << "0" << "\n";
    int k, v;
    vector<vector<double> > Enkv;
    vector<double> sum_output;
    Enkv.resize(num_topic);
    sum_output.resize(num_topic, 0);
    for (k = 0; k < num_topic; k++) {
        Enkv[k].resize(num_vocab);
        for (v = 0; v < num_vocab; v++){
            sum_output[k] += exp(log_sum_qz_di[k][v]);
        }
        for (v = 0; v < num_vocab; v++) {
            Enkv[k][v] = exp(log_sum_qz_di[k][v]) / sum_output[k];
            ofs << to_string(Enkv[k][v]) << " ";
        }
        ofs << "\n";
    }
	for (k = 0; k < num_topic; k++){
	    ofs << alpha[k] << " ";
	}
	ofs << "\n";
    ofs.close();
}

void LDA::show_vlb(){
    calc_vlb();
    cout << "VLB: " << temp_vlb << endl;
    cout << "Improvement point: " << temp_vlb - old_vlb << endl;
}

void LDA::Normalize(vector<double> &vec, int &length){
    double sum = log_sum_exp(vec, length);
    for (int i=0; i < length; i++){
        vec[i] -= sum;
    }
}

double LDA::log_sum_exp(vector<double> &vec, int &length){
    double max_iter = 0;
    for (int i=1; i < length; i++){
        if(vec[i] > vec[max_iter]) max_iter = i;
    }
    double sum = 0;
    for (int i=0; i < length; i++){
        sum += exp(vec[i] -vec[max_iter]);
    }
    double return_value = vec[max_iter] + log(sum);
    return(return_value);
}

void run_VB_LDA(int num_topic, int num_data, int threshold, string cancer_type);

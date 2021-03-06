#include "LDA.h"

void run_VB_LDA(int num_topic, int data_type,
                int threshold, string cancer_type, int experiment);

int main(int argc,char *argv[]){
    if(argc != 6){
        cout << "The number of argument is invalid." << endl;
        return(0);
    }
    int num_topic = atoi(argv[1])+1;
    int data_type = atoi(argv[2]);
    int threshold = atoi(argv[3]);
    string cancer_type = argv[4];
    int experiment = atoi(argv[5]);
    run_VB_LDA(num_topic, data_type, threshold, cancer_type, experiment);
}

void run_VB_LDA(int num_topic, int data_type, 
                int threshold, string cancer_type, int experiment){
    LDA lda(num_topic, data_type, threshold, cancer_type, experiment);
    lda.load_data();
    lda.run_VB();
    lda.write_data();
}

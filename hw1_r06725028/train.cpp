#include "hmm.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

//define parameter
const int max_state = 6;
const int max_observ = 6;
const int max_seq = 10000;
const int SEQ_len = 50;
const int max_model = 5;

//convert allseq[num] char to index
int toInt(const char observ)
{
  switch(observ){
    case 'A': return 0; break;
    case 'B': return 1; break;
    case 'C': return 2; break;
    case 'D': return 3; break;
    case 'E': return 4; break;
    case 'F': return 5; break;
  }
}

void seqUpdate( HMM *hmm, const char *filename)
{   //load data
    char allseq[max_seq][SEQ_len];
    char oneseq[SEQ_len];
    FILE *fp = fopen( filename , "r" );
    for(int i=0; i<max_seq; i++)
    {
        fscanf(fp, "%s", oneseq);
        for(int j=0; j<SEQ_len; j++)
        {
            allseq[i][j] = (char)oneseq[j];
        }     
    }
    fclose(fp);
    printf("%s\n","load data ok" );


    //need to calculate
    double ga_pi[max_state] = {0.};

    double ga_trans_denominator[max_state] = {0.};
    double ep_trans_numerator[max_state][max_state] = {0.};

    double ga_observ_denominator[max_state] = {0.};
    double ga_observ_numerator[max_state][max_observ] = {0.};

    //accumulate through all samples
    for(int num=0;  num<max_seq; num++)
    {
        double alpha[max_state][SEQ_len] = {0.};
        double beta[max_state][SEQ_len] = {0.};
        double gamma_denominator[SEQ_len] = {0.};
        double epsilon_denominator[SEQ_len-1] = {0.};
        //printf("%d\n", num);
        /*for(int nn=0; nn<SEQ_len-1; nn++){
            printf("%c", allseq[num][nn]);
        }
        printf("%c\n", allseq[num][SEQ_len-1]);*/

        for(int i=0; i<max_state; i++)
        {   //initialize alpha & beta
            alpha[i][0] =  hmm->initial[i]*hmm->observation[toInt((char)allseq[num][0])][i];
            beta[i][SEQ_len-1] = 1.;
            //printf("alpha[%d][0] %f\n", i, alpha[i][0]);
        }

        for(int t=1; t<SEQ_len; t++)
        {  //induction alpha
            for(int j=0; j<max_state; j++)
            {
                for(int i=0; i<max_state; i++)
                {   
                    alpha[j][t] += alpha[i][t-1]*hmm->transition[i][j];  
                }
                alpha[j][t] *= hmm->observation[toInt((char)allseq[num][t])][j];
                //printf("%f\n", alpha[j][t]);
            }
        }
        //printf("alpha[%d][%d]: %lf\n", 3,40,alpha[3][40]);


        for(int t=SEQ_len-2; t>=0; t--)
        {  //induction beta
            for(int i=0; i<max_state; i++)
            {
                for(int j=0; j<max_state; j++)
                {
                    beta[i][t] += (hmm->transition[i][j])*(hmm->observation[toInt((char)allseq[num][t+1])][j])*(beta[j][t+1]);
                }
            }
        }

        for( int t=0; t<SEQ_len; t++)
        {  //count gamma_denominator 
            for( int i=0; i<max_state; i++)
            {
                gamma_denominator[t] += (alpha[i][t])*(beta[i][t]);
            }
        }

        for( int t=0; t< (SEQ_len-1); t++)
        {   //count epsilon denominator
            for( int i=0; i<max_state; i++)
            {
                for( int j=0; j<max_state; j++)
                {
                    epsilon_denominator[t] += alpha[i][t]*hmm->transition[i][j]*hmm->observation[toInt((char)allseq[num][t+1])][j]*beta[j][t+1];
                }
            }
        }


        for(int i=0; i<max_state; i++)
        {   //accumulate ga for pi
            ga_pi[i] +=  alpha[i][0]*beta[i][0] / gamma_denominator[0];
                            
            for(int j=0;  j<max_state; j++)
            {   //accumulate epsilon
                for(int t=0; t<(SEQ_len-1); t++)
                {
                    double epsilon = ((alpha[i][t])*(hmm->transition[i][j])*(hmm->observation[toInt((char)allseq[num][t+1])][j])*(beta[j][t+1])) / epsilon_denominator[t];
                    ep_trans_numerator[i][j] += epsilon;    
                }
            }
        }

        for(int i=0; i<max_state; i++)
        {   //accumulate gamma
            for(int t=0; t<SEQ_len; t++)
            {
                double gamma = alpha[i][t]*beta[i][t] / gamma_denominator[t];
                //printf("gamma %f\n", gamma);
                ga_observ_denominator[i] += gamma;

                if(t != SEQ_len-1)
                {
                    ga_trans_denominator[i] += gamma;
                }
                
                for(int k=0; k<max_observ; k++)
                {
                    if(toInt((char)allseq[num][t]) == k)
                    {
                        ga_observ_numerator[i][k] += gamma;
                    }
                }
            }     
        }          
        //printf("epsilon_denominator %f\n", epsilon_denominator[1]);
    }

    /*printf("ep_trans_numerator %e\n", ep_trans_numerator[1][0]);
    printf("ga_trans_denominator %e\n", ga_trans_denominator[1]);
    printf("ga_observ_numerator %e\n", ga_observ_numerator[1][0]);
    printf("ga_observ_denominator %e\n", ga_observ_denominator[1]);*/

    
    for(int i=0; i<max_state; i++)
    {    //re_estimate
        hmm->initial[i] = ga_pi[i] / max_seq;

        for(int j=0; j<max_state; j++)
        {
            hmm->transition[i][j] = ep_trans_numerator[i][j] / ga_trans_denominator[i];
        }

        for(int k=0; k<max_observ; k++)
        {
            hmm->observation[k][i] = ga_observ_numerator[i][k] / ga_observ_denominator[i];
        }
    }

}//update


int main(int argc, char *argv[])
{   
    //read argument
    int iter = atoi(argv[1]);
    char *init_model_name = (char *)malloc( sizeof(char) * (strlen(argv[2])+1));
    char *data_name = (char *)malloc( sizeof(char) * (strlen(argv[3])+1));
    char *output_model_name = (char *)malloc( sizeof(char) * (strlen(argv[4])+1));
    
    strcpy( init_model_name, argv[2] );
    strcpy( data_name, argv[3] );
    strcpy( output_model_name, argv[4] );

    printf("%s\n", init_model_name);
    
    //initial model
    HMM hmm;
    loadHMM( &hmm, init_model_name );
    printf("%s\n", "load init ok");

    //update all
    for(int i=0; i<iter; i++)
    {
        printf("%s", "iter:");
        printf("%d\n", i);
        seqUpdate(&hmm, data_name);
    }
    

    //dump model
    FILE *fp = fopen(output_model_name ,"w");
    dumpHMM( fp, &hmm);
    fclose(fp);

    return 0;
}

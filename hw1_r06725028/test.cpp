#include "hmm.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

//define parameter
const int max_state = 6;
const int max_observ = 6;
const int max_test = 2500;
const int SEQ_len = 50;
const int max_model = 5;

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

void viterbi( HMM *hmms, const char *test_data, const char *result_file )
{	//load test data
	char allseq[max_test][SEQ_len];
    char oneseq[SEQ_len];
    FILE *fp = fopen( test_data , "r" );
    for(int i=0; i<max_test; i++)
    {
        fscanf(fp, "%s", oneseq);
        for(int j=0; j<SEQ_len; j++)
        {
            allseq[i][j] = (char)oneseq[j];
        }     
    }
    fclose(fp);
    printf("%s\n","load data ok" );

    //open result file
    fp = fopen( result_file , "w" );

	for (int n=0; n<max_test; n++)
	{
		printf("test : %d ...\n", n);
		double likelihood = 0.;
		char *max_model_name;

		for(int m=0;  m<max_model; m++)
		{
			double delta[max_state][SEQ_len] = {0.};

			for(int i=0; i<max_state; i++)
			{	//initalization
				delta[i][0] = (hmms[m].initial[i])*(hmms[m].observation[toInt((char)allseq[n][0])][i]);
			}

			for(int t=0; t<SEQ_len-1; t++)
			{	//recursion
				for(int j=0; j<max_state; j++)
				{
					double max_prob = 0.;
					
					for(int i=0; i<max_state; i++)
					{
						double prob = 0.;
						
						prob = (delta[i][t])*(hmms[m].transition[i][j]);
						if(prob > max_prob)
						{
							max_prob = prob;
						}
					}

					delta[j][t+1] = max_prob*(hmms[m].observation[toInt((char)allseq[n][t+1])][j]);
				}
			}
		
			for(int i=0; i<max_state; i++)
			{	//find max prob with model hmm[m]
				if(delta[i][SEQ_len-1] > likelihood)
				{
					likelihood = delta[i][SEQ_len-1];
					max_model_name = hmms[m].model_name;
				}
			}

		}
		//write output
		fprintf( fp, "%s %e\n", max_model_name, likelihood);
	}
	fclose(fp);
}

int main(int argc, char *argv[])
{	
	//read argument
	char *model_list =  (char *)malloc( sizeof(char) * (strlen(argv[1])+1));
	char *test_data =  (char *)malloc( sizeof(char) * (strlen(argv[2])+1));
	char *result_file =  (char *)malloc( sizeof(char) * (strlen(argv[3])+1));

	strcpy( model_list, argv[1] );
    strcpy( test_data, argv[2] );
    strcpy( result_file, argv[3] );

    printf("%s\n", test_data);
	
	//load model
	HMM hmms[5];
	load_models( model_list, hmms, 5);

	//test
	viterbi(hmms, test_data, result_file);

	return 0;
}

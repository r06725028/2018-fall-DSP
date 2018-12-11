//#include "hmm.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <vector>
#include <set>
#include <map>
#include <algorithm>

#include "Prob.h"
#include "LM.h"
#include "Ngram.h"
#include "VocabMap.h"
#include "Vocab.h"
//#include "LMStats.h"
	
/*VocabIndex wid = voc.getIndex("囧");
if(wid == Vocab_None) {
	printf("No word with wid = %d\n", wid);
	printf("where Vocab_None is %d\n", Vocab_None);
}

wid = voc.getIndex("患者");
VocabIndex context[] = {voc.getIndex("癮") , voc.getIndex("毒"), Vocab_None};
printf("log Prob(患者|毒-癮) = %f\n", lm.wordProb(wid, context));

VocabIndex wid = voc.getIndex("囧");
if(wid == Vocab_None) {
	// replace OOV with <unk>
	wid = voc.getIndex(Vocab_Unknown);
}*/


//define parameter
const int max_state = 1024;
const int max_test = 2500;
const int max_len = 256;


int main(int argc, char *argv[])
{
	//讀參數！(output直接print就好！會自動寫入 > 後接的檔名內)
	char *testTxt =  (char *)malloc( sizeof(char) * (strlen(argv[2])+1));
	char *mapFile = (char *)malloc( sizeof(char) * (strlen(argv[4])+1));
	char *modelFile = (char *)malloc( sizeof(char) * (strlen(argv[6])+1));

	strcpy( testTxt, argv[2] );
    strcpy( mapFile, argv[4] );
    strcpy( modelFile, argv[6] );
    int ngramOrder = atoi(argv[8]);

    //宣告變數
	int ngram_order = 2;
	Vocab voc, zy, b5;//word
	//VocabIndex wid;//word's idx
	VocabMap map(zy, b5);//zy map to b5
	Ngram lm( voc, ngram_order );//count word prob
	//LMStats count(voc);//count sentence length
	

	//讀mapping
	//printf("讀mapping...\n", );
	File mapfile(mapFile, "r");
	map.read(mapfile);
	mapfile.close();

	//讀language model
	//printf("讀language model...\n", );
	File lmFile( modelFile, "r" );
	lm.read(lmFile);
	lmFile.close();
	
	//讀seg_xx.txt
	//printf("讀seg.txt start...\n", );	
	File segFile(testTxt, "r");
	char *line;
	int idx = 0;
	while (line = segFile.getline())
	{
		//宣告變數
		static VocabString words[maxWordsPerLine + 1];//放all word
		static VocabIndex wids[maxWordsPerLine + 3];//放all word's idx
		VocabIndex *start = wids;//放句子起始位置
		unsigned int howmany;//放# of words
	
		//從句子轉為詞陣列，再轉為索引陣列
		//從wids[1]開始放起，預留wids[0]放句子開頭，OOV替換成unkIndex
		howmany = voc.parseWords(line, words + 1, maxWordsPerLine + 1);//放all word
		howmany = voc.getIndices(words + 1, wids + 1, maxWordsPerLine + 1, voc.unkIndex());//放all word's idx
		
		//補上句子開頭和結尾
		if (wids[1] != voc.ssIndex())
		{	//補start
			wids[0] = voc.ssIndex();
			words[0] = "<s>"; 
		}
		else
		{
			start = wids + 1;//原本已有的話，把開頭後移一格
			words[1] = "<s>";
			//words = words + 1;
		}
		if (wids[howmany] != voc.seIndex())
		{	//補end
			wids[howmany + 1] = voc.seIndex();
			wids[howmany + 2] = Vocab_None;
			words[howmany + 1] = "</s>";
		}
		else
		{
			words[howmany] = "</s>";
		}

		//宣告viterbi變數
		LogP delta[max_len][max_state] = {0.};//存log機率值
		VocabIndex psi[max_len][max_state]; //存backtrack pointer
		VocabIndex zyforb5[max_len][max_state];
		int numofb5[max_len];//存每個注音可能對應字數

		//初始化viterbi : δ1(qi) = P(W1=qi)
		Prob prob;
		LogP logp;
		numofb5[0] = 0;
		VocabIndex zy_idx, b5_idx;
		VocabIndex check1 = zy.getIndex(words[0]);
		if (check1 == Vocab_None)
		{
			check1 = voc.getIndex(Vocab_Unknown); 
		}
		VocabMapIter iter(map, check1);
		iter.init();//對第一個注音建立iter
		
		for(int i = 0; iter.next(zy_idx, prob); i++)
		{	//依序找出所有能對應到第一個注音的word
			numofb5[0]++;
			b5_idx = voc.getIndex(b5.getWord(zy_idx));	
			if (b5_idx == Vocab_None)
			{	//replace OOV with <unk>
				b5_idx = voc.getIndex(Vocab_Unknown);
			}
			VocabIndex init_context[] = {Vocab_None};
			logp = lm.wordProb(b5_idx, init_context);
			if (logp == LogP_Zero)
			{	//檢查機率是否為零
				logp =  -100;
			}
			delta[0][i] = logp;
			psi[0][i] = 0;
			zyforb5[0][i] = zy_idx;	
		}

		//viterbi inductive
		for(int t = 1; t < howmany+2; t++)
		{
			Prob prob;
			//LogP bi_logp, uni_logp;
			numofb5[t] = 0;
			VocabIndex zy_idx, b5_idx;
			VocabIndex check2 = zy.getIndex(words[t]);
			if (check2 == Vocab_None)
			{
				check2 = voc.getIndex(Vocab_Unknown);
			}
			VocabMapIter iter(map, check2);
			iter.init();
			
			for(int j = 0; iter.next(zy_idx, prob); j++)
			{
				numofb5[t]++;
				b5_idx = voc.getIndex(b5.getWord(zy_idx));
				if (b5_idx == Vocab_None)
				{	//replace OOV with <unk>
					b5_idx = voc.getIndex(Vocab_Unknown);
				}
				LogP max_logp =  -1000.0;//找出有最大機率的bigram
				for(int i = 0; i < numofb5[t-1]; i++)
				{
					VocabIndex tmp = voc.getIndex(b5.getWord(zyforb5[t-1][i]));
					if (tmp == Vocab_None)
					{	//replace OOV with <unk>
						tmp = voc.getIndex(Vocab_Unknown);
					}
					VocabIndex bi_context[] = {tmp, Vocab_None};
					LogP bi_logp = lm.wordProb(b5_idx, bi_context);//bigram prob
					VocabIndex uni_context[] = 	{Vocab_None};
					LogP uni_logp = lm.wordProb(b5_idx, uni_context);//unigram prob
					
					if (bi_logp == LogP_Zero && uni_logp == LogP_Zero)
					{
						bi_logp =  -100;
					}

					bi_logp += delta[t-1][i];		
			
					if(bi_logp > max_logp)
					{
						max_logp = bi_logp;
						psi[t][j] = i;
					}
				}
				delta[t][j] = max_logp;
				zyforb5[t][j] = zy_idx;
			}
		}

		//找出有最大機率的路徑
		LogP max_path =  -1000.0;
		int last_state = -1;
		for(int i = 0; i < numofb5[howmany+1]; i++)
		{
			if(delta[howmany+1][i] > max_path)
			{
				max_path = delta[howmany+1][i];
				last_state = i;
			}
		}

		//backtrack!!!
		VocabString backtrack[max_len];
		backtrack[0] = "<s>";
		backtrack[howmany+1] = "</s>";

		for (int i = howmany+1; i >0; i--) 
		{
		    backtrack[i] = b5.getWord(zyforb5[i][last_state]);
		    last_state = psi[i][last_state];
		}

		//output
		for (int t = 0; t < howmany+2; t++)
		{
			if(t == howmany+1)
			{
				printf("%s%s", backtrack[t], "\n");
			}
			else
			{
				printf("%s%s", backtrack[t], " ");
			}
		}
		
	}
	segFile.close();
	return 0;
}
	

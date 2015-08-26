#include <stdio.h>
#include <stdlib.h>
#include <string.h>


//procedure to get the number of segments of S and T, and the size of \calS and \calT
//using popen and calling shell commands.
void GetNumbers(char Sfile[], char Tfile[], int *Ssegs, int *Tsegs, int *Ssize, int *Tsize){
	FILE *fp;
	char command[200];

	sprintf(command, "head -n+1 %s | wc -c", Sfile);
	fp = popen(command, "r");
	fscanf(fp, "%d", &*Ssize);
	pclose(fp);

	sprintf(command, "head -n+1 %s | wc -c", Tfile);
	fp = popen(command, "r");
	fscanf(fp, "%d", &*Tsize);
	pclose(fp);


	//getting the number of segments from S
	sprintf(command, "tail -n+2 %s | wc -l", Sfile);
	fp = popen(command, "r");
	fscanf(fp, "%d", &*Ssegs);
	pclose(fp);

	//getting the number of segments from T
	sprintf(command, "tail -n+2 %s | wc -l", Tfile);
	fp = popen(command, "r");
	fscanf(fp, "%d", &*Tsegs);
	pclose(fp);

}

/*Get T sequence*/
void GetTSequence(char Tsequence[], char Tfile[], int Tstart[], int Tend[], int Tsize){
	char *sequence;
	int i, j, nSegs;
	FILE *T;

	j = nSegs = 0;
	sequence = malloc(Tsize*sizeof(char));

	T = fopen(Tfile, "r");

	fscanf(T, "%s", sequence);

	while(fscanf(T, "%d", &Tstart[nSegs]) != EOF){
		fscanf(T, "%d", &Tend[nSegs]);

		for (i = Tstart[nSegs]-1; i < Tend[nSegs]; ++i)
		{
			Tsequence[j] = sequence[i];
			j++;
		}
		nSegs++;
	}

	fclose(T);

	if(sequence != NULL) free(sequence);
}

/*Get the consensus sequence from S*/
void GetConsensus(char consensus[], char Sfile[], int Sstart[], int Send[], int Ssize){
	int indexS, indexC, nSegs, max;
	char *sequence;
	FILE *S, *fp;

	S = fopen(Sfile, "r");
	sequence = malloc(Ssize*sizeof(char));

	indexC = indexS = 0;
	nSegs = 0;

	fscanf(S, "%s", sequence);
	fscanf(S, "%d", &Sstart[nSegs]);
	fscanf(S, "%d", &Send[nSegs]);

	indexS = Sstart[nSegs] - 1;
	while(indexS < Send[nSegs]){
		consensus[indexC] = sequence[indexS];
		indexC++;
		indexS++;
	}
	max = Send[nSegs];
	nSegs++;

	while(fscanf(S, "%d", &Sstart[nSegs]) != EOF){
		fscanf(S, "%d", &Send[nSegs]);

		if(Send[nSegs] > max){
			if(Sstart[nSegs] < max){           //overlap
				indexS = max;
				while(indexS < Send[nSegs]){
					consensus[indexC] = sequence[indexS];
					indexC++;
					indexS++;
				}
			}else{
				indexS = Sstart[nSegs] - 1;    //not overlap
				while(indexS < Send[nSegs]){
					consensus[indexC] = sequence[indexS];
					indexC++;
					indexS++;
				}
			}
			max = Send[nSegs];
		}
		nSegs++;

	}
	consensus[indexC] = '\0';


	fclose(S);
	if(sequence != NULL) free(sequence);
}

//read the alignment of consensus and T in the alignment.fa file
void GetAlignment(char consensusAligned[], char TAligned[]){
	char sequence[500], nucl;
	int i;
	FILE *fp;

	fp = fopen("alignment.fa", "r");
	fscanf(fp, "%s", sequence);

	i = 0;

	nucl = fgetc(fp);
	while(nucl != '>'){
		if(nucl > 32 && nucl < 123){
			consensusAligned[i] = nucl;
			i++;
		}
		nucl = fgetc(fp);
	}

	consensusAligned[i] = '\0';

	i = 0;
	fscanf(fp, "%s", sequence);

	nucl = fgetc(fp);
	while(nucl != EOF){
		if(nucl > 32 && nucl < 123){
			TAligned[i] = nucl;
			i++;
		}
		nucl = fgetc(fp);
	}

	TAligned[i] = '\0';

	fclose(fp);
}

/*here, the vector 'gaps' is filled. gaps[i] keeps the number of gaps
  previous the i-th nucleotide from S in the alignment with T*/
void DeterminedGaps(int gaps[], char Conaligned[]){
	int i, j, alignsize, numberGaps;

	alignsize = strlen(Conaligned);

	numberGaps = j = 0;
	for (i = 0; i < alignsize; ++i)
	{
		gaps[j] = 0;
		if(Conaligned[i] == '-'){
			numberGaps++;
		}else{
			gaps[j] = numberGaps;
			j++;
		}
	}
}

//changes the coordinates so the score of each segmente can be determined
void ChangeSCoordinatesAlignment(int Sstart[], int Send[], char SAligned[], int numSSegs){
	int i, max, diff, len, *gaps;

	len = strlen(SAligned);

 //step 1: coordinates from \cal{S} and \cal{T} are converted to coordinates from S and T
	diff = Sstart[0] - 1;
	Sstart[0] = Sstart[0] - diff;
	Send[0] = Send[0] - diff;
	max = Send[0];

	for (i = 1; i < numSSegs; ++i)
	{
		Sstart[i] = Sstart[i] - diff;
		Send[i] = Send[i] - diff;

		if(Sstart[i] > max){							 //do not overlap
			diff = diff + Sstart[i] - max - 1;
			Send[i] = Send[i] - Sstart[i];
			Sstart[i] = max + 1;
			Send[i] = Send[i] + Sstart[i];
			max = Send[i];
		}else{
			if(max < Send[i]) max = Send[i];
		}
	}

 //step 2: consider the gaps in the alignment, shifting the positions.
	gaps = malloc(len*sizeof(int));

	DeterminedGaps(gaps, SAligned);
	
	for (i = 0; i < numSSegs; ++i)
	{
		Sstart[i] = Sstart[i] + gaps[Sstart[i]-1];
		Send[i] = Send[i] + gaps[Send[i]-1];

	}

	if(gaps != NULL) free(gaps);
}

//here is calculated the weights for each vertex and the DAG is created.
void GenerateGraph(char Taligned[], char Conaligned[], int numberReads, int Sstart[], int Send[]){
	int i, j, k, score, *scores;

	scores = malloc(numberReads*sizeof(int));
	k = 0;

		//changing the coordinates acconding to the number of gaps in the alignment
		//start = start + gaps[start-1];
		//end = end + gaps[end-1];
		/*start = start - positions[0] + 1;
		end = end - positions[0] + 1;
		if(start > end0){
			start = end0+1;
		}
		*/

	for (j = 0; j < numberReads; ++j){
		score = 0;

		//determing the weights of each segment
		for (i = Sstart[j]-1; i < Send[j]; ++i)
		{
			if(Taligned[i] == Conaligned[i]){
				score = score + 1;
			}else{
				if(Taligned[i] == '-' || Conaligned[i] == '-'){
					score = score - 2;
				}else{
					score = score - 1;
				}
			}
		}
		scores[k] = score;
		k++;
	}

	
	//first print the number of vertices
	printf("%d\n", numberReads);

	//then print the graph
	for (i = 0; i < numberReads; i++)
	{
		for (j = 0; j < numberReads; j++)
		{
			if(Send[i] < Sstart[j]) printf("1 ");
			else printf("0 ");
		}
		printf("\n");
	}

	//and now the weights of each vertex, in order.
	for (i = 0; i < numberReads; ++i)
	{
		printf("%d\n", scores[i]);
	}

	if(scores != NULL) free(scores);
}

int main(int argc, char *argv[]){
	int i, j, Ssegs, Tsegs, *Sstart, *Send, *Tstart, *Tend, Ssize, Tsize;
	char *consensus, command[500], *Tsequence, nucl;
	FILE *fp;

	GetNumbers(argv[1], argv[2], &Ssegs, &Tsegs, &Ssize, &Tsize);
	//the aligned sequences will be kept in the same variables
	consensus = malloc((Ssize+Tsize)*sizeof(char));
	Tsequence = malloc((Ssize+Tsize)*sizeof(char));
	Sstart = malloc(Ssegs*sizeof(int));
	Tstart = malloc(Tsegs*sizeof(int));
	Send = malloc(Ssegs*sizeof(int));
	Tend = malloc(Tsegs*sizeof(int));

//getting the consensus from S and T sequences
	GetConsensus(consensus, argv[1], Sstart, Send, Ssize);
	GetTSequence(Tsequence, argv[2], Tstart, Tend, Tsize);

//writting sequences
	fp = fopen("consensus.fa", "w");
	fprintf(fp, ">consensus\n%s", consensus);
	fclose(fp);

	fp = fopen("Tconcatenado.fa", "w");
	fprintf(fp, ">Tconcatenado\n%s", Tsequence);
	fclose(fp);


//calling needle from emboss to align the T sequence with the consensus
	sprintf(command, "needle -asequence consensus.fa -bsequence Tconcatenado.fa -gapopen 2 -gapextend 2 -outfile stdout -endweight Y -endopen 2 -endextend 2 -datafile scoring_matrix -aformat3 A2M > alignment.fa");
	system(command);

	GetAlignment(consensus, Tsequence);

	ChangeSCoordinatesAlignment(Sstart, Send, consensus, Ssegs);

	GenerateGraph(Tsequence, consensus, Ssegs, Sstart, Send);

	//memory desallocation
	if(consensus != NULL) free(consensus);
	if(Tsequence != NULL) free(Tsequence);
	if(Sstart != NULL) free(Sstart);
	if(Send != NULL) free(Send);
	if(Tstart != NULL) free(Tstart);
	if(Tend != NULL) free(Tend);

	return 0;
}
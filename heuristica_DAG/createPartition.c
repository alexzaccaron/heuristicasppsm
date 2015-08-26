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

void GetCoordinates(char fileS[], char fileT[], int Sstart[], int Send[], int Tstart[], int Tend[]){
	int coord, i;
	char sequence[500000];
	FILE *fp;

	i = 0;
	fp = fopen(fileS, "r");
	fscanf(fp, "%s", sequence);

	while(fscanf(fp, "%d", &coord) != EOF){
		Sstart[i] = coord;

		fscanf(fp, "%d", &coord);
		Send[i] = coord;

		i++;
	}

	i = 0;
	fp = fopen(fileT, "r");
	fscanf(fp, "%s", sequence);

	while(fscanf(fp, "%d", &coord) != EOF){
		Tstart[i] = coord;

		fscanf(fp, "%d", &coord);
		Tend[i] = coord;

		i++;
	}
}

void DeterminedGaps(int gaps[], char Conaligned[]){
	int i, j, alignsize, numberGaps;

	alignsize = strlen(Conaligned);

	numberGaps = j = 0;
	for (i = 0; i < alignsize; ++i)
	{
		if(Conaligned[i] == '-'){
			numberGaps++;
		}else{
			gaps[j] = numberGaps;
			j++;
		}
	}
}

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

		//printf("%d %d\n", Sstart[i], Send[i]);
	}

	if(gaps != NULL) free(gaps);
}

//procedure to perform the union operation
void Union(int Tsubset[], int value, int *nelements){
	int i, flag;

	flag = 1;

	for (i = 0; i < *nelements; ++i)
	{
		if(Tsubset[i] == value) flag = 0;
	}

	if(flag){
		Tsubset[*nelements] = value;
		*nelements = *nelements + 1;
	}

}

//procedure to print the part si of S and the subset ti of T, where the first number
//is the number of elements in each.
void PrintPair(int Spartition[], int Tsubset[], int Selements, int Telements){
	int i;

	printf("%d", Selements);
	for (i = 0; i < Selements; ++i)
	{
		printf(" %d", Spartition[i]);
	}
	printf("\n%d", Telements);
	for (i = 0; i < Telements; ++i)
	{
		printf(" %d", Tsubset[i]);
	}
	printf("\n");
}
	

void DefineSTPartitions(char Sfile[], int Sstart[], int Send[], int Tstart[], int Tend[], int numSSegs, int numTSegs){
	FILE *fp;
	int i, Telements, Tsubset[1000], Selements, Spartition[1000], Tindex, n, j;

	fp = fopen(Sfile, "r");

	Telements = 0;     //number of elements in the T subset
	Selements = 0;	   //number of elements in the S part

	//reading the paths found in G
	while(fscanf(fp, "%d", &n)!=EOF){

		//n==0 means end of a path, (si, ti) is printed so the next path can be processed
		if(n == 0){
			PrintPair(Spartition, Tsubset, Selements, Telements);

			Telements = 0;
			Selements = 0;
		}else{
			Spartition[Selements] = n;
			Selements++;
			//for each segment from si, check which segments from T it covers.
			for (i = 0; i < numTSegs; ++i)
			{
				//checking the possibilities of overlaping
				if((Tstart[i] <= Sstart[n-1] && Sstart[n-1] <= Tend[i]) ||
				   (Tstart[i] <= Send[n-1] && Send[n-1] <= Tend[i]) ||
				   (Tstart[i] >= Sstart[n-1] && Send[n-1] >= Tend[i])){
				   	//Union operation to avoid repetitions in the subset ti.
					Union(Tsubset, i+1, &Telements);
				}
			}
		}

	}


	fclose(fp);
}

int main(int argc, char *argv[]){
	char *consensusAligned, *TAligned;
	int *Sstart, *Send, *Tstart, *Tend, Ssegs, Tsegs, Ssize, Tsize, i;

	GetNumbers(argv[1], argv[2], &Ssegs, &Tsegs, &Ssize, &Tsize);
	consensusAligned = malloc((Ssize+Tsize)*sizeof(char));
	TAligned = malloc((Ssize+Tsize)*sizeof(char));
	Sstart = malloc(Ssegs*sizeof(int));
	Tstart = malloc(Tsegs*sizeof(int));
	Send = malloc(Ssegs*sizeof(int));
	Tend = malloc(Tsegs*sizeof(int));

	GetAlignment(consensusAligned, TAligned);

	GetCoordinates(argv[1], argv[2], Sstart, Send, Tstart, Tend);
	
	ChangeSCoordinatesAlignment(Tstart, Tend, TAligned, Tsegs);

	ChangeSCoordinatesAlignment(Sstart, Send, consensusAligned, Ssegs);	
	
	DefineSTPartitions(argv[3], Sstart, Send, Tstart, Tend, Ssegs, Tsegs);

	//memory desallocation
	if(consensusAligned != NULL) free(consensusAligned);
	if(TAligned != NULL) free(TAligned);
	if(Sstart != NULL) free(Sstart);
	if(Send != NULL) free(Send);
	if(Tstart != NULL) free(Tstart);
	if(Tend != NULL) free(Tend);


	return 0;
}
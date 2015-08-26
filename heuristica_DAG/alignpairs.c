#include <stdio.h>
#include <stdlib.h>
#include <string.h>

//procedure to get the number of segments of S and T, and the size of \cal{S} and \cal{T}
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

	sprintf(command, "tail -n+2 %s | wc -l", Sfile);
	fp = popen(command, "r");
	fscanf(fp, "%d", &*Ssegs);
	pclose(fp);

	sprintf(command, "tail -n+2 %s | wc -l", Tfile);
	fp = popen(command, "r");
	fscanf(fp, "%d", &*Tsegs);
	pclose(fp);

}

//procedure to get the coordinates of the segments from S and T
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

//procedure the get the \cal{S} and \cal{T} sequences
void GetSequences(char Sfile[], char Tfile[], char Ssequence[], char Tsequence[]){
	FILE *fp;

	fp = fopen(Sfile, "r");
	fscanf(fp, "%s", Ssequence);
	fclose(fp);

	fp = fopen(Tfile, "r");
	fscanf(fp, "%s", Tsequence);
	fclose(fp);
}

//procedure to print a segment from a string 
void PrintSegment(FILE *fp, int startPos, int endPos, char Ssequence[]){
	int i;

	for (i = startPos; i <= endPos; ++i)
	{
		fprintf(fp, "%c", Ssequence[i]);
	}
}


int main(int argc, char *argv[]){
	int *Sstart, *Send, *Tstart, *Tend;
	int Ssegs, Tsegs, Ssize, Tsize, i, nsegs, seg, k;
	char *Ssequence, *Tsequence;
	FILE *fp;

	//memory allocation
	GetNumbers(argv[1], argv[2], &Ssegs, &Tsegs, &Ssize, &Tsize);
	Ssequence = malloc(Ssize*sizeof(char));
	Tsequence = malloc(Tsize*sizeof(char));
	Sstart = malloc(Ssegs*sizeof(int));
	Tstart = malloc(Tsegs*sizeof(int));
	Send = malloc(Ssegs*sizeof(int));
	Tend = malloc(Tsegs*sizeof(int));
	k = 1;

	
	GetCoordinates(argv[1], argv[2], Sstart, Send, Tstart, Tend);
	GetSequences(argv[1], argv[2], Ssequence, Tsequence);

	//each iteration read a pair (si, ti)
	while(scanf("%d", &nsegs) == 1){

		//print si^{\bullet} in Spair file
		fp = fopen("Spair", "w");
		fprintf(fp, ">s%d\n", k);
		for (i = 0; i < nsegs; ++i)
		{
			scanf("%d", &seg);
			PrintSegment(fp, Sstart[seg-1]-1, Send[seg-1]-1, Ssequence);
		}
		fclose(fp);

		//print ti^{\bullet} in Tpair file
		scanf("%d", &nsegs);
		fp = fopen("Tpair", "w");
		fprintf(fp, ">t%d\n", k);
		for (i = 0; i < nsegs; ++i)
		{
			scanf("%d", &seg);
			PrintSegment(fp, Tstart[seg-1]-1, Tend[seg-1]-1, Tsequence);
		}
		fclose(fp);

		//call needle to align si^{\bullet} and ti^{\bullet}
		system("needle -asequence Spair -bsequence Tpair -gapopen 2 -gapextend 2 -outfile stdout -endweight Y -endopen 2 -endextend 2 -datafile scoring_matrix");
		// -aformat3 score | head -n+1 | sed 's/(//g' | sed 's/)//g'

		k++;
	}


	//memory desallocation
	if(Ssequence != NULL) free(Ssequence);
	if(Tsequence != NULL) free(Tsequence);
	if(Sstart != NULL) free(Sstart);
	if(Send != NULL) free(Send);
	if(Tstart != NULL) free(Tstart);
	if(Tend != NULL) free(Tend);

	return 0;
}

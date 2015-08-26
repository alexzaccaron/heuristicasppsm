#include "inout.h"
#include <string.h>

/*Função que conta o número de caracteres presentes no arquivo de entrada, assumindo que o arquivo de entrada é um arquivo FASTA*/
int contaCharFasta (char *file)
{
	int size;
	static char line[FASTA_MAXLINE];
	int read;
	FILE *f;

	f = fopen(file, "r");
	if(!f)
	{
		printf("Error! Cannot Open File %s.", file);
		exit(EXIT_FAILURE);
	}


	size = 0;
	do {
        read = fscanf(f, "%[^\n]", line);
        fscanf(f, "%*c");

        if (read > 0) {
            if(line[0] != '>')
            {
                size += strlen(line);
            }
        }
	}while(!feof(f));

	fclose(f);
	return size;
}

/*Função que recebe um arquivo e seu tamanho, e grava a sequência deste arquivo em seq. Nota: para simplificar, seq conterá um espaço na primeira posição (seq[0]=' '). Seq não termina com \0!*/
void leArquivoFasta(char *arq, int tam, char **seq)
{
	FILE *f;
	static char line[FASTA_MAXLINE];
	int read;
	int  i, pos;

	(*seq) = (char *) malloc (sizeof(char) * (tam));
	eval(*seq);

	f = fopen(arq, "r");
	if(!f)
	{
		printf("Error! Cannot Open File %s.", arq);
		exit(EXIT_FAILURE);
	}

	(*seq)[0] = ' ';
	pos = 1;


	do {
	    read = fscanf(f, "%[^\n]", line);
	    fscanf(f, "%*c");

        if (read > 0) {
            if(line[0] != '>')
            {
                for(i = 0; line[i] !='\0' && i < FASTA_MAXLINE ; i++)
                {
                    (*seq)[pos++] = line[i];
                }
            }
        }
	}while (!feof(f));

	if(pos!=tam)
	{
		printf("Fatal error: %d chars could not be read in %s.\nTerminating program.", tam-pos, arq);
		exit(EXIT_FAILURE);
	}

	fclose(f);
}

/*função que compara dois segmentos para ordenação no qsort()*/
int compare ( const void * elem1, const void * elem2 ) {
    struct segmento *e1 = (struct segmento*) elem1;
    struct segmento *e2 = (struct segmento*) elem2;

    if (e1->first < e2->first)
        return -1;

    if (e1->first > e2->first)
        return 1;

    /*aqui, e1->first == e2->first*/
    if (e1->last < e2->last)
        return -1;

    if (e1->last > e2->last)
        return 1;

    return 0;
}

/*Função que lê os segmentos de um arquivo arq e grava em segs. Retorna a qtde de segmentos lidos*/
int leSegmentos (char *arq, struct segmento **segs)
{
	int amount;
	FILE *in;
	int i, fst, lst;
	int read;

	in = fopen(arq, "r");
	if (!in)
	{
		printf("Error! Cannot Open File %s.", arq);
		exit(EXIT_FAILURE);
	}

	amount = 1;

	do {
	    read = fscanf(in, "%d..%d,", &fst, &lst);

	    if (read>0)
            amount++;
	}while (!feof(in));

    (*segs) = (struct segmento*) malloc(sizeof(struct segmento) * amount);
    eval(*segs);


	(*segs)[0].first = 0;
	(*segs)[0].last = 0;
	(*segs)[0].tam = 1;

	rewind(in);
	for(i = 1; i < amount; i++)
	{
		fscanf(in, "%d..%d,", &fst, &lst);
		(*segs)[i].first = fst;
		(*segs)[i].last = lst;
		(*segs)[i].tam = lst-fst+1;
	}

	qsort(*segs, amount, sizeof(struct segmento), compare);

	for(i = 0; i < amount; i++)
	{
	    (*segs)[i].number = i;
	}

	fclose(in);
	return amount;
}

/*Função que lê as duas sequências*/
void readSeqs (char *arq1, char *arq2) {
    static char temp[512];

    if (INFO) printf("Reading sequence %s.\n", arq1);
    n = contaCharFasta(arq1)+1; /*+1 por causa do ' '*/
    leArquivoFasta(arq1, n, &s1);
    if (INFO) printf("Sequence read successfully!\n");

    if (INFO) printf("Reading sequence %s.\n", arq2);
    m = contaCharFasta(arq2)+1;
    leArquivoFasta(arq2, m, &s2);
    if (INFO) printf("Sequence read successfully!\n");

    strcpy(temp, arq1);
    strcat(temp, ".predex");
    if (INFO) printf("Reading segments %s.\n", temp);
    u = leSegmentos(temp, &B);
    if (INFO) printf("Segments read successfully!\n");

    strcpy(temp, arq2);
    strcat(temp, ".predex");
    if (INFO) printf("Reading segments %s.\n", temp);
    v = leSegmentos(temp, &C);
    if (INFO) printf("Segments read successfully!\n");
}

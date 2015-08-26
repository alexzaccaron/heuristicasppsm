#ifndef INOUT_H_
#define INOUT_H_

#include <stdio.h>
#include <stdlib.h>
#include "basics.h"

#define FASTA_MAXLINE 512

/*Função que conta o número de caracteres presentes na sequência do arquivo passado por parâmetro, assumindo que o arquivo de entrada é um arquivo FASTA*/
int contaCharFasta (char *arq);

/*Função que recebe um arquivo e seu tamanho, e grava a sequência deste arquivo em seq. Nota: para simplificar, seq conterá um espaço na primeira posição (seq[0]=' '). Seq não termina em '\0'!*/
void leArquivoFasta(char *arq, int tam, char **seq);

/*função que compara dois segmentos para ordenação no qsort()*/
int compare ( const void * elem1, const void * elem2 );

/*Função que lê os segmentos de um arquivo arq e grava em segs. Retorna a qtde de segmentos lidos*/
int leSegmentos (char *arq, struct segmento **segs);

/*Função que lê as duas sequências*/
void readSeqs (char *seq1, char *seq2);

#endif /* INOUT_H_ */

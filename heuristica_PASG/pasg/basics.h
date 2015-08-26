#ifndef BASICS_H_
#define BASICS_H_

#include <stdio.h>

#define SPACESCORE -2
#define MATCHSCORE 1
#define MISMATCHSCORE -1

#define INFO 1
#define COMPARE 0
#define SYSTEMINFO 0
#define RELATORIO 0
#define AVALIACAO 0
#define BACKTRACKCELLS 0

char *s1; /*representa a sequência s_1*/
int n;
char *s2; /*representa a sequência s_2*/
int m;


int *gammaB; /*representa a cadeia de segmentos Gamma_B*/
int tamGammaB;
int *gammaC; /*representa a cadeia de segmentos Gamma_C*/
int tamGammaC;

/*struct que representa um segmento*/
struct segmento
{
	int first;
	int last;
	int number;
	int tam;
};

struct segmento *B; /*representa o conjunto B de segmentos de s_1*/
int u;
struct segmento *C; /*representa o conjunto C de segmentos de s_2*/
int v;

/*representa uma célula da matriz de alinhamento M*/
struct cell {
    int value;
    struct cell *before;
    int numSegmentos;
    int k;
    int l;
};

/*representa uma célula da matriz temporária para cálculo da matriz de alinhamento*/
struct tempcell {
    int value;
    struct cell *before;
    int numSegmentos;
};

/*retorna o score de dois caracteres*/
int score (char c1, char c2);

/*retorna o caractere b_k[i]*/
char bki(int k, int i);

/*retorna o caractere c_l[j]*/
char clj(int l, int j);

/*verifica se b_k' precede b_k*/
int precedek (int klinha, int k);

/*verifica se c_l' precede c_l*/
int precedel (int llinha, int l);

/*devolve |b_k'|, usado na equação.*/
int lastLinha(int k);

/*devolve |c_l'|, usado na equação.*/
int lastColuna(int l);

/*verifica se mem é null ou não*/
void eval(void* mem);

#endif /* BASICS_H_ */

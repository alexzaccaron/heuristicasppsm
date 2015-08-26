#include "basics.h"

#include <stdlib.h>

/*retorna o score de dois caracteres*/
int score (char c1, char c2) {
    if (c1!='A' && c1!='C' && c1!='G' && c1!='T' && c1!='-') {
        printf("Error! A char different from ACGT is in input: %c.\n", c1);
        exit(EXIT_FAILURE);
    }
    if (c2!='A' && c2!='C' && c2!='G' && c2!='T' && c2!='-') {
        printf("Error! A char different from ACGT is in input: %c.\n", c2);
        exit(EXIT_FAILURE);
    }

    if (c1=='-' && c2=='-') {
        printf("Error: cannot align 2 spaces.");
        exit(EXIT_FAILURE);
    }

    if (c1=='-' || c2=='-')
        return SPACESCORE;

    if(c1 == c2)
        return MATCHSCORE; /* se iguais retorna valor de match */
    else
       return MISMATCHSCORE; /* se diferentes retorna valor de mismatch */
}

/*retorna o caractere b_k[i]*/
char bki(int k, int i) {
    return s1[B[k].first+i];
}

/*retorna o caractere c_l[j]*/
char clj(int l, int j) {
    return s2[C[l].first+j];
}

/*verifica se b_k' precede b_k*/
int precedek (int klinha, int k) {
    if (B[klinha].last < B[k].first)
        return 1;
    else
        return 0;
}

/*verifica se c_l' precede c_l*/
int precedel (int llinha, int l) {
    if (C[llinha].last < C[l].first)
        return 1;
    else
        return 0;
}

/*devolve |b_k'|, usado na equação.*/
int lastLinha(int k) {
    return B[k].tam-1;
}

/*devolve |c_l'|, usado na equação.*/
int lastColuna(int l) {
    return C[l].tam-1;
}

/*verifica se mem é null ou não*/
void eval(void* mem) {
    if (mem == NULL) {
        printf("Memory error! Terminating program!");
        exit(EXIT_FAILURE);
    }
}

#include "solution.h"

#include "basics.h"
#include <stdlib.h>
#include <limits.h>

/*inicializa a matriz*/
void initMatrix() {
    int i, j, k, l;
    int maxTamB; /*o maior tamanho de um segmento de B*/
    int maxTamC; /*o maior tamanho de um segmento de C*/

    if (INFO) printf("Allocating matrix...\n");

    /*aloca a temp matrix*/
    maxTamB=0;
    for (k=1;k<u;k++) {
        if (B[k].tam > maxTamB)
            maxTamB = B[k].tam;
    }
    maxTamC=0;
    for (l=1;l<v;l++) {
        if (C[l].tam > maxTamC)
            maxTamC = C[l].tam;
    }
    /*aloca as linhas: maxTamB linhas*/
    tempmatrix = (struct tempcell**) malloc(sizeof(struct tempcell*) * maxTamB);
    eval(tempmatrix);

    /*aloca as colunas: maxTamC colunas*/
    for (i=0; i<maxTamB; i++) {
        tempmatrix[i] = (struct tempcell*) malloc(sizeof(struct tempcell) * maxTamC);
        eval(tempmatrix[i]);
    }


    /*aloca a matrix*/
    /*aloca as linhas: |B| linhas*/
    matrix = (struct cell****) malloc (sizeof(struct cell***) * u);
    eval(matrix);

    /*aloca as colunas: |C| colunas*/
    for (k=0; k<u; k++) {
        matrix[k]=(struct cell***) malloc (sizeof(struct cell**) * v);
        eval(matrix[k]);

        /*aloca as submatrizes de alinhamento*/
        for (l=0;l<v;l++) {
            /*aloca as linhas da submatriz*/
            matrix[k][l] = (struct cell**) malloc(sizeof(struct cell*) * B[k].tam);
            eval(matrix[k][l]);

            /*aloca as colunas da submatriz*/
            for (i=0;i<B[k].tam;i++) {
                if (i==B[k].tam-1) {
                    matrix[k][l][i] = (struct cell*) malloc (sizeof(struct cell) * C[l].tam);
                }else {
                    matrix[k][l][i] = (struct cell*) malloc (sizeof(struct cell));
                }

                eval(matrix[k][l][i]);
            }
        }
    }
    if (INFO) printf("Matrix allocated successfully!\n");

    if (INFO) printf("Initializing matrix...\n");
    /*com a matriz alocada, a matriz é inicializada (as três condições base são implementadas)*/
    /*condição base 1*/
    matrix[0][0][0][0].value = 0;
    matrix[0][0][0][0].numSegmentos = 0;
    matrix[0][0][0][0].before = NULL;
    matrix[0][0][0][0].k = 0;
    matrix[0][0][0][0].l = 0;

    /*condição base 2*/
    for (l=1;l<v;l++) {
        matrix[0][l][0][0].value = SPACESCORE;
        matrix[0][l][0][0].numSegmentos = 1;
        matrix[0][l][0][0].before = &matrix[0][0][0][0];
        matrix[0][l][0][0].k = 0;
        matrix[0][l][0][0].l = l;

        for (j=1;j<C[l].tam;j++) {
            matrix[0][l][0][j].value = matrix[0][l][0][j-1].value + SPACESCORE;
            matrix[0][l][0][j].numSegmentos = 1;
            matrix[0][l][0][j].before = &matrix[0][l][0][j-1];
            matrix[0][l][0][j].k = 0;
            matrix[0][l][0][j].l = l;
        }
    }

    /*condição base 3*/
    for (k=1;k<u;k++) {
        matrix[k][0][0][0].value = SPACESCORE;
        matrix[k][0][0][0].numSegmentos = 1;
        matrix[k][0][0][0].before = &matrix[0][0][0][0];
        matrix[k][0][0][0].k = k;
        matrix[k][0][0][0].l = 0;

        for (i=1;i<B[k].tam;i++) {
            matrix[k][0][i][0].value = matrix[k][0][i-1][0].value + SPACESCORE;
            matrix[k][0][i][0].numSegmentos = 1;
            matrix[k][0][i][0].before = &matrix[k][0][i-1][0];
            matrix[k][0][i][0].k = k;
            matrix[k][0][i][0].l = 0;
        }
    }
    if (INFO) printf("Matrix initialized successfully!\n");
}

/*calcula a solução*/
void calculate() {
    int i, j, k, l, klinha, llinha, lastColumn;
    struct tempcell temp;

    if (INFO) printf("Calculating matrix...\n");
    for (k=1;k<u;k++) {
        for (l=1;l<v;l++) {
            for (i=0;i<B[k].tam;i++) {
                for (j=0;j<C[l].tam;j++) {
                    temp.value = INT_MIN;

                    /*equação 1*/
                    if (i>0 && j>0) {
                        /*linha 1*/
                        if (tempmatrix[i-1][j-1].value + score(bki(k,i), clj(l, j)) > temp.value) {
                            temp.value = tempmatrix[i-1][j-1].value + score(bki(k,i), clj(l, j));
                            temp.before = tempmatrix[i-1][j-1].before;
                            temp.numSegmentos = tempmatrix[i-1][j-1].numSegmentos;
                        }

                        /*linha 2*/
                        if (tempmatrix[i][j-1].value + score('-', clj(l, j)) > temp.value) {
                            temp.value = tempmatrix[i][j-1].value + score('-', clj(l, j));
                            temp.before = tempmatrix[i][j-1].before;
                            temp.numSegmentos = tempmatrix[i][j-1].numSegmentos;
                        }

                        /*linha 3*/
                        if (tempmatrix[i-1][j].value + score(bki(k,i), '-') > temp.value) {
                            temp.value = tempmatrix[i-1][j].value + score(bki(k,i), '-');
                            temp.before = tempmatrix[i-1][j].before;
                            temp.numSegmentos = tempmatrix[i-1][j].numSegmentos;
                        }
                    }
                    /*equação 2*/
                    else if (i==0 && j==0) {
                        /*linha 1*/
                        for (klinha=0; klinha<B[k].number; klinha++) {
                            if (precedek(klinha, k)) {
                                for (llinha=0; llinha<C[l].number; llinha++) {
                                    if (precedel(llinha, l)) {
                                        if (matrix[klinha][llinha][lastLinha(klinha)][lastColuna(llinha)].value + score(bki(k, 0), clj(l, 0)) > temp.value
                                            ||
                                            (matrix[klinha][llinha][lastLinha(klinha)][lastColuna(llinha)].value + score(bki(k, 0), clj(l, 0)) == temp.value && matrix[klinha][llinha][lastLinha(klinha)][lastColuna(llinha)].numSegmentos < temp.numSegmentos-2)) {
                                            temp.value = matrix[klinha][llinha][lastLinha(klinha)][lastColuna(llinha)].value + score(bki(k, 0), clj(l, 0));
                                            temp.before = &matrix[klinha][llinha][lastLinha(klinha)]
[lastColuna(llinha)];
                                            temp.numSegmentos = matrix[klinha][llinha][lastLinha(klinha)][lastColuna(llinha)].numSegmentos + 2;
                                        }
                                    }
                                }
                            }
                        }

                        /*linha 2*/
                        for (llinha=0; llinha<C[l].number; llinha++) {
                            if (precedel(llinha, l)) {
                                /*relacionado a otimização...*/
                                lastColumn = B[k].tam > 1 ? 0 : lastColuna(llinha);

                                if (matrix[k][llinha][0][lastColumn].value + score('-', clj(l, 0)) > temp.value ||
                                   (matrix[k][llinha][0][lastColumn].value + score('-', clj(l, 0)) == temp.value && matrix[k][llinha][0][lastColumn].numSegmentos < temp.numSegmentos-1)) {
                                    temp.value = matrix[k][llinha][0][lastColumn].value + score('-', clj(l, 0));
                                    temp.before = &matrix[k][llinha][0][lastColumn];
                                    temp.numSegmentos = matrix[k][llinha][0][lastColumn].numSegmentos + 1;
                                }
                            }
                        }

                        /*linha 3*/
                        for (klinha=0; klinha<B[k].number; klinha++) {
                            if (precedek(klinha, k)) {
                                if (matrix[klinha][l][lastLinha(klinha)][0].value + score(bki(k, 0), '-') > temp.value ||
                                    (matrix[klinha][l][lastLinha(klinha)][0].value + score(bki(k, 0), '-') == temp.value && matrix[klinha][l][lastLinha(klinha)][0].numSegmentos < temp.numSegmentos-1)) {
                                    temp.value = matrix[klinha][l][lastLinha(klinha)][0].value + score(bki(k, 0), '-');
                                    temp.before = &matrix[klinha][l][lastLinha(klinha)][0];
                                    temp.numSegmentos = matrix[klinha][l][lastLinha(klinha)][0].numSegmentos + 1;
                                }
                            }
                        }
                    }
                    /*equação 3*/
                    else if (i==0 && j>0) {
                        /*linha 1*/
                        if (tempmatrix[0][j-1].value + score('-', clj(l, j)) > temp.value) {
                            temp.value = tempmatrix[0][j-1].value + score('-', clj(l, j));
                            temp.before = tempmatrix[0][j-1].before;
                            temp.numSegmentos = tempmatrix[0][j-1].numSegmentos;
                        }

                        /*linhas 2 e 3*/
                        for (klinha=0; klinha<B[k].number; klinha++) {
                            if (precedek(klinha, k)) {
                                /*linha 2*/
                                if (matrix[klinha][l][lastLinha(klinha)][j-1].value + score(bki(k, 0), clj(l, j)) > temp.value ||
    (matrix[klinha][l][lastLinha(klinha)][j-1].value + score(bki(k, 0), clj(l, j)) == temp.value && matrix[klinha][l][lastLinha(klinha)][j-1].numSegmentos < temp.numSegmentos-1)) {
                                    temp.value = matrix[klinha][l][lastLinha(klinha)][j-1].value + score(bki(k, 0), clj(l, j));
                                    temp.before = &matrix[klinha][l][lastLinha(klinha)][j-1];
                                    temp.numSegmentos = matrix[klinha][l][lastLinha(klinha)][j-1].numSegmentos + 1;
                                }

                                /*linha 3*/
                                if (matrix[klinha][l][lastLinha(klinha)][j].value + score(bki(k, 0), '-') > temp.value || (matrix[klinha][l][lastLinha(klinha)][j].value + score(bki(k, 0), '-') == temp.value && matrix[klinha][l][lastLinha(klinha)][j].numSegmentos < temp.numSegmentos-1)) {
                                    temp.value = matrix[klinha][l][lastLinha(klinha)][j].value + score(bki(k, 0), '-');
                                    temp.before = &matrix[klinha][l][lastLinha(klinha)][j];
                                    temp.numSegmentos = matrix[klinha][l][lastLinha(klinha)][j].numSegmentos + 1;
                                }
                            }
                        }
                    }
                    /*equação 4*/
                    else if (i>0 && j==0) {
                        /*linha 1*/
                        if (tempmatrix[i-1][0].value + score(bki(k,i), '-') > temp.value) {
                            temp.value = tempmatrix[i-1][0].value + score(bki(k,i), '-');
                            temp.before = tempmatrix[i-1][0].before;
                            temp.numSegmentos = tempmatrix[i-1][0].numSegmentos;
                        }

                        /*linhas 2 e 3*/
                        for (llinha=0; llinha<C[l].number; llinha++) {
                            if (precedel(llinha, l)) {
                                /*linha 2*/
                                if (matrix[k][llinha][i-1][0].value + score(bki(k, i), clj(l, 0)) > temp.value || (matrix[k][llinha][i-1][0].value + score(bki(k, i), clj(l, 0)) == temp.value && matrix[k][llinha][i-1][0].numSegmentos < temp.numSegmentos-1)) {
                                    temp.value = matrix[k][llinha][i-1][0].value + score(bki(k, i), clj(l, 0));
                                    temp.before = &matrix[k][llinha][i-1][0];
                                    temp.numSegmentos = matrix[k][llinha][i-1][0].numSegmentos + 1;
                                }

                                /*linha 3*/
                                if (matrix[k][llinha][i][i==lastLinha(k) ? lastColuna(llinha) : 0].value + score('-', clj(l, 0)) > temp.value || (matrix[k][llinha][i][i==lastLinha(k) ? lastColuna(llinha) : 0].value + score('-', clj(l, 0)) == temp.value && matrix[k][llinha][i][i==lastLinha(k) ? lastColuna(llinha) : 0].numSegmentos < temp.numSegmentos-1)) {
                                    temp.value = matrix[k][llinha][i][i==lastLinha(k) ? lastColuna(llinha) : 0].value + score('-', clj(l, 0));
                                    temp.before = &matrix[k][llinha][i][i==lastLinha(k) ? lastColuna(llinha) : 0];
                                    temp.numSegmentos = matrix[k][llinha][i][i==lastLinha(k) ? lastColuna(llinha) : 0].numSegmentos + 1;
                                }
                            }
                        }
                    }
                    else {
                        printf("Fatal error! The program tried to calculate a position in matrix which is undefined. Aborting...");
                        exit(EXIT_FAILURE);
                    }
                    tempmatrix[i][j] = temp;
                }
            }

            /*com a tempmatrix calculada, é copiada sua última linha e sua última coluna para a matrix*/
            /*copia última coluna*/
            for (i=0;i<B[k].tam-1;i++) {
                matrix[k][l][i][0].value = tempmatrix[i][C[l].tam-1].value;
                matrix[k][l][i][0].before = tempmatrix[i][C[l].tam-1].before;
                matrix[k][l][i][0].numSegmentos = tempmatrix[i][C[l].tam-1].numSegmentos;
                matrix[k][l][i][0].k=k;
                matrix[k][l][i][0].l=l;
            }

            /*copia última linha*/
            for (j=0;j<C[l].tam;j++) {
                matrix[k][l][B[k].tam-1][j].value = tempmatrix[B[k].tam-1][j].value;
                matrix[k][l][B[k].tam-1][j].before = tempmatrix[B[k].tam-1][j].before;
                matrix[k][l][B[k].tam-1][j].numSegmentos = tempmatrix[B[k].tam-1][j].numSegmentos;
                matrix[k][l][B[k].tam-1][j].k=k;
                matrix[k][l][B[k].tam-1][j].l=l;
            }
        }
    }
    if (INFO) printf("Matrix calculated successfully!\n");
}

/*determina Gamma_B e Gamma_C baseado na matriz calculada*/
void getSolution() {
    int k, l;
    //struct cell *max = &matrix[0][0][0][0]; //with this line, empty solutions are valid.
    struct cell *max = &matrix[1][1][lastLinha(1)][lastColuna(1)]; //with this line, empty solutions are not valid.

    if (INFO) printf("Looking for the maximum value in matrix...\n");
    for (k=1;k<u;k++) {
        for (l=1;l<v;l++) {
            if(matrix[k][l][lastLinha(k)][lastColuna(l)].value > max->value || (matrix[k][l][lastLinha(k)][lastColuna(l)].value == max->value && matrix[k][l][lastLinha(k)][lastColuna(l)].numSegmentos < max->numSegmentos)) {
                max = &matrix[k][l][lastLinha(k)][lastColuna(l)];
            }
        }
    }
    if (INFO) printf("Maximum value in matrix found: %d.\n", max->value);

    if (COMPARE) {
        FILE *saida = fopen("rsolution", "w");
        fprintf(saida, "%d", max->value);
        fclose(saida);
    }


    if (INFO) printf("Building Gamma_B and Gamma_C...\n");
    gammaB = (int*)malloc(u * sizeof(int));
    eval(gammaB);
    tamGammaB = 0;
    gammaC = (int*)malloc(v * sizeof(int));
    eval(gammaC);
    tamGammaC = 0;

    while (max->before != NULL) {
        if (max->k!=0 && (tamGammaB==0 || max->k!=gammaB[tamGammaB-1])) {
            gammaB[tamGammaB] = max->k;
            tamGammaB++;
        }

        if (max->l!=0 && (tamGammaC==0 || max->l!=gammaC[tamGammaC-1])) {
            gammaC[tamGammaC] = max->l;
            tamGammaC++;
        }

        if (BACKTRACKCELLS==1) {
            printf("Going through [%d, %d]\n", max->k, max->l);
        }

        max = max->before;
    }
    if (INFO) printf("Gamma_B and Gamma_C built successfully!\n");
}

/*Desaloca as duas matrizes*/
void DesalocaMatriz(){
    int i, j, k, l;
    int maxTamB; /*o maior tamanho de um segmento de B*/

    if (INFO) printf("Disallocating matrix...\n");

    maxTamB=0;
    for (k=1;k<u;k++) {
        if (B[k].tam > maxTamB)
            maxTamB = B[k].tam;
    }

    /*desaloca a matriz temporária*/
    for (i=0; i<maxTamB; i++) {
        if(tempmatrix[i] != NULL) free(tempmatrix[i]);
    }
    if(tempmatrix != NULL) free(tempmatrix);
    
    /*desaloca a matriz de alinhamento*/
    for (k=0; k<u; k++) {
        
        for (l=0;l<v;l++) {
            
            for (i=0;i<B[k].tam;i++) {

                if(matrix[k][l][i] != NULL) free(matrix[k][l][i]);
            }
            if(matrix[k][l] != NULL) free(matrix[k][l]);
        }
        
        if(matrix[k] != NULL) free(matrix[k]);
    }  
    if(matrix != NULL) free(matrix);

}


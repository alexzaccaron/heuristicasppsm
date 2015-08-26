#ifndef SOLUTION_H_
#define SOLUTION_H_

/*matriz de alinhamento*/
struct cell ****matrix;

/*matriz temporária para cálculo da matriz de alinhamento*/
struct tempcell **tempmatrix;

/*inicializa a matriz*/
void initMatrix();

/*calcula a solução*/
void calculate();

/*determina Gamma_B e Gamma_C baseado na matriz calculada*/
void getSolution();

/*Desaloca as duas matrizes*/
void DesalocaMatriz();

#endif

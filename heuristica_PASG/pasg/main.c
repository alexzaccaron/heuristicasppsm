#include <stdio.h>
#include <stdlib.h>
#include "inout.h"
#include "solution.h"
#include <string.h>

int main(int argc, char *argv[])
{
    int i;
    char c;
    FILE *status, *outBlocks;                  //Changed by Alex
    
    outBlocks = fopen("blocksToRemove", "w");  //Changed by Alex
    
    /*verifica os argumentos passados por parâmetro*/
    if (argc!=3) {
    	printf("Error. Usage: ./pasg <input_file_1> <input_file_2>\n");
    	exit(EXIT_FAILURE);
    }


    /*lê as sequências e seus segmentos*/
    readSeqs(argv[1], argv[2]);

    /*calcula a solução*/
    initMatrix();
    calculate();
    getSolution();
    DesalocaMatriz();

    /*imprime a solução na saída padrão*/
    if (INFO) printf("\n\n--------------------\nPrinting solution:\n");
    if (INFO) printf("Gamma_B = {");
    for (i=tamGammaB-1;i>=0;i--) {
        if (INFO) printf("%d: [%d..%d], ", B[gammaB[i]].number, B[gammaB[i]].first, B[gammaB[i]].last);
    }
    if (INFO) printf("}\nGamma_C = {");
    for (i=tamGammaC-1;i>=0;i--) {
        if (INFO) {printf("%d: [%d..%d], ", C[gammaC[i]].number, C[gammaC[i]].first, C[gammaC[i]].last); fprintf(outBlocks, "%d\n", C[gammaC[i]].number);} //Changed by Alex
    }
    if (INFO) printf("}\n");

    if (SYSTEMINFO) {
        printf("\n\n--------------------\nPrinting process info:\n");
        status = fopen("/proc/self/status", "r");

        do {
            c = fgetc (status);

            if (c!=EOF)
	            putc(c, stdout);
        }while (c != EOF);
        fclose(status);
        printf("[END OF STATUS FILE]\n");
    }

    if (RELATORIO) {
        FILE *saida = fopen("relat", "w");

        fprintf(saida, "%d\n", n-1);
        fprintf(saida, "%d\n", m-1);
        fprintf(saida, "%d\n", u-1);
        fprintf(saida, "%d\n", v-1);
        fprintf(saida, "%d\n", tamGammaB);
        fprintf(saida, "%d\n", tamGammaC);
        for (i=tamGammaB-1;i>=0;i--)
            fprintf(saida, "%d %d ", B[gammaB[i]].first, B[gammaB[i]].last);
        fprintf(saida, "\n");
        for (i=tamGammaC-1;i>=0;i--)
            fprintf(saida, "%d %d ", C[gammaC[i]].first, C[gammaC[i]].last);
        fprintf(saida, "\n");

        fclose(saida);
    }

    if (AVALIACAO) {
        char dirName[1000];
        strncpy(dirName, argv[1], strlen(argv[1])-2);

        char temp[1000];
        strcpy(temp, dirName);
        strcat(temp, "gammab");

        FILE *saida = fopen(temp, "w");
        for (i=tamGammaB-1;i>=0;i--)
            fprintf(saida, "%d..%d,", B[gammaB[i]].first, B[gammaB[i]].last);
        fclose(saida);

        strcpy(temp, dirName);
        strcat(temp, "gammac");
        saida = fopen(temp, "w");
        for (i=tamGammaC-1;i>=0;i--)
            fprintf(saida, "%d..%d,", C[gammaC[i]].first, C[gammaC[i]].last);
        fclose(saida);
    }

    return 0;
}

/*
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

/*estrutura para auxiliar na descoberta do maior caminho*/
typedef struct Vector{ 
  int in;
  int weight;
}Vector;

//Vector vec[8];

/*funcao que aloca a matriz de adjacencia, devolvendo um ponteiro para ela*/
int **AlocaMatriz(int rows, int cols){
  int i;
  int **P;

  P = (int**)malloc(rows * sizeof(int*));
  if(P == NULL){
     printf("\nErro de alocacao\n");
     exit(1);
  }

  for (i=0;i<rows;i++){
     P[i] = (int*)malloc(cols * sizeof(int));
     if(P[i] == NULL){
        printf("\nErro de alocacao\n");
	exit(1);
     }	
  }
  
return P;
}

/*procedimeto que desaloca a matriz de adjacencia*/
void DesalocaMatriz(int *Mat[], int rows, int cols){
  int i;
  
  for(i=0; i<rows; i++){
     if(Mat[i] != NULL) free(Mat[i]);
  }
 
  if(Mat != NULL) free(Mat);
}
/*procedure to read the graph and the vertex weights*/
void ReadData(int *G[], int weights[], int n){
  int i, j;

  for (i = 0; i < n; ++i)
  {
    for (j = 0; j < n; ++j)
    {
      scanf("%d", &G[i][j]);
    }
  }

  for (i = 0; i < n; ++i)
  {
    scanf("%d", &weights[i]);
  }
}



void LongestPath(int *G[], int node, Vector vec[], int n){

  //printf("Determining LongestPath node %d\n", node+1);
  int maxIn, nodeIn, i;
  maxIn = 0;
  nodeIn = -1;
  //max(weight, weight+path)
  for (i = 0; i < n; ++i)
  {
    if (G[i][node] > 0 && vec[i].weight >= maxIn)
    {
      maxIn = vec[i].weight;
      nodeIn = i;
    }
  }

/*save longest path ending at node*/
  vec[node].weight = vec[node].weight + maxIn;
  vec[node].in = nodeIn;
}

void DFS(int i, int n, int *G[], int order[], int *index, int visited[])
{
    int j;
    visited[i]=1;
    for(j=0;j<n;j++)
        if(!visited[j]&&G[i][j]>0)
            DFS(j, n, G, order, &*index, visited);
            order[*index] = i;
            *index = *index-1;
            //printf("%d ",i+1);
              
}

int PathScore(Vector vec[], int n){
  int j, i, max;

  max = vec[0].weight;
  i = 0;
  for (j = 1; j < n; ++j)
  {
    if (vec[j].weight >= max)
    {
      max = vec[j].weight;
      i = j;
    }
  }

  return i;
}

void PrintPath(Vector vec[], int i, int *processed){
  int nodos[1000], count, j;

  count = 0;

  nodos[count] = i+1;
  count++;
  //printf("%d", i+1);
  vec[i].weight = -99999999;
  *processed = *processed + 1;
  while(vec[i].in != -1){
    i = vec[i].in;
    nodos[count] = i+1;
    count++;
    //printf(" %d", i+1);
    vec[i].weight = -99999999;
    *processed = *processed + 1;
  }

  for (j = count-1; j >= 0; --j)
  {
    printf("%d ", nodos[j]);
  }
  printf(" 0\n");
}


int main(int argc, char *argv[]){  
  int **G, *Tstart, *Tend, *Sstart, *Send;
  int i, j, n, processed;
  int *order, *weights;
  int *visited;
  Vector *vec;


  scanf("%d", &n); //number of vertices (matrix order)

  j = n-1;         //index for the order vector
  processed = 0;

  /*memory allocation*/
  order = malloc(n*sizeof(int));
  visited = malloc(n*sizeof(int));
  weights = malloc(n*sizeof(int));
  vec = malloc(n*sizeof(Vector));
  G = AlocaMatriz(n, n);

  /*visited vector set to zero*/
  for (i = 0; i < n; ++i)
  {
    visited[i] = 0;
    vec[i].weight = 0;
  }

  /*reading adjacency matrix*/
  ReadData(G, weights, n);


  for(i=0; i<n; i++){  
    if(!visited[i]) DFS(i, n, G, order, &j, visited);
  }


  do{
    for (i = 0; i < n; ++i)
    {
      if (vec[i].weight > -99999999)
      {
        vec[i].weight = weights[i];
      }
    }

    for (i = 0; i < n; ++i)
    {
      if(vec[order[i]].weight > -9999999){ LongestPath(G, order[i], vec, n); }
    }

    i = PathScore(vec, n);
    
    //printf("\nPath score: %d", vec[i].weight);
    PrintPath(vec, i, &processed);


  }while(processed<n);

  //Rand(2000);

/*memory free*/ 
  DesalocaMatriz(G, n, n);
  if(order != NULL) free(order);
  if(visited != NULL) free (visited);
  if(vec != NULL) free(vec);
  if(weights != NULL) free(weights);

  return 0;
}

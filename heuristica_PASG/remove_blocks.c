/*
Given a file A with intervals (blocks) in the form x..y,w..z,.., and a file B with numbers, one per line:
the blocks which number is in B, will be removed from A.
The first block in A is 1, the second is 2, and so on...*/

#include<stdio.h>
#include<stdlib.h>
#include<string.h>

//while there is no more efficient way, set the maximum number of blocks to remove
#define MAX_TO_REMOVE 9999  

void Usage(){
	printf("Usage: remove_blocks fileA fileB\nfileA: blocks\nfileB: numbers of the blocks to be removed\n");
	exit(1);
}
int main(int argc, char *argv[]){
	
	if(argc < 3 ){
		Usage();
	}
	
	char block[24];
	int countBlocks, j, i, toRemove[MAX_TO_REMOVE], removed, blockNumber;
	FILE *blocks, *blocksToRemove;
	
	//arg1: file with the blocks (A)
	//arg2: file with the number of the blocks to be removed (B)
	blocks = fopen(argv[1], "r");	
	blocksToRemove = fopen(argv[2], "r");

	j = 1;
	/*Reading the blocks to remove*/
	while(fscanf(blocksToRemove, "%d", &blockNumber) > 0) {
        toRemove[j] = blockNumber;
        j++;
    }
    
   /*How many blocks to remove*/ 
	toRemove[0] = j-1;
	
	countBlocks = 0;
	removed = 0;
	
	do{
		i = 0;
		//reading a block
		do{
			block[i] = fgetc(blocks);
			i++;
		}while(block[i-1] != ',' && block[i-1] != EOF);
		block[i] = '\0';
		countBlocks++;

		//if the number match, it will be removed
		for(j=1; j<=toRemove[0]; j++){
			if(countBlocks == toRemove[j]) removed = 1;
		}
		
		//if not removed, then print. It prints only blocks that are are not removed.
		if(!removed && block[i-1] != EOF) printf("%s", block);
	
		removed = 0;
	
	}while(block[i-1] != EOF);
	
	fclose(blocks);
	fclose(blocksToRemove);
return 0;
}

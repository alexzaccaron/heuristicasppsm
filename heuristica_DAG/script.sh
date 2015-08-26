#second MSPP heuristic based on the DAG approach
#requirements: sets S and T, and the needleman-Wunsch algorithm from EMBOSS Package
 
#first argument ($1): set S
#second argument ($2): set T

#first the graph G is created 
/home/said/Copy/heuristica2/theory $1 $2 > graph.in

#then longest paths are found until there is no longer vertices in G
/home/said/Copy/heuristica2/heuristica < graph.in > paths.out

#from the paths, the pairs (si, ti) are defiend 
/home/said/Copy/heuristica2/createPartition $1 $2 paths.out > pairs.out

#lastly, the pairs are alined
/home/said/Copy/heuristica2/alignpairs $1 $2 < pairs.out 

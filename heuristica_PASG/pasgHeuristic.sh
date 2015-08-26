#!/bin/bash

pair_number=1

#copy the files to preserve the original ones
cp $2 $2_temp
cp $2.predex $2_temp.predex


#while the size of read_temp.predex > 0, e. i. while there are segments in the file
while [ -s $2_temp.predex ]
do

#print pair number
echo "Pair "$pair_number

#call pasg
/home/said/Copy/heuristica/pasg/pasg $1 $2_temp
echo "--------------------"$'\n'

pair_number=$((pair_number + 1))


#remove the blocks in the pasg solution from the file
/home/said/Copy/heuristica/remove_blocks $2_temp.predex blocksToRemove > aux.predex
mv aux.predex $2_temp.predex

done

#delete temporary files
rm blocksToRemove $2_temp $2_temp.predex

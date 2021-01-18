#!/bin/bash



for FILE in *_{x,y}.nex; do
    echo $FILE
    sed 's/\(tree\)\([0-9]\)/\1_\2/' $FILE > ${FILE/nex/trees}
done

for FILE in *-0_x.trees; do
    echo $FILE
    bash ~/BEASTv1.10.4/bin/treeloganalyser $FILE
done

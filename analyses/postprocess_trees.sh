#!/bin/bash



for FILE in *_{x,y}.nex; do
    echo $FILE
    sed 's/\(tree\)\([0-9]\)/\1_\2/' $FILE > ${FILE/nex/trees}
done

for FILE in *-0_x.trees; do
    echo $FILE
    bash ~/BEASTv1.10.4/bin/treeloganalyser $FILE
done

for FOLDER in 20210113 20210128 20210129 20210130 20210131 20210132 20210133; do
    echo $FOLDER
    cd $FOLDER/output
    for FILE in *-0.nex; do
        # sed 's/\(tree\)\([0-9]\)/\1_\2/' $FILE > ${FILE/nex/trees}
        # bash ~/BEASTv1.10.4/bin/treeloganalyser -export ${FILE/nex/supp} -burnin 5000 ${FILE/nex/trees}
        sed -n 's/^.*\[\&W \(.*\) \].*$/\1/p' ${FILE/nex/supp} > ${FILE/nex/supps}
    done
    cd ../../
    echo ""
    echo ""
    echo ""
done

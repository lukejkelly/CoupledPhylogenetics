#!/bin/bash

BEAST_LIB="~/BEASTv1.10.4/lib"
java -Xms64m -Xmx1024m -Djava.library.path="$BEAST_LIB" -cp "$BEAST_LIB/beast.jar" dr.app.tools.TreeLogAnalyser $*

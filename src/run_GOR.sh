#!/bin/bash

set -u

# run GOR IV SS prediction from R package DECIPHER from http://www2.decipher.codes/
# --- R package install ---
# if (!requireNamespace("BiocManager", quietly=TRUE))
#    install.packages("BiocManager")
# BiocManager::install("DECIPHER")
# install.packages("stringr")

SEQUENCES_FN=`mktemp`
R_SCRIPT_FN=`mktemp`

# just keep sequences from input file
cut -f1 $1 > $SEQUENCES_FN

cat <<EOF > $R_SCRIPT_FN
library('DECIPHER', quietly = TRUE)
library("stringr", quietly = TRUE)
lines <- readLines('$SEQUENCES_FN')
for (line in lines) {
        aa <- AAStringSet(line)
        hec <- PredictHEC(aa)
        helix_aa <- str_count(hec, "H")
        n <- nchar(hec)
        helicity <- helix_aa / n
        print(sprintf("%s %.2f", hec, helicity))
    }
quit()
EOF

R --vanilla --slave 2>/dev/null < $R_SCRIPT_FN | cut -d '"' -f2

rm -f $SEQUENCES_FN $R_SCRIPT_FN

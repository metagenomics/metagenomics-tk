LOG=$1

grep -q IGNORE ${LOG} ; test $? -eq 1

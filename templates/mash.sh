ls -1 input*.txt > input_list.txt
mash sketch -o reference -l input_list.txt
mash dist -p !{task.cpus} reference.msh  reference.msh | cut -f 1,2,3 > distance.tsv
for i in $(ls input*.txt); do echo "$i\t$(readlink -f $i)";  done > mapping.tsv

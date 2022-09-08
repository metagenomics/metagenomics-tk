# convert simple gfa format to fastg
# the format has S lines: S\t<seq name>\tdp:i:<int>
# and L lines: L\t<seq 1>\t<-/+>\t<seq 2>\t<-/+>\t0M
#
# usage: python metaflye_gfa2fastg.py <infile> <outfile> 

import sys

infile = sys.argv[1]
outfile = sys.argv[2]

seqs_dict = {}
headers = {}

complements = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}
def rc(seq):
  rev = reversed(seq)
  return "".join([complements[i] for i in rev])


with open(infile) as f:
  for line in f:
    splt = line.strip().split()
    if splt[0] == 'S':
       name = splt[1]
       seq = splt[2]
       dp = splt[3]
       depth = float(dp.split(':')[2])
       seqs_dict[name] = (len(seq),depth,seq)
    if splt[0] == 'L':
      headers[splt[1]] = headers.get(splt[1],list())
      values = (splt[2],splt[3],splt[4])
      if values not in headers[splt[1]]:
        headers[splt[1]].append(values)
      values = ('+' if splt[4]=='-' else '-', splt[1],'+' if splt[2]=='-' else '-')
      headers[splt[3]] = headers.get(splt[3],list())
      if values not in headers[splt[3]]:
        headers[splt[3]].append(values)


with open(infile) as f, open(outfile,'w') as o:
  for line in f:
    splt = line.strip().split()
    if splt[0] == 'S':
      name = splt[1]
      fields = seqs_dict[name]
      if name not in headers:
        o.write('>'+name+"_len_"+str(fields[0])+"_cov_"+str(fields[1])+";\n"+fields[2]+'\n')
        o.write('>'+name+"_len_"+str(fields[0])+"_cov_"+str(fields[1])+"';\n"+rc(fields[2])+'\n')
        continue
      pos_neighbours = [x for x in headers[name] if x[0] == "+"]
      neg_neighbours = [x for x in headers[name] if x[0] == "-"]
      if len(pos_neighbours) > 0:
          o.write('>'+name+"_len_"+str(fields[0])+"_cov_"+str(fields[1])+":"+",".join([n[1]+"_len_"+str(seqs_dict[n[1]][0])+"_cov_"+str(seqs_dict[n[1]][1])+("'"if n[2]=="-" else "") for n in pos_neighbours])+ ";\n"+fields[2]+'\n')
      else:
          o.write('>'+name+"_len_"+str(fields[0])+"_cov_"+str(fields[1])+";")
      if len(neg_neighbours) > 0:
          o.write('>'+name+"_len_"+str(fields[0])+"_cov_"+str(fields[1])+"':"+",".join([n[1]+"_len_"+str(seqs_dict[n[1]][0])+"_cov_"+str(seqs_dict[n[1]][1])+("'"if n[2]=="-" else "") for n in neg_neighbours])+ ";\n"+rc(fields[2])+'\n')
      elif len(pos_neighbours) > 0:
        o.write('>'+name+"_len_"+str(fields[0])+"_cov_"+str(fields[1])+";")
        



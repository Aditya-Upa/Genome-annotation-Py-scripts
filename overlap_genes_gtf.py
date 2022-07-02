#syntax: python nameofprogram.py < nameofgenomicfile.gtf 

import re, sys

def overlap(range1, range2):
  return range1[0] == range2[0] and range1[2] >= range2[1] and \
         range1[1] <= range2[2]


def process_gff(gff_file):
  genes = {}

  # read the gff file and store gene name, contig, gene start, gene_end
  for line in open(gff_file, 'r'):
    if line == "" or line.startswith("#"):
      continue

    line = line.rstrip()
    fields = line.split("\t")
    # print(fields)
    seqid = fields[0]
    start = int(fields[3])
    end = int(fields[4])
    attributes = fields[8]

    matches = re.match('gene_id "([^"]+)"', attributes)
    if matches:
      geneid = matches.group(1)
      if geneid not in genes:
        genes[geneid] = [seqid, start, end]
      else:
        assert seqid == genes[geneid][0]
        if start < genes[geneid][1]:
          genes[geneid][1] = start
        if end > genes[geneid][2]:
          genes[geneid][2] = end

  # compare the locations and get the overlap
  list_overlap=[]
  gene_list = list(genes.keys())
  for i in range(len(genes)):
    for j in range(i+1, len(genes)):
      geneid_i = gene_list[i]
      geneid_j = gene_list[j]
      range_i = genes[geneid_i]
      range_j = genes[geneid_j]
      if overlap(range_i, range_j):
        list_overlap.append([geneid_i, geneid_j])
  # print(list_overlap)
  print(f'Gene count: {len(list_overlap)}')
  
  return genes, list_overlap


if __name__ == '__main__':
  gff_file = 'genes.gtf'
  genes, overlap_list = process_gff(gff_file)
  print(overlap_list)

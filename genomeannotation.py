#!usr/bin/python
from __future__ import print_function, division
import sys
from operator import itemgetter
import argparse
from tqdm import tqdm
import gzip
from BCBio import GFF

"""
    Annotate Bed: v0.2.0
    CAUTION: Works currently for only sorted bed files.

    Preferential order of assigning annotation:
        Promoter >> Overlapping >> Intergenic

    Defaults:
        > Promoter distance is 1kb upstream and downstream of TSS.
        > Gene id considered is "gene".

    Output:
        > The column with distance from TSS has the closest distance from TSS from any co-ordinate.
        > Value for upstream is negative and value for downstream is positive
"""

parser = argparse.ArgumentParser(prog="annotate bed", formatter_class=argparse.RawDescriptionHelpFormatter, description=__doc__)
parser.add_argument("-bed", required=True, help="Input bed file.")
parser.add_argument("-gff", required=True, help="Input gff file.")
parser.add_argument("-gid", "--gene-id", dest="gid", default="gene", type=str, help="")
parser.add_argument("-out", dest="out", type=str)
parser.add_argument("-prom", "--prom-distance", dest="prom", nargs='+', type=int, help=""""Define promoter distance around the TSS. 
                                                                                        Upstream distance and downstream distance respectively should be provided.
                                                                                        If given one value same distance will be considered for up and down streams.""")
args = parser.parse_args()

bedFile = args.bed
gffFile = args.gff
geneId = args.gid
if args.out is None:
    args.out = bedFile.split('.')[0] + '_annotation.tsv'
output_file = open(args.out, 'w')
if args.prom is None:
    args.prom = [1000, 500]
if len(args.prom) == 1:
    args.prom = args.prom*2
promUp = args.prom[0]
promDown = args.prom[1]


# Logic for checking overlap-
#   Sort the intervals based on starts.
#   Check if start of second interval is less than stop of first interval.
def checkOverlap(gStart, gStop, rStart, rStop):
    if gStart <= rStart and rStart <= gStop:
        return True
    elif rStart <= gStart and gStart <= rStop:
        return True
    else:
        return False

# Distance returned with respect to TSS-
#   Positive for downstream of TSS
#   Negative for upstream of TSS
def get_TSSdistance(geneStart, geneStop, geneStrand, bedStart, bedStop):
    if geneStrand is None:
        geneStrand = 1
    if geneStrand == 1:
        TSS = geneStart
        if abs(bedStart-TSS) <= abs(bedStop-TSS):
            dist = bedStart - TSS
        else:
            dist = bedStop - TSS
    else:
        TSS = geneStop
        if abs(bedStart-TSS) <= abs(bedStop-TSS):
            dist = TSS - bedStart
        else:
            dist = TSS - bedStop
    return dist

def get_overlapLength(gStart, gStop, rStart, rStop):
    repLength = rStop - rStart
    geneLength = gStart - gStop
    if gStart <= rStart:
        overlap_len = ((gStop - rStart) if gStop <= rStop else repLength) 
    else:
        overlap_len = ((rStop - gStart) if rStop <= gStop else geneLength)
    return round(overlap_len*100/repLength, 2)

def processGFF(gffFile):
    geneObj = {}
    if gffFile.endswith('gz'):
        gffhandle = gzip.open(gffFile, 'rt')
    else:
        gffhandle = open(gffFile, 'r')
    limit_info = dict(gff_type = ['gene', 'exon', 'intron'])
    for i, record in enumerate(GFF.parse(gffhandle, limit_info=limit_info)):
        geneObj[record.id] = []
        for f in record.features:
            geneObj[record.id].append(f)
    gffhandle.close()
    return geneObj

def get_Exons(gffFile):
    exons = {}
    if gffFile.endswith('gz'):
        gffhandle = gzip.open(gffFile, 'rt')
    else:
        gffhandle = open(gffFile, 'r')
    for line in gffhandle:
        if not line.startswith('#'):
            line = line.strip().split('\t')
            description = line[8].split(';')
            if line[2] == 'exon':
                for a in description:
                    if a.startswith('gene='):
                        gene = a
                        if gene in exons:
                            exons[gene].append([line[0], int(line[3]), int(line[4])])
                        else:
                            exons[gene] = [[line[0], int(line[3]), int(line[4])]]
    return exons

def getStrand(strand):
    if strand is None:
        return '+'
    elif strand < 0:
        return '-'
    else:
        return '+'

def get_geneName(gene, geneId):
    try:
        return '%s=%s' %(geneId ,','.join(gene.qualifiers[geneId]))
    except:
        key = gene.qualifiers[geneId].keys()[0]
        return '%s=%s' %(key ,','.join(gene.qualifiers[key]))


def get_closegenes(repStart, repStop, genes):
    closestGene = ['', float('inf')]
    overlapGenes = []
    for gene in genes:
        geneStart = gene.location.start
        geneStop = gene.location.end
        geneStrand = gene.location.strand
        if geneStrand is None:
            geneStrand = 1
        geneTSS = geneStart
        # Considering the region including the upstream promoter distance for overlap check
        regStart = geneStart - promUp
        regStop = geneStop
        if geneStrand == -1:
            geneTSS = geneStop
            regStart = geneStart
            regStop = geneStop + promUp
        # Getting the closest gene w.r.t TSS
        TSSdist = get_TSSdistance(geneStart, geneStop, geneStrand, repStart, repStop)
        # Break the loop-
        #   If the analysing geneStart is at a distance greater than allowed Upstream promoter, from the repeat stop.
        #   AND
        #   Is also at a greater distance than the closest gene. (OR) No gene overlaps are found.
        if (geneStart - repStop > promUp and (abs(TSSdist) > abs(closestGene[1]) or len(overlapGenes) > 0)):
            break
        overlap = checkOverlap(regStart, regStop, repStart, repStop)
        if overlap:
            overlapGenes.append([gene, TSSdist])
        elif abs(TSSdist) < closestGene[1]:
            closestGene = [gene, TSSdist]

    return {'overlap': overlapGenes, 'closest': closestGene}

def get_subfeatureAnno(genes, repStart, repStop):
    annotation_features = {'Intergenic': [], 'Exon': [], 'Intron': [], 'Genic': []}
    exon_overlap_percentage = 0.0
    for gene in genes:
        TSSdist = gene[1]
        gene = gene[0]
        # Intergenic is possible if its a case of upstream promoter
        if TSSdist < 0:
            if len(annotation_features['Intergenic']):
                if abs(TSSdist) < abs(annotation_features['Intergenic'][1]):
                    annotation_features['Intergenic'] = [gene, TSSdist]
            else:
                annotation_features['Intergenic'] = [gene, TSSdist]
        elif TSSdist >= 0:
            sub_features = gene.sub_features
            # Check for sub-annotation when sub features are available
            if len(sub_features) > 0:
                exon_check = 0
                # Check for Exonic overlaps
                for sub_feature in sub_features:
                    exonStart = sub_feature.location.start
                    exonStop = sub_feature.location.end
                    exonOverlap = checkOverlap(exonStart, exonStop, repStart, repStop)
                    if exonOverlap:
                        exon_check = 1
                        new_exon_overlap_percentage = get_overlapLength(exonStart, exonStop, repStart, repStop)
                        if new_exon_overlap_percentage > exon_overlap_percentage:
                            exon_overlap_percentage = new_exon_overlap_percentage
                        if len(annotation_features['Exon']):
                            if abs(TSSdist) < abs(annotation_features['Exon'][1]):
                                annotation_features['Exon'] = [gene, TSSdist]
                        else:
                            annotation_features['Exon'] = [gene, TSSdist]
                # If no exonic overlap is found store sub-annotation as intronic
                if not exon_check:
                    if len(annotation_features['Intron']):
                        if abs(TSSdist) < abs(annotation_features['Intron'][1]):
                            annotation_features['Intron'] = [gene, TSSdist]
                    else:
                        annotation_features['Intron'] = [gene, TSSdist]
            # If no subfeature information available. Store sub-annotation as Genic.
            else:
                if len(annotation_features['Genic']):
                    if abs(TSSdist) < abs(annotation_features['Genic'][1]):
                        annotation_features['Genic'] = [gene, TSSdist]
                else:
                    annotation_features['Genic'] = [gene, TSSdist]
    if len(annotation_features['Exon']):
        return annotation_features['Exon'] + ['Exon', exon_overlap_percentage]
    elif len(annotation_features['Intron']):
        return annotation_features['Intron'] + ['Intron', exon_overlap_percentage]
    elif len(annotation_features['Genic']):
        return annotation_features['Genic'] + ['Genic', exon_overlap_percentage]
    elif len(annotation_features['Intergenic']):
        return annotation_features['Intergenic'] + ['Intergenic', exon_overlap_percentage]

def annotate(bedFile, geneObj, exons):
    lines = 0
    with open(bedFile) as fh:
        for line in fh:
            lines += 1

    with open(bedFile) as fh:
        for line in tqdm(fh, total=lines):
            line = line.strip().split('\t')
            seqname = line[0]
            repStart = int(line[1])
            repStop = int(line[2])
            try:
                genes = geneObj[seqname]
                closeGenes = get_closegenes(repStart, repStop, genes)
                closestGene = closeGenes['closest']
                overlapGenes = closeGenes['overlap']
                if len(overlapGenes) > 0:
                    promoter_check = 0
                    annotations = {'Promoter': [], 'Genic': []}
                    for gene in overlapGenes:
                        TSSdist = gene[1]
                        gene = gene[0]
                        # Check if there is a promoter match
                        if -1*promUp <= TSSdist <= promDown:
                             promoter_check = 1
                             annotations['Promoter'].append([gene, TSSdist])
                             annotations['Genic'] = []
                        # If no promoter match is found yet append to Genic matches 
                        if promoter_check == 0:
                            annotations['Genic'].append([gene, TSSdist])
                    if promoter_check == 1:
                        promoter_anno = 'Promoter'
                        # print(promoter_anno)
                        final_annotation_gene = get_subfeatureAnno(annotations['Promoter'], repStart, repStop)
                    else:
                        promoter_anno = 'Non-Promoter'
                        final_annotation_gene = get_subfeatureAnno(annotations['Genic'], repStart, repStop)
                    if final_annotation_gene[2] == 'Genic':
                        exon_overlap_percentage = 0.0
                        exon_check = 0
                        gene = get_geneName(final_annotation_gene[0], geneId)
                        try:
                            for exon in exons[gene]:
                                gchrom = exon[0]
                                gStart = exon[1]
                                gStop = exon[2]
                                description = exon[2]
                                if gchrom == seqname and checkOverlap(gStart, gStop, repStart, repStop):
                                    exon_check = 1
                                    new_exon_overlap_percentage = get_overlapLength(gStart, gStop, repStart, repStop)
                                    if new_exon_overlap_percentage > exon_overlap_percentage:
                                        exon_overlap_percentage = new_exon_overlap_percentage
                            if exon_check == 1:
                                final_annotation_gene[2] = 'Exon'
                                final_annotation_gene[3] = exon_overlap_percentage
                            else:
                                final_annotation_gene[2] = 'Intron' 
                        except KeyError:
                            pass
                    annotation_fields = [get_geneName(final_annotation_gene[0], geneId), 
                                        str(final_annotation_gene[0].location.start), 
                                        str(final_annotation_gene[0].location.end), 
                                        getStrand(final_annotation_gene[0].location.strand), 
                                        final_annotation_gene[2], promoter_anno, str(final_annotation_gene[1]), str(final_annotation_gene[3])]
                    print('\t'.join(line), '\t'.join(annotation_fields), sep='\t', file=output_file)
                else:
                    annotation = 'Intergenic'
                    promoter_anno = 'Non-Promoter'
                    annotation_fields = [get_geneName(closestGene[0], geneId), 
                                        str(closestGene[0].location.start), 
                                        str(closestGene[0].location.end), 
                                        getStrand(closestGene[0].location.strand), 
                                        annotation, promoter_anno, str(closestGene[1]), '0.0']
                    print('\t'.join(line), '\t'.join(annotation_fields), sep='\t', file=output_file)
            except KeyError:
                print('\t'.join(line), '\t'.join(['NA']*7), sep='\t', file=output_file)
print('Processing "%s" feature file.' %(gffFile), end='\n\n')
geneObj = processGFF(gffFile)
allExons = get_Exons(gffFile)
print('Annotating regions in "%s" bed file' %(bedFile))
annotate(bedFile, geneObj, allExons)
output_file.close()

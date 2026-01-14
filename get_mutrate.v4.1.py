#!/home/users/astar/gis/zhangyu2/.conda/envs/RNA_Struc/bin/python
__doc__="""
    the script is for calculate mutrate and coverage of each line (nucleotide)
    iteratively instead of make a directory

    it's similar but not same with previous get_mutrate.v1.*.py. the previous one used dictionary but
    here the iterator of each line was used. therefore, it would not able to 1) select specified gene,
    2) not all nucleotide would be reserved in the final outpufile,

    changes to release v2.2
    the previous one is to calculated gene_expr normlized coveragee
        normalized_cov = reads_count(nt)/read_count(gene)
    in the v2.2 we add a choice to normlized by total reads counts of the sample
        normalized_cov = reads_count(nt)/total_reads_count(gene, by million)

    changes in release 3.4
    1) don't output each gene mutrate as space limit
    
    changes in release 4.0, many changes
    1) for new bam-readcounts output, there is no = counts for identical nts, so calculate identical counts from base details;
    2) correct an error, the coverage of each nt didn't count insertion ( represented with +), 
       if base_count[+] is not zero, then coverage = coverage + base_count(+)
    3) if the ins or del is too long (>= 3nt by default), then the read was not included into consideration;
    
    change in release 4.1
    1) add an arguement --ref_symbol, the ref symbol is the symbol to represent the ref nucleotide in mutrate.txt file

    change in release 4.1.1
    1) make --gene_readcounts optional; when not provided, normalized coverage calculation is disabled

    Usage: python "$SCRIPT" -i "$f" -e "$fcfile" -o "$out" -w 0 -gc 100
    Usage: python "$SCRIPT" -i "$f" -o "$out" -w 0
    """
__version__="v4.1.1"
__author__="noahzy"
__last_modify__="04-Nov-2025"

import gzip
import io
from sys import stdout
import argparse
import os
import shutil
from copy import deepcopy
from multiprocessing import Pool
from decimal import Decimal

def mk_idx(reffile):
    seqs={}
    genes=[]
    with open(reffile) as ref:
        for l in ref:
            if l[0] == ">":
                genes.append(l[1:].strip("\n"))
                seqs[genes[-1]]=""
            else:
                seqs[genes[-1]]+=l.strip('\n').upper()
    ref_idx={}
    for g in genes:
        ref_idx[g]=[0 for n in range(len(seqs[g]))]
    return ref_idx

def iterator_mutrate(inputfile, fcfile, out_diretory, out_prefix, win_threshold, cov_threshold=-1,
                     gene_output_cov_threshold=10, max_ins_len=5, max_del_len=5, equal_symbol_ref=False):
    ## make output file
    if out_diretory[-1] == "/":
        pass
    else:
        out_diretory = out_diretory+"/"
    if out_prefix != "-":
        if os.path.exists(out_prefix):
            os.remove(out_prefix)
        output = gzip.open(out_prefix, 'wt')
    else:
        output = stdout
    ## make a dictionary for read count in each gene_expr
    ## if fcfile is not provided, disable normalized coverage calculation
    if fcfile is not None:
        gexpr, ref_idx, total_reads_count = make_gene_cov(fcfile)
        total_reads_count = float(total_reads_count)/1000000## million reads
    else:
        gexpr, ref_idx, total_reads_count = {}, {}, None
    ## make variances for the function
    #cov_map = deepcopy(ref_idx)
    #mutrate_map = deepcopy(ref_idx)
    ## input the bam-readcount file (gzipped or not)
    if inputfile.split(".")[-1] in ["gz","gzip"]:
        #myopen = gzip.open
        input = io.TextIOWrapper(io.BufferedReader(gzip.open(inputfile,'rb')))
    else:
        input = open(inputfile)
    ## iterate the bam-readcount file
    ## set up win_cov, for a x nt window for estimate the position

    #    process_line(l,ref_idx, gexpr, output,cov_threshold=1)
    ## make a process_line funtion
    #def process_line(l,ref_idx, gexpr, output,cov_threshold=1):
    def process_line(l, equal_symbol_ref=equal_symbol_ref):
        ## set nonlocal variables
        nonlocal gexpr
        nonlocal total_reads_count ## for the normalization by total reads of the sample (million reads)
        #nonlocal output
        nonlocal cov_threshold
        ## deconstruct each line
        i = l.strip('\n').split()
        gene=i[0]
        # Strip the '>' prefix if present
        if gene.startswith('>'):
            gene = gene[1:]
        pos=int(i[1])-1
        refnt=i[2]
        coverage=int(i[3])
        adj_coverage=int(coverage)
        #print (coverage, total_reads_count, cov_threshold)

        ## calculate normalized coverage
        ## if fcfile is not provided, set normalized_cov and g_readcount to NA
        if total_reads_count is not None:
            g_len, g_readcount = gexpr[gene]
            ## normalized_cov = round(float(coverage)/g_readcount,8)
            normalized_cov = round(Decimal(float(coverage)/total_reads_count),8)
        else:
            g_readcount = "NA"
            normalized_cov = "NA"
        ## if there is any filter for each line, just add here
        # if cov_threshold == -1:
        #     pass
        # elif coverage < cov_threshold:
        #     #continue
        #     return
        ## make detail column to show propotion of each nt(ATCGN)
        ## the format is A:numa;T:numt...
        detail = []
        base_count_df = {}
        ins, delete, large_ins, large_del = 0, 0, 0, 0
        for item in i[4:]:
            ## sum number of  insertions or deleltion together
            if item[0] == "+":
                ins_len = len(item.split(":")[0])-1
                if ins_len <= max_ins_len:
                    ins += int(item.split(":")[1])
                else:
                    large_ins += int(item.split(":")[1])
            elif item[0] == "-":
                del_len = len(item.split(":")[0])-1
                if del_len <= max_del_len:
                    delete += int(item.split(":")[1])
                else:
                    large_del += int(item.split(":")[1])
            else:
                nt = item.split(":")[0]
                num = item.split(":")[1]
                detail.append(nt+":"+num)
                if nt not in base_count_df.keys():
                    base_count_df[nt] = int(num)
                else:
                    base_count_df[nt] += int(num)
        if ins != 0:
            nt = "+"
            num = str(ins)
            detail.append(nt+":"+num)
            if nt not in base_count_df.keys():
                base_count_df[nt] = int(num)
            else:
                base_count_df[nt] += int(num)
            coverage = coverage + ins ## ins is not count in the original coverage
            adj_coverage = adj_coverage + ins ## ins is not count in the original coverage

        if delete != 0:
            nt = "-"
            num = str(delete)
            detail.append(nt+":"+num)
            if nt not in base_count_df.keys():
                base_count_df[nt] = int(num)
            else:
                base_count_df[nt] += int(num)
        if large_ins != 0:
            nt = "+++"
            num = str(large_ins)
            detail.append(nt+":"+num)
            if nt not in base_count_df.keys():
                base_count_df[nt] = int(num)
            else:
                base_count_df[nt] += int(num)
            coverage = coverage + large_ins ## large ins is not count in the original coverage
        if large_del != 0:
            nt = "---"
            num = str(large_del)
            detail.append(nt+":"+num)
            if nt not in base_count_df.keys():
                base_count_df[nt] = int(num)
            else:
                base_count_df[nt] += int(num)
            adj_coverage = adj_coverage - large_del ## remove large del from coverage
        ## the mutant read counts without modification
        ## insertions and deletions are also counted
        #ident=int(i[4].split(":")[1])
        if equal_symbol_ref :
            ident = base_count_df["="]
        else:
            ident = base_count_df[refnt.upper()]
        #mutant=coverage-ident  ## all mutants are counted,
        mutant=adj_coverage-ident  ## remove large_ins/del

        if coverage == 0:
            mutrate = "NA"
        else:
            mutrate = round(Decimal(float(mutant)/coverage),8)
        ## put the coverage and mutrate into pseudo_matrix (map)
        #cov_map[gene][pos]=coverage
        #mutrate_map[gene][pos]=mutrate
        ## output bed-like format
        # it's a simple way to output all lines without filtering by gene length  10-Sep-2024, temp work
        # output.write("\t".join([gene,i[1],str(pos+2),gene+"."+i[1],str(mutrate),".",
        #                         str(coverage), str(adj_coverage), str(mutant),str(normalized_cov),
        #                         str(g_readcount),refnt,";".join(detail)])+"\n")
        return gene,"\t".join([gene,i[1],str(pos+2), gene+"."+i[1],str(mutrate),".",
                    str(coverage),str(adj_coverage),str(mutant),str(normalized_cov),
                    str(g_readcount),refnt,";".join(detail)]), coverage, adj_coverage, normalized_cov, pos


    ## try multiprocessing failed,
    ## single threshod were used
    ## output to files for each gene which reach the threshold
    def parse_single_gene(win_len=100, normalized_cov_threshold=10, len_prop_threshold=None):
        ## notes for thresholds:
        ## normalized_cov_threshold is x reads / total reads, x: cov at position, total reads is measured by million
        ## win_threshold is the threshold for positions number which > normalized_cov_threshold
        nonlocal input

        ## simple run
        for l in input:
            gene, line, cov, adjcov, ncov, pos = process_line(l)
            ## if ncov is "NA", only filter by cov_threshold
            if ncov == "NA":
                if cov >= cov_threshold:
                    output.write(line+'\n')
            else:
                if (cov >= cov_threshold) & (ncov >= normalized_cov_threshold):
                    output.write(line+'\n')
            #print ([gene, line, cov, adjcov, ncov, pos])

        ## complicated run
        """output_temp = []
        gcov_meet_condition = []
        previous_gene = ""
        gene_lens = 0
        for l in input:
            gene, line, cov, adjcov, ncov, pos = process_line(l)
            if gene == previous_gene:
                output_temp.append(line)
                gene_lens += 1
                if ncov >= normalized_cov_threshold:
                    gcov_meet_condition.append(ncov)
                #OUT.write(line)
            else:
                #OUT.write(line)
                ## measure if output the gene
                ## add a option for propotional threshold of win_len
                if len_prop_threshold == None:
                    pass
                else:
                    win_len = float(len_prop_threshold)*gene_lens
                if len(gcov_meet_condition) >= win_len:
                    ## skip output each genes as limit of space
                #    OUT = gzip.open(out_diretory+out_prefix+"."+gene+".mutrate.txt.gz",'wt')
                #    OUT.write("\n".join(output_temp)+"\n")
                #    OUT.close()
                    pass
                else:
                    pass
                ## reset varibles for the next cycle
                previous_gene = gene
                output_temp = [line]
                gcov_meet_condition = []
                gene_lens = 1
                if ncov >= normalized_cov_threshold:
                    gcov_meet_condition.append(ncov)
        ## run the last time for the last one
        ## add a option for propotional threshold of win_len
        if len_prop_threshold == None:
            pass
        else:
            win_len = float(len_prop_threshold)*gene_lens
        ## >>> last one start >>>
        if len(gcov_meet_condition) >= win_len:
        #    OUT = gzip.open(out_diretory+out_prefix+"."+gene+".mutrate.gene.gz",'wt')
        #    OUT.write("\n".join(output_temp)+"\n")
        #    OUT.close()
            pass  ## skip output each genes as limit of storage
        else:
            pass"""
        ## <<< last one end <<<

    # parse_single_gene(len_prop_threshold=0.5)
    if win_threshold < 1:
        parse_single_gene(normalized_cov_threshold=gene_output_cov_threshold, len_prop_threshold=win_threshold)
    else:
        parse_single_gene(win_len=win_threshold,normalized_cov_threshold=gene_output_cov_threshold)
    ## run process_line only
    #for l in input:
    #    gene, line, cov, ncov, pos = process_line(l)

def avg(test,digits=8):
    test = list(map(float,test))
    if len(test) == 0:
        return 0
    else:
        average = round(Decimal(sum(test)/len(test)),digits)
        return average

def make_gene_cov(fcfile):
    d = {}
    ref_idx={}
    total_reads_count = 0
    with open(fcfile) as fc:
        for l in fc:
            if l[0] == "#": continue
            i = l.strip('\n').split('\t')
            if i[0] == "Geneid":    continue
            [gene_len, gene_cov] = i[5:7]
            total_reads_count += int(gene_cov)
            d[i[1]] = (int(gene_len), int(gene_cov))
            ref_idx[i[1]]=[0 for n in range(int(gene_len))]
    return d,ref_idx,total_reads_count

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-i","--input",required=True, help="input bam-readcount output file")
    #parser.add_argument("-f","--ref",required=True, help="input ref fasta file")
    parser.add_argument("-o","--outprefix",default="-")
    parser.add_argument("-d","--directory",default="./temp")
    parser.add_argument("-c","--cov_threshold",default=-1, help="cov threshold for outputing each record")
    parser.add_argument("-gc","--gene_output_cov_threshold",default=10,help="cov threshold for outputing each genes")
    parser.add_argument("-w","--win_threshold",default=50,help="win threshold for outputing each genes")
    parser.add_argument("-e","--gene_readcounts",default=None,help="input the featureCounts file (optional, if not provided, normalized coverage will be set to NA)")
    parser.add_argument("-ins_len", "--ins_len", default=3, help="the maxium length of ins that count as mutation")
    parser.add_argument("-del_len", "--del_len", default=3, help="the maxium length of del that count as mutation")
    parser.add_argument("-ref_symbol", "--ref_symbol", action='store_true', default=False, help="the matched nt was labeled as = True or nt = False")
    args = parser.parse_args()

    #ref_idx=mk_idx(args.ref)
    if os.path.exists(args.directory):
        #os.remove(args.directory)
        shutil.rmtree(args.directory, ignore_errors=True) # remove unempty directory
        #os.makedirs(args.directory)
        folder = args.directory
        folder.strip('/')
    else:
        #os.makedirs(args.directory)
        folder = args.directory
        folder.strip('/')

    #ref_idx = mk_idx(args.ref)
    iterator_mutrate(args.input, args.gene_readcounts, folder,
                     args.outprefix,
                     win_threshold=float(args.win_threshold),
                     cov_threshold=float(args.cov_threshold),
                     gene_output_cov_threshold=float(args.gene_output_cov_threshold),
                     max_ins_len=int(args.ins_len),
                     max_del_len=int(args.del_len),
                     equal_symbol_ref=args.ref_symbol,
                    )
    
    print("Mutation rates calculation completed successfully.")

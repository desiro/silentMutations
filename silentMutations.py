#!/usr/bin/env python3
# script: silentMutations.py
# author: Daniel Desiro'
script_usage="""
usage
    silentMutations.py -f <in_fasta> -p <out_prefix> -s1 
    <name>:<frame>:<start>-<end> -s2 <name>:<frame>:<start>-<end> [options]

version
    silentMutations.py 0.0.1 (alpha)

dependencies
    numpy, ViennaRNA Package 2.4, VARNAv3-93.jar, inkscape

description
    Generates compensatory codon mutations with similar minimum free energy.
    Example call: python3 silentMutations.py -p example -f example.fa -s1 
    seq1:2:20-36 -s2 seq2:0:23-44 -cls ssRNA- -r -c -thr=4

--prefix,-p
    output prefix for result files

--fasta,-f
    fasta file with two sequences; will create all possible silent codon
    combinations for a given range of the two segments

--snip1,-s1
    define the range of the first snip with <name>:<frame>:<start>-<end>

--snip2,-s2
    define the range of the second snip with <name>:<frame>:<start>-<end>

--complement,-c
    creates complements of each strain if set (default: False)

--reverse,-r
    creates reverse of each strain if set (default: False)

--virusclass,-cls
    define virus class; (default: ssRNA+) (choices: dsDNA,ssDNA,dsRNA,ssRNA+,
    ssRNA-,ssRNA-RT,dsDNA-RT)

--filterperc,-prc
    define minimum percentage for the mutant-WT mfe to be below the WT-WT mfe
    (default: 0.5)

--upperdev,-udv
    define upper similarity threshold for the mutant-mutant mfe in contrast to
    the WT-WT mfe (default: 1.05)

--lowerdev,-ldv
    define lower similarity threshold for the mutant-mutant mfe in contrast to
    the WT-WT mfe (default: 0.95)

--mutrange,-mrg
    define percentage of similarity for the mean mfe of both single mutants
    (default: 0.1)

--noncoding1,-nc1
    specify if the first sequence includes a non-coding region (default: False)

--noncoding2,-nc2
    specify if the second sequence includes a non-coding region (default: False)

--mutations,-mut
    define maximum number of allowed mutations per sequence, lower values will
    improve runtime (default: 5)

--dangles,-dng
    use dangling ends for foldings (default: 0) (choices: 0,1,2,3)

--threads,-thr
    number of threads to use for RNAcofold (default: 1)

--varnabin,-var
    use this VARNA binary; example: VARNAv3-93.jar (default: )

--inkscbin,-ink
    use this Inkscape binary; example: inkscape (default: )

reference
    Reference.
"""


import argparse as ap
import sys
import os
import re
import time
from numpy import arange, mean, prod, array
import itertools
try:
    from RNA import fold_compound, cvar, CONSTRAINT_DB, CONSTRAINT_DB_DEFAULT
except:
    print("Error: The mandatory library ViennaRNA 2.4 is missing, please read the README.md for more information!")
    exit()
from subprocess import Popen, PIPE, call
from multiprocessing import Pool, Process, Manager, Lock


################################################################################
## main
################################################################################

def main(var_p  , var_f  , var_s1 , var_s2 , var_c  , var_r  , var_cls, var_prc,
         var_udv, var_ldv, var_mrg, var_nc1, var_nc2, var_mut, var_dng, var_thr,
         var_var, var_ink):
    ############################################################################
    ## read fasta file
    # fastaTest = "/home/vo54saz/projects/dd_influenza_packaging/scripts/codon_mut_test.fa"
    print(f"Status: Read fasta file ...")
    data_dict = readFasta(var_f, var_s1, var_s2)
    if len(data_dict.keys()) > 2:
        print("Error: More than two sequences in fasta file!")
        sys.exit()
    ############################################################################
    ## create snips
    print(f"Status: Create snips ...")
    snip_list = list()
    for header,sequence in data_dict.items():
        snip_list.append(createSnip(header, sequence, var_r, var_c))
    snip1, snip2 = snip_list
    getSnips(snip1, snip2, var_dng)
    ############################################################################
    ## make codon permutations
    print(f"Status: Create permutations ...")
    perm_list, base_list = list(), list()
    for snip,cdn in zip([snip1, snip2],[var_nc1, var_nc2]):
        perm_strings = getPermutations(snip.snip, var_cls, var_mut, cdn)
        filtered_strings, mms_dict = getStats(snip.snip, perm_strings, var_mut)
        perm_list.append(filtered_strings)
        base_list.append(perm_strings)
        mut_str = " ".join([f"{mu}={nu}" for mu,nu in sorted(mms_dict.items())])
        print(f"Status: {snip.name} mutations {mut_str}")
    print(f"Status: Permutations {snip1.name}:{len(base_list[0])} {snip2.name}:{len(base_list[1])}")
    print(f"Status: Max {var_mut} mutations {snip1.name}:{len(perm_list[0])} {snip2.name}:{len(perm_list[1])}")
    ############################################################################
    ## filter permutations
    print(f"Status: Filter permutations ...")
    constraint = snip1.lconst+snip2.rconst
    base_mfe, base_pattern = doCofold(f"{snip1.snip}&{snip2.snip}",constraint, var_dng)
    base_mfe = round(base_mfe,2)
    filtered_list = list()
    for snp,perms in zip([snip2, snip1], perm_list):
        filtered_list.append(filterPermutations(snp, perms, constraint, float(base_mfe), var_prc, var_thr, var_dng))
    if not filtered_list[0] or not filtered_list[1]:
        print("Error: Percentage for filtering too low!")
        sys.exit()
    print(f"Status: Filtered {snip1.name}:{len(filtered_list[0])} {snip2.name}:{len(filtered_list[1])}          ")
    print(f"Status: Range {snip1.name}:{snip1.s_start}-{snip1.s_end} {snip2.name}:{snip2.s_start}-{snip2.s_end}")
    ############################################################################
    ## make cofolds of permutations
    print(f"Status: Fold permutations ...")
    folds = foldPermutations(filtered_list[0], filtered_list[1], float(base_mfe), constraint, var_udv, var_ldv, var_thr, var_dng)
    print(f"Status: ... finished!                     ")
    ############################################################################
    ## get best mutation
    print(f"Status: Get best mutations ...")
    best = getBest(folds, snip1, snip2, constraint, var_dng, var_mrg)
    ############################################################################
    ## save folds
    print(f"Status: Save folds ...")
    var_p = makeDir(var_p)
    saveFolds(var_p, best, snip1, snip2, base_mfe, base_pattern, var_f, var_r, var_c, var_cls, var_prc, var_udv, var_ldv)
    ############################################################################
    ## create varna plots
    print(f"Status: Create structures ...")
    createStructures(var_p, best, var_var, var_ink, var_cls)




################################################################################
## functions
################################################################################

class RNAsnip(object):
    ## RNA object
    def __init__(self, name, seq, frame, start, end):
        self.name = name
        self.seq = seq
        self.frame = frame
        self.start = start
        self.end = end
        self.snip = ""
        self.s_start = 0
        self.s_end = 0
        self.lconst = ""
        self.rconst = ""
        self.extractRNA()
    ## extract RNA snip at the correct frame
    def extractRNA(self):
        sshift = (self.start-1-self.frame)%3
        eshift = 2-(self.end-1-self.frame)%3
        self.s_start, self.s_end = [self.start-sshift,self.end+eshift]
        self.snip = self.seq[self.s_start-1:self.s_end]
        self.lconst = "<"*len(self.snip)
        self.rconst = ">"*len(self.snip)
    ## crop snip left or right
    def lcrop(self):
        self.snip = self.snip[3:]
        self.s_start = self.s_start + 3
        self.lconst = self.lconst[3:]
        self.rconst = self.rconst[3:]
    def rcrop(self):
        self.snip = self.snip[:-3]
        self.s_end = self.s_end - 3
        self.lconst = self.lconst[:-3]
        self.rconst = self.rconst[:-3]

def getSnips(RNA1, RNA2, var_dng):
    ## get minimum snip size for folding from 2 RNAsnip objects
    mfe, pattern = doCofold(f"{RNA1.snip}&{RNA2.snip}",RNA1.lconst+RNA2.rconst, var_dng)
    pattern1, pattern2 = pattern.split("&")
    print(f"RNA: {RNA1.snip}&{RNA2.snip}")
    print(f"Pat: {pattern}")
    #print(pattern1)
    while pattern1[:3] == "...":
        pattern1 = pattern1[3:]
        RNA1.lcrop()
    while pattern1[-3:] == "...":
        pattern1 = pattern1[:-3]
        RNA1.rcrop()
    while pattern2[:3] == "...":
        pattern2 = pattern2[3:]
        RNA2.lcrop()
    while pattern2[-3:] == "...":
        pattern2 = pattern2[:-3]
        RNA2.rcrop()

def createSnip(header, sequence, var_r, var_c):
    ## extract snip position and create snips
    name, frame, srange = header.split(":")
    frame = int(frame)
    start, end = srange.split("-")
    start, end = int(start), int(end)
    if var_r: 
        snip_range = (len(sequence)-end+1, len(sequence)-start+1)
        frame = (len(sequence)-frame) % 3
    else:
        snip_range = (start, end)
        frame = frame
    vRNA = revComp(sequence, var_c, var_r)
    return RNAsnip(name, vRNA, frame, snip_range[0], snip_range[1])

def permutate(RNA, frame, c_dict, a_dict):
    ## get all permutations of an RNA snip
    aa_seq = aaSequence(RNA, frame, c_dict)
    codon_list = []
    for aa in aa_seq:
        codon_list.append(a_dict[aa])
    return codon_list

def getPermutations(snip, var_cls, var_mut, cdn):
    ## create all possible permutation strings 
    co_dict, aa_dict, rc_dict, ra_dict = makeCodons()
    if var_cls == "ssRNA-":
        c_dict = rc_dict
        a_dict = ra_dict
    else:
        c_dict = co_dict
        a_dict = aa_dict
    if not cdn:
        codon_list = permutate(snip,0,c_dict,a_dict)
    else:
        codon_list = [["A","C","G","U"] for i in range(len(snip))]
    nperm = prod(array([len(sub) for sub in codon_list]))
    print(codon_list)
    perm_list = list(itertools.product(*codon_list))
    perm_strings = list()
    for i,ptoup in enumerate(perm_list):
        #print(f"Status: Permutation: {i+1} of {nperm}      ", end="\r")
        perm_strings.append("".join(list(ptoup)))
    return perm_strings

def getStats(snip, perm_strings, var_mut):
    # filter to max mutations and get mutation stats
    filtered_strings, mms_dict = list(), dict()
    for ps in perm_strings:
        mms = sum(c1!=c2 for c1,c2 in zip(ps,snip))
        mms_dict[mms] = mms_dict.get(mms,0) + 1
        if mms <= var_mut:
            filtered_strings.append(ps)
    return filtered_strings, mms_dict

def filterPermutations(snp, perms, constraint, base_mfe, var_prc, var_thr, var_dng):
    ## filter permutations to be less than a percentage of the base mfe
    pdict = dict()
    filtered = list()
    for perm in perms:
        pdict[f"{snp.snip}&{perm}"] = (0.0,"")
    fold_dict = doMulti(pdict, constraint, f"{snp.name} filter", var_thr, var_dng)
    for iRNA,ires in fold_dict.items():
        #print(ires[0])
        if ires[0] >= base_mfe*var_prc:
            filtered.append(iRNA.split("&")[1])
    return filtered

def foldPermutations(perms1, perms2, base_mfe, constraint, var_udv, var_ldv, var_thr, var_dng):
    ## create folds for all permutation combinations
    pdict = dict()
    for prm1 in perms1:
        for prm2 in perms2:
            pdict[f"{prm1}&{prm2}"] = 0.0
    perm_dict = doMulti(pdict, constraint, "permutations", var_thr, var_dng)
    fold_dict = dict()
    lower = var_ldv * base_mfe
    upper = var_udv * base_mfe
    for iRNA,ires in perm_dict.items():
        if ires[0] >= upper and ires[0] <= lower:
            fold_dict[iRNA] = ires[0]
    return fold_dict

def doMulti(ldict, constraint, name, var_thr, var_dng):
    ## multi processing
    res_map = []
    res_items = Manager().dict()
    fold_dict = dict()
    if len(ldict) > var_thr:
        ldsplt = int(len(ldict)/var_thr)
    else:
        ldsplt = int(len(ldict))
    keylist = list(ldict.keys())
    for run_nr,i in enumerate(range(0,len(ldict),ldsplt)):
        if i+ldsplt > len(ldict):
            ldsplt = len(ldict)#-(i+ldsplt)+1
        lsdict = {}
        for key in keylist[i:i+ldsplt]:
            lsdict[key] = ldict[key]
        ## make multi process function tuple
        p = Process(target=multiFold, args=(lsdict, constraint, name, res_items, run_nr+1, var_dng))
        res_map.append(p)
        p.start()
    for res in res_map:
        res.join()
    for iRNA,ires in res_items.items():
        fold_dict[iRNA] = ires
    return fold_dict

def doCofold(RNA, constraint, var_dng):
    ## do Cofold
    cvar.dangles = var_dng
    fc = fold_compound(RNA)
    fc.constraints_add(constraint, CONSTRAINT_DB | CONSTRAINT_DB_DEFAULT)
    pattern, mfe = fc.mfe_dimer()
    ## get pattern
    pattern = pattern[:len(RNA.split("&")[0])]+"&"+pattern[len(RNA.split("&")[0]):]
    return mfe, pattern

def multiFold(lsdict, constraint, name, res_items, run_nr, var_dng):
    ## cofold multi
    #print(f"Status: Do {name} run {run_nr} ...             ")
    total, current = [len(lsdict.keys()), 0]
    for RNA in lsdict.keys():
        percentage = (100/total)*current
        if int(percentage*10) % 100 == 0 and int(percentage) != 0:
            print(f"Status: {name} run {run_nr} ... {int(percentage)} %             ", end="\r")
        mfe, pattern = doCofold(RNA, constraint, var_dng)
        res_items[RNA] = (float(mfe),pattern)
        current += 1
    #if name[-6:] != "filter":
    #    print(f"Status: {name} run {run_nr} ... finished!             ")

def extractRNA(RNA, frame, start, end):
    ## extract RNA snip at the correct frame
    sshift = (start-1-frame)%3
    eshift = 2-(end-1-frame)%3
    start, end = [start-sshift,end+eshift]
    snip = RNA[start-1:end]
    return snip

def makeCodons():
    ## make mutations in IAV strains
    # protein coding (+): 5'-3' AUG 
    # complement     (-): 3'-5' UAC
    # structure      (-): 5'-3' CAU
    co_dict = {"UUU":"F", "UUC":"F", "UUA":"L", "UUG":"L", "UCU":"S", "UCC":"S", "UCA":"S", "UCG":"S",
               "UAU":"Y", "UAC":"Y", "UAA":"*", "UAG":"*", "UGU":"C", "UGC":"C", "UGA":"*", "UGG":"W",
               "CUU":"L", "CUC":"L", "CUA":"L", "CUG":"L", "CCU":"P", "CCC":"P", "CCA":"P", "CCG":"P",
               "CAU":"H", "CAC":"H", "CAA":"Q", "CAG":"Q", "CGU":"R", "CGC":"R", "CGA":"R", "CGG":"R",
               "AUU":"I", "AUC":"I", "AUA":"I", "AUG":"M", "ACU":"T", "ACC":"T", "ACA":"T", "ACG":"T",
               "AAU":"N", "AAC":"N", "AAA":"K", "AAG":"K", "AGU":"S", "AGC":"S", "AGA":"R", "AGG":"R",
               "GUU":"V", "GUC":"V", "GUA":"V", "GUG":"V", "GCU":"A", "GCC":"A", "GCA":"A", "GCG":"A",
               "GAU":"D", "GAC":"D", "GAA":"E", "GAG":"E", "GGU":"G", "GGC":"G", "GGA":"G", "GGG":"G"}
    rc_dict = {} # reverse complement codons
    aa_dict = {} # amino acid dictionary
    ra_dict = {} # reverse complement aa dict
    for codon,aa in sorted(co_dict.items()):
        comp_codon = revComp(codon,True,True)
        rc_dict[comp_codon] = aa
        ra_list, aa_list = [ra_dict.get(aa,list()), aa_dict.get(aa,list())]
        ra_list.append(comp_codon)
        aa_list.append(codon)
        ra_dict[aa] = ra_list
        aa_dict[aa] = aa_list
    return co_dict, aa_dict, rc_dict, ra_dict

def revComp(RNA, var_c, var_r):
    ## complement dictionary, or transform DNA to RNA
    RNA = RNA.upper()
    D2Rc = {"A":"U","T":"A","U":"A","C":"G","G":"C","R":"Y","Y":"R","M":"K",\
            "K":"M","S":"W","W":"S","B":"V","V":"B","D":"H","H":"D"}
    if var_c:
        RNA = "".join(D2Rc[i] for i in RNA)
    else:
        RNA = RNA.replace("T","U")
    if var_r:
        RNA = RNA[::-1]
    return RNA

def aaSequence(RNA, frame, codon_table):
    ## get amino acid sequence
    slen = len(RNA)
    rlen = int((slen-frame)/3)*3
    tseq = RNA[frame:frame+rlen]
    codons = [tseq[i:i+3] for i in range(0, len(tseq), 3)]
    aa_seq = "".join(codon_table[c] for c in codons)
    return aa_seq

def getBest(folds, snip1, snip2, constraint, var_dng, var_mrg):
    ## get best mutations
    worst_1 = -100.0
    worst_2 = -100.0
    best = ()
    for iRNA,imfe in folds.items():
        iRNA1, iRNA2 = iRNA.split("&")
        sequences = [f"{snip1.snip}&{snip2.snip}",        f"{iRNA1}&{iRNA2}", 
                     f"{snip1.snip}&{iRNA2}",             f"{iRNA1}&{snip2.snip}"]
        names =     [f"{snip1.name}_WTx{snip2.name}_WT",  f"{snip1.name}_mutx{snip2.name}_mut",
                     f"{snip1.name}_WTx{snip2.name}_mut", f"{snip1.name}_mutx{snip2.name}_WT"]
        mfes = []
        patterns = []
        mut1pos = [i+1 for i in range(len(snip1.snip)) if snip1.snip[i] != iRNA1[i]]
        mut2pos = [i+1+len(snip1.snip) for i in range(len(snip2.snip)) if snip2.snip[i] != iRNA2[i]]
        mutations = [[], mut1pos+mut2pos, mut2pos, mut1pos]
        for vRNA,vtype in zip(sequences,names):
            mfe, pattern = doCofold(vRNA, constraint, var_dng)
            mfes.append(round(float(mfe),2))
            patterns.append(pattern)
        u_mean = ((mfes[2] + mfes[3]) / 2)*(1+var_mrg)
        l_mean = ((mfes[2] + mfes[3]) / 2)*(1-var_mrg)
        if u_mean <= mfes[2] <= l_mean and u_mean <= mfes[3] <= l_mean and \
           (mfes[2] + mfes[3] > worst_1 + worst_2 or \
           (mfes[2] + mfes[3] == worst_1 + worst_2 and len(mutations[1]) < len(best[4][1]))):
            best = (sequences, names, mfes, patterns, mutations)
            worst_1, worst_2 = mfes[2], mfes[3]
    return best #(sequences, names, mfes, patterns, mutations)

def readFasta(fasta_in, var_s1, var_s2):
    ## read fasta file
    data_dict = dict()
    vRNA = ""
    name = ""
    entry = 0
    with open(fasta_in, "r") as infa:
        for line in infa:
            line = line.strip()
            if re.match(">",line):
                entry += 1
                if entry > 2:
                    break
                if vRNA:
                    data_dict[var_s1] = vRNA
                name = line[1:]
                vRNA = ""
            else:
                vRNA += line
        data_dict[var_s2] = vRNA
    return data_dict

def doVARNA(var_p, RNA, constraint, var_var, title_name, hstring, astring, ystring, pstring, period, algorithm):
    ## calls VARNA
    C_call=["java", "-cp", var_var, "fr.orsay.lri.varna.applications.VARNAcmd",
            "-sequenceDBN", RNA, "-structureDBN", constraint, "-o", f"{var_p}.svg",
            "-algorithm", algorithm, "-periodNum", f"{period}",
            "-spaceBetweenBases", "0.72", "-title", title_name,
            "-highlightRegion" , hstring, "-annotations", astring,
            "-backbone", "#000000", "-bp", "#000000", "-bpStyle", "-lw"]
    if pstring: C_call+=["-basesStyle1", ystring, "-applyBasesStyle1on", pstring,]
    if var_var:
        call(C_call, shell=False)
    with open(f"{var_p}.sh", "w") as cvar:
        cvar.write(" ".join(C_call))
    print(f"Status: VARNA call written to \"{var_p}.sh\"")

def svg2pdf(var_p, var_ink):
    ## calls VARNA
    C_call=[var_ink, "-D", f"{var_p}.svg", "--without-gui", f"--export-pdf={var_p}.pdf"]
    call(C_call, shell=False)

def createStructures(var_p, best, var_var, var_ink, var_cls):
    ## create 2. Structures with varna
    #           dgreen     magenta    yellow      pink       dyellow    lblue      red        lgrey      dgrey
    #crange = ["#0c7226", "#c64fea", "#ffee00", "#ea104e", "#ffb600", "#00edff", "#ff0000", "#babbbc", "#616263"]
    ## get proteins
    co_dict, aa_dict, rc_dict, ra_dict = makeCodons()
    if var_cls == "ssRNA-": c_dict = rc_dict
    else: c_dict = co_dict
    for vRNA,vtype,bmfe,bpattern,bmut in zip(*best):
        Rlen = len(vRNA)-1
        hstring, pstring = "",""
        file_name = f"{var_p}_{vtype}"
        title_name = f"{vtype} ({round(bmfe,2)} kcal/mol)"
        astring = f"5':type=B,anchor=1;3':type=B,anchor={Rlen};"
        aa_seq = aaSequence(vRNA.replace("&",""), 0, c_dict)
        cstring = ["#babbbc","#616263"]
        ## annotate amino acids
        for j,aa in enumerate(aa_seq):
            cs = (j*3)+1
            astring += f"{aa}:type=B,anchor={cs+1};"
            hstring += f"{cs}-{cs+2}:fill={cstring[j%2]},outline={cstring[j%2]},radius=15;"
        ## highlight mutations
        ystring = "fill=#ff0000,outline=#ff0000,label=#ffffff"
        if bmut: pstring += ",".join([str(x) for x in bmut])
        doVARNA(file_name, vRNA, bpattern, var_var, title_name, hstring, astring, ystring, pstring, 10, "naview")
        if var_var and var_ink:
            svg2pdf(file_name, var_ink)

def makeDir(var_p):
    ## create directory
    dir_name, dir_base = var_p, var_p
    i = 1
    while os.path.isdir(dir_name):
        dir_name = f"{dir_base}_{i}"
        i += 1
    os.mkdir(dir_name)
    return os.path.join(dir_name,var_p)

def saveFolds(var_p, best, snip1, snip2, base_mfe, base_pattern, var_f, var_r, var_c, var_cls, var_prc, var_udv, var_ldv):
    ## save folds to file
    if not best:
        print("Error: No mutations found, try to adjust the parameters!")
        sys.exit()
    file_name = f"{var_p}.cmut"
    with open(file_name, "w") as outfile:
        co_dict, aa_dict, rc_dict, ra_dict = makeCodons()
        outfile.write(f"Settings: -f {var_f} -r {var_r} -c {var_c} -cls {var_cls} -prc {var_prc} -udv {var_udv} -ldv {var_ldv}\n")
        if var_cls == "ssRNA-": 
            outfile.write(f"Positive snip 1:   {snip1.name} {len(snip1.seq)-int(snip1.s_end)+1}-{len(snip1.seq)-int(snip1.s_start)+1} {revComp(snip1.snip,1,1)} {aaSequence(revComp(snip1.snip,1,1),0,co_dict)}\n")
            outfile.write(f"Positive snip 2:   {snip2.name} {len(snip2.seq)-int(snip2.s_end)+1}-{len(snip2.seq)-int(snip2.s_start)+1} {revComp(snip2.snip,1,1)} {aaSequence(revComp(snip2.snip,1,1),0,co_dict)}\n")
            c_dict = rc_dict
            sniptype = "Negative"
        else:
            c_dict = co_dict
            sniptype = "Positive"
        outfile.write(f"{sniptype} snip 1:   {snip1.name} {snip1.s_start}-{snip1.s_end} {snip1.snip} {aaSequence(snip1.snip,0,c_dict)}\n")
        outfile.write(f"{sniptype} snip 2:   {snip2.name} {snip2.s_start}-{snip2.s_end} {snip2.snip} {aaSequence(snip2.snip,0,c_dict)}\n")
        outfile.write(f"Sequence: {snip1.snip}&{snip2.snip}\n")
        outfile.write(f"Fold:     {base_pattern} {base_mfe} kcal/mol\n")
        outfile.write(f"Codon Mutations:\n")
        for vRNA,vtype,bmfe,bpattern,bmut in zip(*best):
            vRNA = "".join([s2 if s1==s2 else s2.lower() for s1,s2 in zip(f"{snip1.snip}&{snip2.snip}",vRNA)])
            outfile.write(f"{vtype}: {vRNA}\n")
            outfile.write(f"{vtype}: {bpattern} {round(bmfe,2)} kcal/mol\n")
        # get aa sequences
        if var_cls == "ssRNA-":
            c_dict = rc_dict
        else:
            c_dict = co_dict
        best_iRNA1, best_iRNA2 = best[0][1].split("&")
        best_name1, best_name2 = best[1][1].split("x")
        for vRNA,vtype in zip([best_iRNA1,best_iRNA2],[best_name1,best_name2]):
            outfile.write(f"Positive {vtype} aa: {revComp(vRNA,1,1)} {aaSequence(vRNA,0,c_dict)[::-1]}\n")
            outfile.write(f"Negative {vtype} aa: {vRNA} {aaSequence(vRNA,0,c_dict)}\n")




################################################################################
## parser
################################################################################

if __name__ == "__main__":

    ############################################################################
    ## get time and save call
    sscript = sys.argv[0]
    start_time = time.time()
    current_time = time.strftime('%x %X')
    scall = " ".join(sys.argv)
    with open(f"{sscript}.log", "a") as calllog:
        calllog.write(f"{sscript} started at {current_time}\n")
        calllog.write(f"Call: {scall}\n")
    print(f"Call: {scall}")
    print(f"Status: Started at {current_time}")
    ############################################################################
    ## transform string into int, float, bool if possible
    def trans(s):
        if isinstance(s, str):
            try: return int(s)
            except ValueError:
                try: return float(s)
                except ValueError:
                    if s in ["True", "False"]: return s == "True"
                    else: return s
        else: return s
    ############################################################################
    ## save documentation
    rx_text = re.compile(r"\n^(.+?)\n((?:.+\n)+)",re.MULTILINE)
    rx_oneline = re.compile(r"\n+")
    rx_options = re.compile(r"\((.+?)\:(.+?)\)")
    help_dict, type_dict, text_dict, mand_list = {}, {}, {}, []
    for match in rx_text.finditer(script_usage):
        argument = match.groups()[0].strip()
        text = " ".join(rx_oneline.sub("",match.groups()[1].strip()).split())
        argopts = {"action":"store", "help":None, "default":None, "choices":None}
        for option in rx_options.finditer(text):
            key = option.group(1).strip()
            var = option.group(2).strip()
            if var == "False": argopts["action"] = "store_true"
            if var == "True": argopts["action"] = "store_false"
            if key == "choices": var = [vs.strip() for vs in var.split(",")]
            if key == "default": var = trans(var)
            argopts[key] = var
        if argopts["default"]: add_default = f" (default: {str(argopts['default'])})"
        else: add_default = ""
        argopts["help"] = rx_options.sub("",text).strip()+add_default
        argnames = argument.split(",")
        if len(argnames) > 1:
            if argopts["default"] == None:
                mand_list.append(f"var_{argnames[1][1:]}")
            type_dict[f"var_{argnames[1][1:]}"] = argopts["default"]
            argopts["argshort"] = argnames[1]
            help_dict[argnames[0]] = argopts
        else:
            text_dict[argnames[0]] = argopts["help"]
    ############################################################################
    ## get arguments
    if text_dict["dependencies"]:
        desc = f"{text_dict['description']} (dependencies: {text_dict['dependencies']})"
    else:
        desc = text_dict['description']
    p = ap.ArgumentParser(prog=sscript, prefix_chars="-", usage=text_dict["usage"],
                          description=desc, epilog=text_dict["reference"])
    p.add_argument("-v", "--version", action="version", version=text_dict["version"])
    for argname,argopts in help_dict.items():
        argshort = argopts["argshort"]
        if argopts["choices"]:
            p.add_argument(argshort, argname,            dest=f"var_{argshort[1:]}",\
                           action=argopts["action"],     help=argopts["help"],\
                           default=argopts["default"],   choices=argopts["choices"])
        else:
            p.add_argument(argopts["argshort"], argname, dest=f"var_{argshort[1:]}",\
                           action=argopts["action"],     help=argopts["help"],\
                           default=argopts["default"])
    p._optionals.title = "arguments"
    opt = vars(p.parse_args())
    ############################################################################
    ## validate arguments
    if None in [opt[mand] for mand in mand_list]:
        print("Error: Mandatory arguments missing!")
        print(f"Usage: {text_dict['usage']} use -h or --help for more information.")
        sys.exit()
    for key,var in opt.items():
        if key not in mand_list:
            arg_req, arg_in = type_dict[key], trans(var)
            if type(arg_req) == type(arg_in):
                opt[key] = arg_in
            else:
                print(f"Error: Argument {key} is not of type {type(arg_req)}!")
                sys.exit()
    ############################################################################
    ## call main function
    try:
        main(**opt)
    except KeyboardInterrupt:
        print("Error: Interrupted by user!")
        sys.exit()
    except SystemExit:
        print("Error: System exit!")
        sys.exit()
    except Exception:
        print("Error: Script exception!")
        traceback.print_exc(file=sys.stderr)
        sys.exit()
    ############################################################################
    ## finish
    elapsed_time = time.time()-start_time
    elapsed_time = time.strftime("%H:%M:%S", time.gmtime(elapsed_time))
    current_time = time.strftime('%x %X')
    print(f"Status: Finished at {current_time} in {elapsed_time}")
    sys.exit(0)

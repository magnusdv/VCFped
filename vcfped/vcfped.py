#!/usr/bin/env python

"""
Module for identifying colse relationships (e.g. trios) in a multisample VCF file.

Author: Magnus Dehli Vigeland
Email: magnusdv at medisin.uio.no

This module exports a single function: vcfped.

Example usage as module:
>>>import vcfped
>>>res = vcfped(path_to_file)
>>>res

Example of command line usage:
>python vcfped.py path_to_vcf

To show help info:
>python vcfped.py -h
"""

__version__ = "1.0.2"
__author__ = "Magnus Dehli Vigeland"

import collections
import operator
import time
import itertools
import random
import sys
import os
import argparse
    

# Hard coding coordinates of pseudo-autosomal regions on X (hg19)
_PAR1_ = (60000, 2699521)
_PAR2_ = (154931043, 155270561)

def XminusPAR(linesplit, chromcol):
    ''' Determine if a variant is on X outside PARs.'''
    if not 'X' in linesplit[chromcol]:
        return False
    pos = float(linesplit[chromcol+1])
    par = _PAR1_[0] < pos < _PAR1_[1] or _PAR2_[0] < pos < _PAR2_[1]
    return not par
  
def _writeout(txt, fout, quiet=False):
    if not quiet:
        print >> fout, txt
        
def inspectFile(filename, fout=sys.stdout, quiet=False):
    _writeout("\n====FILE INFO====", fout, quiet)
    _writeout("File path: " + os.path.abspath(filename), fout, quiet)

    with open(filename, 'r') as f:
        # Identify VCF format if stated in first line
        preamble_firstline = f.readline().rstrip()
        line = preamble_firstline
        # Skip preamble lines
        while line.startswith("##"): 
            line = f.readline()
        
        headers = line.rstrip().split('\t')
        startpos = f.tell()
        firstvar = f.readline().rstrip().split('\t')
        
        if headers[0] == "#CHROM" and headers[5] == "QUAL" and headers[8] == 'FORMAT' and firstvar[8].startswith("GT"):
            if 'fileformat=' in preamble_firstline:
                fileformat = preamble_firstline[preamble_firstline.index('=')+1:]
            else: 
                fileformat = "Not stated in the first line, but looks like VCF format"
            chromcol, qualcol, formatcol = 0,5,8
        else:
            try:
                chromcol = next(i for i,h in enumerate(headers) if h.lower() in ("chr", "chrom", "chromosome"))
                formatcol = next(i for i,field in enumerate(firstvar) if field.startswith("GT:") or field=="GT")
            except:
                raise RuntimeError("Could not identify chromosome and/or format columns")
            try:
                qualcol = next(i for i,field in enumerate(headers) if field.lower() in ['qual', 'vcf_qual'])
                float(firstvar[qualcol])
            except:
                qualcol = None
            fileformat = "Non-standard VCF"
        
        _writeout("File format: " + fileformat, fout, quiet)
        _writeout("Chromosome column (index): %s (%d)" %(headers[chromcol], chromcol), fout, quiet)
        _writeout("QUAL column (index): %s (%d)" %(headers[qualcol], qualcol), fout, quiet)
        _writeout("FORMAT column (index): %s (%d)" %(headers[formatcol], formatcol), fout, quiet)

        format = firstvar[formatcol]
        formatfields = format.split(':')
        ADind = formatfields.index('AD') if 'AD' in formatfields else None
        DPind = formatfields.index('DP') if 'DP' in formatfields else None
        GQind = formatfields.index('GQ') if 'GQ' in formatfields else None
        _writeout("Format fields: " + ", ".join(formatfields), fout, quiet)
        
        samples = headers[(formatcol + 1):]
        nSamples = len(samples)
        
        _writeout("Number of samples: %d" % nSamples, fout, quiet)
        _writeout("Samples names:", fout, quiet)
        for i, s in enumerate(samples):
            _writeout('  %d: %s' %(i+1, s), fout, quiet)
        
        # Estimating the number of variants: If less than 2000, report exactly
        linecount = None
        skipsize = 2000
        for i in xrange(skipsize-1): 
            if f.readline() == '':
                linecount = i
                break
            
        # If > 2000: Estimate variant count using mean line width of first 2000 
        if linecount is None:
            after_skip = f.tell()
            f.seek(0, 2) # end of file
            endpos = float(f.tell())
            linecount = int((endpos-startpos)/(after_skip-startpos)*skipsize)
            counttype = 'Approximate'
        else:
            counttype = "Exact"
        _writeout("%s variant count: %d" %(counttype,linecount), fout, quiet)
        
    return dict(filename=filename, fileformat=fileformat, chromcol=chromcol, formatcol=formatcol, qualcol=qualcol, formatfields=formatfields, DPind=DPind, GQind=GQind, ADind=ADind, 
                counttype=counttype, linecount=linecount, samples=samples, nSamples=nSamples, startpos=startpos, endpos=endpos)
            
def getGTtuple(linesplit, firstSampleCol):
    # Mapping VCF genotypes to {AA, AB, BB}
    GT2AB = {'0/0':'AA', '0|0':'AA', '0/1':'AB', '0|1':'AB', '1|0':'AB', '1/1':'BB', '1|1':'BB'}
    try:
        return tuple(GT2AB[gt[:3]] for gt in linesplit[firstSampleCol:])
    except KeyError:
        return None

        
class Filter(object):
    def __init__(self, formatinfo, PASS=None, DPmin=None, GQmin=None, QUALmin=None, ADmin=None):
        self.nSamples=formatinfo['nSamples']
        self.qualcol=formatinfo['qualcol']
        self.ADind=formatinfo['ADind']
        self.DPind=formatinfo['DPind']
        self.GQind=formatinfo['GQind']
        self.PASS=PASS
        self.DPmin=DPmin
        self.GQmin=GQmin
        self.QUALmin=QUALmin
        self.ADmin=ADmin
        
    def __str__(self):
        pss =  'PASS: ' + str(self.PASS)
        qual = 'QUAL >= %d' %self.QUALmin if self.QUALmin is not None else 'QUAL: NA'
        dp = 'DP >= %d' %self.DPmin if self.DPmin is not None else 'DP: NA'
        gq = 'GQ >= %d' %self.GQmin if self.GQmin is not None else 'GQ: NA'
        ad = 'ALT ratio: %d - %d' %(self.ADmin, 100-self.ADmin) if self.ADmin is not None else 'ALT ratio: NA'
        return '\n'.join([pss, qual, dp, gq, ad])
        
    def filter(self):
        nSamples, qualcol, ADind, DPind, GQind = self.nSamples, self.qualcol, self.ADind, self.DPind, self.GQind
        PASS, QUALmin, DPmin, GQmin, ADmin = self.PASS, self.QUALmin, self.DPmin, self.GQmin, self.ADmin
        PASSind = qualcol + 1 if qualcol is not None else None
        
        PASStest = PASS and PASSind is not None
        QUALtest = QUALmin>0 and qualcol is not None
        ADtest = ADmin > 0 and ADind is not None
        DPtest = DPmin > 0 and DPind is not None
        GQtest = GQmin > 0 and GQind is not None

        if ADmin is not None:
            ADlow = ADmin/100.0
            ADhigh = 1 - ADlow
        
        split_it = ADtest or DPtest or GQtest
        pass_strings = ('PASS', '.')
        
        def _tester(linesplit):
            """Define the function to be returned"""
            try:
                if PASStest and linesplit[PASSind] not in pass_strings: return 0
                if QUALtest and float(linesplit[qualcol]) < QUALmin: return 0
                if split_it:
                    for sample in linesplit[-nSamples:]:
                        spl = sample.split(':')
                        if DPtest and int(spl[DPind]) < DPmin: return 0
                        if GQtest and int(spl[GQind]) < GQmin: return 0
                        if ADtest and spl[0][0] != spl[0][2]: # if heterozygot
                            ad = map(float, spl[ADind].split(','))
                            if not ADlow <= ad[1]/(ad[0]+ad[1]) <= ADhigh:
                                return 0
                return 1
            except Exception as e:
                return 0
        
        return _tester
 
        
def qualityDistrib(variantlines, formatinfo, fields=['QUAL', 'DP', 'GQ', 'AD']):
    if isinstance(variantlines, dict):
        variantlines = itertools.chain(*variantlines.values())
    ADind=formatinfo['ADind']
    DPind=formatinfo['DPind']
    GQind=formatinfo['GQind']
    qualcol=formatinfo['qualcol'] 
    firstSampleCol = formatinfo['formatcol'] + 1
    useAD = 'AD' in fields and ADind is not None
    useGQ = 'GQ' in fields and GQind is not None
    useDP = 'DP' in fields and DPind is not None
    useQUAL = 'QUAL' in fields and qualcol is not None
    distrib_data = collections.defaultdict(list)
    
    def _ALTratio(ADfield):
        ad0, ad1 = map(float, ADfield.split(','))
        rat = ad1/(ad0+ad1)
        return round(min(rat, 1-rat) * 100, 1)
    
    hets = ('0/1', '0|1', '1/0', '1|0')
    for s in variantlines:
        try:
            if useQUAL: distrib_data['QUAL'].append(float(s[qualcol]))
            gtcols = [a.split(':') for a in s[firstSampleCol:]]
            if useDP: distrib_data['DP'].append(min(int(b[DPind]) for b in gtcols))
            if useGQ: distrib_data['GQ'].append(min(int(b[GQind]) for b in gtcols))
            if useAD: 
                ad = [_ALTratio(b[ADind]) for b in gtcols if b[0] in hets]
                if ad:
                    distrib_data['AD'].append(min(ad))
        except Exception as e:
            continue
    
    for key in distrib_data:
        distrib_data[key].sort()
    return distrib_data
    
def qualityPercentile(distrib_data, p, fout, quiet):
    if isinstance(p, list):
        res = {key : [data[int(pp/100.0*len(data))] for pp in p] for key,data in distrib_data.iteritems()}
    else:
        res = {key : data[int(p/100.0*len(data))] for key,data in distrib_data.iteritems()}
        
    if not quiet:
        _writeout("\n====SCORE PERCENTILES (%s)====" % ", ".join(map(str, p)), fout)
        rows = ["  %s: %s" % (key, res[key]) for key in ['QUAL', 'DP', 'GQ', 'AD'] if key in res]
        _writeout('\n'.join(rows), fout)
    
    return res
 
    
# Utility function for the trio tests
def _conditionalDist(GTtriple, GT1, GT2, order, decimals=2):
    """Find the observed genotype freqs of third sample, conditional on given genotypes in the first two."""
    _perm = operator.itemgetter(*order)
    if GT1 == GT2:
        cond_dist = [GTtriple[_perm([GT1, GT2, gt])] for gt in ('AA', 'AB', 'BB')]
    else:
        cond_dist = [GTtriple[_perm([GT1, GT2, gt])] + GTtriple[_perm([GT2, GT1, gt])] for gt in ('AA', 'AB', 'BB')]
    
    cond_tot = sum(cond_dist) # total variants in this class 
    cond_dist_perc = [round(100.0 * x / cond_tot, decimals) for x in cond_dist] if cond_tot else [0,0,0] # distribution in percent
    return cond_tot, cond_dist_perc 


def checkTriple(AUTOSOMfilt, triple, percentile, threshold1, threshold2):
    # AUTOSOMfilt: dictionary. Keys are tuples (GT1, GT2, ...) whose value is a list of all variants with that gt combo.
    # percentile is just for book-keeping in output
    
    distr = collections.Counter()
    g = operator.itemgetter(*triple)
    for k,v in AUTOSOMfilt.iteritems():
        distr[g(k)] += len(v)
    SURV = sum(distr.values()) 
    
    testresults = []
    # test each permutation:
    for order in [(2,0,1), (1,2,0), (0,1,2)]:
        pivot = triple[order.index(2)]
        _perm = operator.itemgetter(*order)
        
        # test 1: AA + BB = AB
        tot_AA_BB = sum(distr[_perm(['AA', 'BB', _])] for _ in ('AA','AB','BB')) + \
                    sum(distr[_perm(['BB', 'AA', _])] for _ in ('AA','AB','BB'))
        AA_BB_AB = distr[_perm(['AA','BB','AB'])]+distr[_perm(['BB','AA','AB'])]
        AA_BB_ABp = round(100.0 * AA_BB_AB / tot_AA_BB, 2) if tot_AA_BB > 0 else 0
        
        # test 2: BB + BB = BB
        tot_BB_BB = sum(distr[_perm(['BB', 'BB', _])] for _ in ('AA','AB','BB'))
        BB_BB_BB = distr[('BB','BB','BB')] 
        BB_BB_BBp = round(100.0 * BB_BB_BB/tot_BB_BB, 2) if tot_BB_BB > 0 else 0
                
        if tot_AA_BB>=100: # TODO: should this be user argument?
            if AA_BB_ABp >= threshold1:
                verdict = 'Regular trio' if BB_BB_BBp >= threshold2 else 'Inverted/generational trio'
            else:
                verdict = 'Not trio'
        else:
            verdict = 'na'
        result_dict = dict(triple=triple, order=order, pivot=pivot, percentile=percentile, verdict=verdict, autos=SURV, tot_AA_BB=tot_AA_BB, 
                           AA_BB_AB=AA_BB_AB, AA_BB_ABp=AA_BB_ABp, tot_BB_BB=tot_BB_BB, BB_BB_BB=BB_BB_BB, BB_BB_BBp=BB_BB_BBp)
        testresults.append(result_dict)
    return testresults

def pretty_trio_table(TRIORES):
    table = ['\t'.join(['Triple', 'Pivot', 'Perc', 'Autos', 'AA+BB', '=AB', 'Test1', 'BB+BB', '=BB', 'Test2', 'Verdict'])]
    for d in TRIORES:
        triple_txt = '%d,%d,%d' % tuple(p+1 for p in d['triple'])
        pivot_out = d['pivot'] + 1
        test1 = '%.1f%%' % d['AA_BB_ABp']
        test2 = '%.1f%%' % d['BB_BB_BBp']
        printdat = [triple_txt, pivot_out, d['percentile'], d['autos'], d['tot_AA_BB'], d['AA_BB_AB'], \
                               test1, d['tot_BB_BB'], d['BB_BB_BB'], test2, d['verdict']]
        table.append('\t'.join(map(str,printdat)))
    return '\n'.join(table)  
    
def checkPair(AUTOSOMfilt, pair, percentile):
    testresults = []
    distr = collections.Counter()
    g = operator.itemgetter(*pair)
    for k,v in AUTOSOMfilt.iteritems():
        distr[g(k)] += len(v)
    SURV = sum(distr.values()) 
    
    # test for MZ twins: Percentage of all variants with identical genotypes
    MZ = sum(distr[(gt,gt)] for gt in ('AA','AB','BB'))
    MZp = round(100.0 * MZ/SURV, 2)
        
    # test for parent-offspring: Percentage of AA+BB out of all with BB. Should be (close to) 0.
    AA_BB = distr[('AA','BB')] + distr[('BB', 'AA')]
    tot_BB = AA_BB + distr[('AB','BB')] + distr[('BB','AB')] + distr[('BB','BB')]
    plenty_BB = tot_BB >= 100
    AA_BBp = round(100.0 * AA_BB/tot_BB, 2) if tot_BB>0 else 0
    
    if SURV >= 100 and MZp > 95:
        verdict = 'MZ twins'
    elif plenty_BB:
        verdict = 'Parent-child' if AA_BBp < 1 else 'Other/unrelated'
    else:
        verdict ='na'
    
    return dict(pair=pair, percentile=percentile, autos=SURV, MZ=MZ, MZp=MZp, tot_BB=tot_BB, AA_BB=AA_BB, AA_BBp=AA_BBp, verdict=verdict)

def pretty_pairwise_table(PAIRRES):
    table = ['\t'.join(['Pair', 'Perc', 'Autos', 'IBS2', 'MZtest', 'anyBB', 'IBS0', 'POtest', 'Verdict'])]
    for d in PAIRRES:
        pair_txt = '%d,%d' % tuple(p+1 for p in d['pair'])
        MZtest = '%.1f%%' % d['MZp']
        POtest = '%.1f%%' % d['AA_BBp']
        printdat = [pair_txt, d['percentile'], d['autos'], d['MZ'], MZtest, d['tot_BB'], d['AA_BB'], POtest, d['verdict']]
        table.append('\t'.join(map(str,printdat)))
    return '\n'.join(table) 
    
def all_variants(file, info, fout, quiet):
    _writeout('\n====VARIANT SAMPLE (BEFORE FILTERING)====', fout, quiet)
    chromcol = info['chromcol']
    firstSampleCol = info['formatcol'] + 1
    nSamples = info['nSamples']
    
    AUTOSOM = collections.defaultdict(list)
    XCHR = collections.defaultdict(list)
    AABB = 0
    with open(file, 'r') as f:
        f.seek(info['startpos']) # in files with LF, this ends in the middle of the 1st var!
        f.readline()
        for line in f:
            if not 'PASS' in line and not '\t.\t' in line: 
                continue
            linesplit = line.split('\t')
            if 'Y' in linesplit[chromcol]:
                continue
            GTtuple = getGTtuple(linesplit, firstSampleCol)
            if not GTtuple:
                continue
            if XminusPAR(linesplit, chromcol):
                XCHR[GTtuple].append(linesplit)
            else:
                AUTOSOM[GTtuple].append(linesplit)
                if 'AA' in GTtuple and 'BB' in GTtuple:
                    AABB += 1
    
    if not quiet:
        AUTant = sum(len(v) for v in AUTOSOM.values())
        Xant = sum(len(v) for v in XCHR.values())
        summary = 'Using all PASS variants:\n %d autosomal (including %d with AA+BB)\n %d X-linked.'%(AUTant, AABB, Xant)
        _writeout(summary, fout, quiet)
    
    return AUTOSOM, XCHR
    
def sample_variants(file, maxsize, AABBsize, Xsize, info, fout, quiet, chunks=4):
    _writeout('\n====VARIANT SAMPLE (BEFORE FILTERING)====', fout, quiet)
    chromcol = info['chromcol']
    firstSampleCol = info['formatcol'] + 1
    startpos, endpos = info['startpos'], info['endpos']
    randint = random.randint
    def _randpos():
        return randint(startpos, endpos)
    
    AUTOSOM = collections.defaultdict(list)
    XCHR = collections.defaultdict(list)
    total, AUTant, Xant, AABB = 0, 0, 0, 0
    with open(file, 'r') as f:
        while total < maxsize and (AABB < AABBsize or Xant < Xsize):
            total += 1
            f.seek(_randpos())
            f.readline() # Skip current line (because we might be in the middle of a line)
        
            for _ in xrange(chunks):
                line = f.readline()
                while line and 'PASS' not in line:  # avoiding EOF infinite loop
                    line = f.readline()
                if not line: #EOF!
                    break
                linesplit = line.split('\t')
                if 'Y' in linesplit[chromcol]:
                    continue
                GTtuple = getGTtuple(linesplit, firstSampleCol)
                if not GTtuple:
                    continue
                if XminusPAR(linesplit, chromcol):
                    XCHR[GTtuple].append(linesplit)
                    Xant += 1
                else:
                    AUTOSOM[GTtuple].append(linesplit)
                    AUTant += 1
                    if 'AA' in GTtuple and 'BB' in GTtuple:
                        AABB += 1
    _writeout('Sampled random variants:\n %d autosomal variants (including %d with AA+BB)\n %d on X.'%(AUTant, AABB, Xant), fout, quiet)
     
    return AUTOSOM, XCHR

def inferGenders(XCHRfilt, formatinfo, threshMale, percentile):
    nSamples = formatinfo['nSamples']
    Xdistribs = [collections.Counter() for _ in range(nSamples)]
    Xtot = 0
    for gtcomb, lines in XCHRfilt.iteritems():
        nvars = len(lines)
        Xtot += nvars
        for i, gt in enumerate(gtcomb):
            Xdistribs[i][gt] += nvars 
    
    common = dict(percentile=percentile, Xtot=Xtot)
    def dist2gender(Xdist):
        denom = Xdist['AB']+Xdist['BB']
        Xhetp = round(100.0 * Xdist['AB']/denom, 1) if denom > 0 else 0
        if denom >=25:
            if Xhetp < threshMale:
                gender = 'Male' 
            elif Xhetp < 2*threshMale:
                gender = '?'
            else:
                gender = 'Female' 
        else: gender = 'na'
        return dict(Xhetp=Xhetp, gender=gender)
    
    res = []
    for i, Xdist in enumerate(Xdistribs):
        r = dist2gender(Xdist)
        ## adding other info. Note AA,AB,BB: if absent, they default to 0.
        r.update(dict(sample=i, AA=Xdist['AA'], AB=Xdist['AB'], BB=Xdist['BB'], **common))
        res.append(r)
    
    return res
    
def pretty_gender_table(GENDERRES):
    table = ['\t'.join(['Sample', 'Perc', 'Xvars', 'AA', 'AB', 'BB', 'Xhet', 'Gender'])]
    for d in GENDERRES:
        xhet_perc = '%.1f%%' % d['Xhetp']
        printdat = [d['sample']+1, d['percentile'], d['Xtot'], d['AA'], d['AB'], d['BB'], xhet_perc, d['gender']]
        table.append('\t'.join(map(str,printdat)))
    return '\n'.join(table)  
    
def bestTrios(triodata, fout, quiet, report_all=False):
    best = []
    for k,v in itertools.groupby(triodata, key=lambda x: x['triple']):
        use = [r for r in v if r['verdict'] != 'na']
        if not use:
            continue
        b = sorted(use, key=lambda x: (-round(x['AA_BB_ABp'], 1), x['percentile']))
        best.append(b[0])
    
    if not quiet:
        _writeout('\nMost confident trio calls:', fout)
        printtable = pretty_trio_table(best) if best else "No trios detected"
        _writeout(printtable, fout)
    return best

def bestPairs(pairdata, fout, quiet):
    best = []
    for k,v in itertools.groupby(pairdata, key=lambda x: x['pair']):
        use = [r for r in v if r['verdict'] != 'na']
        if not use:
            continue
        verdict = use[-1]['verdict'] 
        if verdict == 'Parent-child':
            use.sort(key=lambda x: (round(x['AA_BBp'], 1), x['percentile']))
        elif verdict == 'MZ twins':
            use.sort(key=lambda x: (-round(x['MZp'], 1), x['percentile']))
        else:
            continue
        best.append(use[0])
    
    if not quiet:
        _writeout('\nMost confident pairwise calls:', fout)
        printtable = pretty_pairwise_table(best) if best else "No parent-child pairs or MZ twins detected"
        _writeout(printtable, fout)
    return best

def bestGenders(genderdata, fout, quiet):
    best = []
    sortfun = lambda x: (round(x['Xhetp'], 1), x['percentile'])
    for k,v in itertools.groupby(genderdata, key=lambda x: x['sample']):
        allcalls = list(v)
        noNA = [r for r in allcalls if r['gender'] != 'na']
        if not noNA: 
            best.append(allcalls[0])
            continue
        verdicts = {r['gender'] for r in noNA}
        if 'Male' in verdicts and 'Female' in verdicts:
            # if both genders are called, choose the last one
            noQ = [r for r in noNA if r['gender'] != '?']
            best.append(noQ[-1])
            continue
        if verdicts == {'?'}:
            verd = '?'
        else:
            verd = list(verdicts.difference({'?'}))[0]
        use = sorted([r for r in noNA if r['gender'] == verd], key=sortfun)
        best.append(use[0])
    
    if not quiet:
        _writeout('\nMost confident gender calls:', fout)
        _writeout(pretty_gender_table(best), fout)
    return best
    
def vcfped(file, threshTrio1=90, threshTrio2=95, threshMale=5, exactmax=100000, samplesize=10000, samplesizeAABB=1000, percentiles=[10,30,50], nogender=False, nopairwise=False, notrio=False, out=sys.stdout, quiet=False):
    """Detect close relationships (e.g. trios and parent-child pairs) in a multisample VCF file.
    
    The input file should contain jointly called variants from two or more individuals.
    
    Args:
        file (str): The path to a variant file with exactly 3 samples. The file format must be VCF or
            VCF-like. More precisely, the format requirements are:
                * Initial preamble lines start with '##'
                * Columns (after preamble) are tab-separated.
                * The final columns of the file must contain genotype data, in the form of a (VCF-like) format column 
                  (e.g. GT:AD:DP:GQ:PL) followed by one column per sample.
        thresh1 (int): Threshold for test 1. To pass test1, the percentage of "AA + BB" variants 
            resulting in AB must be at least thresh1. (Default = 90%)
        thresh2 (int): Threshold for test 1. To pass test1, the percentage of "BB + BB" variants 
            resulting in BB must be at least thresh2. (Default = 95%)
        threshMale (int): Threshold for gender estimation. If heterozygosity (in percent) on X (minus 
            pseudoautosomal regions) is lower than this, the sample is set to 'male'. (Default = 10%)
        
        quiet (bool): If True, print only conclusion (Default = False)
    """
    
    fout = sys.stdout if out in ("", "STDOUT") else open(out, 'w')
    
    callTrios, callPairs, callGenders = not notrio, not nopairwise, not nogender 
    st = time.time()
    
    # inspect file and guess various format details.
    info = inspectFile(file, fout=fout, quiet=quiet)
    nSamples = info['nSamples']
    samples = info['samples']
    linecount = info['linecount']
    chromcol = info['chromcol']
    firstSampleCol = info['formatcol'] + 1
    
    if linecount > exactmax:
        AUTOSOM, XCHR = sample_variants(file, maxsize=samplesize*2, AABBsize=samplesizeAABB, Xsize=100, info=info, fout=fout, quiet=quiet)
    else:
        AUTOSOM, XCHR = all_variants(file, info=info, fout=fout, quiet=quiet)
    
    qualdist = qualityDistrib(AUTOSOM, info) # fields = ['DP']
    qualityPercentile(qualdist, percentiles, fout, quiet)

    autosom_filt = AUTOSOM.copy()
    x_filt = XCHR.copy()
    triples = list(itertools.combinations(range(nSamples), 3))
    pairs = list(itertools.combinations(range(nSamples), 2))
     
    TRIORES, PAIRRES, GENDERRES = [], [], []
    for p in percentiles:
        qualperc = qualityPercentile(qualdist, p, fout=None, quiet=True)
        filterargs = {key+'min':val for key,val in qualperc.iteritems()}
        filterObject = Filter(formatinfo=info, PASS=True, **filterargs)
        filterFUN = filterObject.filter()
        
        autosom_filt = {k:filter(filterFUN, v) for k,v in autosom_filt.iteritems()}
        x_filt = {k:filter(filterFUN, v) for k,v in x_filt.iteritems()}
        
        N = sum(len(v) for v in autosom_filt.values()) 
        
        if N == 0: break #raise RuntimeError("No variants surviving filter.")

        #### GENDERS ###
        if callGenders:
            genders = inferGenders(x_filt, info, threshMale, percentile=p)
            GENDERRES.extend(genders)
          
        #### PAIRWISE ###
        if callPairs:
            for pair in pairs:
                testres = checkPair(autosom_filt, pair, percentile=p)
                PAIRRES.append(testres)
                    
        #### TRIOS ###
        if callTrios:
            for triple in triples:
                testres = checkTriple(autosom_filt, triple, percentile=p, threshold1=threshTrio1, threshold2=threshTrio2)
                TRIORES.extend(testres)
    
    _writeout('\nAll analysis finished.\nTime used: %.2f seconds' % (time.time()-st), fout, quiet)
        
    if callGenders:
        GENDERRES.sort(key=lambda x: (x['sample'], x['percentile']))
        if not quiet:
            _writeout('\n====GENDERS====', fout)
            _writeout(pretty_gender_table(GENDERRES), fout)
        
    if callPairs:    
        PAIRRES.sort(key=lambda x: (x['pair'], x['percentile']))
        if not quiet:
            _writeout('\n====PAIRWISE RELATIONS====', fout)
            printtable = pretty_pairwise_table(PAIRRES) if PAIRRES else "Not enough individuals in the file"
            _writeout(printtable, fout)
            
    if callTrios:
        TRIORES.sort(key=lambda x: (x['triple'], x['pivot'], x['percentile']))
        if not quiet:
            _writeout('\n====TRIO RELATIONS====', fout)
            printtable = pretty_trio_table(TRIORES) if TRIORES else "Not enough individuals in the file"
            _writeout(printtable, fout)
        
    if callGenders:
        bestG = bestGenders(GENDERRES, fout=fout, quiet=False)
    if callPairs:
        best2 = bestPairs(PAIRRES, fout=fout, quiet=False)
    if callTrios:
        best3 = bestTrios(TRIORES, report_all=True, fout=fout, quiet=False)
          
    if fout is not sys.stdout:
        fout.close()
    
def main():    
    parser = argparse.ArgumentParser()
    parser.add_argument("file", help="Path to VCF file")
    parser.add_argument("-o", "--out", type=str, help="Output file name", default="STDOUT")
    parser.add_argument("--no-gender", help="Skip gender analysis", dest='nogender', action='store_true')
    parser.add_argument("--no-pairwise", help="Skip pairwise analysis", dest='nopairwise', action='store_true')
    parser.add_argument("--no-trio", help="Skip trio analysis", dest='notrio', action='store_true')
    parser.add_argument("-e", "--exact-max", dest='exactmax', type=int, help="If approx. line count exceeds this, apply random sampling", default=100000)
    parser.add_argument("-s", "--sample-size", dest='samplesize', type=int, help="Sample at least this many variant lines (if sampling)", default=10000)
    parser.add_argument("-sAABB", "--sample-size-AABB", dest='samplesizeAABB', type=int, help="Sample at least this many lines where both 0/0 and 1/1 occur as genotypes (if sampling)", default=1000)
    parser.add_argument("-p", "--percentiles", nargs='+', type=int, help="Filtering percentiles", default=[10,20,30,40,50])
    parser.add_argument("-t1", "--threshTrio1", type=int, help="Lower limit for passing trio test 1 (AA + BB = AB)", default=90)
    parser.add_argument("-t2", "--threshTrio2", type=int, help="Lower limit for passing trio test 2 (BB + BB = BB)", default=95)
    parser.add_argument("-male", "--threshMale", type=int, help="Maximum male heterozygosity (percent) on X", default=10)
    parser.add_argument("--quiet", help="Print less information", action='store_true')
    parser.add_argument('--version', action='version', version='VCFped %s'%__version__)
    
    args = parser.parse_args(sys.argv[1:])  
    for pp in args.percentiles:
        if not 0 <= pp < 100: parser.error("argument -p/--percentiles: invalid percentile: %d" %pp)
    thresh_error = "argument -%s/--%s: invalid threshold (must be in [0, 100]): %d"
    if not 0<= args.threshTrio1 <=100: parser.error(thresh_error % ('t1','threshTrio1', args.threshTrio1))
    if not 0<= args.threshTrio2 <=100: parser.error(thresh_error % ('t2','threshTrio2', args.threshTrio2))
    if not 0<= args.threshMale <=100: parser.error(thresh_error % ('male','threshMale', args.threshMale))
    
    if not args.quiet:
        form = '{:>50} : {}'
        print 'VCFped %s\n'%__version__
        print 'Options in effect:'
        print form.format('Input file', args.file)
        print form.format('Output to', args.out)
        print form.format('Sample if approximate count exceeds', args.exactmax)
        print form.format('Minimum sample size', args.samplesize)
        print form.format('Minimum # samples with both AA and BB', args.samplesizeAABB)
        print form.format('Percentile ranks used in filtering procedure', ' '.join(map(str, args.percentiles)))
        print form.format('Threshold (%) for trio test 1 (AA+BB=AB)', args.threshTrio1)
        print form.format('Threshold (%) for trio test 2 (BB+BB=BB)', args.threshTrio2)
        print form.format('Max male heterozygosity (%) on X', args.threshMale)
    vcfped(**vars(args))
    
if __name__ == '__main__':
    main()
    
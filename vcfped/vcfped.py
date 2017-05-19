#!/usr/bin/env python

"""
Module for identifying trios and close pairwise relationships in a multisample VCF file.

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

__version__ = "1.2.0"
__author__ = "Magnus Dehli Vigeland"

import collections
import operator
import time
import itertools
import random
import sys
import os
import argparse
from tabulate import tabulate
    

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
  
def _writeout(txt, log, quiet=False):
    if log:
        print >> log, txt
    if not quiet:
        print txt
    
def inspectFile(filename, log, quiet=False):
    _writeout("\n====FILE INFO====", log, quiet)
    _writeout("File path: " + os.path.abspath(filename), log, quiet)

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
        
        _writeout("File format: " + fileformat, log, quiet)
        _writeout("Chromosome column [index]: %s [%d]" %(headers[chromcol], chromcol), log, quiet)
        _writeout("QUAL column [index]: %s [%d]" %(headers[qualcol], qualcol), log, quiet)
        _writeout("FORMAT column [index]: %s [%d]" %(headers[formatcol], formatcol), log, quiet)

        format = firstvar[formatcol]
        formatfields = format.split(':')
        ADind = formatfields.index('AD') if 'AD' in formatfields else None
        DPind = formatfields.index('DP') if 'DP' in formatfields else None
        GQind = formatfields.index('GQ') if 'GQ' in formatfields else None
        _writeout("Format fields: " + ", ".join(formatfields), log, quiet)
        
        available_vars = [v for v,ind in zip(['QUAL','DP','GQ','AD'], [qualcol, DPind, GQind, ADind]) if ind is not None]
        _writeout("Available filtering variables: " + ", ".join(available_vars), log, quiet)
        
        samples = headers[(formatcol + 1):]
        nSamples = len(samples)
        
        _writeout("Number of samples: %d" % nSamples, log, quiet)
        _writeout("Sample names:", log, quiet)
        for i, s in enumerate(samples):
            _writeout('  %d: %s' %(i+1, s), log, quiet)
        
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
            endpos = None
        _writeout("%s variant count: %d" %(counttype,linecount), log, quiet)
        
    return dict(filename=filename, fileformat=fileformat, chromcol=chromcol, formatcol=formatcol, qualcol=qualcol, 
                formatfields=formatfields, DPind=DPind, GQind=GQind, ADind=ADind, counttype=counttype, available_vars=available_vars,
                linecount=linecount, samples=samples, nSamples=nSamples, startpos=startpos, endpos=endpos)

_GT2AB_ = {'0/0':'AA', '0|0':'AA', '0/1':'AB', '0|1':'AB', '1|0':'AB', '1/1':'BB', '1|1':'BB'}                
def getGTtuple(linesplit, firstSampleCol):
    # Mapping VCF genotypes to {AA, AB, BB}, anything else to '-'
    return tuple(_GT2AB_.get(gt[:3], '-') for gt in linesplit[firstSampleCol:])
    
        
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

def all_variants(file, info, log, quiet):
    _writeout('\n====VARIANT SAMPLE (BEFORE FILTERING)====', log, quiet)
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
    
    if log or not quiet:
        AUTant = sum(len(v) for v in AUTOSOM.values())
        Xant = sum(len(v) for v in XCHR.values())
        summary = 'Using all PASS variants:\n %d autosomal (including %d with AA+BB)\n %d X-linked.'%(AUTant, AABB, Xant)
        _writeout(summary, log, quiet)
    
    return AUTOSOM, XCHR
    
def sample_variants(file, maxsize, AABBsize, Xsize, info, log, quiet, chunks=4):
    _writeout('\n====VARIANT SAMPLE (BEFORE FILTERING)====', log, quiet)
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
    _writeout('Sampled random variants:\n %d autosomal variants (including %d with AA+BB)\n %d on X.'%(AUTant, AABB, Xant), log, quiet)
     
    return AUTOSOM, XCHR
        
def qualityDistrib(variantlines, formatinfo, variables=['QUAL', 'DP', 'GQ', 'AD']):
    if isinstance(variantlines, dict):
        variantlines = itertools.chain(*variantlines.values())
    useVariables = [v for v in formatinfo['available_vars'] if v in variables]
    ADind = formatinfo['ADind'] if 'AD' in useVariables else None
    DPind = formatinfo['DPind'] if 'DP' in useVariables else None
    GQind = formatinfo['GQind'] if 'GQ' in useVariables else None
    qualcol = formatinfo['qualcol'] if 'QUAL' in useVariables else None
    firstSampleCol = formatinfo['formatcol'] + 1
    distrib_data = collections.defaultdict(list)
    
    def _ALTratio(ADfield):
        ad0, ad1 = map(float, ADfield.split(','))
        rat = ad1/(ad0+ad1)
        return round(min(rat, 1-rat) * 100, 1)
    
    hets = ('0/1', '0|1', '1/0', '1|0')
    for s in variantlines:
        try:
            if qualcol: distrib_data['QUAL'].append(float(s[qualcol]))
            gtcols = [a.split(':') for a in s[firstSampleCol:]]
            if DPind: distrib_data['DP'].append(min(int(b[DPind]) for b in gtcols))
            if GQind: distrib_data['GQ'].append(min(int(b[GQind]) for b in gtcols))
            if ADind: 
                ad = [_ALTratio(b[ADind]) for b in gtcols if b[0] in hets]
                if ad:
                    distrib_data['AD'].append(min(ad))
            # TODO change to indiv distrib:
            # for b in gtcols:
            #    if not b[0] in _GT2AB_: # skip e.g. 1/2 or ./.
            #        continue
            #    if DPind: distrib_data['DP'].append(int(b[DPind]))
            #    if GQind: distrib_data['GQ'].append(int(b[GQind]))
            #    if ADind and b[0] in hets: distrib_data['AD'].append(_ALTratio(b[ADind]))
        except Exception as e:
            continue
    
    for key in distrib_data:
        distrib_data[key].sort()
    return distrib_data
    
def qualityPercentile(distrib_data, p, log, quiet):
    if isinstance(p, list):
        res = {key : [data[int(pp/100.0*len(data))] for pp in p] for key,data in distrib_data.iteritems()}
    else:
        res = {key : data[int(p/100.0*len(data))] for key,data in distrib_data.iteritems()}
        
    if log or not quiet:
        _writeout("\n====SCORE PERCENTILES (%s)====" % ", ".join(map(str, p)), log, quiet)
        rows = ["  %s: %s" % (key, res[key]) for key in ['QUAL', 'DP', 'GQ', 'AD'] if key in res]
        _writeout('\n'.join(rows), log, quiet)
    
    return res
    
def checkTriple(AUTOSOMfilt, triple, percentile, T1_thresh, T2_thresh):
    # AUTOSOMfilt: dictionary. Keys are tuples (GT1, GT2, ...) whose value is a list of all variants with that gt combo.
    # percentile is just for book-keeping in output
    
    distr = collections.Counter()
    g = operator.itemgetter(*triple)
    for k,v in AUTOSOMfilt.iteritems():
        if '-' in g(k): continue
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
            if AA_BB_ABp >= T1_thresh:
                verdict = 'Regular trio' if BB_BB_BBp >= T2_thresh else 'Inverted/generational trio'
            else:
                verdict = 'Not trio'
        else:
            verdict = 'na'
        result_dict = dict(triple=triple, order=order, pivot=pivot, percentile=percentile, verdict=verdict, autos=SURV, tot_AA_BB=tot_AA_BB, 
                           AA_BB_AB=AA_BB_AB, AA_BB_ABp=AA_BB_ABp, tot_BB_BB=tot_BB_BB, BB_BB_BB=BB_BB_BB, BB_BB_BBp=BB_BB_BBp)
        testresults.append(result_dict)
    return testresults

def checkPair(AUTOSOMfilt, pair, percentile, MZ_thresh=95, PO_thresh=99):
    testresults = []
    distr = collections.Counter()
    g = operator.itemgetter(*pair)
    for k,v in AUTOSOMfilt.iteritems():
        if '-' in g(k): continue
        distr[g(k)] += len(v)
    autos_total = sum(distr.values()) 
    
    # MZ score = freq(IBS=2 | none are AA)
    neitherAA = distr[('AB','AB')] + distr[('AB', 'BB')] + distr[('BB','AB')] + distr[('BB', 'BB')]
    IBS2_neitherAA = distr[('AB','AB')] + distr[('BB', 'BB')]
    MZp = 100.0 * IBS2_neitherAA/neitherAA if neitherAA>0 else 0
        
    # Parent-offspring test: PO score = freq(IBS>0 | either is BB)
    eitherBB = distr[('AA','BB')] + distr[('BB', 'AA')] + distr[('AB','BB')] + distr[('BB','AB')] + distr[('BB','BB')]
    IBS12 = eitherBB - distr[('AA','BB')] - distr[('BB', 'AA')]
    POp = 100.0 * IBS12/eitherBB if eitherBB>0 else 0
    
    if neitherAA >= 100 and MZp >= MZ_thresh:
        verdict = 'MZ twins'
    elif eitherBB >= 100 and POp >= PO_thresh:
        verdict = 'Parent-child' 
    elif eitherBB >= 100 and POp < PO_thresh:
        verdict = 'Other/unrelated'
    else:
        verdict ='na'
    
    return dict(pair=pair, percentile=percentile, autos=autos_total, neitherAA=neitherAA, IBS2=IBS2_neitherAA, MZp=MZp, eitherBB=eitherBB, IBS12=IBS12, POp=POp, verdict=verdict)

def inferGenders(XCHRfilt, formatinfo, percentile, FEMALE_thresh, MALE_thresh):
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
        Xhetp = 100.0 * Xdist['AB']/denom if denom > 0 else 0
        if denom >=25:
            if Xhetp <= MALE_thresh:
                gender = 'Male' 
            elif Xhetp < FEMALE_thresh:
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
    
def trio_table(TRIORES):
    table = [['Triple', 'Pivot', 'Perc', 'Autos', 'AA+BB', '=AB', 'Test1', 'BB+BB', '=BB', 'Test2', 'Verdict']]
    for d in TRIORES:
        triple_txt = '%d,%d,%d' % tuple(p+1 for p in d['triple'])
        pivot_out = d['pivot'] + 1
        test1 = '%.1f' % d['AA_BB_ABp']
        test2 = '%.1f' % d['BB_BB_BBp']
        printdat = [triple_txt, pivot_out, d['percentile'], d['autos'], d['tot_AA_BB'], d['AA_BB_AB'], \
                               test1, d['tot_BB_BB'], d['BB_BB_BB'], test2, d['verdict']]
        table.append(printdat)
    return table
    
def pairwise_table(PAIRRES):
    table = [['Pair', 'Perc', 'Autos', 'NeitherAA', 'IBS2', 'MZscore', 'EitherBB', 'IBS>0', 'POscore', 'Verdict']]
    for d in PAIRRES:
        pair_txt = '%d,%d' % tuple(p+1 for p in d['pair'])
        MZscore = '%.1f' % d['MZp']
        POscore = '%.1f' % d['POp']
        printdat = [pair_txt, d['percentile'], d['autos'], d['neitherAA'], d['IBS2'], MZscore, d['eitherBB'], d['IBS12'], POscore, d['verdict']]
        table.append(printdat)
    return table
    
def gender_table(GENDERRES):
    table = [['Sample', 'Perc', 'Xvars', 'AA', 'AB', 'BB', 'Xhet', 'Gender']]
    for d in GENDERRES:
        XHET = '%.1f' % d['Xhetp']
        printdat = [d['sample']+1, d['percentile'], d['Xtot'], d['AA'], d['AB'], d['BB'], XHET, d['gender']]
        table.append(printdat)
    return table
    
def bestTrios(triodata, reportall):
    best = []
    for k,v in itertools.groupby(triodata, key=lambda x: x['triple']):
        allcalls = list(v)
        nonNA = [r for r in allcalls if r['verdict'] != 'na']
        if not nonNA:
            if reportall: 
                best.append(allcalls[0])
            continue
        nonNA.sort(key=lambda x: (-round(x['AA_BB_ABp'], 1), x['percentile']))
        b = nonNA[0]
        if b['verdict'][:3] in ('Reg', 'Inv')  or reportall: 
            best.append(b)
    
    return trio_table(best)
    
def bestPairs(pairdata, reportall):
    best = []
    for k,v in itertools.groupby(pairdata, key=lambda x: x['pair']):
        allcalls = list(v)
        nonNA = [r for r in allcalls if r['verdict'] != 'na']
        if not nonNA:
            if reportall: 
                best.append(allcalls[0])
            continue
        verdict = nonNA[-1]['verdict'] 
        if verdict == 'MZ twins':
            nonNA.sort(key=lambda x: (-round(x['MZp'], 1), x['percentile']))
        elif verdict == 'Parent-child' or reportall:
            nonNA.sort(key=lambda x: (round(x['POp'], 1), x['percentile']))
        else:
            continue
        best.append(nonNA[0])
    
    return pairwise_table(best)

def bestGenders(genderdata):
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
    
    return gender_table(best)
    
    
def vcfped(file, quiet=True, reportall=False, prefix=None, variables=['QUAL','DP','GQ','AD'], percentiles=[10,30,50], 
           T1_thresh=90, T2_thresh=95, MZ_thresh=95, PO_thresh=99, MALE_thresh=5, FEMALE_thresh=25, 
           nogender=False, nopairwise=False, notrio=False, exactmax=150000, samplesize=10000, samplesizeAABB=500):
    """Detect close relationships (e.g. trios and parent-child pairs) in a multisample VCF file.
    
    The input file should contain jointly called variants from two or more individuals.
    
    Args:
        file (str): The path to a variant file with exactly 3 samples. The file format must be VCF or
            VCF-like. More precisely, the format requirements are:
                * Initial preamble lines start with '##'
                * Columns (after preamble) are tab-separated.
                * The final columns of the file must contain genotype data, in the form of a (VCF-like) format column 
                  (e.g. GT:AD:DP:GQ:PL) followed by one column per sample.
        T1_thresh (int): Threshold for trio test 1. To pass, the percentage of "AA + BB" variants 
            resulting in AB must be at least T1_thresh. (Default = 90)
        T2_thresh (int): Threshold for trio test 2. To pass, the percentage of "BB + BB" variants 
            resulting in BB must be at least T2_thresh. (Default = 95)
        MZ_thresh (int): Threshold for MZ score. To pass, the percentage of variants with IBS=2 among 
            those where neither is AA, must be at least MZ_thresh. (Default = 95)
        PO_thresh (int): Threshold for parent-offspring. To pass, the percentage of variants with IBS>0 among 
            those where either is BB, must be at least PO_thresh. (Default = 99)
        MALE_thresh (int): If heterozygosity (in percent) on X (minus 
            pseudoautosomal regions) is lower than this, the sample is set to 'male'. (Default = 5%)
        FEMALE_thresh (int): If heterozygosity (in percent) on X (minus 
            pseudoautosomal regions) is higher than this, the sample is set to 'female'. (Default = 25%)
        quiet (bool): If True, print only conclusion. (Default = False)
        nogender (bool): If True, skip gender estimation. (Default = False)
        nopairwise (bool): If True, skip pairwise analysis. (Default = False)
        notrio (bool): If True, skip trio analysis. (Default = False)
        percentiles (list): List of (integral) percentile ranks to be used in filtering. (Default = [10,30,50])
    """
    
    log = None if prefix is None or prefix in ("", "STDOUT") else open('%s.log'%prefix, 'w')
    callTrios = 'YES' if not notrio else 'NO'
    callPairs = 'YES' if not nopairwise else 'NO'
    callGenders = 'YES' if not nogender else 'NO'
    
    form = '{:>40} : {}'
    _writeout('VCFped %s\n'%__version__, log, quiet)
    _writeout('Options in effect:', log, quiet)
    _writeout(form.format('Input file', file), log, quiet)
    _writeout(form.format('Output file prefix', prefix), log, quiet)
    _writeout(form.format('Analysis:', 'Genders [%s], pairwise [%s], trios [%s]' %(callGenders, callPairs, callTrios)), log, quiet)
    _writeout(form.format('Include in output:', 'All pairs/triples' if reportall else 'Only inferred pairs/triples'), log, quiet)
    _writeout(form.format('Trio test thresholds', 'Test 1 (AA+BB=AB) [%d%%], test 2 (BB+BB=BB) [%d%%]'%(T1_thresh,T2_thresh)), log, quiet)
    _writeout(form.format('Gender thresholds (X heterozygosity)', 'Male [<%d%%], female [>%d%%]' %(MALE_thresh, FEMALE_thresh)), log, quiet)
    _writeout(form.format('Sample if approx. count exceeds', exactmax), log, quiet)
    _writeout(form.format('If sampling: Minimum # variants', samplesize), log, quiet)
    _writeout(form.format('If sampling: Minimum with AA and BB', samplesizeAABB), log, quiet)
    _writeout(form.format('Filtering variables', ' '.join(variables)), log, quiet)
    _writeout(form.format('Filtering percentile ranks', ' '.join(map(str, percentiles))), log, quiet)
    st = time.time()
    
    # inspect file and guess various format details.
    info = inspectFile(file, log=log, quiet=quiet)
    nSamples = info['nSamples']
    samples = info['samples']
    linecount = info['linecount']
    chromcol = info['chromcol']
    firstSampleCol = info['formatcol'] + 1
    
    if linecount > exactmax:
        AUTOSOM, XCHR = sample_variants(file, maxsize=samplesize*2, AABBsize=samplesizeAABB, Xsize=100, info=info, log=log, quiet=quiet)
    else:
        AUTOSOM, XCHR = all_variants(file, info=info, log=log, quiet=quiet)
    
    qualdist = qualityDistrib(AUTOSOM, info, variables) # fields = ['DP']
    qualityPercentile(qualdist, percentiles, log=log, quiet=quiet)

    autosom_filt = AUTOSOM.copy()
    x_filt = XCHR.copy()
    triples = list(itertools.combinations(range(nSamples), 3))
    pairs = list(itertools.combinations(range(nSamples), 2))
     
    TRIORES, PAIRRES, GENDERRES = [], [], []
    for p in percentiles:
        qualperc = qualityPercentile(qualdist, p, log=None, quiet=True)
        filterargs = {key+'min':val for key,val in qualperc.iteritems()}
        filterObject = Filter(formatinfo=info, PASS=True, **filterargs)
        filterFUN = filterObject.filter()
        
        autosom_filt = {k:filter(filterFUN, v) for k,v in autosom_filt.iteritems()}
        x_filt = {k:filter(filterFUN, v) for k,v in x_filt.iteritems()}
        
        N = sum(len(v) for v in autosom_filt.values()) 
        
        if N == 0: break #raise RuntimeError("No variants surviving filter.")

        #### GENDERS ###
        if not nogender:
            genders = inferGenders(x_filt, info, percentile=p, FEMALE_thresh=FEMALE_thresh, MALE_thresh=MALE_thresh)
            GENDERRES.extend(genders)
          
        #### PAIRWISE ###
        if not nopairwise:
            for pair in pairs:
                testres = checkPair(autosom_filt, pair, percentile=p, MZ_thresh=MZ_thresh, PO_thresh=PO_thresh)
                PAIRRES.append(testres)
                    
        #### TRIOS ###
        if not notrio:
            for triple in triples:
                testres = checkTriple(autosom_filt, triple, percentile=p, T1_thresh=T1_thresh, T2_thresh=T2_thresh)
                TRIORES.extend(testres)
    
    _writeout('\nAll analysis finished.\nTime used: %.2f seconds' % (time.time()-st), log, quiet)
     
    TEMP = "\n===={heading}====\n{table}\n\nMost confident {analysis} calls written to {outfile}"     
    
    if not nogender:
        GENDERRES.sort(key=lambda x: (x['sample'], x['percentile']))
        outfile = "%s.gender" % prefix
        
        # write all results to log (pretty-printed)
        complete = tabulate(gender_table(GENDERRES), tablefmt="plain") 
        _writeout(TEMP.format(heading="GENDERS", table=complete, analysis="gender", outfile=outfile), log, quiet=True)
        
        # most confident calls
        best = bestGenders(GENDERRES)
        
        # write most confident calls to screen, if not quiet (pretty-printed)
        best_pretty = tabulate(best, tablefmt="plain") 
        _writeout(TEMP.format(heading="GENDERS", table=best_pretty, analysis="gender", outfile=outfile), log=None, quiet=quiet)
        
        # write most confident calls to separate file (tab separated!)
        best_tsv = '\n'.join('\t'.join(map(str, line)) for line in best)
        with open(outfile, 'w') as out:
            out.write(best_tsv)   
        
    if not nopairwise:    
        PAIRRES.sort(key=lambda x: (x['pair'], x['percentile']))
        outfile = "%s.pair" % prefix
        # write all results to log (pretty-printed)
        complete = tabulate(pairwise_table(PAIRRES), tablefmt="plain") if PAIRRES else "Not enough individuals in the file"
        _writeout(TEMP.format(heading="PAIRWISE RELATIONS", table=complete, analysis="pairwise", outfile=outfile), log, quiet=True)
        
        # most confident calls
        best = bestPairs(PAIRRES, reportall=reportall)
        
        # write most confident calls to screen, if not quiet (pretty-printed)
        best_pretty = tabulate(best, tablefmt="plain") if PAIRRES else "Not enough individuals in the file"
        _writeout(TEMP.format(heading="PAIRWISE RELATIONS", table=best_pretty, analysis="pairwise", outfile=outfile), log=None, quiet=quiet)
        
        # write most confident calls to separate file (tab separated!)
        best_tsv = '\n'.join('\t'.join(map(str, line)) for line in best)
        with open(outfile, 'w') as out:
            out.write(best_tsv)
                   
    if not notrio:
        TRIORES.sort(key=lambda x: (x['triple'], x['pivot'], x['percentile']))
        outfile = "%s.trio" % prefix
        
        # write all results to log (pretty-printed)
        complete = tabulate(trio_table(TRIORES), tablefmt="plain") if TRIORES else "Not enough individuals in the file"
        _writeout(TEMP.format(heading="TRIO RELATIONS", table=complete, analysis="trio", outfile=outfile), log, quiet=True)
        
        # most confident calls
        best = bestTrios(TRIORES, reportall=reportall)
        
        # write most confident calls to screen, if not quiet (pretty-printed)
        best_pretty = tabulate(best, tablefmt="plain") if TRIORES else "Not enough individuals in the file"
        _writeout(TEMP.format(heading="TRIO RELATIONS", table=best_pretty, analysis="trio", outfile=outfile), log=None, quiet=quiet)
        
        # write most confident calls to separate file (tab separated!)
        best_tsv = '\n'.join('\t'.join(map(str, line)) for line in best)
        with open(outfile, 'w') as out:
            out.write(best_tsv)        
    
    if log:
        log.close()
    
def main():    
    parser = argparse.ArgumentParser()
    parser.add_argument("file", help="path to VCF file")
    parser.add_argument("--quiet", help="do not print information on the screen", action='store_true')
    parser.add_argument('--version', action='version', version='VCFped %s'%__version__)
    parser.add_argument("--no-gender", help="skip gender analysis", dest='nogender', action='store_true')
    parser.add_argument("--no-pairwise", help="skip pairwise analysis", dest='nopairwise', action='store_true')
    parser.add_argument("--no-trio", help="skip trio analysis", dest='notrio', action='store_true')
    parser.add_argument("--all", help="show results for all pairs/triples (not only the inferred)", dest="reportall", action='store_true')
    parser.add_argument("-o", help="prefix for output files", dest="prefix")
    parser.add_argument("-v", dest="variables", nargs='+', help="quality variables to be used for filtering", default=['QUAL', 'DP', 'GQ', 'AD'])
    parser.add_argument("-p", dest="percentiles", nargs='+', type=int, help="filtering percentile ranks", default=[10,30,50])
    parser.add_argument("-e", dest='exactmax', type=int, help="if approx. line count exceeds this, apply random sampling", default=150000)
    parser.add_argument("-s", dest='samplesize', type=int, help="sample at least this many variant lines (if sampling)", default=10000)
    parser.add_argument("-d", dest='samplesizeAABB', type=int, help="sample at least this many lines where both 0/0 and 1/1 occur as genotypes (if sampling)", default=500)
    parser.add_argument("-t1", dest="T1_thresh", type=int, help="threshold (%%) for T1 score (AA + BB = AB)", default=90)
    parser.add_argument("-t2", dest="T2_thresh", type=int, help="threshold (%%) for T2 score (BB + BB = BB)", default=95)
    parser.add_argument("-mz", dest="MZ_thresh", type=int, help="threshold (%%) for MZ score (IBS=2 | neither is AA)", default=95)
    parser.add_argument("-po", dest="PO_thresh", type=int, help="threshold (%%) for PO score (IBS>0 | either is BB)", default=99)
    parser.add_argument("-female", dest="FEMALE_thresh", type=int, help="lower limit (%%) for female heterozygosity on X", default=25)
    parser.add_argument("-male", dest="MALE_thresh", type=int, help="upper limit (%%) for male heterozygosity on X", default=5)
    
    args = parser.parse_args(sys.argv[1:])  
    for pp in args.percentiles:
        if not 0 <= pp < 100: parser.error("argument -p: invalid percentile: %d" %pp)
    thresh_error = "argument -%s/--%s: invalid threshold (must be in [0, 100]): %d"
    if not 0<= args.T1_thresh <=100: parser.error(thresh_error % ('t1','T1_thresh', args.T1_thresh))
    if not 0<= args.T2_thresh <=100: parser.error(thresh_error % ('t2','T2_thresh', args.T2_thresh))
    if not 0<= args.MZ_thresh <=100: parser.error(thresh_error % ('mz','MZ_thresh', args.MZ_thresh))
    if not 0<= args.PO_thresh <=100: parser.error(thresh_error % ('po','PO_thresh', args.PO_thresh))
    if not 0<= args.FEMALE_thresh <=100: parser.error(thresh_error % ('female','FEMALE_thresh', args.FEMALE_thresh))
    if not 0<= args.MALE_thresh <=100: parser.error(thresh_error % ('male','MALE_thresh', args.MALE_thresh))
    if args.MALE_thresh > args.FEMALE_thresh: parser.error("argument -male/--MALE_thresh: invalid threshold (must be below female threshold which is %d): %d" %(args.FEMALE_thresh, args.MALE_thresh))
    
    if args.prefix is None: args.prefix = os.path.splitext(os.path.basename(args.file))[0]
    vcfped(**vars(args))
    
if __name__ == '__main__':
    main()
    
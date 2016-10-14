"""Module for identifying trios and parent-offspring relationships in a multisample VCF file.

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


import collections
import operator
import time
import itertools
import random
import sys

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
    
# Mapping VCF genotypes to {AA, AB, BB}
GT2AB = {'0/0':'AA', '0|0':'AA', '0/1':'AB', '0|1':'AB', '1|0':'AB', '1/1':'BB', '1|1':'BB'}

def inspectFile(filename, verbose=1):
    vfc = 0
    with open(filename, 'r') as f:
        
        # Identify VCF format if stated in first line
        line = f.readline().rstrip()
        if 'fileformat=' in line:
            vcf = line[line.index('=')+1:]

        # Skip preamble lines
        while line.startswith("##"): 
            line = f.readline()
        
        headers = line.rstrip().split('\t')
        startpos = f.tell()
        firstvar = f.readline().rstrip().split('\t')
        
        # Identify chromosome column
        if headers[0] == "#CHROM" or "chrom" in headers[0].lower():
            chromcol = 0
        else:
            vcf = 0
            chromcol = next(i for i,h in enumerate(headers) if h.lower() in ("chr", "chrom", "chromosome"))

        # Identify QUAL column
        try:
            QUALind = next(i for i,field in enumerate(headers) if field.lower() in ['qual', 'vcf_qual'])
            float(firstvar[QUALind])
        except:
            QUALind = None
            
        # Identify format column
        formatcol = next(i for i,field in enumerate(firstvar) if field.startswith("GT:") or field=="GT")
        if formatcol != 8:
            vcf = 0

        chromheader = headers[chromcol]
        formatheader = headers[formatcol]
        format = firstvar[formatcol]
        formatfields = format.split(':')
        ADind = formatfields.index('AD') if 'AD' in formatfields else None
        DPind = formatfields.index('DP') if 'DP' in formatfields else None
        GQind = formatfields.index('GQ') if 'GQ' in formatfields else None
        
        samples = headers[(formatcol + 1):]
        nSamples = len(samples)
        if not vcf:
            vcf = "Non-standard vcf"
            
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
        
    return dict(filename=filename, fileformat=vcf, chromheader=chromheader, chromcol=chromcol, formatheader=formatheader, formatcol=formatcol, QUALind=QUALind, formatfields=formatfields, DPind=DPind, GQind=GQind, ADind=ADind, 
                counttype=counttype, linecount=linecount, samples=samples, nSamples=nSamples, startpos=startpos, endpos=endpos)
            
def getGTtuple(linesplit, firstSampleCol):
    try:
        return tuple(GT2AB[gt[:3]] for gt in linesplit[firstSampleCol:])
    except KeyError:
        return None

        
class Filter(object):
    def __init__(self, formatinfo, PASS=None, DPmin=None, GQmin=None, QUALmin=None, ADmin=None):
        self.nSamples=formatinfo['nSamples']
        self.QUALind=formatinfo['QUALind']
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
        nSamples, QUALind, ADind, DPind, GQind = self.nSamples, self.QUALind, self.ADind, self.DPind, self.GQind
        PASS, QUALmin, DPmin, GQmin, ADmin = self.PASS, self.QUALmin, self.DPmin, self.GQmin, self.ADmin
        PASSind = QUALind + 1 if QUALind is not None else None
        
        PASStest = PASS and PASSind is not None
        QUALtest = QUALmin>0 and QUALind is not None
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
                if QUALtest and float(linesplit[QUALind]) < QUALmin: return 0
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
    QUALind=formatinfo['QUALind'] 
    firstSampleCol = formatinfo['formatcol'] + 1
    useAD = 'AD' in fields and ADind is not None
    useGQ = 'GQ' in fields and GQind is not None
    useDP = 'DP' in fields and DPind is not None
    useQUAL = 'QUAL' in fields and QUALind is not None
    distrib_data = collections.defaultdict(list)
    
    def _ALTratio(ADfield):
        ad0, ad1 = map(float, ADfield.split(','))
        rat = ad1/(ad0+ad1)
        return round(min(rat, 1-rat) * 100, 1)
    
    hets = ('0/1', '0|1', '1/0', '1|0')
    for s in variantlines:
        try:
            if useQUAL: distrib_data['QUAL'].append(float(s[QUALind]))
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
    
def qualityPercentile(distrib_data, p, fout=None):
    if isinstance(p, list):
        res = {key : [data[int(pp/100.0*len(data))] for pp in p] for key,data in distrib_data.iteritems()}
    else:
        res = {key : data[int(p/100.0*len(data))] for key,data in distrib_data.iteritems()}
        
    if fout:
        print >>fout, "====SCORE PERCENTILES (%s)====" % p #", ".join(str(i) for i in p)
        for key in ['QUAL', 'DP', 'GQ', 'AD']:
            if key in res: 
                print >>fout, "  %s: %s" % (key, res[key])
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
    
    triple_txt = '%d,%d,%d' % tuple(p+1 for p in triple)
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
                verdict = 'RegularTrio' if BB_BB_BBp >= threshold2 else 'OtherTrio'
            else:
                verdict = 'NotTrio'
        else:
            verdict = 'na'
        
        rr = (triple_txt, pivot+1, percentile, SURV, tot_AA_BB, AA_BB_AB, tot_BB_BB, BB_BB_BB, '%.1f%%'%AA_BB_ABp, '%.1f%%'%BB_BB_BBp, verdict)
        testresults.append(rr)
    
    return testresults
            
def checkPair(AUTOSOMfilt, pair, percentile):
    pair_txt = '%d,%d' % tuple(p+1 for p in pair)
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
    
    return (pair_txt, percentile, SURV, MZ, '%.1f%%'%MZp, tot_BB, AA_BB, '%.1f%%'%AA_BBp, verdict)

  
def all_variants(file, info):
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
    
    AUTant = sum(len(v) for v in AUTOSOM.values())
    Xant = sum(len(v) for v in XCHR.values())
    text = 'Using all PASS variants:\n %d autosomal (including %d with AA+BB)\n %d X-linked.'%(AUTant, AABB, Xant)
    
    return AUTOSOM, XCHR, text
    
def sample_variants(file, maxsize, AABBsize, Xsize, info, chunks=4):
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
    text = 'Sampled variants:\n %d autosomal variants (including %d with AA+BB)\n %d on X.'%(AUTant, AABB, Xant)
        
    return AUTOSOM, XCHR, text

def inferGenders(XCHRfilt, formatinfo, threshMale, percentile):
    nSamples = formatinfo['nSamples']
    Xdistribs = [collections.Counter() for _ in range(nSamples)]
    Xtot = 0
    for gtcomb, lines in XCHRfilt.iteritems():
        nvars = len(lines)
        Xtot += nvars
        for i, gt in enumerate(gtcomb):
            Xdistribs[i][gt] += nvars 
    
    def _dist2gender(Xdist):
        denom = Xdist['AB']+Xdist['BB']
        Xhetp = round(100.0 * Xdist['AB']/denom, 1) if denom > 0 else 'na'
        if denom >=25:
            if Xhetp < threshMale:
                gender = 'male' 
            elif Xhetp < 2*threshMale:
                gender = '?'
            else:
                gender = 'female' 
        else: gender = 'na'
        return Xdist['AA'], Xdist['AB'], Xdist['BB'], Xhetp, gender
    
    res = [(i+1, percentile, Xtot) + _dist2gender(Xdist) for i, Xdist in enumerate(Xdistribs)]
    return res
    
def interpretTrioResult(triodata, report_all=False):
    best = []
    for k,v in itertools.groupby(triodata, key=operator.itemgetter(0)):
        use = [r for r in v if r[10] != 'na']
        if not use:
            if report_all:
                use = v
            else:
                continue
        b = sorted(use, key=lambda x: (-round(100.0*x[5]/x[4], 1), -x[2]))
        best.append(b[0])
    return best


def vcfped(file, threshTrio1=90, threshTrio2=95, threshMale=5, exactmax=100000, samplesize=10000, samplesizeAABB=1000, percentiles=[10,30,50], nogender=False, nopairwise=False, notrio=False, out=sys.stdout, verbose=True):
    """Determine the pedigree structure (order and genders) of a trio.
    
    The input file should contain jointly called variants from three individuals.
    
    Args:
        file (str): The path to a variant file with exactly 3 samples. The file format must be VCF or
            VCF-like. More precisely, the format requirements are:
                * Initial preamble lines start with '##'
                * Columns (after preamble) are tab-separated.
                * The final columns of the file must contain genotype data, in the form of a (VCF-like) format column 
                  (e.g. GT:AD:DP:GQ:PL) followed by one column per sample.
        threshMale (int): Threshold for gender estimation. If heterozygosity (in percent) on X (minus 
            pseudoautosomal regions) is lower than this, the sample is set to 'male'. (Default = 10%)
        thresh1 (int): Threshold for test 1. To pass test1, the percentage of "AA + BB" variants 
            resulting in AB must be at least thresh1. (Default = 95%)
        verbose (bool): If True, print conclusion and evidence to the screen. (Default = True)
        
    Returns:
        A tuple with the following entries in order:
            - (string) Trio type, either "regular trio", "inverted trio", "generational trio" or "not trio".
            - (int) The index (0, 1 or 2) of the proband, relative to the order in the input file. 
                The proband is defined as the child in a regular trio, the parent in an inverted trio, 
                or the middle individual in a generational trio. If not a trio, the proband is set to -1.
            - (int) Total number of variants in the input file
            - (int) The number of autosomal variants passing all filters. Note that in addition to the 
                user-specified filters, multiallelic variants and variants with missing data are always ignored.
            - (int) The number of X-linked variants passing all filters. 
            - (list) The heterozygosity on X (-PAR) for each sample, in the same order as in the input file.
            - (dict) Test results, given as entries key:value where 'key' is an ordering of the input samples,
                and value is a list with 4 elements:
                    * (string) A description of the test, on the form GT1 + GT2 = GT3.
                    * (float) Test score: Percent of autosomal variants with GT3 in the third sample, 
                        among those with GT1, GT2 in the first two.
                    * (int) The number of autosomal variants with the combination GT1 + GT2 in the first two samples.
                    * (float) The previous count in percent of the total number of autosomal variants (after filtering). 
                If a trio is recognized, only the test results for the correct order is given, otherwise
                results for all orderings are included.
    """
    
    st = time.time()
    fout = sys.stdout if out in ("", "STDOUT") else open(out, 'w')
    
    # inspect file and guess various format details.
    info = inspectFile(file, verbose=verbose)
    nSamples = info['nSamples']
    samples = info['samples']
    linecount = info['linecount']
    chromcol = info['chromcol']
    firstSampleCol = info['formatcol'] + 1
    
    if verbose:
        info['formatheads_joined'] = ', '.join(info['formatfields'])
        info['samples_joined'] = '\n'.join('  %d: %s' %(i+1, s) for i,s in enumerate(info['samples']))
        print >> fout, '''\
====FILE INFO====
File name: {filename}
File format: {fileformat}
{counttype} variant count: {linecount}
Chromosome column (index): {chromheader} ({chromcol})
Format column (index): {formatheader} ({formatcol})
Format fields: {formatheads_joined}
Samples:
{samples_joined}'''.format(**info)
        
        
    if linecount > exactmax:
        AUTOSOM, XCHR, text = sample_variants(file, maxsize=samplesize*2, AABBsize=samplesizeAABB, Xsize=100, info=info)
    else:
        AUTOSOM, XCHR, text = all_variants(file, info=info)
        
    if verbose:
        print >> fout, '\n====VARIANT SAMPLE (BEFORE FILTERING)===='
        print >> fout, text
        print >> fout, 'Time used so far:', time.time()-st, '\n'
    
    qualdist = qualityDistrib(AUTOSOM, info) # fields = ['DP']
    qualityPercentile(qualdist, percentiles, fout)

    autosom_filt = AUTOSOM.copy()
    x_filt = XCHR.copy()
    triples = list(itertools.combinations(range(nSamples), 3))
    pairs = list(itertools.combinations(range(nSamples), 2))
     
    TRIORES, PAIRRES, GENDERRES = [], [], []
    for p in percentiles:
        qualperc = qualityPercentile(qualdist, p)
        filterargs = {key+'min':val for key,val in qualperc.iteritems()}
        filterObject = Filter(formatinfo=info, PASS=True, **filterargs)
        filterFUN = filterObject.filter()
        
        autosom_filt = {k:filter(filterFUN, v) for k,v in autosom_filt.iteritems()}
        x_filt = {k:filter(filterFUN, v) for k,v in x_filt.iteritems()}
        
        N = sum(len(v) for v in autosom_filt.values()) 
        
        if N == 0: break #raise RuntimeError("No variants surviving filter.")

        #### GENDERS ###
        if not nogender:
            genders = inferGenders(x_filt, info, threshMale, percentile=p)
            GENDERRES.extend(genders)
          
        #### PAIRWISE ###
        if not nopairwise:
            for pair in pairs:
                testres = checkPair(autosom_filt, pair, percentile=p)
                PAIRRES.append(testres)
                    
        #### TRIOS ###
        if not notrio:
            for triple in triples:
                testres = checkTriple(autosom_filt, triple, percentile=p, threshold1=threshTrio1, threshold2=threshTrio2)
                TRIORES.extend(testres)
    
    if verbose:
        print >> fout, '\nAll analysis finished.\nTime used: %.2f seconds' % (time.time()-st)
        
    if not nogender:
        GENDERRES.sort(key=operator.itemgetter(0,1))
        print >> fout, '\n====GENDERS===='
        print >> fout, '\t'.join(['Sample', 'Perc', 'Xvars', 'AA', 'AB', 'BB', 'Xhet', 'Verdict'])
        for r in GENDERRES:
            print >> fout, '\t'.join(map(str,r))
    
    if not nopairwise:    
        PAIRRES.sort(key=operator.itemgetter(0,1))
        print >> fout, '\n====PAIRWISE RELATIONS===='
        if not pairs:
            print >> fout, "Not enough individuals in the file (only %d present)." % nSamples
        else:
            print >> fout, '\t'.join(['Pair', 'Perc', 'Autos', '==', 'MZtest', 'anyBB', '=AA', 'POtest', 'Verdict'])
            for r in PAIRRES:
                print >> fout, '\t'.join(map(str,r))
    
    if not notrio:
        TRIORES.sort(key=operator.itemgetter(0,1,2))
        print >> fout, '\n====TRIO RELATIONS===='
        if not triples:
            print >> fout, "Not enough individuals in the file (only %d present)." % nSamples
        else:
            print >> fout, '\t'.join(['Triple', 'Pivot', 'Perc', 'Autos', 'AA+BB', '=AB', 'BB+BB', '=BB', 'T1', 'T2', 'Verdict'])
            for r in TRIORES:
                print >> fout, '\t'.join(map(str,r))
        
        best = interpretTrioResult(TRIORES, report_all=True)
        if best:
            print >> fout, '\nMost confident trio calls:'
            print >> fout, '\t'.join(['Triple', 'Pivot', 'Perc', 'Autos', 'AA+BB', '=AB', 'BB+BB', '=BB', 'T1', 'T2', 'Verdict'])
            for r in best:
                print >> fout, '\t'.join(map(str,r))
        else:
            print >> fout, "No trios detected"
            
    if fout is not sys.stdout:
        fout.close()
    
    
    
if __name__ == '__main__':
    import sys
    import argparse
    
    # For testing purposes: Defining argument parsing as a function
    def parse_args(args):
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
        return parser.parse_args(args)  
      
    args = parse_args(sys.argv[1:])
    vcfped(**vars(args))
    
    
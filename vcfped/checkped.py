def readped(pedfile):
    # read file in pre-makeped format
    # output: dictionary with samplenames of father, mother, child, and child gender
    with open(pedfile, 'r') as inf:
        lines = [line.strip() for line in inf if line.strip()]
    
    # remove header row if present
    if all(y in lines[0].upper() for y in ['FID', 'MID', 'SEX']):
        lines = lines[1:]
    
    if len(lines) != 3:
        raise RuntimeError("Input pedigree file does not have exactly 3 lines. Pedigree checking is only implemented for trios.")
    
    # separate on tab or space
    sep = '\t' if '\t' in lines[0] else ' '
    ped = [line.split(sep) for line in lines]
    
    try:
        father = next(ped[i][1] for i in range(3) if ped[i][2]==ped[i][3]=='0' and ped[i][4]=='1')
        mother = next(ped[i][1] for i in range(3) if ped[i][2]==ped[i][3]=='0' and ped[i][4]=='2')
        child_index = next(i for i in range(3) if ped[i][2] == father and ped[i][3] == mother)
    except StopIteration:
        raise RuntimeError("Pedfile does not describe a trio.")
    
    child = ped[child_index][1]
    child_gender_int = int(ped[child_index][4])
    child_gender = ('', 'Male', 'Female')[child_gender_int]
    
    res = dict(father=father, mother=mother, child=child, child_gender=child_gender, sampleset=[father, mother, child])
    return res

def checkped(pedfile, bestG, bestP, bestT, info):
    ped = readped(pedfile)
    samples = info['samples']
    if len(samples) != 3: 
        raise RuntimeError("File does not contain exactly 3 samples.")
    if set(samples) != set(ped['sampleset']): 
        raise RuntimeError("Sample names don't match pedfile ID's.\nVariant file: %s\nPedfile: %s" %(samples, ped['sampleset']))
    
    bestTdict = dict(zip(*bestT))
    pivot = bestTdict['Pivot']
    if not bestTdict['Verdict'] == "Regular trio":
        raise RuntimeError("This is not recognized as a regular trio. Check output files for details.")
    if not samples[pivot-1] == ped['child']:
        raise RuntimeError("Regular trio, but wrong labels. \nPedfile: child = %s\nVariant file: child = %s" %(ped['child'], samples[pivot]))
    
    # Sample names
    fa, mo, ch = ped['father'], ped['mother'], ped['child']
    
    # genders predicted by vcfped
    chi, fai, moi = samples.index(ch), samples.index(fa), samples.index(mo)
    cgend, fgend, mgend = bestG[chi+1][7], bestG[fai+1][7], bestG[moi+1][7]
    
    # Check gender of child
    if not cgend == ped['child_gender']:
        raise RuntimeError("Regular trio, but wrong child gender.\nPedfile: child = %s\nVariant file: child = %s" %(ped['child_gender'], cgend))
    
    if (fgend,mgend) == ('Female', 'Male'):
        raise RuntimeError("Regular trio, but parents are switched.")
    if fgend != 'Male':
        raise RuntimeError("Regular trio, but father does not appear to be male: \nPedfile: %s = %s\nVariant file: %s = %s" %(fa, 'Male', fa, fgend))
    if mgend != 'Female':
        raise RuntimeError("Regular trio, but mother does not appear to be female: \nPedfile: %s = %s\nVariant file: %s = %s" %(mo, 'Female', mo, mgend))
    
    
    
    
if __name__ == '__main__':
    print readped("testfiles/trio.ped")

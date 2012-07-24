# -*- coding: utf-8 -*-
# <nbformat>2</nbformat>

# <codecell>

num_haplo = 100
pct_var = 5
target_cov = 10
outdir = os.path.expanduser('~/code/AFS/sim-%shaplos_exp-dist_theta%0.2f_%sxCov/' % (num_haplo,pct_var/100.0,target_cov))
print outdir

# <codecell>

import re,os
from glob import glob

# <codecell>

f = 'agouti-mc1r.fa'
seqd = dict([(seqid,sequence.replace('\n','').upper()) for seqid,sequence in re.findall('>(.+?)\n([\w\n]+)',open(f).read())])

# <codecell>

[(k,len(v)) for k,v in seqd.items()]

# <codecell>

e_dist = numpy.random.exponential(scale=num_haplo/10,size=100*num_haplo)
h,b,p = hist(e_dist,bins=range(num_haplo))

# <codecell>

common_alleles = h[int(round(num_haplo*0.1))]
pad_h = h+numpy.random.normal(common_alleles,common_alleles/20,size=len(h))
invar = (100.0-pct_var)/pct_var
sfs = [sum(pad_h)*invar] + list(pad_h)
sfs = array(sfs)/sum(sfs)
csfs = cumsum(sfs)
print sum(sfs[1:])/float(sfs[0]), len(sfs)
plot(sfs)
ax = figure(1).axes[0]
ax.set_yscale('log')
#print sfs
figure(2)
plot(sfs[1:])

# <codecell>

simseq = {}
for i in range(len(csfs)):
    simseq[i] = {}
for sn,ss in seqd.items():
    for i in range(len(csfs)):
        simseq[i][sn] = []
    for c,nvar in zip(ss,[[i for i in range(len(csfs)) if n<csfs[i]][0] for n in numpy.random.random_sample(size=len(ss))]):
        if nvar == 0:
            for i in range(len(csfs)):
                simseq[i][sn].append(c)
        else:
            alt = numpy.random.permutation(list(set('ACGT')-set(c)))[0]
            this_ord = numpy.random.permutation(len(csfs))
            for i in this_ord[:nvar]:
                simseq[i][sn].append(alt)
            for i in this_ord[nvar:]:
                simseq[i][sn].append(c)

for i in range(len(csfs)):
    for sn in simseq[i].keys():
        simseq[i][sn] = ''.join(simseq[i][sn])

# <codecell>

try:
    os.makedirs(outdir)
except:
    pass
for ind,simseqd in simseq.items():
    outf = open(outdir+'hap%s_agouti-mc1r.fa' % ind,'w')
    for sn,ss in simseqd.items():
        outf.write('>%s\n%s\n' % (sn,ss))
    outf.close()

# <codecell>

num_ind = num_haplo/2
tot_seq = sum(map(len,seqd.values()))
# number of bases generated per haplo / bp per read-pair
nreads = int((tot_seq*target_cov/2.0) / 200)
print 'target coverage %s; reads per haplotype %s' % (target_cov,nreads)
try:
    os.makedirs(outdir+'sim-reads/')
except:
    pass
ind_haplos = reduce(lambda x,y:x+y, [[outdir+'sim-reads/sim%sa' % i,outdir+'sim-reads/sim%sb'% i] for i in range(1,num_ind+1)])
cmds = []
map_reference = os.path.expanduser('~/code/AFS/agouti-mc1r.masked.fasta')
for read_prefix, source_reference in zip(ind_haplos,numpy.random.permutation(glob(outdir+'hap*_agouti-mc1r.fa'))):
    cmds.append('simulate_and_map.py %s %s %s %s %s %s.vcf' % (read_prefix,os.path.basename(read_prefix)[:-1],source_reference,map_reference,nreads,read_prefix))
len(cmds)

# <codecell>

for cmd in cmds:
    !$cmd

# <codecell>

bams = glob(outdir+'sim-reads/*.rg.bam')
Istr = ' -I '.join(bams)
vcf = outdir+'all.vcf'
java_ram = '1g'
jar = os.path.expanduser('~/jars/gatk/GenomeAnalysisTK.jar')
gatkcmd = 'java -Xmx%s -jar %s -T UnifiedGenotyper -I %s -R %s -o %s -maxAlleles 1 -out_mode EMIT_ALL_SITES' % (java_ram, jar, Istr, map_reference, vcf)
!$gatkcmd

# <codecell>

afli = map(float,re.findall('AF=([\d\.]+)',open(vcf).read()))
len(afli)

# <codecell>

hist(array(afli)*num_haplo,bins=range(1,num_haplo+1))
figure(2)
h,b = histogram(array(afli)*num_haplo,bins=range(0,num_haplo+1))
plot(b[1:-1],(array(h,dtype=float)/sum(h))[1:],'--r')
plot(b[1:-1],sfs[1:],'g')

# <codecell>

h,b = histogram(array(afli)*num_haplo,bins=range(0,num_haplo+1))
plot(array(h,dtype=float)/sum(h),'--r')
plot(sfs,'g')

ax = figure(1).axes[0]
ax.set_yscale('log')

# <codecell>

rare = round(num_haplo*0.05)
common = round(num_haplo*0.1)
print rare,common
print sum(sfs[1:rare+1])/sum(sfs[common:])
print sum((array(h,dtype=float)/sum(h))[1:rare+1])/sum((array(h,dtype=float)/sum(h))[common:])

# <codecell>



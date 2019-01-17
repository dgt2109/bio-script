
import sys
import statsmodels.sandbox.stats.multicomp as mc

with open(sys.argv[1],"r") as f:
    lines=f.readlines()
header=int(sys.argv[2])
if header == 1:
    plist=[float(x) for x in lines[1:]]
elif header == 0:
    plist=[float(x) for x in lines[0:]]
    
res = mc.multipletests(plist,alpha=0.05,method='fdr_bh',is_sorted=False,returnsorted=False)

with open(sys.argv[3],'a') as g:
    g.write("Sig\tFDR_q\n")
    for i in range(0,len(res[0])):
        g.write(str(res[0][i]) + '\t'+ str(res[1][i]) + '\n')

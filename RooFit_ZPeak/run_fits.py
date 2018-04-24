import itertools
from multiprocessing import Process,Pool
from os import system

def fit(args):#sig="",bkg="",isData=0,qmult=(0,1)):
    sig=args[0]
    bkg=args[1]
    isData=args[2]
    qmult=args[3]
    command = 'root -l -n -b -q "fit_prel.C({2},{3},{4},model::{0},model::{1})" > results_{0}_{1}_qmult_{3}_{4}_isData_{2}.log'.format(sig,bkg,isData,qmult[0],qmult[1])
    print '[STARTING]',command
    system(command)
    print '[DONE]',command

sigModels=["kBWCrystBall","kBWGauss"]
bkgModels=["kPoly1","kPoly2","kPoly3","kExpPoly1"]
qmult = [(0,1),(2,9999)]
isData = [0,1]

permutations = itertools.product(sigModels,bkgModels,isData,qmult)
pool = Pool(processes=4)
pool.map(fit,permutations)

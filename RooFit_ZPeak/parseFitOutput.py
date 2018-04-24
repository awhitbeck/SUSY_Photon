import numpy as np
import itertools
import pandas as pd

sigModels=["kBWCrystBall","kBWGauss"]
bkgModels=["kPoly1","kPoly2","kPoly3","kExpPoly1"]
qmult = [(0,1),(2,9999)]
isData = [0,1]

permutations = itertools.product(sigModels,bkgModels,isData,qmult)

results={'sig_model':[],
         'bkg_model':[],
         'qmult_low':[],
         'qmult_high':[],
         'isData':[],
         'fake-rate':[],
         'fake-rate-error':[],
         # 'fake-rate-data':[],
         # 'fake-rate-error-data':[],
         'chi-square-pho':[],
         'ndof-pho':[],
         'chi-square-ele':[],
         'ndof-ele':[],
         # 'chi-square-pho-data':[],
         # 'ndof-pho-data':[],
         # 'chi-square-ele-data':[],
         # 'ndof-ele-data':[]
         }

for p in permutations:
    inputFile = open('results_{0}_{1}_qmult_{3}_{4}_isData_{2}.log'.format(p[0],p[1],p[2],p[3][0],p[3][1]),'r')

    results['sig_model'].append(p[0])
    results['bkg_model'].append(p[1])
    results['qmult_low'].append(p[3][0])
    results['qmult_high'].append(p[3][1])
    results['isData'].append(p[2])
    
    suffix=''

    for line in inputFile : 
        if "fake rate" in line : 
            #fake rate: 0.0170182 +/- 1.50343e-05
            results['fake-rate'+suffix].append(float(line.split(' ')[2]))
            results['fake-rate-error'+suffix].append(float(line.split(' ')[4]))
        if "my chi-squared (pho)" in line : 
            #my chi-squared (pho): 64.9886 ndof: 99 and model params: 7
            results['chi-square-pho'+suffix].append(float(line.split(' ')[3]))
            results['ndof-pho'+suffix].append(int(line.split(' ')[5]) - int(line.split(' ')[9]))
        if "my chi-squared (ele)" in line : 
            results['chi-square-ele'+suffix].append(float(line.split(' ')[3]))
            results['ndof-ele'+suffix].append(int(line.split(' ')[5]) - int(line.split(' ')[9]))
            
for key in results : 
    print key,len(results[key])

df = pd.DataFrame.from_dict(results)
#print df.head(100)

#df['scale-factor'] = df['fake-rate-data']/df['fake-rate']
#df['scale-factor-error'] = sqrt(df['fake-rate-error']*df['fake-rate-error']/df['fake-rate']/df['fake-rate']+df['fake-rate-error-data']*df['fake-rate-error-data']/df['fake-rate-data']/df['fake-rate-data'])*df['fake-rate-data']/df['fake-rate']
#print df[(df['qmult_low']==6)&(df['isData']==0)][['fake-rate','fake-rate-error','sig_model','bkg_model']]#,'scale-factor','scale-factor-error']]

df_MC = df[df['isData']==0]
df_data = df[df['isData']==1]
df_MC = df_MC.sort(['sig_model','bkg_model','qmult_low'])
df_data = df_data.sort(['sig_model','bkg_model','qmult_low'])
df_MC.reset_index(inplace=True)
df_data.reset_index(inplace=True)

df_merge = df_MC.merge(df_data,how='outer',on=['sig_model','bkg_model','qmult_low','qmult_high'],suffixes=['_MC','_data'])

df_merge['scale-factor'] = df_merge['fake-rate_data']/df_merge['fake-rate_MC']
df_merge['scale-factor-error'] = np.sqrt(df_merge['fake-rate-error_MC']*df_merge['fake-rate-error_MC']/df_merge['fake-rate_MC']/df_merge['fake-rate_MC']+df_merge['fake-rate-error_data']*df_merge['fake-rate-error_data']/df_merge['fake-rate_data']/df_merge['fake-rate_data'])*df_merge['scale-factor']

print df_merge[['sig_model','bkg_model','qmult_low','qmult_high','fake-rate_MC','fake-rate-error_MC','fake-rate_data','fake-rate-error_data','scale-factor','scale-factor-error']]

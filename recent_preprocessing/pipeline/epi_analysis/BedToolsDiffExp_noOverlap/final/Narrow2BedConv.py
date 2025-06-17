import pandas as pd
import os

#path= 'narrowPeak/10pct/'
#outdir= 'bed5/10pct/'
path= 'narrowPeak/all/'
outdir= 'bed5/all/'
os.makedirs(outdir, exist_ok= True)
for filename in os.listdir(path):
    if not filename.startswith('.'):
        df= pd.read_csv(path+filename, sep= '\t', header= None)
        df.drop([4,5,6,7,8], axis= 1, inplace= True)
        df.to_csv(outdir+filename.split('.')[0]+'.bed', sep= '\t', index= False, header= False)

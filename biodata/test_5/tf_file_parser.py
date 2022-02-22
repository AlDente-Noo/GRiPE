#%%
import random as rnd
import pandas as pd
import numpy as np

#%%
#pd.set_option('display.max_columns', None)

tfs = "cad hkb kni gt tll bcd hb Kr".split()
# copy_numbers = [0] * len(tfs)

# for i, tf in enumerate(tfs):
# 	if tf == 'kni':
# 		copy_numbers[i] = 10
# 	if tf == 'gt':
# 		copy_numbers[i] = 1


with open('tf_2.csv') as f:
	df = pd.read_csv(f, sep=';')

df

#df.to_csv(path_or_buf='tf_1.csv', sep=';', index=False)

#%% remove old columns and add the new ones
#df = df.drop(['REPRESSOR', 'REPLENLEFT', 'REPLENRIGHT', 'PWMREPTHRESHOLD', 
#              'REPRESSIONPROBABILITY'], axis=1)

repr_rate = [0.0, ] + [1.0,]*3 + [0.0]*4
repr_size = [0, ] + [100,]*3 + [0]*4
unrepr_rate = [0.0, ] + [0.1,]*3 + [0.0]*4

#repr_rate = [0.0, ]*8
#repr_size = [0, ]*8
#unrepr_rate = [0.0, ]*8

df['REPRESSIONRATE'] = repr_rate
df['UNREPRESSIONRATE'] = unrepr_rate
df['REPRLENLEFT'] = repr_size
df['REPRLENRIGHT'] = repr_size

df[['PREBOUNDTOHIGHESTAFFINITY','TFISIMMOBILE','ISBIASEDRANDOMWALK',
    'ISTWOSTATERANDOMWALK']] = False

df[['SIZELEFT','SIZERIGHT']] = 0

#df = df.drop(['REPRESSIONRATE', 'UNREPRESSIONRATE', 'REPRLENLEFT',
#              'REPRLENRIGHT'], axis=1)
df.to_csv(path_or_buf='tf_100_repr_2.csv', sep=';', index=False)

# %%
rnd.seed(9672391)
dna_len = 100
bps = ['a', 'c', 'g', 't']
freqs = [0.25] * 4
seq = ''
for bp, freq in zip(bps, freqs):
    seq += bp * int(freq * dna_len)
l = list(seq)
rnd.shuffle(l)
seq = ''.join(l)
print(seq)

with open('seq.fasta', 'w') as f:
    f.write(seq)
f.close()
# %%
with open('seq.btrack', 'w') as f:
    f.write('1\n' * dna_len)
f.close()
# %%

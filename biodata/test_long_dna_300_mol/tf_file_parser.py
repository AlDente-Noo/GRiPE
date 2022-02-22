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


with open('tf_1.csv') as f:
	df = pd.read_csv(f, sep=';')

df

#df.to_csv(path_or_buf='tf_1.csv', sep=';', index=False)

#%%
df[['SLIDELEFTPROBABILITY', 'SLIDERIGHTPROBABILITY']] = 0.4999
df['UNBINDINGPROBABILITY'] = 1 - 0.4999*2
df['JUMPINGPROBABILITY'] = 0.1
df['COPYNUMBER'] = 1
df.to_csv(path_or_buf='tf_1.csv', sep=';', index=False)

# %%
rnd.seed(9672391)
dna_len = 16
bps = ['a', 'c', 'g', 't']
freqs = [0.25] * 4
seq = ''
for bp, freq in zip(bps, freqs):
    seq += bp * int(freq * dna_len)
l = list(seq)
rnd.shuffle(l)
seq = ''.join(l)
print(seq)

#with open('seq.fasta', 'w') as f:
#    f.write(seq)
#f.close()
# %%
with open('seq.btrack', 'w') as f:
    f.write('1\n' * dna_len)
f.close()
# %%

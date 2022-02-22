#%%
import random as rnd
import pandas as pd
import numpy as np

with open('tf_100.csv') as f:
	df = pd.read_csv(f, sep=';')

#%%
df['ES'] = 1
#df.at[:, 'COPYNUMBER'] = 100
df.at[:, 'COPYNUMBER'] = 1
df.drop(['REPRESSOR'], axis='columns', inplace=True)
df.drop(['REPLENLEFT','REPLENRIGHT'], axis='columns', inplace=True)
df.drop(['PWMREPTHRESHOLD','REPRESSIONPROBABILITY'], axis='columns', inplace=True)
#df['REPRESSIONRATE'] = np.array([0., 0., 1., 1., 0., 0., 1., 1.]) * 0.2
df['REPRESSIONRATE'] = [0., ] * 8
df['UNREPRESSIONRATE'] = df['REPRESSIONRATE'] * 0.1
df['REPRLENLEFT']  = [100,] * 8
df['REPRLENRIGHT'] = [100,] * 8
df['SIZELEFT'] = [0,] * 8
df['SIZERIGHT'] = [0,] * 8
#df

#%%
#df.to_csv(path_or_buf='tf_100_new.csv', sep=';', index=False)
df.to_csv(path_or_buf='tf_1_no_repr.csv', sep=';', index=False)

#%%
df['COPYNUMBER'] = [0,] * 8
tf1 = 'hb'
df.at[df['name'].to_list().index(tf1), 'COPYNUMBER'] = 1
#df.to_csv(path_or_buf='tf_1_'+tf1+'.csv', sep=';', index=False)

fn = 'kr_ts_new_dir'
sites = open(fn + '.txt')
sites_tf1 = open(fn + '_' + tf1 + '.txt', 'w')
for site in sites:
	if site.startswith(tf1):
		sites_tf1.write(site)
sites.close()
sites_tf1.close()

# %%
df['ASSOCRATE'] = [10.,] * 8
df['UNBINDINGPROBABILITY'] = [1.,] * 8
df['JUMPINGPROBABILITY'] = [1.,] * 8
df['SLIDELEFTPROBABILITY'] = [0.,] * 8
df['SLIDERIGHTPROBABILITY'] = [0.,] * 8
#%%
df.to_csv(path_or_buf='tf_1_'+tf1+'_no_facilit.csv', sep=';', index=False)

# %%
df['ASSOCRATE'] = [10.,] * 8
df['UNBINDINGPROBABILITY'] = [0.,] * 8
df['JUMPINGPROBABILITY'] = [0.,] * 8
df['SLIDELEFTPROBABILITY'] = [0.5,] * 8
df['SLIDERIGHTPROBABILITY'] = [0.5,] * 8
# %%
df.to_csv(path_or_buf='tf_1_'+tf1+'_no_3d_diff.csv', sep=';', index=False)
# %%

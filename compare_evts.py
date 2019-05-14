from ROOT import TFile, gDirectory, TH1F
import root_pandas as rp
import numpy as np
import matplotlib.pyplot as plt

path1 = "/afs/hephy.at/data/higgs01/sync_2018/VBF_2018_v6/mt-NOMINAL_sync_2018.root"
path2 = "/afs/cern.ch/work/c/ccaillol/public/cecile_mt2018.root"

name1 = "Vienna" 
name2 = "Cecile"

tkey1 = "TauCheck"
tkey2 = "mutau_tree"

col = ['evt','mjj','jpt_1','jpt_2']

df1 = rp.read_root(path1, key= tkey1,columns=col)
df2 = rp.read_root(path2, key= tkey2,ignore='addlepton_p4',columns=col)

evt_list1 = np.array(df1['evt'])
evt_list2 = np.array(df2['evt'])


# for e in evt_list1[:100] :
# 	if e in evt_list2 :
# 		mask1 = df1['evt']==e
# 		mask2 = df2['evt']==e

# 		print "====Event: {0}====".format(e)
# 		print df1[col][mask1]
# 		print df2[col][mask2]

# 	else :
# 		print "event {0} not in both samples".format(e)

mjj_fixed = []
for e in evt_list2[:-1] :
	mask2 = df2['evt']==e
	if float(df2["jpt_1"][mask2]) < 0 or float(df2["jpt_2"][mask2]) < 0 :
		mjj_fixed.append(-10.)
	else :
		mjj_fixed.append(float(df2["mjj"][mask2]))

mjj_unchanged = np.array(df1["mjj"])[:-1]
#plt.hist(np.array(df1["mjj"]),bins=100)
fig, (ax1, ax2) = plt.subplots(nrows=2)

ax1.hist(np.array(mjj_fixed),bins=100,range=[0.,5000.],histtype="step")
ax1.hist(mjj_unchanged,bins=100,range=[0.,5000.],histtype="step")
ax1.legend([name1,name2])
ax1.set_yscale("log")

numerator, be1 = np.histogram(mjj_unchanged,bins=100,range=[0.,5000.])
denominator, be2 = np.histogram(np.array(mjj_fixed),bins=100,range=[0.,5000.])


ratio = [float(x)/(denominator[i]+0.1) for i,x in enumerate(numerator)]
print ratio

ax2.plot(np.linspace(0.,5000.,100),ratio,'k*')
ax2.set_ylabel("{0} / {1}".format(name1,name2))
ax2.set_ylim([0.5,1.5])
plt.savefig("test_ratio.png")

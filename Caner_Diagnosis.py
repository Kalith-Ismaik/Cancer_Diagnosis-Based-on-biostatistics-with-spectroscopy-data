#!/usr/bin/env python
# coding: utf-8

# # Mathematical Methods Of Biopectroscopy Invivo: Assignment

# ## Reading the spectral data  from CSV file and visualizing for finding out location of peaks.

# In[2]:


import pandas as pd
import matplotlib.pyplot as plt
import math

#reading csv file.
data = pd.read_csv('spectra1.csv')
#print(data.info())

#preparing sets of different samples and taking a mean of their absorbance at different wavelenghts.
Tumor1_before_pdt = data.iloc[0:,1:9].mean(axis=1)
Tumor2_before_pdt = data.iloc[0:,9:13].mean(axis=1)
Tumor_before_pdt = data.iloc[0:,1:13].mean(axis=1)
Tumoredge_before_pdt = data.iloc[0:,13:18].mean(axis=1)
Tumor_after_pdt = data.iloc[0:,21:33].mean(axis=1)
Tumoredge_after_pdt = data.iloc[0:,[18,19,20,36,37,38]].mean(axis=1)
Normal = data.iloc[0:,18:21].mean(axis=1)

wave = data.iloc[0:,0]
#print(wave)

#plotting graph with respect to wavelength and the test values.
plt.title("Tumor1_before_pdt vs Tumor2_before_pdt")

plt.plot(wave,Tumor1_before_pdt,label='Tumor1_before_pdt')
plt.plot(wave,Tumor2_before_pdt,label='Tumor2_before_pdt')
plt.legend()
plt.xlabel("wavelength")
plt.ylabel("absorbance")
plt.show()

plt.title("Tumor_before_pdt vs Normal")
plt.plot(wave,Tumor_before_pdt,label='Tumor_before_pdt')
plt.plot(wave,Normal,label='Normal')
plt.legend()
plt.xlabel("wavelength")
plt.ylabel("absorbance")
plt.show()

plt.title("Tumor_after_pdt vs Normal")
plt.plot(wave,Tumor_after_pdt,label='Tumor_after_pdt')
plt.plot(wave,Normal,label='Normal')
plt.legend()
plt.xlabel("wavelength")
plt.ylabel("absorbance")
plt.show()

plt.title("Tumor_before_pdt vs Tumor_after_pdt")
plt.plot(wave,Tumor_before_pdt,label='Tumor_before_pdtt')
plt.plot(wave,Tumor_after_pdt,label='Tumor_after_pdt')
plt.legend()
plt.xlabel("wavelength")
plt.ylabel("absorbance")
plt.show()

plt.title("Tumoredge_before_pdt vs Tumoredge_after_pdt vs normal")
plt.plot(wave,Tumoredge_before_pdt,label='Tumoredge_before_pdt')
plt.plot(wave,Tumoredge_after_pdt,label='Tumoredge_after_pdt')
plt.plot(wave,Normal,label='Normal')
plt.legend()
plt.xlabel("wavelength")
plt.ylabel("absorbance")
plt.show()

plt.title("Tumor_before_pdt vs Tumoredge_before_pdt vs normal")
plt.plot(wave,Tumor_before_pdt,label='Tumor_before_pdt')
plt.plot(wave,Tumoredge_before_pdt,label='Tumoredge_before_pdt')
plt.plot(wave,Normal,label='Normal')
plt.legend()
plt.xlabel("wavelength")
plt.ylabel("absorbance")
plt.show()


# ## Localizing peaks and calculation of peak area, amplitude and width

# In[1]:


import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import find_peaks,peak_widths
from scipy import signal
from scipy.integrate import trapz
import pandas as pd
from tabulate import tabulate

#Read spectral data 
data = pd.read_csv('spectra1.csv')
peak_I = []
peak_F = []
peakw_I = []
peakw_F = []
area_I = []
area_F = []

#Estimating the peak location
for i in range(1,len(data.columns)):
    x = data.iloc[675:738,i]
    y = data.iloc[970:1110,i]
    
#calcualtion of area of peak
    arex = trapz(x)
    area_I.append(arex)
    arey = trapz(y)
    area_F.append(arey)
    
#Finding peaks and estimating repective peaks amplitude and half width
#Peak data acquisition for incident laser light
    peak1,v1 = find_peaks(x,height= -250)
    peak2,v2 = find_peaks(y,height= -250)
    l1 = []
    a = str(v1.values())[20:-5].split(',')
    for ii1 in a:
        l1.append(float(ii1))
    peak_I.append(max(l1))
    jj = []
    loc = l1.index(max(l1))
    peak = list(peak1)
    jj.append(peak[loc])
    width = peak_widths(x,jj,rel_height = 0.5)
    w=width[0]
    peakw_I.append(w[0])

#Peak data acquisition for flourescent light
    l2 = []
    b = str(v2.values())[20:-5].split(',')
    for ii2 in b:
        l2.append(float(ii2))
    peak_F.append(max(l2))
    jj = []
    loc = l2.index(max(l2))
    peak = list(peak2)
    jj.append(peak[loc])
    width = peak_widths(y,jj,rel_height = 0.5)
    w=width[0]
    peakw_F.append(w[0])
    
title=list(data.columns.values.tolist())
title.pop(0)  

#Calculating Ratios with respect to peak area, peak amplitude and peak half width
ratio = {}
for j in range(0,len(peak_I)):
    rat_amp = abs(peak_F[j]/peak_I[j])*100
    rat_width = abs(peakw_F[j]/peakw_I[j])*100
    rat_area = abs(area_F[j]/area_I[j])*100
    ratio[title[j]] = rat_area,rat_amp,rat_width
    
#Visualizing data calculated from ratios in both graphical and tabulated forms
val = [*range(1,39,1)]
headers = ['Sample','Area_ratio','Amplitude_ratio','Width_ratio']
print(tabulate([(k,) + v for k,v in ratio.items()],headers = headers))
are = [a_tuple[0] for a_tuple in ratio.values()]
amp = [a_tuple[1] for a_tuple in ratio.values()]
wid = [a_tuple[2] for a_tuple in ratio.values()]

plt.bar(val,are)
plt.ylabel("Peak area Ratio")
plt.xlabel("Spectra")
plt.show()

plt.bar(val,amp)
plt.ylabel("Peak Amplitude Ratio")
plt.xlabel("Spectra")
plt.show()

plt.bar(val,wid)
plt.ylabel("Peak Width Ratio")
plt.xlabel("Spectra")
plt.show()


# ## Processing datas from calculated ratios and estimate statistical distribution parameter

# In[3]:


from scipy.stats import skew, kurtosis, mode
import math

Ratio_of_Area = np.array(are)
Ratio_of_Amplitude = np.array(amp)
Ratio_of_width = np.array(wid)
n = 0
sample = [Ratio_of_Area,Ratio_of_Amplitude]
stats = []
for zoy in sample:
    if  n==0:print("Distribution data for the Ratio of Area:")
    elif n==1:print("Distribution data for the Ratio of Width:")
    nor = list(zoy[17:20])+list(zoy[36:39])
    norm = np.array(nor)
    data_dict={
        'Tumor before PDT 1':zoy[0:8],
        'Tumor before PDT 2':zoy[8:12],
        'Tumor before PDT':zoy[0:12],
        'Tumor edge before PDT':zoy[12:17],
        'Normal [Doctor Hand]':norm,
        'Tumor after PDT':zoy[20:32],
        'Tumor edge after PDT':zoy[32:36]
    }
    n+=1
    for k in data_dict:
        print("Distribution data for",k,":")
        x = data_dict[k]
        stat_dict={
            'mode': mode(x),
            'skewness': skew(x),
            'kurtosis': kurtosis(x),
            'max': x.max(),
            'min': x.min(),
            'sum': x.sum(),
            'mean': x.mean(),
            'std': x.std(),
            'median': np.median(x)
        }
        
        guass = []
        xy = np.sort(x)
        sam = list(xy)
        for i in sam:
            ex = np.exp(-0.5*np.square((i-x.mean())/x.std()))/(x.std()*np.sqrt(2*math.pi))
            guass.append(ex)
        xx = np.array(guass)
        plt.plot(xy,xx)
        plt.show()
        
        print(tabulate(zip(stat_dict.keys(),stat_dict.values())))
        stats.append(stat_dict)


# ## Testing multiple comparison hypothesis testing of Bornferroni correction

# In[5]:


import numpy as np
import pandas as pd
from scipy.stats import ttest_rel, f_oneway
from statsmodels.sandbox.stats.multicomp import MultiComparison
from scikit_posthocs import posthoc_tukey as tukey

Ratio_of_Area = np.array(are)
Ratio_of_Amplitude = np.array(amp)
Ratio_of_width = np.array(wid)
n = 0
sample = [Ratio_of_Area,Ratio_of_Amplitude]
stats = []
for zoy in sample:
    nor = list(zoy[17:20])+list(zoy[36:39])
    norm = np.array(nor)
    
    stat,pvalue = f_oneway(zoy[0:8],zoy[8:12],zoy[0:12],zoy[12:17],norm,zoy[20:32],zoy[32:36])
    conc = "in" if pvalue>0.05 else ""
    lu = "the data from peak area ratio." if n==0 else "the data from peak amplitude ratio."
    print(f'Results are statically {conc}significant for {lu}'
          f'\nWe can {"not" if pvalue>0.05 else ""} apply pairwise comparison.')
    n+=1
   
    if pvalue<0.05:
        dat = pd.DataFrame(dict([ (k,pd.Series(v)) for k,v in data_dict.items() ]))
        dat = dat.melt(var_name='groups', value_name='values')
    # Multiple comparison
    # Student criteria with Bonferroni correction
        mcls = MultiComparison(dat['values'],dat['groups'])
        print(mcls.allpairtest(ttest_rel, method='bonf')[0])
        #print(mcls.allpairtest(ttest_rel, method='holm')[0])
        #print(mcls.tukeyhsd().summary())
    


# In[ ]:





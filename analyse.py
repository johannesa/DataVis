import numpy as np
import glob 
import matplotlib.pyplot as plt
plt.rcParams.update({"figure.max_open_warning": 0}) #plt won't complain about opening more than 20 figures
from matplotlib.ticker import NullFormatter
from astropy.stats import LombScargle

#Acquire the data files 
filename = sorted(glob.glob("lmc*.dat"))
def find_period(filename):
    """
    Apply the Lomb-Scargle periodogram to compute the period if any exists. 
    Compute the L-S periodogram using an automatically-determined frequency 
    grid and find the frequency of maximum power. L-S power is always a unitless
    quantity, because it is related to the chi-square of the best-fit periodic 
    model at each frequency.
    """
    file_base = filename.split('/')[-1].split('.')[0]
    [hjd, mag, magerr] = np.genfromtxt(fname=filename, usecols=(0,1,2), unpack=True)    
    freq, ls_power = LombScargle(t=hjd, y=mag, dy=magerr).autopower(nyquist_factor=5)
    T = 1/freq[np.argmax(ls_power)] #period is inverse of frequency
    """
    Split the data over several intervals, and fold each interval about the determined period. Binning and
    folding are done here. Refer to http://www.astrobetter.com/wiki/Periodicity+Analysis and
    https:www.southampton.ac.uk/~sdc1g08/BinningFolding.html for more information. 
    """
    #compute the phase
    phase = (hjd/T) % 1 #use fractional phase
    #define bins and binwidth for data binning and folding over each bin
    nbins = 35; binwidth = 1.0/float(nbins)
    #create arrays for bin values and weights
    bins = np.zeros(nbins)
    weights = np.zeros(nbins)   
    for i in range(len(mag)):
        #calculate bin number for this value-fold each interval about the calculated period
        n = int(phase[i]/binwidth)
        #calculate weight (w) - inverse square of magerr
        w = magerr[i]**-2.0
        #add weighted value to bin value (value times weight)
        bins[n] += mag[i] * w
        #add weight to bin weights (weights)
        weights[n] += w
    bins /= weights #normalise weighted values using sum of weights 
    binErr = np.sqrt(1.0/weights) #calulate bin errors from squared weights
    binEdges = np.arange(nbins)*binwidth #create array of bin edge for plotting
    #Find the area under the curve (auc) for each star
    auc = np.cumsum(ls_power)
    """
    Find the 1st derivative of the cumulative sum in steps of 100. We would prefer the sample distance
    to be large enough in order to avoid wrongful identifcation of steps over short distance.
    """
    delta = np.gradient(auc[np.min(np.where(freq > 0.02)):-1,], 100)
    max_del = np.max(delta)
    return file_base,T,auc,hjd,mag,freq,ls_power,phase,binEdges,bins,delta,max_del

def make_plot(file_base,T,auc,hjd,mag,freq,ls_power,phase,binEdges,bins,delta,max_del,rank):
    """
    Lets make plots 
    """
    plt.figure(figsize=[16,16])
    plt.suptitle("Variation of %s with T = %.2f days and auc = %.2f and max_del = %.4f and rank = %.f" %(file_base,T,auc[-1],max_del,rank))

    #mag vs time (hjd)
    plt.subplot(321)
    plt.plot(hjd, mag, linestyle="None", marker=".")
    plt.xlabel("Time [HJD]", fontsize='x-large'), plt.ylabel("Magnitude", fontsize='x-large')
    
    #Lomb-Scargle power vs period
    plt.subplot(322)
    plt.plot(freq, ls_power)
    #plt.xlim(min(1/freq),T+0.8) #view the spectra over whole range of periods
    plt.xlabel("Frequency [per day]", fontsize='x-large'), plt.ylabel("Lomb-Scargle Power", fontsize='x-large')

    #mag vs phase (folded time)
    plt.subplot(323)
    plt.errorbar(phase, mag, linestyle="None", marker=".")
    plt.ylim(0.02+max(mag), 0.02+min(mag))
    plt.xlabel("Folded Time", fontsize='x-large'), plt.ylabel("Magnitude", fontsize='x-large')
    plt.subplot(324)
    plt.errorbar(binEdges, bins, linestyle="None", marker=".")
    plt.ylim(max(mag), min(mag))
    plt.xlabel("Folded Time", fontsize='x-large'), plt.ylabel("Weighted Average Magnitude", fontsize='x-large')
        
    #cumulative sum of ls_power vs frequency
    plt.subplot(325)
    plt.plot(freq, auc)
    plt.xlabel("Frequency [per day]", fontsize='x-large'), plt.ylabel("$\Sigma$(Lomb-Scargle Power)", fontsize='x-large')
    
    #differential cumulative sum vs frequency
    freq_new = freq[np.min(np.where(freq > 0.02)):-1,]
    plt.subplot(326)
    plt.plot(freq_new, delta, label="Data")
    plt.axvline(x=freq_new[np.argmax(delta)], ymin=0.05, ymax=0.95, color="r", label="max($\delta$)")
    plt.xlabel("Frequency [per day]", fontsize='x-large'), plt.ylabel("$\delta$($\Sigma$Lomb-Scargle Power)", fontsize='x-large')
    plt.legend()
    
    #figure settings 
    plt.gca().yaxis.set_minor_formatter(NullFormatter())
    #plt.subplots_adjust(top=0.92, bottom=0.08, left=0.10, right=0.95, hspace=0.25, wspace=0.35)
    plt.savefig("LS_vs_Period/rank%.f_%s.png" %(rank,file_base))

#run the find_period and make_plots functions respectively.
rank_table = np.zeros((1,2))
for f in filename:
    #We will only want to work with data files with more than 50 data points. Anything below this shall not be used for our analysis.
    NumLines = sum(1 for line in open(f)) #the no. of lines represent the no. of data points since the data doesn't have any header
    if NumLines <= 100:
        print ("%s has few data points so it won't be used for the analysis" %f)
        print ("=========================")
    else:
        #analyse the data since it should have adequate points
        print (f)
        file_base,T,auc,hjd,mag,freq,ls_power,phase,binEdges,bins,delta,max_del = find_period(f)
        print ("The period of %s is %.2f days" %(file_base,T))
        print ("The area under the curve is %.2f" %auc[-1])
        print ("======================================")
        if max_del > 0.003 and auc[-1] < 120:
            rank_table = np.append(rank_table, [[f,max_del]], axis=0)

"""
for the classification of the star systems assending order by the maximum of the gradient of the cumulative LombSargle power
"""
rank_table = rank_table[1:,]
rank_table = rank_table[np.argsort(rank_table[:, 1])]

info_table = np.zeros((1,3)) # table for doing PCA and k-means clustering. Copy and paste this table into excel and export into Orange interactive data visualisation tool for classification of the stars 
for f in rank_table[:,0]:
    rank = len(rank_table[:,0])-np.where(rank_table[:,0] == f)[0][0]
    file_base,T,auc,hjd,mag,freq,ls_power,phase,binEdges,bins,delta,max_del = find_period(f)
    for i in np.arange(0,len(bins)): 
        info_table = np.append(info_table, [[f,binEdges[i],bins[i]]], axis=0)
    make_plot(file_base,T,auc,hjd,mag,freq,ls_power,phase,binEdges,bins,delta,max_del,rank)
    
info_table = info_table[1:,] # remove unwanted 0,0,0 row


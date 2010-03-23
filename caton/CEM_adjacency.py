# This file has a mixture-of-gaussian clustering algorithm. It optionally uses sparse matrices along
# the functions in estimate_cov.py to estimate sparse concentration matrices, instead of using the
# full covariance matrices.
# This is currently not being used, because I found that splitting and recombining the probe worked better than clustering
# on the whole probe at once (though the whole-probe method did give reasonable results on 32 channels), and
# KlustaKwik works fine for this purpose.

# Requires compiled CEM_extensions

from __future__ import division
from CEM_extensions import *
import numpy as np
if False: import matplotlib.pyplot as plt
if False: import matplotlib
import time
import os
import estimate_cov
import scipy.sparse as sparse
from scipy.stats import scoreatpercentile
import multiprocessing as mp
import itertools as it


PlottedX = None
UnitsPerChannel = 3
NoiseFrac = .01
MinIter = 2
MaxIterDefault = 200
FeaturesPerChannelDefault = 3
TAB_LEVEL = -1
ChangedFrac = .001
UseMyPenalty = True
HUGE = 9.e20

VERBOSE = 1
N_NOISE = 0
BestSplitClu = 0
BestDelClu = 0

def get_InvSqrtCov(Cov_mff,GraphInfo):
    "Calculates concentration matrices given sample covariance matrices and GraphInfo object"
    InvSqrtCov_mff = [None]*len(Cov_mff)
    ### we want rows contiguous for left multiplication
    for m in range(len(Cov_mff)):
        try:
            if GraphInfo is None or GraphInfo.IsComplete:
                InvSqrtCov_mff[m] = np.matrix(np.linalg.cholesky(np.linalg.inv(Cov_mff[m])).T.copy())
                ### copy() because we want rows stored contiguously
            else:
                InvSqrtCov_mff[m] = sparse.csr_matrix(np.linalg.cholesky(GraphInfo.concentration(Cov_mff[m])).T)
        except np.linalg.LinAlgError:
            InvSqrtCov_mff[m] = np.matrix(np.zeros_like(Cov_mff[0]))
            tprint("singular covariance matrix encountered with cluster %i/%i"%(m+1,len(Cov_mff)),colorstring=bcolors.WARNING)
    return InvSqrtCov_mff


def e_step(X_nf,Mean_mf,Cov_mff,Weight_m,F,M,N,GraphInfo):
    "Finds responsibilities of data in clusters"    
    LogP_nm = np.zeros((N,M),dtype=np.float32)        
    InvSqrtCov_bff = get_InvSqrtCov(Cov_mff[1:],GraphInfo)                    
    LogP_nm[:,1:] = compute_LogP(X_nf,Mean_mf[1:],Weight_m[1:],InvSqrtCov_bff)
    lowscore = scoreatpercentile(LogP_nm[:,1:].max(axis=1),NoiseFrac*100)
    global Outliers
    Outliers = np.argsort(LogP_nm[:,1:].max(axis=1))[:N_NOISE]
    noise_cutoff = a_little_under(lowscore) if not np.isnan(lowscore) else -HUGE
    LogP_nm[:,0] = noise_cutoff
    # are we ever actually using LogP_nm[:,0]? Or could we set it to whatever we want?
    return LogP_nm
    
    
def c_step(LogP,Class_n,F,M,N):
    "Finds the best cluster for each point"
    OldClass_n = Class_n.copy()
    Class_n,Class2_n = argmax_two_per_row(LogP)
    NChanged = (Class_n != OldClass_n).sum()
    return Class_n,Class2_n,NChanged


def m_step(Class_n,X_nf,F,M,N):
    "Finds sample means and covariance"
    Mean_mf = class_means(X_nf,Class_n,M)
    Cov_mff = class_covs(X_nf,Mean_mf,Class_n,M)
    Weight_m = class_wts(Class_n,M)
    return Mean_mf, Cov_mff, Weight_m

def try_delete(LogP,Class_n,Class2_n,GraphInfo,N_m,F,M,N):
    "Check if deletion of some cluster will reduce score"
    global BestDelClu
    
    tprint("Testing if deletion helps")
    ###only run this if you're not killing any other clusters
    if M <= 2: return None
    dScore_m = deltapen(LogP,Class_n,Class2_n,N_m,F,GraphInfo)
    
    argmax_dScore,max_dScore = dScore_m.argmax(),dScore_m.max()
    BestDelClu = argmax_dScore
    
    if max_dScore > 0:
        Clu2Del = argmax_dScore
        tprint("Yes: score changes by %.1f"%max_dScore,colorstring=bcolors.OKBLUE)
        return Clu2Del
    else:
        tprint("No: score changes by %.1f"%max_dScore,colorstring=bcolors.OKGREEN)      
        return None

def deltapen(LogP,Class_n,Class2_n,N_m,F,GraphInfo):
    "Returns change in penalty if each of M clusters is removed, and points are assigned to next best cluster"
    N,M = LogP.shape
    dLikelihood_m = class_sums(LogP[xrange(N),Class2_n] - LogP[xrange(N),Class_n],Class_n,M) #negative
    dLikelihood_m[0] = -np.inf ###so we won't delete it
    T_mm = transitions(Class_n,Class2_n,M)
    dPenalty_m = params_per_cluster(F,GraphInfo)/2 * (np.log(N_m)-np.dot(T_mm,1/N_m)) #positive
    return dPenalty_m + dLikelihood_m
    
def channel_init(X_nf,NChannels):
    "Initialize by clustering separately on channels 0:FPC, 1:FPC+1, ..."
    N,F = X_nf.shape
    X_ncq = X_nf.reshape(N,NChannels,-1)
    BestChannel_n = ((X_ncq**2).sum(axis=2)).argmax(axis=1)
    PeakChannelIs = subset_inds(BestChannel_n)
    Class_n = np.zeros(len(X_nf),np.int32)
    TotalClasses = 0
    for ch in range(NChannels):
        if len(PeakChannelIs[ch]) != 0:
            tprint("initialization: sorting on channel %i"%ch)        
            ClassHere,Score,M = cem(X_ncq[PeakChannelIs[ch],ch,:],Recurse=True,init=('random',UnitsPerChannel),GraphInfo=None)
            Class_n[PeakChannelIs[ch]] = (ClassHere + TotalClasses)*(ClassHere != 0)   
            TotalClasses += M-1
    return Class_n
    
    
def params_per_cluster(F,GraphInfo):
    return GraphInfo.NParams + F + 1 if GraphInfo is not None else F*(F+1)/2+F+1
    
def classic_penalty(Class_n,GraphInfo,F,M,N):
    "Usual BIC penalty"
    ParamsPerCluster = GraphInfo.NParams + F + 1 if GraphInfo is not None else F*(F-1)/2+F+1
    NParams = ParamsPerCluster*M
    return NParams*np.log(N)/2 #- M*np.log(M) # from permutations

def my_penalty(N_m,GraphInfo,F,M,N):
    "Modified penalty, which should be more accurate in clustering: Sum_m NParams/2 * log(N_m)/2"
    return params_per_cluster(F,GraphInfo)/2 * np.log(N_m[1:]+1).sum() # +1 just so its not zero

penalty = my_penalty if UseMyPenalty else classic_penalty

def compute_score(LogP,N_m,Class_n,GraphInfo,F,M,N):
    "Compute likelihood of data given parameters plus penalty"
    LogL = np.delete(LogP[xrange(N),Class_n],Outliers).sum()
    #LogL = LogP[xrange(N),Class_n].sum()    
    BIC = penalty(N_m,GraphInfo,F,M,N)
    return LogL,BIC,LogL-BIC


def cem(X_nf,Recurse,init, MaxIter = MaxIterDefault,GraphInfo=None,Plotting=False):
    """
    If GraphInfo = None, a complete graph is assumed.
    Cluster the data, in X_nf
    Return class assignments, score, and number of clusters.
    
    X_nf : ndarray. Rows (1..N) are observations, cols (1..F) are features
    Recurse : boolean. If yes, then try to split, delete, and tunnel clusters
    init :  ("random",number_of_starting_clusters)
          | ("channel", number_of_channels)
          | ("given", integer_ndarray)
    """
    global TAB_LEVEL,N_NOISE
    TAB_LEVEL += 1        
    tprint("__________",verb_level=2)
    
    N,F = X_nf.shape
    R = GraphInfo.Rank if GraphInfo is not None else F+1
    R = max(R,10)
    tprint("Note: min(R) = 10",bcolors.WARNING)
    
    N_NOISE = N*NoiseFrac
    
    if init[0] == 'random': ###init = ('random',Number of Starting Clusters (other than noise))
        Class_n = 1+np.random.randint(init[1],size = N).astype(np.int32)
    elif init[0] == 'channel': ###init = ('channel',NChannels)
        NChannels = init[1]
        assert (NChannels is not None)
        Class_n = channel_init(X_nf,NChannels)
    elif init[0] == 'given':
        Class_n = init[1]
    else: raise Exception
    M = Class_n.max()+1           

    IsTooSmall = (bincount(Class_n,M)< R) & (np.arange(M) > 0)    
    if IsTooSmall.sum() > 0:            
        Class_n,M = remove_clusters(IsTooSmall,Class_n, np.zeros(N,np.int32),"too few points",M)
        DidSomething = True
    
    LogL = -HUGE; Score = -HUGE  
    Iteration = 0        
    tprint("Starting EM iterations")
    while(M > 1):
        Mean_mf, Cov_mff, Weight_m = m_step(Class_n,X_nf,F,M,N)
        N_m = bincount(Class_n,M)

        LogP = e_step(X_nf,Mean_mf,Cov_mff,Weight_m,F,M,N,GraphInfo)  
        OldLogL = LogL        
        LogL,BIC,Score = compute_score(LogP,N_m,Class_n,GraphInfo,F,M,N)
        tprint("After M. Iteration: %i, Score: %.1f-%.1f=%.1f, N_m: %s"%(Iteration, LogL,BIC,Score,N_m),verb_level=2)                       

        if Iteration > MinIter:
            if (NotMuchChanged and not DidSomething) or (Iteration > MaxIter):
                break
            if (LogL < OldLogL):
                tprint("LogL decreased",bcolors.FAIL)

        Class_n,Class2_n,NChanged = c_step(LogP,Class_n,F,M,N)        
        if Plotting:
            cluster_plot(X_nf,Class_n,Mean_mf,Cov_mff)        
        Iteration += 1        

        LogP = e_step(X_nf,Mean_mf,Cov_mff,Weight_m,F,M,N,GraphInfo)  
        OldLogL = LogL        
        LogL,BIC,Score = compute_score(LogP,N_m,Class_n,GraphInfo,F,M,N)
        tprint("After C. Iteration: %i, Score: %.1f-%.1f=%.1f, NChanged: %i, N_m: %s"%(Iteration, LogL,BIC,Score,NChanged,N_m),verb_level=2)                       
        
        
        DidSomething = False
        NotMuchChanged = NChanged < max(2,N*ChangedFrac)        
        IsTooSmall = (Weight_m*N < R) & (np.arange(M) > 0)        
        if IsTooSmall.sum() > 0:            
            Class_n,M = remove_clusters(IsTooSmall,Class_n, Class2_n,"too few points",M)
            DidSomething = True
        
        LetsTryOtherStuff = NotMuchChanged and not DidSomething and Recurse and Iteration >= MinIter
            
        if LetsTryOtherStuff:
            ToKill = try_delete(LogP,Class_n,Class2_n,GraphInfo,N_m,F,M,N)
            if ToKill is not None:
                KillIndicator = np.zeros(M,dtype=bool); KillIndicator[ToKill] = True
                DidSomething = True
                Class_n,M = remove_clusters(KillIndicator,Class_n,Class2_n,"for net score increase",M)
    
        if LetsTryOtherStuff and not DidSomething:
            Class_n,Score,M,DidSplit = try_splits(X_nf,Class_n,Score,F,M,N,R,GraphInfo)                 
            DidSomething = DidSomething or DidSplit
            
        if LetsTryOtherStuff and not DidSomething:
            Class_n,Score,M,DidTunnel = try_tunnel(X_nf,Class_n,Class2_n,N_m,Score,F,M,N,R,GraphInfo)                                     
            DidSomething = DidSomething or DidTunnel
            
        if DidSomething:
            Iteration = 0            
                 
        if Iteration >= MinIter:
            if (NotMuchChanged and not DidSomething) or (Iteration > MaxIter):
                break
            if (LogL < OldLogL):
                tprint("LogL decreased",bcolors.FAIL,1)
                #break
                        
    tprint("Done searching",1)
    TAB_LEVEL -= 1
    return Class_n,Score,M

def remove_clusters(KillIndicator,Class_n,Class2_n,Message,M):
    """Remove clusters where KillIndicator is True. Switch points from removed cluster to
    next best cluster"""
    for clu in np.flatnonzero(KillIndicator):
        tprint("Removing cluster %i: %s"%(clu,Message))

    
    np.putmask(Class_n,KillIndicator[Class_n],Class2_n)
    RelabelArr = np.zeros(M,np.int32)
    RelabelArr[np.flatnonzero(~KillIndicator)] = np.arange(M,dtype=np.int32)
    Class_n = RelabelArr[Class_n]

    M = M-KillIndicator.sum()
    tprint("%i clusters remain"%M)    
    return Class_n,M
    



def try_splits(X_nf,Class_n,Score,F,M,N,R,GraphInfo):
    "Try to split clusters"
    global BestSplitClu
    BestSplitScore,BestSplitClu = -HUGE,0
    
    
    ClassInds = subset_inds(Class_n,M)
    for m in range(1,M):
        if len(ClassInds[m]) <= 2.5*R: continue        
        tprint("Trying split on cluster %i"%m)
        tprint("Iterating with unsplit cluster")            
        SubsetUnsplitClass,SubsetUnsplitScore,SubsetUnsplitM = cem(X_nf[ClassInds[m]],Recurse=False,init=('random',1),GraphInfo=GraphInfo)
        tprint("Iterating with split cluster")                        
        SubsetSplitClass,SubsetSplitScore,SubsetSplitM = cem(X_nf[ClassInds[m]],Recurse=False,init=('random',2),GraphInfo=GraphInfo)
        BestSplitScore,BestSplitClu = max_update(SubsetSplitScore,m,BestSplitScore,BestSplitClu) #!!!!#
        if SubsetSplitScore > SubsetUnsplitScore:
            tprint("Split cluster has %i members. Unsplit cluster has %i members"%(SubsetSplitM,SubsetUnsplitM))
            if SubsetSplitM == SubsetUnsplitM:
                tprint("WEIRD. Split and Unsplit clusters have same number of members.")
            SplitClass = Class_n.copy()
            RelabelArr = np.array([0,m,M],np.int32)
            SplitClass[ClassInds[m]] = RelabelArr[SubsetSplitClass]
            tprint("Seeing if splitting cluster reduces TOTAL score")            
            SplitClass,SplitScore,SplitM = cem(X_nf,Recurse=False,init=('given',SplitClass),GraphInfo=GraphInfo)
            if SplitScore > Score:
                tprint("Splitting cluster %i"%m)
                return SplitClass,SplitScore,SplitM,True
    tprint("No splits")
    return Class_n,Score,M,False

def try_tunnel(X_nf,Class_n,Class2_n,N_m,Score,F,M,N,R,GraphInfo):
    "Try to simultaneously delete the best deletion candidate and split the best split candidate"
    if BestDelClu == BestSplitClu:
        tprint("Best split candidate is best deletion candidate: skipping tunnelling")
        return Class_n,Score,M,False
    tprint("Try tunneling. Delete %i and split %i"%(BestDelClu,BestSplitClu))
    Class_n_copy = Class_n.copy()
    np.putmask(Class_n_copy,Class_n_copy == BestDelClu,Class2_n)
    split_inds = np.flatnonzero(Class_n_copy == BestDelClu)
    Class_n_copy[split_inds] = np.array([BestSplitClu,BestDelClu],dtype=np.int32)[np.random.randint(2,size=split_inds.size)]
    NewClass,NewScore,NewM = cem(X_nf,Recurse=False,init=('given',Class_n_copy),GraphInfo=GraphInfo)
    if NewScore > Score:
        tprint("Tunneling successful. Deleting %i and spliting %i"%(BestDelClu,BestSplitClu),colorstring=bcolors.HEADER)
        return NewClass,NewScore,NewM,True    
    else:
        tprint("Tunneling is no good")
        return Class_n,Score,M,False

def compute_LogP(X_nf, Mean_mf, Weight_m, InvSqrtCov_mff):
    "Find likelihood of each point in each cluster"
    N = X_nf.shape[0]
    M,F = Mean_mf.shape
    
    LogP = np.zeros((N,M),dtype=np.float32)
    for m in xrange(M):
        Vec2Mean_nf = Mean_mf[m] - X_nf
        LogInvSqrtDet = np.log(InvSqrtCov_mff[m].diagonal()).sum()        
        Mahal_n = sqnorm(InvSqrtCov_mff[m]*Vec2Mean_nf.T)
        LogP[:,m] = - Mahal_n/2 + LogInvSqrtDet + np.log(Weight_m[m]) - np.log(2*np.pi)*F/2 ###+ logrootdet because we're using the inverse
    return LogP     

#def logp_helper(tup):
    #(X_nf,Mean_f,Weight,InvSqrtCov_ff) = tup
    #F = X_nf.shape[1]
    #Vec2Mean_nf = Mean_f - X_nf
    #LogInvSqrtDet = np.log(InvSqrtCov_ff.diagonal()).sum()
    #Mahal_n = sqnorm(InvSqrtCov_ff*Vec2Mean_nf.T)
    #return -Mahal_n/2 + LogInvSqrtDet + np.log(Weight) - np.log(2*np.pi)*F/2

#def compute_LogP_par(X_nf, Mean_mf, Weight_m, InvSqrtCov_mff,pool):
    #N = X_nf.shape[0]
    #M,F = Mean_mf.shape
    #return np.vstack(pool.map(logp_helper,zip(it.repeat(X_nf),Mean_mf,Weight_m,InvSqrtCov_mff))).T

    
def cluster_plot(X_nf,Class_n, Mean_mf=None,Cov_mff=None):
    "Plot different clusters as different colors, along the first two dimensions"
    global ScatterCollection,PlottedX
        
    plt.ioff()            
    ax = plt.gca()
    M = Class_n.max()+1
    class_colors = plt.cm.spectral(np.linspace(.1,.9,M))
    point_colors = class_colors[Class_n]

    
    if PlottedX is None or (PlottedX != X_nf.shape):
        plt.cla()
        PlottedX = X_nf.shape
        print "replotting everything"
        ScatterCollection = plt.scatter(X_nf[:,0],X_nf[:,1],c=point_colors)
    else:
        ScatterCollection.set_facecolor(point_colors)
    for artist in ax.artists:
        artist.remove()
        
    for m in xrange(1,M):
        if Mean_mf is not None:
            x0,y0 = Mean_mf[m,0],Mean_mf[m,1]
        if Cov_mff is not None:
            eigvals,eigvecs = np.linalg.eig(Cov_mff[m])
            len0,len1 = np.sqrt(eigvals[0]),np.sqrt(eigvals[1])
            plt.arrow(x0,y0,len0*eigvecs[0,0],len0*eigvecs[1,0],ec=class_colors[m],lw=5)
            plt.arrow(x0,y0,len1*eigvecs[0,1],len1*eigvecs[1,1],ec=class_colors[m],lw=5)
            
            plt.arrow  
    plt.ion()
    plt.draw()
    #time.sleep(.5)

def mybasic_cluster(Fet_nc3,n_starts =5):
    "Cluster with usual EM, i.e. no graph"
    global VERBOSE
    VERBOSE = 1
    n_spikes,n_ch,fpc = Fet_nc3.shape
    ScoreList,Class_nList = [],[]
    print("clustering using EM...")
    for i in xrange(n_starts):
        print("try %i"%i)
        Class_n,Score,M = cem(Fet_nc3.reshape(n_spikes,n_ch*fpc).astype(np.float32),Recurse=True,init=('channel',n_ch),GraphInfo=None)
        Class_nList.append(Class_n)
        ScoreList.append(Score)
        print("score: %.1f"%Score)
        
    return Class_nList[np.argmax(ScoreList)]
    
def cluster_withgraph(X_nf,Graph):
    "EM with a graph"
    GraphInfo = estimate_cov.graph_info(ChannelGraph,1)
    return cem(X_nf,Recurse=True,init='random',GraphInfo=GraphInfo)

def write_clu(Class_n,CluFilePath):
    "Write cluster file"
    CluFile = open( CluFilePath,'w')
    CluFile.write( '1000\n')
    np.savetxt(CluFile,Class_n,fmt="%i")

def min_update(new_val, new_arg, best_val, best_arg):
    "update arg if new_val is less than best_val"
    return (new_val, new_arg) if new_val < best_val else (best_val, best_arg)
def max_update(new_val, new_arg, best_val, best_arg):
    return (new_val, new_arg) if new_val > best_val else (best_val, best_arg)
    
    
    
class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'

    def disable(self):
        self.HEADER = ''
        self.OKBLUE = ''
        self.OKGREEN = ''
        self.WARNING = ''
        self.FAIL = ''
        self.ENDC = ''    
    
def tprint(String,verb_level = 1,colorstring = ''):
    "Print string with appropriate number of tabs"
    if verb_level > VERBOSE: return
    if colorstring != '':
        print('\t'*TAB_LEVEL+colorstring+str(String)+bcolors.ENDC)
    else:
        print('\t'*TAB_LEVEL+str(String))
      
def a_little_under(x):
    return x*1.00001 if x<0 else x*.99999    
    
    

    
def test_fake():
    X1 = 30*np.random.random((1000,2)).astype(np.float32)
    X2 = 30*np.random.random((1000,2)).astype(np.float32) + 40
    X3 = 30*np.random.random((1000,2)).astype(np.float32) -40
    X = np.concatenate((X1,X2,X3),axis=0)

    Sd_f = np.sqrt(np.var(X,axis=0))
    X /= Sd_f    
    Class_n = (1+np.random.randint(9,size=3000)).astype(np.int32)
    cem(X,Recurse=True,init=('given',Class_n))    
        
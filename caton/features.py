from __future__ import division, with_statement
import numpy as np
import os
try:
    import json
except Exception:
    import simplejson as json
import output

feature_path = os.path.join(os.path.split(__file__)[0],"data/features.txt")
ts_path = os.path.join(os.path.split(__file__)[0],"data/timeseries.txt")


def compute_pcs(X_ns):
    Cov_ss = np.cov(X_ns.T.astype(np.float64))
    Vals,Vecs = np.linalg.eigh(Cov_ss)

    return Vecs.astype(np.float32).T[np.argsort(Vals)[::-1]]

def save_feature_info(X_ns,s_before,s_after,sample_rate):
    PC_fs = compute_pcs(X_ns)
    TS_s = np.arange(-s_before,s_after)/sample_rate
    np.savetxt(feature_path,PC_fs)
    np.savetxt(ts_path,TS_s)

def make_features_from_spk(SpkFileName,n_ch=None,s_before=None,s_after=None,sample_rate=None):
    if n_ch == None:
        SpkDir = os.path.dirname(SpkFileName)
        with open(os.path.join(SpkDir,"params.json"),"r") as fd: params = json.load(fd)
        sample_rate = float(params["SAMPLE_RATE"])
        s_before = params["T_BEFORE"]*sample_rate
        s_after = params["T_AFTER"]*sample_rate
        n_ch = params["N_CH"]
        s_before = params["T_BEFORE"]*sample_rate
        
    X_ns = np.fromfile(SpkFileName,dtype=np.int16).reshape(-1,s_before+s_after,n_ch)[:,:,0]
    save_feature_info(X_ns,s_before,s_after,sample_rate)

def matnorm2(X_ns):
    return np.dot(X_ns.T,X_ns)

def sumofsquarediffs(X_ns):
    return 3*matnorm2(X_ns) - matnorm2(X_ns.sum(axis=0).reshape(1,-1))
    
def make_features_from_spk_and_clu(SpkFileName,CluFileName,F,n_ch = None):
    if n_ch == None:
        SpkDir = os.path.dirname(SpkFileName)
        print os.path.join(SpkDir,"params.json")
        with open(os.path.join(SpkDir,"params.json"),"r") as fd: params = json.load(fd)
        sample_rate = float(params["SAMPLE_RATE"])
        s_before = params["T_BEFORE"]*sample_rate
        s_after = params["T_AFTER"]*sample_rate
        n_ch = params["N_CH"]
        s_before = params["T_BEFORE"]*sample_rate
    
    X_ns = np.fromfile(SpkFileName,dtype=np.int16).reshape(-1,s_before+s_after,n_ch)[:,:,0]
    CluArr = output.read_clu(CluFileName)
    n_clu = CluArr.max()+1

    X_ns -= X_ns.mean(axis=0)
    Q_ss = sumofsquarediffs(X_ns[CluArr > 0])
    for i_clu in xrange(1,n_clu):
        Q_ss -= sumofsquarediffs(X_ns[CluArr == i_clu])
    
    Vals,Vecs = np.linalg.eigh(Q_ss)
    print Vals,Vecs

    PC_fs = Vecs.astype(np.float32).T[np.argsort(Vals)[::-1]]
    TS_s = np.arange(-s_before,s_after)/sample_rate
    np.savetxt(feature_path,PC_fs)
    np.savetxt(ts_path,TS_s)
        
    
def get_features(s_before,s_after,sample_rate,F):
    Feats_fs = np.loadtxt(feature_path)
    TS_s = np.loadtxt(ts_path)
    return np.array([np.interp(np.arange(-s_before,s_after,dtype=np.float32)/sample_rate,TS_s,Feat_s) for Feat_s in Feats_fs]).astype(np.float32)[:F]

def plot_features(F=None,savefig=False,output_dir=None):
    import matplotlib.pyplot as plt
    TS_s = np.loadtxt(ts_path)
    Feats_fs = np.loadtxt(feature_path)[:F]
    for Feat_s in Feats_fs:
        plt.plot(TS_s,Feat_s)
    plt.legend([str(i) for i in range(F)])
    if savefig:
        outpath = os.path.join(output_dir or ".","features.png")
        plt.savefig(outpath)
    else:
        plt.show()

if __name__ == "__main__":
    make_features_from_spk_and_clu("../doc/generated_data/d11221/d11221.002_just_tet_batch/d11221.002_just_tet_batch.spk.1",
                                   "../doc/generated_data/d11221/d11221.002_just_tet_batch/d11221.002_just_tet_batch.clu.1",
                                   3)
    #make_features_from_spk("../doc/generated_data/d11221/d11221.002_just_tet_batch/d11221.002_just_tet_batch.spk.1")
    #print np.linalg.eigvals(sumofsquarediffs(np.array([[1,0],[0,1],[1,1],[6,5]])))
    

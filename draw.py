from matplotlib import pyplot as plt
import numpy as np
import re
import pandas as pd

class ZeroDenomError(RuntimeError):
  pass

class EstimatorNotFoundError(RuntimeError):
  pass

# USAGE:
# >>> draw.draw()     draw.draw(fname,esimator,rolling) 
#       *Args:
#          fname: str, filename of FCIMCStats, fciqmc_stats
#          estimator: str, "projE" or "trial"
#          rolling: rolling window to smooth the curve
#
# >>> draw.diagnostic()

# ver 1.1
# --- fixed bug: # of projE(1) is not fixed in fciqmc_stats

# ver 1.2
# --- added draw_all() for mneci

# ver 1.3
# --- using Trial WF
# --- in draw: plot HF/Tiral num and denom at the same fig.

# ver 1.3.1
# --- draw.draw(1), and fix 0 denominator
# --- draw.draw_all() can be plotted while allowing python operations

# ver 2.0
# --- draw proj. E with "projE columns"
# --- adding rolling to smooth energy plot
# --- draw_all() is depricated, use draw() instead
# TODO: distinguish continue calc

def draw(fname=None,estimator="trial",rolling=500,cutoff=0):
  # Args: 
  #   fname
  # read in data
  if fname is None:
    for fname in ["fciqmc_stats","FCIMCStats","FCIMCStats2"]:
      try:   
        step, data = get_data(fname, estimator=estimator)
        break
      except FileNotFoundError as e:
        continue
  else:
    try:   
      step, data = get_data(fname, estimator=estimator)
    except FileNotFoundError as e:
      print("please check the input file")

  # draw
  fig, ax = plt.subplots(1,1)
  for s,d,i in zip(step,data,1+np.arange(len(data),dtype=int) ):
    if cutoff > s[0]:
      d = d[s > cutoff]
      s = s[s > cutoff]
    d = pd.Series(d).rolling(window=rolling).mean()
    ax.plot(s,d,label="state"+str(i))
  ax.legend()
  ax.set_title("step-E: "+estimator)
  ax.set_xlabel("STEP")
  ax.set_ylabel("ENERGY (t)")
  plt.show()

draw_all=draw



def diagnostic(fname=None,index=None):
  def get_denom_and_numer(fname):
    with open(fname,'r') as f:
      header = f.readline()
      if len(re.findall("(Step|Iter)",header)) == 0:
        header = f.readline()
    str_ind_denom = re.findall('(\d+)\.\s*Trial[^0-9]*Denom',header)
    str_ind_numer = re.findall('(\d+)\.\s*Trial[^0-9]*Num',header)
    ind_denom = np.array([int(i)-1 for i in str_ind_denom],dtype=int)
    ind_numer = np.array([int(i)-1 for i in str_ind_numer],dtype=int)
    trial_denom_collection = np.loadtxt(fname,usecols=ind_denom).T if len(ind_denom)>0 else []
    trial_numer_collection = np.loadtxt(fname,usecols=ind_numer).T if len(ind_numer)>0 else []
    
    str_ind_denom = re.findall('(\d+)\.\s*ProjE[^0-9]*Denom',header)
    str_ind_numer = re.findall('(\d+)\.\s*ProjE[^0-9]*Num',header)
    ind_denom = np.array([int(i)-1 for i in str_ind_denom],dtype=int)
    ind_numer = np.array([int(i)-1 for i in str_ind_numer],dtype=int)
    proj_denom_collection = np.loadtxt(fname,usecols=ind_denom).T if len(ind_denom)>0 else []
    proj_numer_collection = np.loadtxt(fname,usecols=ind_numer).T if len(ind_numer)>0 else []
    return trial_denom_collection, trial_numer_collection, proj_denom_collection, proj_numer_collection
  

  if fname is None:
    for fname in ["fciqmc_stats","FCIMCStats","FCIMCStats2"]:
      try:
        trial_denom_collection, trial_numer_collection, proj_denom_collection, proj_numer_collection = get_denom_and_numer(fname)
        break
      except FileNotFoundError as e:
        continue
  else:
    try:   
      trial_denom_collection, trial_numer_collection, proj_denom_collection, proj_numer_collection = get_denom_and_numer(fname)
    except FileNotFoundError as e:
      print("please check the input file")

  if index is not None:
    trial_denom_collection = [ trial_denom_collection[index] ]
    trial_numer_collection = [ trial_numer_collection[index] ]
    proj_denom_collection = [ proj_denom_collection[index] ]
    proj_numer_collection = [ proj_numer_collection[index] ]


  
  step = np.loadtxt(fname,usecols=0)
  fig = plt.figure()
  ax1, ax2 = fig.subplots(2,1)
  ax1p = ax1.twinx()
  ax2p = ax2.twinx()

  for i, trial_numer in enumerate(trial_numer_collection):
    ax1.plot(step,trial_numer,color='C'+str(i),label='trial: '+str(i) )
  for i, proj_numer in enumerate(proj_numer_collection):
    ax1p.plot(step,proj_numer,'--',color='C'+str(i),label='proj: '+str(i) )
  ax1.set_title('numer',loc='left')
  ax1.legend(loc='upper left')
  ax1p.legend(loc='upper right')
  #ax1.set_xlabel("STEP")
  ax1.set_xticklabels([])
  ax1p.set_xticklabels([])
  ax1.set_ylabel("(solid) <Psi_T|H|Psi>",labelpad=8.0)
  ax1p.set_ylabel("(dashed) <HF|H|Psi>",labelpad=8.0)

  for i, trial_denom in enumerate(trial_denom_collection):
    ax2.plot(step,trial_denom,color='C'+str(i),label='trial: '+str(i) )
  for i, proj_denom in enumerate(proj_denom_collection):
    ax2p.plot(step,proj_denom,'--',color='C'+str(i),label='proj: '+str(i) )
  ax2.set_title('denom',loc='left')
  ax2.legend(loc='upper left')
  ax2p.legend(loc='upper right')
  ax2.set_xlabel("STEP")
  ax2.set_ylabel("(solid) <Psi_T|Psi>",labelpad=8.0)
  ax2p.set_ylabel("(dashed) <HF|Psi>",labelpad=8.0)

  plt.show()
 
   
def get_data(fname, estimator='trial'):
  # --- get the column index ---
  # Returns:
  #   step: step where energy is not NaN
  #   data: energy, shaped (nstep,nreplica)
  #         for neci, shaped (nstep,1)

  # NECI: FCIMCStats
  #   projE (HF estimator): Tot-Proj.E.ThisCyc
  #   trialE (trial estimator): TrialNumerator/TrialDenom 
  # MNECI: fciqmc_stats
  #   projE (HF): Tot ProjE
  #               this should equal to 
  #                 ProjE Num/ProjE Denom + Shift
  #   trialE (trial): TrialE Num/TrialE Denom
  # DNECI: FCIMCStats & FCIMCStats2
  #   same as NECI
  
  
  def _get_neci_data(fname, estimator):
    with open(fname,'r') as f:
      header = f.readline()
      if len(re.findall("(Step|Iter)",header)) == 0:
        header = f.readline()
    if estimator in ["trial"]:
      str_ind_denom = re.findall('(\d+)\.\s*Trial.*Denom',header)
      str_ind_numer = re.findall('(\d+)\.\s*Trial.*Num',header)
      if len(str_ind_denom) == 0:
        raise EstimatorNotFoundError("No TrialE Denom found. Please check the FCIMCStats file and pyblock manually")
      ind_denom = int(str_ind_denom[0])-1
      ind_numer = int(str_ind_numer[0])-1
      denom = np.loadtxt(fname,usecols=ind_denom)
      numer = np.loadtxt(fname,usecols=ind_numer)
      # exclude 0 denom ...
      step = np.loadtxt(fname,usecols=0)
      if np.allclose(denom,0):
        raise ZeroDenomError("All Trial Denom is 0. Please try draw(estimator='projE')")
      step = step[np.isclose(denom,0) == False]
      numer = numer[np.isclose(denom,0) == False]
      denom = denom[np.isclose(denom,0) == False]
      data = numer/denom
      return step, data
    elif estimator in ["projE",'proj E','proj','HF']:
      str_ind_projE = re.findall('(\d+)\.\s*Tot[^0-9]*Proj.*ThisCyc',header)
      print("Please be careful, now using Proj.E (HF) estimator\n"+re.findall('\d+\.Tot[^0-9]*Proj.*ThisCyc',header)[0])
      ind_projE = int(str_ind_projE[0])-1
      data = np.loadtxt(fname,usecols=ind_projE)
      step = np.loadtxt(fname,usecols=0)
      return step,data

  def _get_mneci_data(fname, estimator):
    with open(fname,'r') as f:
      header = f.readline()
      if len(re.findall("(Step|Iter)",header)) == 0:
        header = f.readline()
    if estimator in ["trial"]:
      str_ind_denom = re.findall('(\d+)\. TrialE Denom',header)
      str_ind_numer = re.findall('(\d+)\. TrialE Num',header)
      if len(str_ind_denom) == 0: 
        raise EstimatorNotFoundError("No TrialE Denom found. Please use draw(estimator='projE') or check the fciqmc_stats file and pyblock manually")
      ind_denom = np.array([int(i)-1 for i in str_ind_denom],dtype=int)
      ind_numer = np.array([int(i)-1 for i in str_ind_numer],dtype=int)
      denom = np.loadtxt(fname,usecols=ind_denom)
      numer = np.loadtxt(fname,usecols=ind_numer)
      # exclude 0 denom ... can be different for replicas
      step = np.loadtxt(fname,usecols=0)
      if np.allclose(denom,0):
        raise ZeroDenomError("All Trial Denom is 0. Please try dray(estimator='projE')")
      nreplica = len(denom)
      step_collection = []
      data_collection = []
      for i in range(nreplica):
        step_i,numer_i,denom_i = step[:,i], numer[:,i], denom[:,i]
        step_i = step_i[np.isclose(denom_i,0) == False] 
        numer_i = numer_i[np.isclose(denom_i,0) == False]
        denom_i = denom_i[np.isclose(denom_i,0) == False]
        data_i = numer_i/denom_i
        step_collection = step_collection.append(step_i)
        data_collection = data_collection.append(data_i)
      return step_collection, data_collection
    elif estimator in ["projE",'proj E','proj','HF']:
      str_ind_projE = re.findall('(\d+)\.\s*Tot[^0-9]*ProjE',header)
      print("Please be careful, now using Proj.E (HF) estimator\n"+re.search('\d+\.\s*Tot[^0-9]*ProjE',header)[0])
      ind_projE = np.array([int(i)-1 for i in str_ind_projE],dtype=int)
      nreplica = len(ind_projE)
      data = np.loadtxt(fname,usecols=ind_projE).T
      step = np.loadtxt(fname,usecols=0)
      step = np.tile(step,(nreplica,1))
      return step,data



  try:
    if fname.split('/')[-1] in ["FCIMCStats","FCIMCStats2"]:
      step,data = _get_neci_data(fname, estimator=estimator)
      return [step],[data]
    elif fname.split('/')[-1] in ["fciqmc_stats"]:
      return _get_mneci_data(fname, estimator=estimator)
  except ZeroDenomError as e:
    print("All Trial Denom is 0. Please try dray(estimator='projE')")
    exit()
  except EstimatorNotFoundError as e:
    print("No TrialE Denom found. Please check the FCIMCStats or fciqmc_stats file and pyblock manually")
    exit()


        
 
def main():
  try: 
    draw_all()
  except:
    draw()
    input()

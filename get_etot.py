import re
import numpy as np
import pandas as pd
import pyblock

class ZeroDenomError(RuntimeError):
  pass

class EstimatorNotFoundError(RuntimeError):
  pass


# --- ver 1.0 ------

# --- ver 1.1 -----
# blocking with trial wavefunction 

# --- ver 1.2 ----
# adding etot_kernel for reading FCIMCStats2

# --- ver 2.0 ---
# dealing with continued NECI fciqmc_stas
# using ndarray instead of pandas
# formated output in etot(), for data, use get_etot()

# --- ver 2.0.1 ---
# fixed support for FCIMCStats

# --- ver 2.0.2 ---
# adding *arg "proj" for HF energy estimator

# --- ver 2.1.0 ---
# apply re to numer and denom

# --- ver 3.0 ---
# fixed projection energy, with offset
# re-structured the codes
# using pandas
# cutoff should be chosen manually

def etot(fname=None,cutoff=0,estimator="trial"):
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

  block_data = pd.DataFrame(columns=['start','energy','err','converged'])
  for s,d in zip(step,data):
    if cutoff > s[0]:
      d = d[s>cutoff]
      s = s[s>cutoff]
    start = max(s[0],cutoff)

    _,r,_ = pyblock.pd_utils.reblock(pd.Series(d))
    b = pyblock.pd_utils.reblock_summary(r)
    if not b.empty:
      e = b['mean'].data
      err = b['standard error'].data
      conv = True
    else:
      e = r.data['mean'].loc[0]
      err = r.data['standard error'].loc[0]
      conv = False
    block_data.loc[len(block_data)] = [start, e, err, conv]

  print('================')
  print('{:<12}'.format('start:')
        +' '.join('{:<20d}'.format(int(s)) for s in block_data['start'] ))
  print('{:<12}'.format('energy:')
        +' '.join('{:<20}'.format(e) for e in block_data['energy']))
  print('{:<12}'.format('std_err:')
        +' '.join('{:<20}'.format(err) for err in block_data['err']))
  print('{:<12}'.format('converged:')
        +' '.join('{:<20s}'.format(str(c)) for c in block_data['converged']))
  print('\n')

       
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
      step = np.loadtxt(fname,usecols=0,dtype=int)
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

   

def blocking(data):
    # input: ndarray, [# of step, # of mneci states]

    # pyblock handle 1d array and 2d array differenetly 
    # 1D: [...] -> E (float), ...
    # 2D: [...,0] [...,1] ...
    #     and (...,0) shaped ndarray is not supported :(
    #     -> E ndarry 
    E = np.array([])
    Err = np.array([])
    Conv = np.array([],dtype=bool)
    nstep, nstate = data.shape
    if nstate != 1:
        block_data = pyblock.blocking.reblock(data.T)
        ind_opt = pyblock.blocking.find_optimal_block(nstep, block_data)
        for i,opt in enumerate(ind_opt): 
            # rslt is a list (blocking1, blocking2, ...)
            # each contains blocking results
            if np.isnan( opt ):
                E = np.append(E,block_data[-1].mean[i])
                Err = np.append(Err,block_data[-1].std_err[i])
                Conv = np.append(Conv,False)
            else:
                E = np.append(E,block_data[ind_opt[i]].mean[i])
                Err = np.append(Err,block_data[ind_opt[i]].std_err[i])
                Conv = np.append(Conv,True)
    else: # nstate = 1, for FCIMCstats input
        block_data = pyblock.blocking.reblock(data[:,0])
        ind_opt = pyblock.blocking.find_optimal_block(nstep, block_data)
        for i,opt in enumerate(ind_opt): 
            # rslt is a list (blocking1, blocking2, ...)
            # each contains blocking results
            if np.isnan( opt ):
                E = np.append(E,block_data[-1].mean)
                Err = np.append(Err,block_data[-1].std_err)
                Conv = np.append(Conv,False)
            else:
                E = np.append(E,block_data[ind_opt[i]].mean)
                Err = np.append(Err,block_data[ind_opt[i]].std_err)
                Conv = np.append(Conv,True)

    return {'energy':E, 'err':Err, 'converged':Conv}



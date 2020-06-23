import numpy as np
import pandas as pd
import inspect
import scipy
import scipy.stats
from pathlib import Path
import olga3
import olga3.load_model
import olga3.generation_probability
import mpmath


translate_header = {
    'sequenceStatus': 'frame_type',
    'productive': 'frame_type',
    'aminoAcid': 'amino_acid',
    'vGeneName': 'v_gene',
    'v_call': 'v_gene',
    'jGeneName': 'j_gene',
    'j_call': 'j_gene',
    'rearrangement': 'rearrangement',
    'aminoAcid': 'cdr3_aa',
    'cdr3_aa': 'cdr3_aa',
    'cdr3s_aa': 'cdr3_aa',
    'cdr3': 'cdr3_nt',
    'nSeqCDR3': 'cdr3_nt',
    'aaSeqCDR3': 'cdr3_aa',
    'cdr3s_nt': 'cdr3_nt',
    'junction_aa': 'cdr3_aa',
    'CDR3.amino.acid.sequence': 'cdr3_aa',
    'CDR3.nucleotide.sequence': 'cdr3_nt',
    'bestVGene': 'v_gene',
    'bestJGene': 'j_gene',
    'count': 'read_count',
    'consensus_count': 'read_count',
    'Read.count': 'read_count',
    'reads': 'read_count'}


def rename_cols(df):
    """
    Standardize the column names of the dataframe
    """
    columns = df.columns.values
    new_cols = []
    for col in columns:
        if col in translate_header:
            new_cols.append(translate_header[col])
        else:
            new_cols.append(col)
    df.columns = new_cols
    return df


def load_olga_model(model_directory, D=True):
    """ Fully load an olga model.
        model: directory where the model is stored
        D: True if the model is VDJ (default)
    """
    params_file_name = model_directory / 'model_params.txt'
    marginals_file_name =  model_directory / 'model_marginals.txt'
    V_anchor_pos_file =  model_directory / 'V_gene_CDR3_anchors.csv'
    J_anchor_pos_file =  model_directory / 'J_gene_CDR3_anchors.csv'

    if D:
        genomic_data = olga3.load_model.GenomicDataVDJ()
        generative_model = olga3.load_model.GenerativeModelVDJ()
    else:
        genomic_data = olga3.load_model.GenomicDataVJ()
        generative_model = olga3.load_model.GenerativeModelVJ()
        
    genomic_data.load_igor_genomic_data(params_file_name, V_anchor_pos_file, J_anchor_pos_file)
    generative_model.load_and_process_igor_model(marginals_file_name)
    
    if D:
        pgen_model = olga3.generation_probability.GenerationProbabilityVDJ(generative_model, genomic_data)
    else:
        pgen_model = olga3.generation_probability.GenerationProbabilityVJ(generative_model, genomic_data)

    return genomic_data, generative_model, pgen_model


# load pgen model for TRA & TRB
local_directory = Path(inspect.getfile(olga3)).parent 
pgen = {}
_, _, pgen['TRA'] = load_olga_model(
    model_directory=local_directory/'default_models'/'human_T_alpha',
    D=False)
_, _, pgen['TRB'] = load_olga_model(
    model_directory=local_directory/'default_models'/'human_T_beta',
    D=True)


def compute_pgen(seq, full=False):
    """
    Compute the pgen of a CDR3 nucleotide sequence.
    If full, the sequence must be written in the format: "TRA:TGT...TT;TRB:TGT...GC"
    """
    if full:
        seqa = [u[4:] for u in seq.split(";") if u.startswith("TRA:")][0]
        seqb = [u[4:] for u in seq.split(";") if u.startswith("TRB:")][0]
        return pgen['TRA'].compute_nt_CDR3_pgen(seqa) * pgen['TRB'].compute_nt_CDR3_pgen(seqb)
    return pgen['TRB'].compute_nt_CDR3_pgen(seq)


def boundI1(r, meanS1, µp, σp, maxM=None):
    """ return an upper bound for P[I(A, B) < r | A&B autologous] """
    λ = meanS1 # parameter of the poisson law
    if maxM is None:
        maxM = int(5 * meanS1)
    ks = range(1, maxM)
    poiss = [scipy.stats.poisson.pmf(k, λ) for k in ks]
    norms =  [1/2*(1 + scipy.special.erf((r - μp * k)/(np.sqrt(2 * k) * σp))) for k in ks]
    if r <= 0:
        return 0.
    else:
        return sum([p * n for p, n in zip(poiss, norms)]) + np.exp(-λ)


def boundI2(r, meanS2, µp, σp, maxM=None):
    """ return an upper bound for P[I(A, B) >= r | A&B not autologous] """
    λ = meanS2 # parameter of the poisson law
    if maxM is None:
        maxM = int(5*meanS2)
    ks = range(0, maxM)
    poiss =  [scipy.stats.poisson.pmf(k, λ) for k in ks]
    norms = [1/2*(1 + scipy.special.erf((r - μp * k)/(np.sqrt(2 * k) * σp))) 
             if k > 0 else 1 for k in ks]

    if r > 0:
        return sum([p * (1-n) for p, n in zip(poiss, norms)])
    else:
        return 1.


def regularized_gamma_Q(n, mu):
    return float(mpmath.gammainc(n, a=mu, regularized=True))

    

def S_graph(m):
    """
    Return two array x, y representing the graph of S (with mean value m)
    """
    rg = np.arange(0, max(int(3*m), 200))
    return (np.array(rg),
            np.array([scipy.stats.poisson.pmf(k, m) for k in rg]))


def I_graph(m, µ, σ):
    """
    Return two array x, y
    representing the graph of I
    (with mean value m, mean
    log_pgen µ, var log_pgen σ)
    """
    rg = np.linspace(0, 3 * m * µ, 1000)
    S = [(k, scipy.stats.poisson.pmf(k, m)) for k in  np.arange(1, max(int(3*m), 200))]
    Is = [float(mpmath.fsum(1/mpmath.sqrt(2 * np.pi * σ**2 * k) *
                            mpmath.exp(-m) * mpmath.power(m, k) / mpmath.factorial(k) *
                            mpmath.exp(-(x - μ*k)**2/(2 * σ**2 * k))
                            for k, s in S if s > 1e-6))
          + (scipy.stats.poisson.pmf(0, m) if x == 0. else 0)
          for x in rg]
    return np.array(rg), np.array(Is)


def estimate_S_from_counts(dfA, dfB):
    """
    Give an unbiased estimate for the value of S given the counts. 
    dfA & dfB need to have a column "read_count" & "cdr3_nt"
    """
    cAs = dfA.groupby("cdr3_nt").read_count.sum().to_frame()
    cBs = dfB.groupby("cdr3_nt").read_count.sum().to_frame()
    cs = cAs.merge(cBs, left_index=True, right_index=True,
                   how='outer', suffixes={'_A', '_B'})
    cs = cs.fillna(0)
    cs["read_count"] = cs["read_count_A"] + cs["read_count_B"]
    M2 = len(cBs)
    M1 = len(cAs)
    # The formula below remove the bias associated with sampling
    Sest = sum([1 + (mpmath.binomial(0, int(c)) 
                  - mpmath.binomial(M2, int(c))
                  - mpmath.binomial(M1, int(c)))
             /mpmath.binomial(M1+M2, int(c))  
             for c in cs["read_count"].values if int(c) != 0])
    return float(Sest)
    

def immprint(dfA, dfB, full=False, use_counts=True, onlyS=False, max_shared=np.inf, γ=12):
    """Main immprint function
    Arguments:
    - dfA, dfB: dataframe containing a column `cdr3_nt` (and optionally a column `read_count`)
    - full: if True, use both chains alpha & beta (rather than only beta) (default False)
    - use_counts: if True, use the `read_count` column in the
      dataframe dfA & dfB to give an estimate of S
    - onlyS: if True, do not compute the statistic I
    - max_shared: If the number of sequences shared between dfA & dfB is
      larger than max_shared, do not compute I.
    Return:
    tuple: S, I, shared_sequences, parms
    parms is a dict that contains all the parameters of the theoretical distribution:
    - µS1 (mean of S when autologous), µS2 (mean of S when not autologous)
    - µ/σ_logpgen(_shared): mean and variance of the logpgen distribution in both scenarios
    - pv1/pv2: the p-values associated with the autologous and non-autologous case respectively
    - rS / rI: the threshold values for S & I
    """
    MA = dfA.cdr3_nt.nunique()
    MB = dfB.cdr3_nt.nunique()

    if full: # check that the format matches
        if (np.sum(
                dfA.cdr3_nt.apply(
                    lambda x: (len(x.split(";")) != 2) or ("TRA:" not in x) or ("TRB:" not in x)
                )) > 0 or 
            np.sum(
                dfB.cdr3_nt.apply(
                    lambda x: (len(x.split(";")) != 2) or ("TRA:" not in x) or ("TRB:" not in x)
                )) > 0):
            raise Exception("CDR3 sequence format doesn't match a full receptor. "
                            "Correct format: TRA:XXXXX;TRB:XXXXXX")

    
    shared_sequences = set(dfA.cdr3_nt.values) & set(dfB.cdr3_nt.values)
    S = len(set(dfA.cdr3_nt.values) & set(dfB.cdr3_nt.values)) 
    I = None
    if S < max_shared and (not onlyS):
        # if "full", the cdr3_nt column should have the standard form:
        # TRA:..,TRB:...
        pgens = [compute_pgen(s, full=full) for s in shared_sequences]
        # too low pgens are ignored (sequencing error / bad model)
        I = sum(-np.log(p) - γ for p in pgens if p > 1e-50)
    
    # the mean of pgen is used to estimate I, it depends on
    # the choice of model and on the scenario (homologous or
    # not)
    if full:
        µ_logpgen_shared = 33
        σ_logpgen_shared = 3
        µ_logpgen = 47
        σ_logpgen = 8
    else: 
        µ_logpgen_shared = 18 - γ
        σ_logpgen_shared = 2
        µ_logpgen = 26 - γ
        σ_logpgen = 6


    µS2 = (1e-13 if full else 3.0e-07) * MA * MB  # default value (not autologous case)
    µS1 = µS1_default = 7.96e-07 * MA * MB  # default values (autologous)
    if "read_count" in dfA.keys() and "read_count" in dfB.keys() and use_counts:
        µS1 = estimate_S_from_counts(dfA, dfB)
    elif use_counts:
        print("Frequency counts information not found, falling back to default values")


    # Probability that the number of sequences shared is lower than S
    # if both samples come from the same individual.
    pv1 = regularized_gamma_Q(S+1, µS1)
    # Probability that the number of sequences shared is higher than S
    # if both samples come from different individuals.
    pv2 = 1 - regularized_gamma_Q(S+1, µS2)
    if I is not None:
        # Probability that a random Immprit, pgen-based, estimate is lower than I
        # if both samples come from the same individual.
        pvI1 = boundI1(I, µS1, µ_logpgen, σ_logpgen)
        # Probability that a random Immprit, pgen-based, estimate is higher than I
        # if both samples come from different individuals
        pvI2 = boundI2(I, µS2, µ_logpgen_shared, σ_logpgen_shared)
        if pvI1 < pv1:
            pv1 = pvI1
        if pvI2 < pv2:
            pv2 = pvI2

    # the "threshold" choices for the "full" case are mostly arbitrary.
    # their exact value does not matter, in general.
    return (S, I,
            shared_sequences,
            {"µS1": µS1, "µS2": µS2,
             "µ_logpgen_shared": µ_logpgen_shared,
             "σ_logpgen_shared": σ_logpgen_shared,
             "µ_logpgen": µ_logpgen,
             "σ_logpgen": σ_logpgen,
             "pv1": pv1, "pv2": pv2,
             "rS": (µS1 + µS2)/2,
             "rI": (µS1 * µ_logpgen + µS2 * µ_logpgen_shared)/2})



    

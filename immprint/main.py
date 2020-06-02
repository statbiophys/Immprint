#!/usr/bin/env python3
import immprint
import pandas as pd
import argparse
import matplotlib.pyplot as plt
import seaborn as sns



def parse_arguments():
    parser = argparse.ArgumentParser(description='ImmPrint')
    parser.add_argument('sampleA', type=str,
                        help="csv file containing the sequences of the first sample")
    parser.add_argument('sampleB', type=str,
                        help="csv file containing the sequences of the second sample")
    parser.add_argument('-f', '--full', dest='full',
                        help=("If paired CDR3a / CDR3b are provided, use the full receptor "
                              "to discriminate between the two scenarios. "
                              "The CDR3 sequences should be provided in the column `cdr3_nt` (or `cdr3s_nt`)"
                              "with the standard form: \"TRA:C...F;TRB:C...F\""),
                        action='store_true')
    parser.add_argument('-S', '--no-I', dest='onlyS',
                        help=("Don't compute the \"I\" statistics (potentially more precise but slower)"),
                        action='store_true')
    parser.add_argument('-n', '--no-counts', dest='nocounts',
                        help="If specified, don't use read counts, even if provided",
                        action='store_true')

    args = parser.parse_args()

    try:
        dfA = pd.read_csv(args.sampleA, sep=None, engine='python')
    except IOError as e:
        print("Error with opening the file corresponding to sample A.")
        print(e)
        return

    try:
        dfB = pd.read_csv(args.sampleB, sep=None, engine='python')
    except IOError as e:
        print("Error with opening the file corresponding to sample B.")
        print(e)
        return

    dfA = immprint.rename_cols(dfA)
    dfB = immprint.rename_cols(dfB)

    S, I, shared_sequences, parms = immprint.immprint(dfA, dfB,
                                                      full=args.full,
                                                      use_counts=(not args.nocounts),
                                                      onlyS=args.onlyS)
    print(shared_sequences)

    print(f"Number of sequences shared between the two datasets: S = {S}")
    if not args.onlyS:
        print(f"Immprint value between the two datasets: I = {float(I):.1}")

            
    under_threshold = (S < parms['rS'] if args.onlyS else I < parms['rI'])
    ambiguous = (under_threshold and parms['pv1'] > 1e-4) or (not under_threshold and pv2 > 1e-4)

    if ambiguous:
        print("Ambiguous case: maybe not enough sequences.\n"
              f"- Probability the two samples come from the same patient <= {parms['pv1']:.2}\n"
              f"- Probability the two samples come from different patients <= {pv2:.2}\n")
    else:
        if under_threshold: # different patients
            print("The samples come from two different patients.\n"
                  f"The probability that they are autologous is lower than {parms['pv1']:.2e}")
        if not under_threshold: # same patient
            print(f"The samples come from the same patient.\n"
                  f"The probability that they were extracted in two "
                  f" different individuals is lower than {parms['pv1']:.2e}")

    ## XXXX for "full", the number of sequences needed to associate two samples to the same patient is ridicully low, so just sequencing error + chimera pcr + sequencing multiple samples together => risk of identification error increasing

    ## plotting can be quite sloooooow, deal with that
        
    fig, axes = plt.subplots((2 if I is not None else 1), 1, figsize=(5, 3))
    if I is None:
        axes = [axes]
    x1, y1 = immprint.S_graph(parms["µS1"])
     axes[0].plot(x1, y1, color=sns.color_palette()[0],
            label="Same patient")
    axes[0].axvline(S, color=sns.color_palette()[1], label="Measured")
    axes[0].axvline(parms['rS'], color='k', ls='--', label="Threshold")
    x2, y2 = immprint.S_graph(parms["µS2"])
    axes[0].plot(x2, y2, color=sns.color_palette()[1],
            label="Different patients")
    axes[0].set_ylim((0, 1.5*min(max(y2), max(y1))))
    axes[0].set_xlabel(r"$\mathcal{S}$")
    axes[0].legend()

    if not args.onlyS:
        x1, y1 = immprint.I_graph(parms["µS1"], parms["µ_logpgen"], parms["σ_logpgen"])
        axes[1].plot(x1, y1, color=sns.color_palette()[0],
                label="Same patients")
        axes[1].axvline(I, color=sns.color_palette()[1], label="Measured")
        axes[1].axvline(parms['rI'], color='k', ls='--', label="Threshold")
        x2, y2 = immprint.I_graph(parms["µS2"], parms["µ_logpgen_shared"], parms["σ_logpgen_shared"])
        axes[1].plot(x2, y2, color=sns.color_palette()[1],
                label="Different patient")
        axes[1].set_ylim((0, 1.5*min(max(y2), max(y1))))
        axes[1].set_xlabel(r"$\mathcal{I}$")
        axes[1].legend()
    plt.show()

if __name__ == "__main__":
    parse_arguments()

import os
def save_pwm_heatmap_from_probs(probs, outpath):
    try:
        import matplotlib.pyplot as plt
        import numpy as np
    except Exception as e:
        print('visualize: matplotlib required for heatmap:', e)
        return False
    order = ['A','C','G','T']
    arr = [[probs[i][j] for j in range(len(probs[0]))] for i in range(4)]
    arr = np.array(arr)
    fig, ax = plt.subplots(figsize=(max(4, arr.shape[1]*0.6), 3))
    im = ax.imshow(arr, aspect='auto', cmap='Blues', vmin=0, vmax=1)
    ax.set_yticks(range(len(order)))
    ax.set_yticklabels(order)
    ax.set_xticks(range(arr.shape[1]))
    ax.set_xticklabels([str(i+1) for i in range(arr.shape[1])])
    ax.set_xlabel('Position')
    ax.set_title('PWM (probabilities)')
    fig.colorbar(im, ax=ax, orientation='vertical', fraction=0.02)
    fig.tight_layout()
    fig.savefig(outpath, dpi=150)
    plt.close(fig)
    return True

def save_sequence_logo_from_counts(counts, outpath):
    try:
        import pandas as pd
        import matplotlib.pyplot as plt
        import logomaker
    except Exception as e:
        print('visualize: logomaker/matplotlib/pandas required for logo:', e)
        return False
    # counts: list of 4 lists (A,C,G,T) length k
    k = len(counts[0])
    df = pd.DataFrame({ 'A': counts[0], 'C': counts[1], 'G': counts[2], 'T': counts[3] })
    # convert to frequencies
    df = df.div(df.sum(axis=1), axis=0)
    plt.figure(figsize=(max(6, k*0.6), 3))
    logo = logomaker.Logo(df)
    plt.title('Sequence logo')
    plt.xlabel('Position')
    plt.ylabel('Frequency')
    plt.tight_layout()
    plt.savefig(outpath, dpi=150)
    plt.close()
    return True

from matplotlib import pyplot as plt
import numpy as np

def show_atommap(atommap):
    """
    show_atommap
    """
    coverage_per_sample = np.sum(atommap,axis=0)
    coverage_per_feature = np.sum(atommap,axis=1)
    fig = plt.figure(figsize=(12,4))
    ax = fig.add_subplot(1,3,1)
    ax.set_title('Atom map')
    ax.set_xlabel('list of CIF files')
    ax.set_ylabel('list of atoms')
    plt.imshow(atommap, aspect='auto', cmap='Greys')
    ax = fig.add_subplot(1,3,2)
    ax.set_title('coverage per atom')
    ax.set_xlabel('list of atoms')
    plt.plot(coverage_per_feature, '.')
    ax = fig.add_subplot(1,3,3)
    ax.set_title('coverage per structure')
    ax.set_xlabel('list of CIF files')
    plt.plot(coverage_per_sample, '.')
    plt.tight_layout()

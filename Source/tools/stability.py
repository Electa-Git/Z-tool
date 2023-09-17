__all__ = ['nyquist']
import numpy as np


# Ongoing

class results:
    def __init__(self, omegas):
        self.eigenvalues = {i: 3.1415 for i in range(len(omegas))}
        self.eigenvectors = {i: 3.1415 for i in range(len(omegas))}

def nyquist(L, omegas):
    1 ==1


def oscillation_modes(G, omegas):
    1 == 1


def participations(G, omegas):
    output = results(omegas)
    for idx, w in enumerate(omegas):
        eigenvalues, eigenvectors = np.linalg.eig(G[idx])
        output.eigenvalues[idx] = eigenvalues
        output.eigenvectors[idx] = eigenvectors


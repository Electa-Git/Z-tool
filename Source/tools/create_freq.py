__all__ = ['loglist']
import numpy as np

def loglist(f_min = 1, f_max = 1000, f_points = 50, f_base = 1):
    # Create the frequency perturbation vector, float by deafult, use dtype='int16' for integers
    freq = np.logspace(np.log10(f_min), np.log10(f_max), num=int(f_points))
    # Modify the list so the values are multiples of a base frequency
    for j in range(freq.size): freq[j] = freq[j] - (freq[j] % f_base)
    freq = np.unique(freq)  # Delete the repeated values
    multiples = np.arange(f_min, f_max + f_base, f_base)  # Compute the multiples in the given range
    if len(freq) < int(f_points) and len(multiples) > len(freq):  # Not enough values
        scope = np.setxor1d(np.round(freq,5), np.round(multiples,5))  # XOR operator + round to avoid float errors
        to_be_added = min(int(f_points) - len(freq), len(multiples) - len(freq))
        idx = np.floor(len(scope) / to_be_added) * np.arange(to_be_added)  # Add the values indexed uniformly
        for i in idx: freq = np.append(freq, scope[int(i)])
        freq.sort()
    return np.around(freq,8)  # Round to get rid of floating-point inaccuracies
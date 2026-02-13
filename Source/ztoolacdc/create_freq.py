"""
Function to creat a frequency vector between the minimum and the maximum which are all multiples of the base frequency
"""

__all__ = ['loglist','linlist']
import numpy as np
from itertools import combinations

def loglist(f_min = 1, f_max = 1000, f_points = 80, f_base = 1, f_exclude = None):
    # Create the frequency vector, float by deafult, use dtype='int16' for integers
    freq = np.logspace(np.log10(f_min), np.log10(f_max), num=int(f_points))
    # Modify the list so the values are multiples of a base frequency
    for j in range(freq.size): freq[j] = np.rint(freq[j] / f_base) * f_base
    freq = freq[~np.isin(freq, f_exclude)] # Remove the excluded frequencies
    freq = np.unique(freq)  # Delete the repeated values
    multiples = np.arange(f_min, f_max + f_base, f_base)  # Compute the multiples in the given range
    multiples = multiples[~np.isin(multiples, f_exclude)] # Remove the excluded frequencies
    multiples = np.unique(multiples)  # Delete the repeated values
    if len(freq) < int(f_points) and len(multiples) > len(freq):  # Not enough values
        scope = np.setxor1d(np.round(freq,5), np.round(multiples,5))  # XOR operator + round to avoid float errors
        to_be_added = min(int(f_points) - len(freq), len(multiples) - len(freq))
        idx = np.floor(len(scope) / to_be_added) * np.arange(to_be_added)  # Add the values indexed uniformly
        for i in idx: freq = np.append(freq, scope[int(i)])
        freq.sort()
    return np.round(freq,8)  # Round to get rid of floating-point inaccuracies


def linlist(f_min = 1, f_max = 1000, f_points = 80, f_base = 1, f_exclude = None):
    # Create the frequency vector, float by deafult, use dtype='int16' for integers
    freq = np.linspace(np.log10(f_min), np.log10(f_max), num=int(f_points))
    # Modify the list so the values are multiples of a base frequency
    for j in range(freq.size): freq[j] = np.rint(freq[j] / f_base) * f_base
    freq = freq[~np.isin(freq, f_exclude)] # Remove the excluded frequencies
    freq = np.unique(freq)  # Delete the repeated values
    multiples = np.arange(f_min, f_max + f_base, f_base)  # Compute the multiples in the given range
    multiples = multiples[~np.isin(multiples, f_exclude)] # Remove the excluded frequencies
    multiples = np.unique(multiples)  # Delete the repeated values
    if len(freq) < int(f_points) and len(multiples) > len(freq):  # Not enough values
        scope = np.setxor1d(np.round(freq,5), np.round(multiples,5))  # XOR operator + round to avoid float errors
        to_be_added = min(int(f_points) - len(freq), len(multiples) - len(freq))
        idx = np.floor(len(scope) / to_be_added) * np.arange(to_be_added)  # Add the values indexed uniformly
        for i in idx: freq = np.append(freq, scope[int(i)])
        freq.sort()
    return np.round(freq,8)  # Round to get rid of floating-point inaccuracies

loglist.__doc__ = """
Function to create a list of log-spaced frequencies between the minimum and the maximum which are all multiples of the base frequency.
The function can enforce the exclusion of specific frequencies from the generated list. Parameters are self-explanatory
"""

linlist.__doc__ = """
Function to create a list of linearly-spaced frequencies between the minimum and the maximum which are all multiples of the base frequency.
The function can enforce the exclusion of specific frequencies from the generated list. Parameters are self-explanatory
"""
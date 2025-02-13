"""
Function to creat a list of log-spaced frequencies between the minimum and the maximum which are all multiples of the base frequency
"""

__all__ = ['loglist']
import numpy as np
from itertools import combinations

def loglist(f_min = 1, f_max = 1000, f_points = 80, f_base = 1):
    # Create the frequency vector, float by deafult, use dtype='int16' for integers
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
    return np.round(freq,8)  # Round to get rid of floating-point inaccuracies

def minPAPR_8freqs(f_min = 1, f_max = 1000, f_points = 80, f_base = 1):
    # Compute sampling frequency and total available frequencies
    fs = 0.5 / f_max
    available_frequencies = np.round(np.arange(f_min, f_max + f_base, f_base),5)
    used_frequencies = set()  # Track globally used frequencies
    freq_multi = []

    # Generate combinations of two frequencies for f1 and f2
    for f1, f2 in combinations(available_frequencies, 2):
        # Compute the frequencies for minimum PAPR
        f1p = fs / 2 - f1
        f2p = fs / 2 - f2
        f3 = f1 + fs / 4
        f4 = f2 + fs / 4
        f3p = fs / 2 - f3
        f4p = fs / 2 - f4
        
        # Check if all derived frequencies are distinct, within range, and unused
        group = [f1, f2, f1p, f2p, f3, f4, f3p, f4p]
        if (all(f_min <= freq <= f_max for freq in group) and
                len(set(group)) == 8 and
                not used_frequencies.intersection(group)):
            freq_multi.append(group)  # Add them to the group list
            used_frequencies.update(group)  # Mark them as used

        # Stop if M groups are found
        if len(freq_multi) >= int(round(f_points/8)):
            break

    return np.round(freq_multi,8)  # Round to get rid of floating-point inaccuracies
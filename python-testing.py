# This file is used to check stuff with julia code in notebook

from scipy import signal
import numpy as np

kernel = np.array([
    [0.05, 0.1, 0.05],
    [0.1, 0.5, 0.1],
    [0.05, 0.1, 0.05]
])

matrix = np.array([
    [1, 1, 1, 1, 1],
    [1, 2, 3, 2, 1],
    [1, 3, 5, 3, 1],
    [1, 2, 3, 2, 1],
    [1, 1, 1, 1, 1],
])



# print(signal.fftconvolve(kernel, matrix))
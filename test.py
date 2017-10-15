import numpy as np 
import matplotlib.pyplot as plt  
import scipy.interpolate
from scipy.ndimage.interpolation import rotate
from insert import insert

ph = np.load("phantom_small.npy")

ft_ph = np.fft.fftshift(np.fft.fft2(ph))

r = ft_ph.real
i = ft_ph.imag
a = ft_ph.angle

plt.plot(ft_ph[ft_ph.shape[0]//2, :])
plt.show()
import numpy as np 
import matplotlib.pyplot as plt  
import scipy.interpolate
from scipy.ndimage.interpolation import rotate
from insert import insert

origin = np.load("phantom_small.npy")

ft_ph = np.fft.fft2(origin)

noise = np.random.normal(loc = 1, scale = 10, size = origin.shape)
# noise = np.fft.fft2(noise)

ft_ph = np.abs(ft_ph) * np.exp( 1j* np.pi) 

# r = ft_ph.real
# i = ft_ph.imag
# a = ft_ph.angle

# plt.plot(ft_ph[ft_ph.shape[0]//2, :])
# plt.show()

recon = np.fft.ifft2(ft_ph).real

imgs = [origin, recon]

fig, axes = plt.subplots(ncols = len(imgs))

for i, ax in enumerate(axes):
    ax.imshow(imgs[i])

plt.show()


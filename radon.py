import numpy as np 
import matplotlib.pyplot as plt  
import scipy.interpolate
from scipy.ndimage.interpolation import rotate
from insert import insert

phantom = np.load("phantom_small.npy")
# X, Y = np.mgrid[-5*np.pi:5*np.pi:64*1j, -5*np.pi:5*np.pi:64*1j]

# y1 = np.sin(X)
# y2 = np.sin(X/np.sqrt(2) + Y/np.sqrt(2))

# phantom = y1 + y2

# phantom = np.pad(phantom, pad_width = [[100,101],[100,101]], mode = "constant", constant_values = 0)


def low_pass_filter(img):
    n0, n1 = img.shape

    x, y = np.mgrid[0:n0, 0:n1]
    x -= n0//2
    y -= n1//2
    r = np.sqrt(x*x + y*y)

    fil = np.fft.ifftshift(np.where(r < min(n0, n1)/3, 1, 0))

    filtered_img = np.fft.ifft2( np.fft.fft2(img) * fil ).real
    return filtered_img


phantom = low_pass_filter(phantom)


def rotate_img(img, theta):
    theta *= np.pi/180
    n0, n1 = img.shape
    y_1d = np.arange(0, n0)
    x_1d = np.arange(0, n1)
    f = scipy.interpolate.interp2d(x_1d, y_1d, img, fill_value = 0)

    rotate = lambda x,y : [x*np.cos(theta) - y*np.sin(theta), +x*np.sin(theta) + y*np.cos(theta)]
    xx, yy = np.mgrid[0:n0,0:n1]

    tmp = xx[:,:]
    xx = yy[:,:] - n1//2
    yy = tmp - n0//2
    rotated_frame = np.array([rotate(x, y) for x,y in zip(xx.flatten(),yy.flatten())])

    rotated_frame[:,0] += n1//2
    rotated_frame[:,1] += n0//2
    rotated_img = np.array([f(v[0], v[1]) for v in rotated_frame])
    
    return rotated_img.reshape([n0,n1])

def slices_gen(img, thetas):
    l_max = np.max(img.shape)
    slices = np.zeros([len(thetas), l_max], dtype = np.float64)
    for i, theta in enumerate(thetas):
        rotated_img = rotate(img, -theta, reshape = False, prefilter = False)
        mid = rotated_img[rotated_img.shape[0]//2, :]
        pad_width = [int(np.floor(l_max - len(mid))/2.), int(np.ceil((l_max - len(mid))/2.))]
        mid = np.pad(mid, pad_width, mode = "constant", constant_values = 0)
        slices[i] = mid
    return slices



def radon(img, thetas):
    l_max = np.max(img.shape)
    sinogram = np.zeros([len(thetas), l_max], dtype = np.float64)
    for i, theta in enumerate(thetas):
        rotated_img = rotate(img, -theta, reshape = False, prefilter = False)
        proj = np.sum(rotated_img, axis = 0)
        pad_width = [int(np.floor(l_max - len(proj))/2.), int(np.ceil((l_max - len(proj))/2.))]
        proj = np.pad(proj, pad_width, mode = "constant", constant_values = 0)
        sinogram[i] = proj
    return sinogram.T

def iradon(sinogram, thetas):
    n0, n1 = sinogram.shape
    assert(n1 == len(thetas))
    thetas = thetas / 180 * np.pi
    fil = np.abs(np.fft.fftfreq(n0)).reshape((-1, 1))
    filtered_sino = np.fft.ifft(np.fft.fft(sinogram, axis = 0) * fil, axis = 0).real

    x, y = np.mgrid[0:n0,0:n0]
    x -= n0//2
    y -= n0//2
    t = np.arange(n0) - n0//2

    recon = np.zeros([n0,n0], dtype = np.float64)
    for i,theta in enumerate(thetas):
        tt = -x*np.sin(theta) + y*np.cos(theta)
        back_proj = np.interp(tt, t, filtered_sino[:, i])
        recon += back_proj
    return recon

def recon_f(sinogram, thetas):
    n0, n1 = sinogram.shape

    x, y = np.mgrid[0:n0, 0:n0]
    x -= n0//2
    y -= n1//2
    r = np.sqrt(x*x + y*y)

    fil = np.fft.ifftshift(np.where(r < n0/3, 1, 0))
    # fil = 1 #np.abs(np.fft.fftfreq(n0)).reshape((-1, 1))
    ft_sino = np.fft.fftshift( np.fft.fft(sinogram, axis = 0), axes = 0 )
    slices_real = np.real(ft_sino).T
    slices_imag = np.imag(ft_sino).T

    model_real = insert(slices_real, thetas)
    model_imag = insert(slices_imag, thetas)

    model = model_real + 1j*model_imag
    # model *= fil

    recon = np.fft.ifft2(np.fft.ifftshift(model)).real

    return model, recon


ft_origin = np.fft.fftshift(np.fft.fft2(phantom))

thetas = np.arange(0, 180, 1)
thetas = np.asarray(thetas, dtype = np.float64)

sinogram = radon(phantom, thetas)
recon1 = iradon(sinogram, thetas)
model, recon2 = recon_f(sinogram, thetas)

sinogram2 = radon(recon1, thetas)
sinogram3 = radon(recon2, thetas)

# a, b = np.real(ft_origin), np.real(model)

# fig, (ax1, ax2) = plt.subplots(ncols = 2)
# ax1.plot(a[a.shape[0]//2, :])
# ax2.plot(b[a.shape[0]//2, :])
# ax2.axis([0,128+400, -600, 600])
# ax1.axis([0,128+400, -600, 600])


# imgs = [phantom, recon1, recon2]
# imgs = [phantom, np.abs(ft_origin), np.abs(model), recon2]
# imgs = [sinogram, sinogram2, sinogram3]
# # titles = ["origin", "FT origin",  "inserted", "recon by insertion"]
# # titles = ["origin", "recon by FBP", "recon by Insertion"]
# titles = ["a", "b", "c"]
# fig, axes = plt.subplots(ncols = len(imgs))
# for i, ax in enumerate(axes):
#     ax.imshow(imgs[i])
#     ax.set_yticks([])
#     ax.set_title(titles[i])
# ax.set_xlabel("angles = 0:180:1")

plt.plot(sinogram3[:,1], 'b')
plt.plot(sinogram2[:,1], "r")
plt.tight_layout()
# plt.savefig("imgs/phan3_pad.png")

plt.show()
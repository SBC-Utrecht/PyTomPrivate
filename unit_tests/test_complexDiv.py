from pytom.tompy.io import read, write
import  numpy as np

import matplotlib
matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt
from pytom.tompy.tools import create_sphere
from numpy.fft import fftn, fftshift, ifftn

start_array = np.zeros((10,10,10), dtype=np.float32)
sphere = create_sphere(start_array.shape, 3)





start_array[5,5,5] = 1


start_array = ifftn(fftshift(fftshift(fftn(start_array))*sphere)).real


fname = '/data2/ctsanchez/graphs/templates/complexDiv/singlepix.mrc'
write('/data2/ctsanchez/graphs/templates/complexDiv/singlepix.mrc', start_array)


from pytom.basic.files import read as readC
from pytom.basic.fourier import fft, ifft
from pytom_volume import complexDiv, abs, real, sum, imag, limit
from pytom_numpy import vol2npy
from pytom.tompy.transform import fourier_reduced2full as red

# Cpu version

v = readC(fname)
print(vol2npy(v).sum())
v_np = vol2npy(v)

ft = fft(v)
amp = real(abs(ft))
aa = red(vol2npy(amp))
threshold = 9.99999999e-8
amp.write('/data2/ctsanchez/graphs/templates/complexDiv/amp.em')
#limit(amp, threshold, 0, 0, 0, True, False)
v2 = complexDiv(ft,amp)
divided = red(vol2npy(v2))
v2 = ifft(v2)
v2.shiftscale(shiftV=0, scaleV=1/(start_array.size))
v2 = vol2npy(v2).copy()

print(v2.sum(), v2.max())

#Gpu version

import cupy as xp
template = xp.array(start_array,dtype=xp.float32)
ftemplate = xp.fft.fftn(template)

size = ftemplate.shape
a = xp.abs(ftemplate)
aaa = a.copy()
#ftemplate[a == 0]= -1
#a[a == 0] = 1

print(a.max(), a.min(), aa.max(), aa.min())

#
# for i in range(10):
#     for j in range(10):
#         for k in range(10):
#             if a[i,j,k] >= threshold:
#                 ftemplate[i,j,k] /= a[i,j,k]
#                 if ftemplate[i,j,k].real < -0.05: print(i,j,k, ftemplate[i,j,k], divided[i,j,k])
#             else:
#                 ftemplate[i,j,k] = ftemplate[i,j,k]

a[a<=threshold] = 0
ftemplate /= a
ftemplate[xp.isnan(ftemplate)] = 0

# print(new_a.sum(), new_a.size)

# phase = xp.angle(ftemplate)
#
# ftemplateftemplate = new_a * xp.exp(1j * phase)

#ftemplate /= xp.abs(ftemplate)
print('max is', ftemplate.max())
#ftemplate[ftemplate ] = 0-

template2 = (xp.fft.ifftn(ftemplate).real)
print(template2.sum())
#template2[0:2,0,0] = 0
#template2[-2:,0,0] = 0
diff = np.abs(np.abs(ftemplate).get()-np.abs(divided))

x,y,z = (np.unravel_index(diff.argmax(),diff.shape))

print(template2[x,y,z], v2[x,y,z], divided[x,y,z], ftemplate[x,y,z], a[x,y,z])


diff = np.abs(template2.get()-v2)


print(diff.max(), diff.min(), diff.mean(), diff.std())
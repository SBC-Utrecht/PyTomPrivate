import numpy as np
from pytom.agnostic.transform import fourier_full2reduced, fourier_reduced2full
from pytom.agnostic.tools import subvolume, putSubVolume
from pytom.gpu.initialize import xp, device
import numpy as cp
from numba import njit, prange

class vol(np.ndarray):

    def __new__(subtype, data, sy=0, sz=0, dtype=None, buffer=None, offset=0,
                strides=None, order=None, info=None):
        if isinstance(data,(vol, vol_comp, np.ndarray)):
            sx,sy,sz = data.shape
        else:
            sx = data
            data = np.zeros((sx,sy,sz),dtype=np.float32)

        obj = super().__new__(subtype, (sx, sy, sz), np.float32,
                              buffer, offset, strides, order)

        # set the new 'info' attribute to the value passed
        obj.info = info
        obj.FtSize = list(data.shape)
        obj[:] = data
        # Finally, we must return the newly created object:
        return obj

    def __array_finalize__(self, obj):
        if obj is None: return
        self.info = getattr(obj, 'info', None)



    def sizeX(self):
        return self.shape[0]

    def sizeY(self):
        return self.shape[1]

    def sizeZ(self):
        return self.shape[2]

    def setFtSizeX(self, sizeX):
        self.FtSize[0] = sizeX

    def setFtSizeY(self, sizeY):
        self.FtSize[1] = sizeY

    def setFtSizeZ(self, sizeZ):
        self.FtSize[2] = sizeZ

    def getFtSizeX(self):
        return self.FtSize[0]

    def getFtSizeY(self):
        return self.FtSize[1]

    def getFtSizeZ(self):
        return self.FtSize[2]

    def strideX(self):
        return self.strides[2]

    def strideY(self):
        return self.strides[1]

    def strideZ(self):
        return self.strides[0]

    def write(self, fileName, fileType='mrc'):
        """
        write(vol self, std::string fileName)
        write(vol self, std::string fileName, std::string fileType)
        """
        from pytom.agnostic.io import write

        if not fileName.endswith(fileType):
            fileName += f'.{fileType}'

        return write(fileName, np.array(self[:,:,:]))

    def info(self, name):
        """info(vol self, std::string const & name)"""

        string = f'''Volume {name}:
  size[zyx] = [{self.sizeZ()},{self.sizeY()},{self.sizeX()}]; of tom_io_type 5
  stride[zyx] = [{self.strideZ()},{self.strideY()},{self.strideX()}]; gap[zyx] = [0,0,0]; 
  range = [{self.min()} .. {self.mean()} .. {self.ma()}]
  stddev = {self.std()} (stddev_sample = 0)
'''

        print(string)

    def numelem(self):
        """numelem(vol self) -> std::size_t"""
        return self.size

    def equalsTo(self, v):
        """equalsTo(vol self, vol v) -> bool"""
        return np.array_equal(self, v)

    def getV(self, x, y, z):
        """getV(vol self, std::size_t x, std::size_t y, std::size_t z) -> float"""
        return self[x, y, z]

    def setV(self, val, x, y, z):
        """setV(vol self, float val, std::size_t x, std::size_t y, std::size_t z)"""
        self[x, y, z] = val

    def setAll(self, val):
        """setAll(vol self, float val)"""
        self[:,:,:] = val

    def copyVolume(self, v2):
        """copyVolume(vol self, vol v2)"""
        self[:,:,:] = v2[:,:,:]
        return v2

    def shiftscale(self, shiftV, scaleV):
        """shiftscale(vol self, float const & shiftV, float const & scaleV)"""
        self[:,:,:] = self[:,:,:] + shiftV
        self[:,:,:] = self[:,:,:] * scaleV

    def dims(self):
        """dims(vol self) -> std::size_t"""
        return len(self.shape)

    def duplicate(self):
        v2 = vol(*self.shape)
        v2[:] = self[:,:,:]
        return v2

    def __call__(self, val, x, y, z=None):
        if z is None:
            return self[val,x,y]
        else:
            self[x,y,z] = val

class vol_comp(np.ndarray):
    def __new__(subtype, data, sy=0, sz=0, buffer=None, offset=0,
                strides=None, order=None, info=None):
        # Create the ndarray instance of our type, given the usual
        # ndarray input arguments.

        if isinstance(data, (vol, vol_comp, np.ndarray)):
            sx,sy,sz = data.shape
        else:
            sx = data
            data = np.zeros((sx,sy,sz),dtype=np.complex64)


        obj = super().__new__(subtype, (sx, sy, sz), np.complex64,
                              buffer, offset, strides, order)
        # set the new 'info' attribute to the value passed
        obj.info = info
        obj.FtSize = list(data.shape)
        obj[:] = data

        # Finally, we must return the newly created object:
        return obj

    def __array_finalize__(self, obj):
        if obj is None: return
        self.info = getattr(obj, 'info', None)

    def sizeX(self):
        return self.shape[0]

    def sizeY(self):
        return self.shape[1]

    def sizeZ(self):
        return self.shape[2]

    def setFtSizeX(self, sizeX):
        self.FtSize[0] = sizeX

    def setFtSizeY(self, sizeY):
        self.FtSize[1] = sizeY

    def setFtSizeZ(self, sizeZ):
        self.FtSize[2] = sizeZ

    def getFtSizeX(self):
        return self.FtSize[0]

    def getFtSizeY(self):
        return self.FtSize[1]

    def getFtSizeZ(self):
        return self.FtSize[2]

    def strideX(self):
        return self.strides[2]

    def strideY(self):
        return self.strides[1]

    def strideZ(self):
        return self.strides[0]

    def write(self, fileName, fileType=None):
        """
        write(vol self, std::string fileName)
        write(vol self, std::string fileName, std::string fileType)
        """
        from pytom.agnostic.io import write
        if not fileName.endswith(fileType):
            fileName += f'.{fileType}'

        return write(fileName, np.array(self[:,:,:]))

    def info(self, name):
        """info(vol self, std::string const & name)"""

        string = f'''Volume {name}:
  size[zyx] = [{self.sizeZ()},{self.sizeY()},{self.sizeX()}]; of tom_io_type 5
  stride[zyx] = [{self.strideZ()},{self.strideY()},{self.strideX()}]; gap[zyx] = [0,0,0]; 
  range = [{self.min()} .. {self.mean()} .. {self.ma()}]
  stddev = {self.std()} (stddev_sample = 0)
'''

        print(string)

    def numelem(self):
        """numelem(vol self) -> std::size_t"""
        return self.size

    def equalsTo(self, v):
        """equalsTo(vol self, vol v) -> bool"""
        return np.array_equal(self, v)

    def getV(self, x, y, z):
        """getV(vol self, std::size_t x, std::size_t y, std::size_t z) -> float"""
        return self[x, y, z]

    def setV(self, val, x, y, z):
        """setV(vol self, float val, std::size_t x, std::size_t y, std::size_t z)"""
        self[x, y, z] = val

    def setAll(self, val):
        """setAll(vol self, float val)"""
        self[:, :, :] = val

    def copyVolume(self, v2):
        """copyVolume(vol self, vol v2)"""
        self[:, :, :] = v2[:, :, :]
        return v2

    def shiftscale(self, shiftV, scaleV):
        """shiftscale(vol self, float const & shiftV, float const & scaleV)"""
        self[:, :, :] = self[:, :, :] + shiftV
        self[:, :, :] = self[:, :, :] * scaleV

    def dims(self):
        """dims(vol self) -> std::size_t"""
        return len(self.shape)

    def duplicate(self):
        v2 = vol_comp(*self.shape)
        v2[:] = self[:,:,:]
        return v2

    def __call__(self, val, x, y, z=None):
        if z is None:
            return self[val,x,y]
        else:
            self[x,y,z] = val
# TESTED

def read(fileName, subregion1=None, subregion2=None, subregion3=None, subregion4=None, subregion5=None, subregion6=None,
         sampling1=1, sampling2=1, sampling3=1, binning1=1, binning2=1, binning3=1):
    """
    read(std::string fileName) -> vol
    read(std::string fileName, std::size_t sampling1, std::size_t sampling2, std::size_t sampling3, std::size_t binning1, std::size_t binning2, std::size_t binning3) -> vol
    """
    from pytom.agnostic.io import read
    from pytom.agnostic.transform import resize
    from pytom.gui.mrcOperations import downsample

    data = read(fileName)
    v = vol(*data.shape)
    v[:] = data


    if not subregion1 is None and not (subregion1 == 0 and subregion4 == 0):
        v = v[subregion1:subregion4+1, subregion2:subregion5+1, subregion3:subregion6+1]
    if not sampling1 is None and sampling1 != 0:
        v = v[::sampling1,::sampling2, ::sampling3]

    #TODO resize is uniform over all dimensions
    binning1 = 1 if binning1 == 0 else binning1
    binning2 = 1 if binning2 == 0 else binning2
    binning3 = 1 if binning3 == 0 else binning3

    if not (binning1 == 1 and binning2 == 1 and binning3 == 1):
        # Non-integer binning is allowed for read, while integer binning shows same behavior as pytomc
        if True in [np.abs(bin//1 - bin) > 1E-6 for bin in (binning1, binning2, binning3)]:
            v = resize(v, 1/binning1)
        else:
            v = downsample(v,[binning1,binning2,binning3])

    return v

def sum(volume):
    """
    sum(vol volume) -> float
    sum(vol_comp volume) -> std::complex< float >
    """
    return volume.sum()

def mean(volume):
    """mean(vol volume) -> double"""
    return volume.mean()

def variance(volume, use_sample_standard_deviation):
    """variance(vol volume, bool use_sample_standard_deviation) -> double"""
    return volume.std()**2

def min(volume):
    """min(vol volume) -> float"""
    return volume.min()

def max(volume):
    """max(vol volume) -> float"""
    return volume.max()

def power(volume, exponent):
    """
    power(vol volume, float const & exponent)
    power(vol_comp volume, float const & exponent)
    """
    volume[:,:,:] = volume**exponent

def conjugate(volume):
    """conjugate(vol_comp volume)"""
    volume[:,:,:] = volume.conj()

def limit(volume, lowerBound, lowerReplacement, upperBound, upperReplacement, doLower, doUpper):
    """limit(vol volume, float lowerBound, float lowerReplacement, float upperBound, float upperReplacement, bool doLower, bool doUpper)"""
    if doLower:
        volume[volume < lowerBound] = lowerReplacement
    if doUpper:
        volume[volume> upperBound] = upperReplacement

def abs(volume):
    return npy2vol(np.abs(volume))

# Complex functions
def real(volume):
    """real(vol_comp vol) -> vol"""
    return volume.real

def imag(volume):
    """imag(vol_comp vol) -> vol"""
    im = (volume.imag) * 1

    return im

def conj_mult(v1, v2):
    """conj_mult(vol_comp v1, vol_comp v2)"""
    v1[:,:,:]=  v1.conj() * v2

def complexDiv(volume, div):
    """complexDiv(vol_comp volume, vol div) -> vol_comp"""
    return volume / div

def complexRealMult(vol, otherVol):
    """complexRealMult(vol_comp vol, vol otherVol) -> vol_comp"""
    return vol * otherVol

def mergeRealImag(real, imag):
    """mergeRealImag(vol real, vol imag) -> vol_comp"""
    v = vol_comp(*real.shape)
    v.real[:,:,:] = real
    v.imag[:,:,:] = imag
    return v


def maskNorm(volume, mask, is_boolean_mask=False):
    from pytom.tools.maths import epsilon

    if is_boolean_mask:
        mask = mask.astype(xp.float32)

    volume -= (volume[mask>0].mean())

    volume *= mask
    volume -= volume[mask>0].mean() * (mask>0)
    volumeStd = volume.std()
    volumeStd = 1 if volumeStd < epsilon else volumeStd

    volume /= volumeStd

def initSphere(vol, radius, sigma, max_radius, centerx=None, centery=None, centerz=None):
    """initSphere(vol vol, float radius, float sigma, float max_radius, float centerx, float centery, float centerz)"""

    size = vol.shape

    if size.__class__ == float or len(size) == 1:
        size = (size, size, size)
    assert len(size) == 3


    vol.setAll(0)

    if centerx is None:
        center = [size[0] // 2 -1, size[1] // 2 -1, size[2] // 2 -1]
    else:
        center = [centerx, centery, centerz]

    if radius == -1:
        radius = 0
    if (max_radius < radius):
        max_radius = radius + sqrt(2. * sigma)

    [x, y, z] = np.mgrid[0:size[0], 0:size[1], 0:size[2]]
    r = np.sqrt((x - center[0]) ** 2 + (y - center[1]) ** 2 + (z - center[2]) ** 2)
    vol[r <= radius] = 1

    if sigma > 0:
        ind = np.logical_and(r > radius, r <= max_radius)
        vol[ind] = np.exp(-((r[ind] - radius) / sigma) ** 2 )

def mirrorVolume(src, des):
    """mirrorVolume(vol src, vol des)"""
    des[:] = src[::-1,:,:]

def paste(volume, destination, x, y, z):
    """paste(vol volume, vol destination, long x, long y, long z)"""
    #from pytom.gpu.gpuFunctions import general_transform_crop as transvt

    sx, sy, sz = volume.shape
    destination[x:x+sx, y:y+sy, z:z+sz] = volume#transvt(a, destination, translation=(dx//2-sx, dy//2-sy, dz//2-sz))

    #return _pytom_volume.paste(volume, destination, x, y, z)

def pasteCenter(volume, destination):
    """pasteCenter(vol volume, vol destination)"""
    from pytom.agnostic.tools import paste_in_center

    destination[:] = paste_in_center(volume, destination)

def numberSetVoxels(volume):
    """numberSetVoxels(vol volume) -> std::size_t"""
    return (1*(np.abs(volume) >1E-46)).sum()

def projectSum(volume, axis):
    """projectSum(vol volume, int axis) -> vol"""
    return volume.sum(axis=axis)

def gaussianNoise(volume, mean, sigma):
    """gaussianNoise(vol volume, double mean, double sigma)"""
    volume[:] = np.random.normal(mean, sigma, volume.shape)

def writeSubregion(source, filename, positionX, positionY, positionZ):
    """writeSubregion(vol source, std::string filename, std::size_t positionX, std::size_t positionY, std::size_t positionZ)"""
    v = read(filename)
    v[positionX:positionX+source.shape[0], positionY:positionY+source.shape[1], positionZ:positionZ+source.shape[2]] = source
    v.write(filename)

def updateResFromIdx(resV, newV, orientV, orientIdx):
    """updateResFromIdx(vol resV, vol newV, vol orientV, std::size_t const orientIdx)"""
    orientV[resV<newV] = orientIdx
    resV[:] = np.maximum.reduce([resV, newV])

def updateResFromVol(resV, newV, orientV, neworientV):
    """updateResFromVol(vol resV, vol newV, vol orientV, vol neworientV)"""
    orientV[resV<newV] = 0
    orientV += neworientV*(resV<newV)
    resV[:] = np.maximum.reduce([resV, newV])


# Interpolation functions
def interpolate(source, x, y, z):
    """interpolate(vol source, float const x, float const y, float const z) -> float"""
    return linearInterpolation(source, x, y, z)

def interpolateCubic(source, x, y, z):
    """interpolateCubic(vol source, float const x, float const y, float const z) -> float"""
    return cubicInterpolation(source, x, y, z)

def interpolateSpline(source, x, y, z):
    """interpolateSpline(vol source, float const x, float const y, float const z) -> float"""
    return splineInterpolation(source, x, y, z)

def interpolateFourierSpline(source, x, y, z):
    """interpolateFourierSpline(vol_comp source, float const x, float const y, float const z) -> std::complex< float >"""
    return splineInterpolation(source, x, y, z)





# NON TESTED

# Scaling functions
def rescale(source, destination):
    """rescale(vol source, vol destination)"""
    from pytom.voltools.utils.matrices import scale_matrix
    from  pytom.agnostic.transform import resize
    sx, sy, sz = source.shape
    dx, dy, dz = destination.shape

    # destination[:] = resize(source, dx/sx)
    mtx = scale_matrix(coefficients=(sx/dy, sy/dy, sz/dz))

    general_transform(source, destination, mtx, 'linear')

def rescaleCubic(source, destination):
    """rescaleCubic(vol source, vol destination)"""
    from pytom.voltools.utils.matrices import scale_matrix

    sx, sy, sz = source.shape
    dx, dy, dz = destination.shape

    mtx = scale_matrix(coefficients=(sx / dy, sy / dy, sz / dz))

    general_transform(source, destination, mtx, 'cubic')

def rescaleSpline(source, destination):
    """rescaleSpline(vol source, vol destination)"""
    from pytom.voltools.utils.matrices import scale_matrix

    sx, sy, sz = source.shape
    dx, dy, dz = destination.shape

    mtx = scale_matrix(coefficients=(sx / dy, sy / dy, sz / dz))

    general_transform(source, destination, mtx, 'spline')

# Rotation functions
def rotate(src, dst, phi, psi, theta):
    """rotate(vol src, vol dst, double phi, double psi, double theta)"""
    from pytom.voltools.utils.matrices import transform_matrix

    mtx = transform_matrix(rotation=(phi,theta,psi), rotation_order='rzxz', center=[src.sizeX()//2, src.sizeY()//2, src.sizeZ()//2])
    general_transform(src, dst, mtx, 'linear')

def rotateCubic(src, dst, phi, psi, theta):
    """rotateCubic(vol src, vol dst, double phi, double psi, double theta)"""
    from pytom.voltools.utils.matrices import transform_matrix

    mtx = transform_matrix(rotation=(phi,theta,psi), rotation_order='rzxz', center=[src.sizeX()//2, src.sizeY()//2, src.sizeZ()//2])
    general_transform(src, dst, mtx, 'cubic')

def rotateSpline(src, dst, phi, psi, theta):
    """rotateSpline(vol src, vol dst, double phi, double psi, double theta)"""
    from pytom.voltools.utils.matrices import transform_matrix

    mtx = transform_matrix(rotation=(phi,theta,psi), rotation_order='rzxz', center=[src.sizeX()//2, src.sizeY()//2, src.sizeZ()//2])
    general_transform(src, dst, mtx, 'spline')

def rotateSplineInFourier(src, phi, psi, theta):
    """rotateSplineInFourier(vol_comp src, double phi, double psi, double theta) -> vol_comp"""
    from pytom.voltools.utils.matrices import transform_matrix

    full = xp.fft.fftshift(reducedToFull(src, isodd=(False if (src.shape[0] %2 ==0) else False)))
    dst = vol_comp(*full.shape)

    dst[:] = 0

    mtx = transform_matrix(rotation=(phi,theta,psi), rotation_order='rzxz', center=[full.shape[0]//2, full.shape[1]//2, full.shape[2]//2])
    general_transform(full, dst, mtx, 'complex')

    return np.fft.fftshift(dst)[:,:,:src.sizeZ()]

#Translation functions
def shift(src, dst, preX, preY, preZ):
    """shift(vol src, vol dst, double preX, double preY, double preZ)"""
    return _pytom_volume.shift(src, dst, preX, preY, preZ)

def shiftFourier(src, dst, shiftX, shiftY, shiftZ):
    """shiftFourier(vol_comp src, vol_comp dst, double shiftX, double shiftY, double shiftZ)"""
    sizex, sizey, sizez = src.shape
    dst[:] = src[:] * create_fourier_shift(sizex, sizey, sizez, shiftX, shiftY, shiftZ)

def create_fourier_shift(sizex, sizey, sizez, shiftX, shiftY, shiftZ, grids=None):
    M_PI = np.pi
    if grids is None:
        x,y,z = np.meshgrid(np.fft.fftshift(np.arange(sizex)-sizex//2), np.fft.fftshift(np.arange(sizey)-sizey//2), np.arange(sizez), indexing='ij')
    else:
        x,y,z = grids

    shift = mergeRealImag(xp.cos(2 * M_PI / sizex * x * shiftX), -xp.sin(2 * M_PI / sizex * x * shiftX)) * \
            mergeRealImag(xp.cos(2 * M_PI / sizey * y * shiftY), -xp.sin(2 * M_PI / sizey * y * shiftY)) * \
            mergeRealImag(xp.cos(2 * M_PI / sizez * z * shiftZ), -xp.sin(2 * M_PI / sizez * z * shiftZ))

    return shift

# General transformation functions
def transform(src, dst, phi, psi, theta, centerX, centerY, centerZ, preX=0, preY=0, preZ=0, postX=0, postY=0, postZ=0):
    """transform(vol src, vol dst, double phi, double psi, double theta, double centerX, double centerY, double centerZ, double preX, double preY, double preZ, double postX, double postY, double postZ)"""
    from pytom.voltools.utils.matrices import rotation_matrix, translation_matrix, transform_matrix

    tpre = translation_matrix(translation=(preX,preY,preZ))
    tpost = translation_matrix(translation=(postX, postY, postZ))
    rot = transform_matrix(rotation=(phi,theta,psi), rotation_order='rzxz', center=(centerX, centerY, centerZ))
    mtx = tpre.dot(rot.dot(tpost))
    general_transform(src, dst, mtx, 'linear')

def is_identity(mtx):
    sum = []
    [sum.append(np.abs(value-1) < 1E-6) for value in mtx.flatten()]
    if False in sum:
        return False
    else:
        return True

def general_transform(src, dst, mtx, interpolation_func='spline', func=None, shift=None):
    """general_transform(vol src, vol dst, vol mtx)"""

    from pytom.agnostic.interpolation import fill_values_complex, fill_values_real_cubic, fill_values_real_fourierSpline, fill_values_real_linear, fill_values_real_spline

    interpolation = {'spline': fill_values_real_spline, 'cubic': fill_values_real_cubic, 'linear': fill_values_real_linear,
                     'fourierSpline': fill_values_real_fourierSpline, 'complex': fill_values_complex}

    dims = dst.shape
    a = xp.identity(4)
    if mtx.size == 9:
        mtx = mtx.flatten().reshape((3,3))
        a[:3,:3] = mtx
        mtx = a

    mtx = np.linalg.inv(mtx)
    mtx = mtx.flatten()

    if is_identity(mtx):
        try:
            dst[:src.shape[0],:src.shape[1],:src.shape[2]] = src[:]
        except:
            dst[:] = src[:dims[0], :dims[1], :dims[2]]
    else:
        if shift is None:
            interpolation[interpolation_func](src, dst, mtx, src.shape, dims)

        else: interpolation[interpolation_func](src, dst, mtx, src.shape, dims, shift[0], shift[1], shift[2])

def transformCubic(src, dst, phi, psi, theta, centerX, centerY, centerZ, preX, preY, preZ, postX, postY, postZ):
    """transformCubic(vol src, vol dst, double phi, double psi, double theta, double centerX, double centerY, double centerZ, double preX, double preY, double preZ, double postX, double postY, double postZ)"""
    from pytom.voltools.utils.matrices import rotation_matrix, translation_matrix, transform_matrix

    tpre = translation_matrix(translation=(preX,preY,preZ))
    tpost = translation_matrix(translation=(postX, postY, postZ))
    rot = transform_matrix(rotation=(phi,theta,psi), rotation_order='rzxz', center=(centerX, centerY, centerZ))
    mtx = tpre.dot(rot.dot(tpost))

    general_transform(src, dst, mtx, 'cubic')

def transformSpline(src, dst, phi, psi, theta, centerX, centerY, centerZ, preX, preY, preZ, postX, postY, postZ):
    """transformSpline(vol src, vol dst, double phi, double psi, double theta, double centerX, double centerY, double centerZ, double preX, double preY, double preZ, double postX, double postY, double postZ)"""
    from pytom.voltools.utils.matrices import rotation_matrix, translation_matrix, transform_matrix

    tpre = translation_matrix(translation=(preX,preY,preZ))
    tpost = translation_matrix(translation=(postX, postY, postZ))
    rot = transform_matrix(rotation=(phi,theta,psi), rotation_order='rzxz', center=(centerX, centerY, centerZ))
    mtx = tpre.dot(rot.dot(tpost))

    general_transform(src, dst, mtx, 'spline')

def transformFourierSpline(src, phi, psi, theta, shiftX, shiftY, shiftZ):
    """transformFourierSpline(vol_comp src, double phi, double psi, double theta, double shiftX, double shiftY, double shiftZ) -> vol_comp"""
    from pytom.voltools.utils.matrices import rotation_matrix
    full = np.fft.fftshift(reducedToFull(src))
    dst = vol_comp(*full.shape)

    mtx = rotation_matrix(rotation=(phi,theta,psi), rotation_order='rzxz')

    general_transform(full, dst, mtx, 'complex', shift=[shiftX, shiftY, shiftZ])

    dst[:] = np.fft.fftshift(dst[:])
    return fullToReduced(dst)

def peak(volume, mask=1):
    """
    peak(vol volume) -> tom::st_idx
    peak(vol volume, vol mask) -> tom::st_idx
    """
    return np.unravel_index((volume*mask).argmax(), volume.shape)

def backProject(*args):
    """
    backProject(vol src, vol dst, vol phi, vol theta, vol offsetPaticle, vol offsetProjections)
    backProject(vol src, vol dst, vol phi, vol theta, vol psi, vol offsetPaticle, vol offsetProjections)
    """
    from pytom.agnostic.backprojectionCpu import backProjectCPU
    if len(args) == 6:
        src, dst, phi, theta, offsetParticle, offsetProjections = args
    else:
        src, dst, phi, theta, psi, offsetParticle, offsetProjections = args

    rec = backProjectCPU(src, dst.copy(), phi, theta, offsetParticle, offsetProjections)
    dst[:] = rec

def vectorize(source):
    """vectorize(vol source) -> vol"""
    out = vol(source.size,1,1)
    out[:] = np.expand_dims(np.expand_dims(source.T.flatten(),2),1)
    return out

def vol2npy(volume):
    return volume[:]

def npy2vol(volume):
    if 'complex' in str(volume.dtype):
        v = vol_comp(*volume.shape)
    else:
        v = vol(*volume.shape)

    v[:] = volume
    return v

def reducedToFull(vol_comp):

    full = fourier_reduced2full(vol_comp[:],isodd=(vol_comp.sizeX()%2 ==1))
    v2 = vol_comp(*full.shape)
    v2[:] = full
    return v2

def fullToReduced(vol_comp):

    red = fourier_full2reduced(vol_comp[:])
    v2 = vol_comp(*red.shape)
    v2[:] = red

    return v2

def ftshift(volume):
    v2 = volume.duplicate()
    v2[:] = np.fft.fftshift(volume)
    return v2

def iftshift(volume):
    v2 = volume.duplicate()
    v2[:] = np.fft.ifftshift(volume)
    return v2

def fft(data,scaling=''):
    return vol_comp(np.atleast_3d( np.fft.rfftn(data.squeeze())))

def ifft(data, scaling=''):
    return vol(np.atleast_3d(np.fft.irfftn(data.squeeze()).real))

def volTOsf(data, res, r, b, m_x, m_y, m_z):
    """volTOsf(vol data, double const r, int const b, float const m_x, float const m_y, float const m_z) -> vol"""
    from pytom.agnostic.frm import vol2sf
    res = xp.zeros((4 * b * b, 1, 1),dtype=np.float32)

    vol2sf(data,res,r,b,m_x,m_y,m_z)

    return npy2vol(res)

def fvolTOsf(data, r, b):
    """fvolTOsf(vol_comp data, double const r, int const b) -> vol_comp"""

    from pytom.agnostic.frm import fvol2sf

    res = xp.zeros((4 * b * b, 1, 1),dtype=np.complex64)

    fvol2sf(data,res,r,b)

    return npy2vol(res)

def sfFourierShift(data, r, b, sizex, sizey, sizez, shiftX, shiftY, shiftZ):
    """sfFourierShift(vol_comp vol, double const r, int const b, int sizex, int sizey, int sizez, double shiftX, double shiftY, double shiftZ) -> vol_comp"""

    from pytom.agnostic.interpolation import splineInterpolation
    from pytom.agnostic.frm import apply_shift_sf

    if (4*b*b,1,1) != data.shape:
        raise Exception('Wrong size of data')

    res = np.zeros((4 * b * b, 1, 1),dtype=np.float32)
    M_PI = xp.pi

    apply_shift_sf(data, res, b,r,M_PI, sizex, sizey, sizez, shiftX, shiftY, shiftZ )

    res = npy2vol(res)
    return res

def powerspectrum(volume):
    """
    compute power spectrum of a volume

    @param volume: input volume
    @type volume: L{pytom_volume.vol}
    @return: power spectrum of vol
    @rtype: L{pytom_volume.vol}

    @author: FF
    """

    fvol = ftshift(reducedToFull(fft(volume)), inplace=False)

    ps = vol(*fvol.shape)
    ps[:] = fvol*fvol.conj()

    return ps



# This file is compatible with both classic and new-style classes.

if __name__ == '__main__':
    a = vol(10,10,10)
    a.setAll(2.)

    print(a)

    a.shiftscale(1,2)

    print(a)

    a.write('return.em')


    print(a.shape)


    print(a.sizeY())
    import cupy

    b = cupy.array(a,dtype=cupy.float32)

    print(b.shape)


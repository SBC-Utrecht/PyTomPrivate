# This file was automatically generated by SWIG (https://www.swig.org).
# Version 4.1.1
#
# Do not make changes to this file unless you know what you are doing - modify
# the SWIG interface file instead.

from sys import version_info as _swig_python_version_info
# Import the low-level C/C++ module
if __package__ or "." in __name__:
    from . import _pytom_volume
else:
    import _pytom_volume

try:
    import builtins as __builtin__
except ImportError:
    import __builtin__

def _swig_repr(self):
    try:
        strthis = "proxy of " + self.this.__repr__()
    except __builtin__.Exception:
        strthis = ""
    return "<%s.%s; %s >" % (self.__class__.__module__, self.__class__.__name__, strthis,)


def _swig_setattr_nondynamic_instance_variable(set):
    def set_instance_attr(self, name, value):
        if name == "this":
            set(self, name, value)
        elif name == "thisown":
            self.this.own(value)
        elif hasattr(self, name) and isinstance(getattr(type(self), name), property):
            set(self, name, value)
        else:
            raise AttributeError("You cannot add instance attributes to %s" % self)
    return set_instance_attr


def _swig_setattr_nondynamic_class_variable(set):
    def set_class_attr(cls, name, value):
        if hasattr(cls, name) and not isinstance(getattr(cls, name), property):
            set(cls, name, value)
        else:
            raise AttributeError("You cannot add class attributes to %s" % cls)
    return set_class_attr


def _swig_add_metaclass(metaclass):
    """Class decorator for adding a metaclass to a SWIG wrapped class - a slimmed down version of six.add_metaclass"""
    def wrapper(cls):
        return metaclass(cls.__name__, cls.__bases__, cls.__dict__.copy())
    return wrapper


class _SwigNonDynamicMeta(type):
    """Meta class to enforce nondynamic attributes (no new attributes) for a class"""
    __setattr__ = _swig_setattr_nondynamic_class_variable(type.__setattr__)


class vol(object):
    r"""Proxy of C++ swigTom::swigVolume< float,float > class."""

    thisown = property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc="The membership flag")
    __repr__ = _swig_repr

    def __init__(self, *args):
        r"""
        __init__(vol self, std::size_t sizex, std::size_t sizey, std::size_t sizez) -> vol
        __init__(vol self, tom::Volume< float > const & v) -> vol
        __init__(vol self, vol v) -> vol
        __init__(vol self, float * data, std::size_t sizex, std::size_t sizey, std::size_t sizez, std::size_t stridex, std::size_t stridey, std::size_t stridez) -> vol
        """
        _pytom_volume.vol_swiginit(self, _pytom_volume.new_vol(*args))
    __swig_destroy__ = _pytom_volume.delete_vol

    def size_x(self):
        r"""size_x(vol self) -> std::size_t"""
        return _pytom_volume.vol_size_x(self)

    def size_y(self):
        r"""size_y(vol self) -> std::size_t"""
        return _pytom_volume.vol_size_y(self)

    def size_z(self):
        r"""size_z(vol self) -> std::size_t"""
        return _pytom_volume.vol_size_z(self)
    shape = property(_pytom_volume.vol_shape_get, _pytom_volume.vol_shape_set, doc=r"""shape : std::tuple<(std::size_t,std::size_t,std::size_t)>""")

    def getFtSizeX(self):
        r"""getFtSizeX(vol self) -> float"""
        return _pytom_volume.vol_getFtSizeX(self)

    def getFtSizeY(self):
        r"""getFtSizeY(vol self) -> float"""
        return _pytom_volume.vol_getFtSizeY(self)

    def getFtSizeZ(self):
        r"""getFtSizeZ(vol self) -> float"""
        return _pytom_volume.vol_getFtSizeZ(self)

    def setFtSizeX(self, size_x):
        r"""setFtSizeX(vol self, float size_x)"""
        return _pytom_volume.vol_setFtSizeX(self, size_x)

    def setFtSizeY(self, size_y):
        r"""setFtSizeY(vol self, float size_y)"""
        return _pytom_volume.vol_setFtSizeY(self, size_y)

    def setFtSizeZ(self, size_z):
        r"""setFtSizeZ(vol self, float size_z)"""
        return _pytom_volume.vol_setFtSizeZ(self, size_z)

    def strideX(self):
        r"""strideX(vol self) -> std::size_t"""
        return _pytom_volume.vol_strideX(self)

    def strideY(self):
        r"""strideY(vol self) -> std::size_t"""
        return _pytom_volume.vol_strideY(self)

    def strideZ(self):
        r"""strideZ(vol self) -> std::size_t"""
        return _pytom_volume.vol_strideZ(self)

    def write(self, *args):
        r"""
        write(vol self, std::string fileName)
        write(vol self, std::string fileName, std::string fileType)
        """
        return _pytom_volume.vol_write(self, *args)

    def info(self, name):
        r"""info(vol self, std::string const & name)"""
        return _pytom_volume.vol_info(self, name)

    def numelem(self):
        r"""numelem(vol self) -> std::size_t"""
        return _pytom_volume.vol_numelem(self)

    def equalsTo(self, v):
        r"""equalsTo(vol self, vol v) -> bool"""
        return _pytom_volume.vol_equalsTo(self, v)

    def getV(self, x, y, z):
        r"""getV(vol self, std::size_t x, std::size_t y, std::size_t z) -> float"""
        return _pytom_volume.vol_getV(self, x, y, z)

    def setV(self, val, x, y, z):
        r"""setV(vol self, float val, std::size_t x, std::size_t y, std::size_t z)"""
        return _pytom_volume.vol_setV(self, val, x, y, z)

    def setAll(self, val):
        r"""setAll(vol self, float val)"""
        return _pytom_volume.vol_setAll(self, val)

    def copyVolume(self, v2):
        r"""copyVolume(vol self, vol v2)"""
        return _pytom_volume.vol_copyVolume(self, v2)

    def shiftscale(self, shiftV, scaleV):
        r"""shiftscale(vol self, float const & shiftV, float const & scaleV)"""
        return _pytom_volume.vol_shiftscale(self, shiftV, scaleV)

    def dims(self):
        r"""dims(vol self) -> std::size_t"""
        return _pytom_volume.vol_dims(self)

    def __add__(self, *args):
        r"""
        __add__(vol self, vol otherVol) -> vol
        __add__(vol self, float const & value) -> vol
        """
        return _pytom_volume.vol___add__(self, *args)

    def __mul__(self, *args):
        r"""
        __mul__(vol self, vol otherVol) -> vol
        __mul__(vol self, float const & value) -> vol
        """
        return _pytom_volume.vol___mul__(self, *args)

    def __sub__(self, *args):
        r"""
        __sub__(vol self, vol otherVol) -> vol
        __sub__(vol self, float const & value) -> vol
        """
        return _pytom_volume.vol___sub__(self, *args)

    def __truediv__(self, *args):
        return _pytom_volume.vol___truediv__(self, *args)
    __div__ = __truediv__



    def __call__(self, *args):
        r"""
        __call__(vol self, std::size_t const & x, std::size_t const & y, std::size_t const & z) -> float const
        __call__(vol self, float const & value, std::size_t const & x, std::size_t const & y, std::size_t const & z)
        """
        return _pytom_volume.vol___call__(self, *args)

# Register vol in _pytom_volume:
_pytom_volume.vol_swigregister(vol)
class vol_comp(object):
    r"""Proxy of C++ swigTom::swigVolume< std::complex< float >,float > class."""

    thisown = property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc="The membership flag")
    __repr__ = _swig_repr

    def __init__(self, *args):
        r"""
        __init__(vol_comp self, std::size_t sizex, std::size_t sizey, std::size_t sizez) -> vol_comp
        __init__(vol_comp self, tom::Volume< std::complex< float > > const & v) -> vol_comp
        __init__(vol_comp self, vol_comp v) -> vol_comp
        __init__(vol_comp self, std::complex< float > * data, std::size_t sizex, std::size_t sizey, std::size_t sizez, std::size_t stridex, std::size_t stridey, std::size_t stridez) -> vol_comp
        """
        _pytom_volume.vol_comp_swiginit(self, _pytom_volume.new_vol_comp(*args))
    __swig_destroy__ = _pytom_volume.delete_vol_comp

    def size_x(self):
        r"""size_x(vol_comp self) -> std::size_t"""
        return _pytom_volume.vol_comp_size_x(self)

    def size_y(self):
        r"""size_y(vol_comp self) -> std::size_t"""
        return _pytom_volume.vol_comp_size_y(self)

    def size_z(self):
        r"""size_z(vol_comp self) -> std::size_t"""
        return _pytom_volume.vol_comp_size_z(self)
    shape = property(_pytom_volume.vol_comp_shape_get, _pytom_volume.vol_comp_shape_set, doc=r"""shape : std::tuple<(std::size_t,std::size_t,std::size_t)>""")

    def getFtSizeX(self):
        r"""getFtSizeX(vol_comp self) -> float"""
        return _pytom_volume.vol_comp_getFtSizeX(self)

    def getFtSizeY(self):
        r"""getFtSizeY(vol_comp self) -> float"""
        return _pytom_volume.vol_comp_getFtSizeY(self)

    def getFtSizeZ(self):
        r"""getFtSizeZ(vol_comp self) -> float"""
        return _pytom_volume.vol_comp_getFtSizeZ(self)

    def setFtSizeX(self, size_x):
        r"""setFtSizeX(vol_comp self, float size_x)"""
        return _pytom_volume.vol_comp_setFtSizeX(self, size_x)

    def setFtSizeY(self, size_y):
        r"""setFtSizeY(vol_comp self, float size_y)"""
        return _pytom_volume.vol_comp_setFtSizeY(self, size_y)

    def setFtSizeZ(self, size_z):
        r"""setFtSizeZ(vol_comp self, float size_z)"""
        return _pytom_volume.vol_comp_setFtSizeZ(self, size_z)

    def strideX(self):
        r"""strideX(vol_comp self) -> std::size_t"""
        return _pytom_volume.vol_comp_strideX(self)

    def strideY(self):
        r"""strideY(vol_comp self) -> std::size_t"""
        return _pytom_volume.vol_comp_strideY(self)

    def strideZ(self):
        r"""strideZ(vol_comp self) -> std::size_t"""
        return _pytom_volume.vol_comp_strideZ(self)

    def write(self, *args):
        r"""
        write(vol_comp self, std::string fileName)
        write(vol_comp self, std::string fileName, std::string fileType)
        """
        return _pytom_volume.vol_comp_write(self, *args)

    def info(self, name):
        r"""info(vol_comp self, std::string const & name)"""
        return _pytom_volume.vol_comp_info(self, name)

    def numelem(self):
        r"""numelem(vol_comp self) -> std::size_t"""
        return _pytom_volume.vol_comp_numelem(self)

    def equalsTo(self, v):
        r"""equalsTo(vol_comp self, vol_comp v) -> bool"""
        return _pytom_volume.vol_comp_equalsTo(self, v)

    def getV(self, x, y, z):
        r"""getV(vol_comp self, std::size_t x, std::size_t y, std::size_t z) -> std::complex< float >"""
        return _pytom_volume.vol_comp_getV(self, x, y, z)

    def setV(self, val, x, y, z):
        r"""setV(vol_comp self, std::complex< float > val, std::size_t x, std::size_t y, std::size_t z)"""
        return _pytom_volume.vol_comp_setV(self, val, x, y, z)

    def setAll(self, val):
        r"""setAll(vol_comp self, std::complex< float > val)"""
        return _pytom_volume.vol_comp_setAll(self, val)

    def copyVolume(self, v2):
        r"""copyVolume(vol_comp self, vol_comp v2)"""
        return _pytom_volume.vol_comp_copyVolume(self, v2)

    def shiftscale(self, shiftV, scaleV):
        r"""shiftscale(vol_comp self, float const & shiftV, float const & scaleV)"""
        return _pytom_volume.vol_comp_shiftscale(self, shiftV, scaleV)

    def dims(self):
        r"""dims(vol_comp self) -> std::size_t"""
        return _pytom_volume.vol_comp_dims(self)

    def __add__(self, *args):
        r"""
        __add__(vol_comp self, vol_comp otherVol) -> vol_comp
        __add__(vol_comp self, float const & value) -> vol_comp
        """
        return _pytom_volume.vol_comp___add__(self, *args)

    def __mul__(self, *args):
        r"""
        __mul__(vol_comp self, vol_comp otherVol) -> vol_comp
        __mul__(vol_comp self, float const & value) -> vol_comp
        """
        return _pytom_volume.vol_comp___mul__(self, *args)

    def __sub__(self, *args):
        r"""
        __sub__(vol_comp self, vol_comp otherVol) -> vol_comp
        __sub__(vol_comp self, float const & value) -> vol_comp
        """
        return _pytom_volume.vol_comp___sub__(self, *args)

    def __truediv__(self, *args):
        return _pytom_volume.vol_comp___truediv__(self, *args)
    __div__ = __truediv__



    def __call__(self, *args):
        r"""
        __call__(vol_comp self, std::size_t const & x, std::size_t const & y, std::size_t const & z) -> std::complex< float > const
        __call__(vol_comp self, std::complex< float > const & value, std::size_t const & x, std::size_t const & y, std::size_t const & z)
        """
        return _pytom_volume.vol_comp___call__(self, *args)

# Register vol_comp in _pytom_volume:
_pytom_volume.vol_comp_swigregister(vol_comp)

def read(*args):
    r"""
    read(std::string fileName) -> vol
    read(std::string fileName, std::size_t subregion1, std::size_t subregion2, std::size_t subregion3, std::size_t subregion4, std::size_t subregion5, std::size_t subregion6, std::size_t sampling1, std::size_t sampling2, std::size_t sampling3, std::size_t binning1, std::size_t binning2, std::size_t binning3) -> vol
    """
    return _pytom_volume.read(*args)

def conj_mult(v1, v2):
    r"""conj_mult(vol_comp v1, vol_comp v2)"""
    return _pytom_volume.conj_mult(v1, v2)

def initSphere(vol, radius, sigma, max_radius, centerx, centery, centerz):
    r"""initSphere(vol vol, float radius, float sigma, float max_radius, float centerx, float centery, float centerz)"""
    return _pytom_volume.initSphere(vol, radius, sigma, max_radius, centerx, centery, centerz)

def maskNorm(vol, mask, is_boolean_mask):
    r"""maskNorm(vol vol, vol mask, bool is_boolean_mask)"""
    return _pytom_volume.maskNorm(vol, mask, is_boolean_mask)

def rotate(src, dst, phi, psi, theta):
    r"""rotate(vol src, vol dst, double phi, double psi, double theta)"""
    return _pytom_volume.rotate(src, dst, phi, psi, theta)

def rotateCubic(src, dst, phi, psi, theta):
    r"""rotateCubic(vol src, vol dst, double phi, double psi, double theta)"""
    return _pytom_volume.rotateCubic(src, dst, phi, psi, theta)

def rotateSpline(src, dst, phi, psi, theta):
    r"""rotateSpline(vol src, vol dst, double phi, double psi, double theta)"""
    return _pytom_volume.rotateSpline(src, dst, phi, psi, theta)

def rotateSplineInFourier(src, phi, psi, theta):
    r"""rotateSplineInFourier(vol_comp src, double phi, double psi, double theta) -> vol_comp"""
    return _pytom_volume.rotateSplineInFourier(src, phi, psi, theta)

def shift(src, dst, preX, preY, preZ):
    r"""shift(vol src, vol dst, double preX, double preY, double preZ)"""
    return _pytom_volume.shift(src, dst, preX, preY, preZ)

def shiftFourier(src, dst, shiftX, shiftY, shiftZ):
    r"""shiftFourier(vol_comp src, vol_comp dst, double shiftX, double shiftY, double shiftZ)"""
    return _pytom_volume.shiftFourier(src, dst, shiftX, shiftY, shiftZ)

def transform(src, dst, phi, psi, theta, centerX, centerY, centerZ, preX, preY, preZ, postX, postY, postZ):
    r"""transform(vol src, vol dst, double phi, double psi, double theta, double centerX, double centerY, double centerZ, double preX, double preY, double preZ, double postX, double postY, double postZ)"""
    return _pytom_volume.transform(src, dst, phi, psi, theta, centerX, centerY, centerZ, preX, preY, preZ, postX, postY, postZ)

def general_transform(src, dst, mtx):
    r"""general_transform(vol src, vol dst, vol mtx)"""
    return _pytom_volume.general_transform(src, dst, mtx)

def transformCubic(src, dst, phi, psi, theta, centerX, centerY, centerZ, preX, preY, preZ, postX, postY, postZ):
    r"""transformCubic(vol src, vol dst, double phi, double psi, double theta, double centerX, double centerY, double centerZ, double preX, double preY, double preZ, double postX, double postY, double postZ)"""
    return _pytom_volume.transformCubic(src, dst, phi, psi, theta, centerX, centerY, centerZ, preX, preY, preZ, postX, postY, postZ)

def transformSpline(src, dst, phi, psi, theta, centerX, centerY, centerZ, preX, preY, preZ, postX, postY, postZ):
    r"""transformSpline(vol src, vol dst, double phi, double psi, double theta, double centerX, double centerY, double centerZ, double preX, double preY, double preZ, double postX, double postY, double postZ)"""
    return _pytom_volume.transformSpline(src, dst, phi, psi, theta, centerX, centerY, centerZ, preX, preY, preZ, postX, postY, postZ)

def transformFourierSpline(src, phi, psi, theta, shiftX, shiftY, shiftZ):
    r"""transformFourierSpline(vol_comp src, double phi, double psi, double theta, double shiftX, double shiftY, double shiftZ) -> vol_comp"""
    return _pytom_volume.transformFourierSpline(src, phi, psi, theta, shiftX, shiftY, shiftZ)

def peak(*args):
    r"""
    peak(vol volume) -> tom::st_idx
    peak(vol volume, vol mask) -> tom::st_idx
    """
    return _pytom_volume.peak(*args)

def conjugate(volume):
    r"""conjugate(vol_comp volume)"""
    return _pytom_volume.conjugate(volume)

def paste(volume, destination, x, y, z):
    r"""paste(vol volume, vol destination, long x, long y, long z)"""
    return _pytom_volume.paste(volume, destination, x, y, z)

def pasteCenter(volume, destination):
    r"""pasteCenter(vol volume, vol destination)"""
    return _pytom_volume.pasteCenter(volume, destination)

def power(*args):
    r"""
    power(vol volume, float const & exponent)
    power(vol_comp volume, float const & exponent)
    """
    return _pytom_volume.power(*args)

def numberSetVoxels(volume):
    r"""numberSetVoxels(vol volume) -> std::size_t"""
    return _pytom_volume.numberSetVoxels(volume)

def sum(*args):
    r"""
    sum(vol volume) -> float
    sum(vol_comp volume) -> std::complex< float >
    """
    return _pytom_volume.sum(*args)

def projectSum(volume, axis):
    r"""projectSum(vol volume, int axis) -> vol"""
    return _pytom_volume.projectSum(volume, axis)

def mean(volume):
    r"""mean(vol volume) -> double"""
    return _pytom_volume.mean(volume)

def variance(volume, use_sample_standard_deviation):
    r"""variance(vol volume, bool use_sample_standard_deviation) -> double"""
    return _pytom_volume.variance(volume, use_sample_standard_deviation)

def min(volume):
    r"""min(vol volume) -> float"""
    return _pytom_volume.min(volume)

def max(volume):
    r"""max(vol volume) -> float"""
    return _pytom_volume.max(volume)

def limit(volume, lowerBound, lowerReplacement, upperBound, upperReplacement, doLower, doUpper):
    r"""limit(vol volume, float lowerBound, float lowerReplacement, float upperBound, float upperReplacement, bool doLower, bool doUpper)"""
    return _pytom_volume.limit(volume, lowerBound, lowerReplacement, upperBound, upperReplacement, doLower, doUpper)

def abs(*args):
    r"""
    abs(vol volume) -> vol
    abs(vol_comp volume) -> vol_comp
    """
    return _pytom_volume.abs(*args)

def reducedToFull(*args):
    r"""
    reducedToFull(vol volume) -> vol
    reducedToFull(vol_comp volume) -> vol_comp
    """
    return _pytom_volume.reducedToFull(*args)

def gaussianNoise(volume, mean, sigma):
    r"""gaussianNoise(vol volume, double mean, double sigma)"""
    return _pytom_volume.gaussianNoise(volume, mean, sigma)

def complexDiv(volume, div):
    r"""complexDiv(vol_comp volume, vol div) -> vol_comp"""
    return _pytom_volume.complexDiv(volume, div)

def vectorize(source):
    r"""vectorize(vol source) -> vol"""
    return _pytom_volume.vectorize(source)

def subvolume(*args):
    r"""
    subvolume(vol source, std::size_t startX, std::size_t startY, std::size_t startZ, std::size_t endX, std::size_t endY, std::size_t endZ) -> vol
    subvolume(vol_comp source, std::size_t startX, std::size_t startY, std::size_t startZ, std::size_t endX, std::size_t endY, std::size_t endZ) -> vol_comp
    """
    return _pytom_volume.subvolume(*args)

def putSubVolume(*args):
    r"""
    putSubVolume(vol source, vol destination, std::size_t positionX, std::size_t positionY, std::size_t positionZ)
    putSubVolume(vol_comp source, vol_comp destination, std::size_t positionX, std::size_t positionY, std::size_t positionZ)
    """
    return _pytom_volume.putSubVolume(*args)

def writeSubregion(source, filename, positionX, positionY, positionZ):
    r"""writeSubregion(vol source, std::string filename, std::size_t positionX, std::size_t positionY, std::size_t positionZ)"""
    return _pytom_volume.writeSubregion(source, filename, positionX, positionY, positionZ)

def updateResFromIdx(resV, newV, orientV, orientIdx):
    r"""updateResFromIdx(vol resV, vol newV, vol orientV, std::size_t const orientIdx)"""
    return _pytom_volume.updateResFromIdx(resV, newV, orientV, orientIdx)

def updateResFromVol(resV, newV, orientV, neworientV):
    r"""updateResFromVol(vol resV, vol newV, vol orientV, vol neworientV)"""
    return _pytom_volume.updateResFromVol(resV, newV, orientV, neworientV)

def mirrorVolume(src, des):
    r"""mirrorVolume(vol src, vol des)"""
    return _pytom_volume.mirrorVolume(src, des)

def rescale(source, destination):
    r"""rescale(vol source, vol destination)"""
    return _pytom_volume.rescale(source, destination)

def rescaleCubic(source, destination):
    r"""rescaleCubic(vol source, vol destination)"""
    return _pytom_volume.rescaleCubic(source, destination)

def rescaleSpline(source, destination):
    r"""rescaleSpline(vol source, vol destination)"""
    return _pytom_volume.rescaleSpline(source, destination)

def interpolate(source, x, y, z):
    r"""interpolate(vol source, float const x, float const y, float const z) -> float"""
    return _pytom_volume.interpolate(source, x, y, z)

def interpolateCubic(source, x, y, z):
    r"""interpolateCubic(vol source, float const x, float const y, float const z) -> float"""
    return _pytom_volume.interpolateCubic(source, x, y, z)

def interpolateSpline(source, x, y, z):
    r"""interpolateSpline(vol source, float const x, float const y, float const z) -> float"""
    return _pytom_volume.interpolateSpline(source, x, y, z)

def interpolateFourierSpline(source, x, y, z):
    r"""interpolateFourierSpline(vol_comp source, float const x, float const y, float const z) -> std::complex< float >"""
    return _pytom_volume.interpolateFourierSpline(source, x, y, z)

def backProject(*args):
    r"""
    backProject(vol src, vol dst, vol phi, vol theta, vol offsetPaticle, vol offsetProjections)
    backProject(vol src, vol dst, vol phi, vol theta, vol psi, vol offsetPaticle, vol offsetProjections)
    """
    return _pytom_volume.backProject(*args)

def complexRealMult(vol, otherVol):
    r"""complexRealMult(vol_comp vol, vol otherVol) -> vol_comp"""
    return _pytom_volume.complexRealMult(vol, otherVol)

def real(vol):
    r"""real(vol_comp vol) -> vol"""
    return _pytom_volume.real(vol)

def imag(vol):
    r"""imag(vol_comp vol) -> vol"""
    return _pytom_volume.imag(vol)

def mergeRealImag(real, imag):
    r"""mergeRealImag(vol real, vol imag) -> vol_comp"""
    return _pytom_volume.mergeRealImag(real, imag)

def volTOsf(vol, r, b, m_x, m_y, m_z):
    r"""volTOsf(vol vol, double const r, int const b, float const m_x, float const m_y, float const m_z) -> vol"""
    return _pytom_volume.volTOsf(vol, r, b, m_x, m_y, m_z)

def fvolTOsf(vol, r, b):
    r"""fvolTOsf(vol_comp vol, double const r, int const b) -> vol_comp"""
    return _pytom_volume.fvolTOsf(vol, r, b)

def sfFourierShift(vol, r, b, sizex, sizey, sizez, shiftX, shiftY, shiftZ):
    r"""sfFourierShift(vol_comp vol, double const r, int const b, int sizex, int sizey, int sizez, double shiftX, double shiftY, double shiftZ) -> vol_comp"""
    return _pytom_volume.sfFourierShift(vol, r, b, sizex, sizey, sizez, shiftX, shiftY, shiftZ)

def fullToReduced(*args):
    r"""
    fullToReduced(vol volume) -> vol
    fullToReduced(vol_comp volume) -> vol_comp
    """
    return _pytom_volume.fullToReduced(*args)


# This file was automatically generated by SWIG (http://www.swig.org).
# Version 3.0.12
#
# Do not make changes to this file unless you know what you are doing--modify
# the SWIG interface file instead.

from sys import version_info as _swig_python_version_info
if _swig_python_version_info >= (2, 7, 0):
    def swig_import_helper():
        import importlib
        pkg = __name__.rpartition('.')[0]
        mname = '.'.join((pkg, '_pytom_volume')).lstrip('.')
        try:
            return importlib.import_module(mname)
        except ImportError:
            return importlib.import_module('_pytom_volume')
    _pytom_volume = swig_import_helper()
    del swig_import_helper
elif _swig_python_version_info >= (2, 6, 0):
    def swig_import_helper():
        from os.path import dirname
        import imp
        fp = None
        try:
            fp, pathname, description = imp.find_module('_pytom_volume', [dirname(__file__)])
        except ImportError:
            import _pytom_volume
            return _pytom_volume
        try:
            _mod = imp.load_module('_pytom_volume', fp, pathname, description)
        finally:
            if fp is not None:
                fp.close()
        return _mod
    _pytom_volume = swig_import_helper()
    del swig_import_helper
else:
    import _pytom_volume
del _swig_python_version_info

try:
    _swig_property = property
except NameError:
    pass  # Python < 2.2 doesn't have 'property'.

try:
    import builtins as __builtin__
except ImportError:
    import __builtin__

def _swig_setattr_nondynamic(self, class_type, name, value, static=1):
    if (name == "thisown"):
        return self.this.own(value)
    if (name == "this"):
        if type(value).__name__ == 'SwigPyObject':
            self.__dict__[name] = value
            return
    method = class_type.__swig_setmethods__.get(name, None)
    if method:
        return method(self, value)
    if (not static):
        if _newclass:
            object.__setattr__(self, name, value)
        else:
            self.__dict__[name] = value
    else:
        raise AttributeError("You cannot add attributes to %s" % self)


def _swig_setattr(self, class_type, name, value):
    return _swig_setattr_nondynamic(self, class_type, name, value, 0)


def _swig_getattr(self, class_type, name):
    if (name == "thisown"):
        return self.this.own()
    method = class_type.__swig_getmethods__.get(name, None)
    if method:
        return method(self)
    raise AttributeError("'%s' object has no attribute '%s'" % (class_type.__name__, name))


def _swig_repr(self):
    try:
        strthis = "proxy of " + self.this.__repr__()
    except __builtin__.Exception:
        strthis = ""
    return "<%s.%s; %s >" % (self.__class__.__module__, self.__class__.__name__, strthis,)

try:
    _object = object
    _newclass = 1
except __builtin__.Exception:
    class _object:
        pass
    _newclass = 0

class vol(_object):
    """Proxy of C++ swigTom::swigVolume<(float,float)> class."""

    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, vol, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, vol, name)
    __repr__ = _swig_repr

    def __init__(self, *args):
        """
        __init__(swigTom::swigVolume<(float,float)> self, std::size_t sizex, std::size_t sizey, std::size_t sizez) -> vol
        __init__(swigTom::swigVolume<(float,float)> self, tom::Volume< float > const & v) -> vol
        __init__(swigTom::swigVolume<(float,float)> self, vol v) -> vol
        __init__(swigTom::swigVolume<(float,float)> self, float * data, std::size_t sizex, std::size_t sizey, std::size_t sizez, std::size_t stridex, std::size_t stridey, std::size_t stridez) -> vol
        """
        this = _pytom_volume.new_vol(*args)
        try:
            self.this.append(this)
        except __builtin__.Exception:
            self.this = this
    __swig_destroy__ = _pytom_volume.delete_vol
    __del__ = lambda self: None

    def size_x(self):
        """size_x(vol self) -> std::size_t"""
        return _pytom_volume.vol_size_x(self)


    def size_y(self):
        """size_y(vol self) -> std::size_t"""
        return _pytom_volume.vol_size_y(self)


    def size_z(self):
        """size_z(vol self) -> std::size_t"""
        return _pytom_volume.vol_size_z(self)


    def getFtSizeX(self):
        """getFtSizeX(vol self) -> float"""
        return _pytom_volume.vol_getFtSizeX(self)


    def getFtSizeY(self):
        """getFtSizeY(vol self) -> float"""
        return _pytom_volume.vol_getFtSizeY(self)


    def getFtSizeZ(self):
        """getFtSizeZ(vol self) -> float"""
        return _pytom_volume.vol_getFtSizeZ(self)


    def setFtSizeX(self, size_x):
        """setFtSizeX(vol self, float size_x)"""
        return _pytom_volume.vol_setFtSizeX(self, size_x)


    def setFtSizeY(self, size_y):
        """setFtSizeY(vol self, float size_y)"""
        return _pytom_volume.vol_setFtSizeY(self, size_y)


    def setFtSizeZ(self, size_z):
        """setFtSizeZ(vol self, float size_z)"""
        return _pytom_volume.vol_setFtSizeZ(self, size_z)


    def strideX(self):
        """strideX(vol self) -> std::size_t"""
        return _pytom_volume.vol_strideX(self)


    def strideY(self):
        """strideY(vol self) -> std::size_t"""
        return _pytom_volume.vol_strideY(self)


    def strideZ(self):
        """strideZ(vol self) -> std::size_t"""
        return _pytom_volume.vol_strideZ(self)


    def write(self, *args):
        """
        write(vol self, std::string fileName)
        write(vol self, std::string fileName, std::string fileType)
        """
        if len(args) == 1:
            try:
                args= [args[0], args[0].split('.')[-1]]
            except:
                pass
        return _pytom_volume.vol_write(self, *args)


    def info(self, name):
        """info(vol self, std::string const & name)"""
        return _pytom_volume.vol_info(self, name)


    def numelem(self):
        """numelem(vol self) -> std::size_t"""
        return _pytom_volume.vol_numelem(self)


    def equalsTo(self, v):
        """equalsTo(vol self, vol v) -> bool"""
        return _pytom_volume.vol_equalsTo(self, v)


    def getV(self, x, y, z):
        """getV(vol self, std::size_t x, std::size_t y, std::size_t z) -> float"""
        return _pytom_volume.vol_getV(self, x, y, z)


    def setV(self, val, x, y, z):
        """setV(vol self, float val, std::size_t x, std::size_t y, std::size_t z)"""
        return _pytom_volume.vol_setV(self, val, x, y, z)


    def setAll(self, val):
        """setAll(vol self, float val)"""
        return _pytom_volume.vol_setAll(self, val)


    def copyVolume(self, v2):
        """copyVolume(vol self, vol v2)"""
        return _pytom_volume.vol_copyVolume(self, v2)


    def shiftscale(self, shiftV, scaleV):
        """shiftscale(vol self, float const & shiftV, float const & scaleV)"""
        return _pytom_volume.vol_shiftscale(self, shiftV, scaleV)


    def dims(self):
        """dims(vol self) -> std::size_t"""
        return _pytom_volume.vol_dims(self)


    def __add__(self, *args):
        """
        __add__(vol self, vol otherVol) -> vol
        __add__(vol self, float const & value) -> vol
        """
        return _pytom_volume.vol___add__(self, *args)


    def __mul__(self, *args):
        """
        __mul__(vol self, vol otherVol) -> vol
        __mul__(vol self, float const & value) -> vol
        """
        return _pytom_volume.vol___mul__(self, *args)


    def __sub__(self, *args):
        """
        __sub__(vol self, vol otherVol) -> vol
        __sub__(vol self, float const & value) -> vol
        """
        return _pytom_volume.vol___sub__(self, *args)


    def __truediv__(self, *args):
        return _pytom_volume.vol___truediv__(self, *args)
    __div__ = __truediv__



    def __call__(self, *args):
        """
        __call__(vol self, std::size_t const & x, std::size_t const & y, std::size_t const & z) -> float const
        __call__(vol self, float const & value, std::size_t const & x, std::size_t const & y, std::size_t const & z)
        """
        return _pytom_volume.vol___call__(self, *args)

vol_swigregister = _pytom_volume.vol_swigregister
vol_swigregister(vol)

class vol_comp(_object):
    """Proxy of C++ swigTom::swigVolume<(std::complex<(float)>,float)> class."""

    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, vol_comp, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, vol_comp, name)
    __repr__ = _swig_repr

    def __init__(self, *args):
        """
        __init__(swigTom::swigVolume<(std::complex<(float)>,float)> self, std::size_t sizex, std::size_t sizey, std::size_t sizez) -> vol_comp
        __init__(swigTom::swigVolume<(std::complex<(float)>,float)> self, tom::Volume< std::complex< float > > const & v) -> vol_comp
        __init__(swigTom::swigVolume<(std::complex<(float)>,float)> self, vol_comp v) -> vol_comp
        __init__(swigTom::swigVolume<(std::complex<(float)>,float)> self, std::complex< float > * data, std::size_t sizex, std::size_t sizey, std::size_t sizez, std::size_t stridex, std::size_t stridey, std::size_t stridez) -> vol_comp
        """
        this = _pytom_volume.new_vol_comp(*args)
        try:
            self.this.append(this)
        except __builtin__.Exception:
            self.this = this
    __swig_destroy__ = _pytom_volume.delete_vol_comp
    __del__ = lambda self: None

    def size_x(self):
        """size_x(vol_comp self) -> std::size_t"""
        return _pytom_volume.vol_comp_size_x(self)


    def size_y(self):
        """size_y(vol_comp self) -> std::size_t"""
        return _pytom_volume.vol_comp_size_y(self)


    def size_z(self):
        """size_z(vol_comp self) -> std::size_t"""
        return _pytom_volume.vol_comp_size_z(self)


    def getFtSizeX(self):
        """getFtSizeX(vol_comp self) -> float"""
        return _pytom_volume.vol_comp_getFtSizeX(self)


    def getFtSizeY(self):
        """getFtSizeY(vol_comp self) -> float"""
        return _pytom_volume.vol_comp_getFtSizeY(self)


    def getFtSizeZ(self):
        """getFtSizeZ(vol_comp self) -> float"""
        return _pytom_volume.vol_comp_getFtSizeZ(self)


    def setFtSizeX(self, size_x):
        """setFtSizeX(vol_comp self, float size_x)"""
        return _pytom_volume.vol_comp_setFtSizeX(self, size_x)


    def setFtSizeY(self, size_y):
        """setFtSizeY(vol_comp self, float size_y)"""
        return _pytom_volume.vol_comp_setFtSizeY(self, size_y)


    def setFtSizeZ(self, size_z):
        """setFtSizeZ(vol_comp self, float size_z)"""
        return _pytom_volume.vol_comp_setFtSizeZ(self, size_z)


    def strideX(self):
        """strideX(vol_comp self) -> std::size_t"""
        return _pytom_volume.vol_comp_strideX(self)


    def strideY(self):
        """strideY(vol_comp self) -> std::size_t"""
        return _pytom_volume.vol_comp_strideY(self)


    def strideZ(self):
        """strideZ(vol_comp self) -> std::size_t"""
        return _pytom_volume.vol_comp_strideZ(self)


    def write(self, *args):
        """
        write(vol_comp self, std::string fileName)
        write(vol_comp self, std::string fileName, std::string fileType)
        """
        return _pytom_volume.vol_comp_write(self, *args)


    def info(self, name):
        """info(vol_comp self, std::string const & name)"""
        return _pytom_volume.vol_comp_info(self, name)


    def numelem(self):
        """numelem(vol_comp self) -> std::size_t"""
        return _pytom_volume.vol_comp_numelem(self)


    def equalsTo(self, v):
        """equalsTo(vol_comp self, vol_comp v) -> bool"""
        return _pytom_volume.vol_comp_equalsTo(self, v)


    def getV(self, x, y, z):
        """getV(vol_comp self, std::size_t x, std::size_t y, std::size_t z) -> std::complex< float >"""
        return _pytom_volume.vol_comp_getV(self, x, y, z)


    def setV(self, val, x, y, z):
        """setV(vol_comp self, std::complex< float > val, std::size_t x, std::size_t y, std::size_t z)"""
        return _pytom_volume.vol_comp_setV(self, val, x, y, z)


    def setAll(self, val):
        """setAll(vol_comp self, std::complex< float > val)"""
        return _pytom_volume.vol_comp_setAll(self, val)


    def copyVolume(self, v2):
        """copyVolume(vol_comp self, vol_comp v2)"""
        return _pytom_volume.vol_comp_copyVolume(self, v2)


    def shiftscale(self, shiftV, scaleV):
        """shiftscale(vol_comp self, float const & shiftV, float const & scaleV)"""
        return _pytom_volume.vol_comp_shiftscale(self, shiftV, scaleV)


    def dims(self):
        """dims(vol_comp self) -> std::size_t"""
        return _pytom_volume.vol_comp_dims(self)


    def __add__(self, *args):
        """
        __add__(vol_comp self, vol_comp otherVol) -> vol_comp
        __add__(vol_comp self, float const & value) -> vol_comp
        """
        return _pytom_volume.vol_comp___add__(self, *args)


    def __mul__(self, *args):
        """
        __mul__(vol_comp self, vol_comp otherVol) -> vol_comp
        __mul__(vol_comp self, float const & value) -> vol_comp
        """
        return _pytom_volume.vol_comp___mul__(self, *args)


    def __sub__(self, *args):
        """
        __sub__(vol_comp self, vol_comp otherVol) -> vol_comp
        __sub__(vol_comp self, float const & value) -> vol_comp
        """
        return _pytom_volume.vol_comp___sub__(self, *args)


    def __truediv__(self, *args):
        return _pytom_volume.vol_comp___truediv__(self, *args)
    __div__ = __truediv__



    def __call__(self, *args):
        """
        __call__(vol_comp self, std::size_t const & x, std::size_t const & y, std::size_t const & z) -> std::complex< float > const
        __call__(vol_comp self, std::complex< float > const & value, std::size_t const & x, std::size_t const & y, std::size_t const & z)
        """
        return _pytom_volume.vol_comp___call__(self, *args)

vol_comp_swigregister = _pytom_volume.vol_comp_swigregister
vol_comp_swigregister(vol_comp)


def read(*args):
    """
    read(std::string fileName) -> vol
    read(std::string fileName, std::size_t subregion1, std::size_t subregion2, std::size_t subregion3, std::size_t subregion4, std::size_t subregion5, std::size_t subregion6, std::size_t sampling1, std::size_t sampling2, std::size_t sampling3, std::size_t binning1, std::size_t binning2, std::size_t binning3) -> vol
    """
    return _pytom_volume.read(*args)

def conj_mult(v1, v2):
    """conj_mult(vol_comp v1, vol_comp v2)"""
    return _pytom_volume.conj_mult(v1, v2)

def initSphere(vol, radius, sigma, max_radius, centerx, centery, centerz):
    """initSphere(vol vol, float radius, float sigma, float max_radius, float centerx, float centery, float centerz)"""
    return _pytom_volume.initSphere(vol, radius, sigma, max_radius, centerx, centery, centerz)

def maskNorm(vol, mask, is_boolean_mask):
    """maskNorm(vol vol, vol mask, bool is_boolean_mask)"""
    return _pytom_volume.maskNorm(vol, mask, is_boolean_mask)

def rotate(src, dst, phi, psi, theta):
    """rotate(vol src, vol dst, double phi, double psi, double theta)"""
    return _pytom_volume.rotate(src, dst, phi, psi, theta)

def rotateCubic(src, dst, phi, psi, theta):
    """rotateCubic(vol src, vol dst, double phi, double psi, double theta)"""
    return _pytom_volume.rotateCubic(src, dst, phi, psi, theta)

def rotateSpline(src, dst, phi, psi, theta):
    """rotateSpline(vol src, vol dst, double phi, double psi, double theta)"""
    return _pytom_volume.rotateSpline(src, dst, phi, psi, theta)

def rotateSplineInFourier(src, phi, psi, theta):
    """rotateSplineInFourier(vol_comp src, double phi, double psi, double theta) -> vol_comp"""
    return _pytom_volume.rotateSplineInFourier(src, phi, psi, theta)

def shift(src, dst, preX, preY, preZ):
    """shift(vol src, vol dst, double preX, double preY, double preZ)"""
    return _pytom_volume.shift(src, dst, preX, preY, preZ)

def shiftFourier(src, dst, shiftX, shiftY, shiftZ):
    """shiftFourier(vol_comp src, vol_comp dst, double shiftX, double shiftY, double shiftZ)"""
    return _pytom_volume.shiftFourier(src, dst, shiftX, shiftY, shiftZ)

def transform(src, dst, phi, psi, theta, centerX, centerY, centerZ, preX, preY, preZ, postX, postY, postZ):
    """transform(vol src, vol dst, double phi, double psi, double theta, double centerX, double centerY, double centerZ, double preX, double preY, double preZ, double postX, double postY, double postZ)"""
    return _pytom_volume.transform(src, dst, phi, psi, theta, centerX, centerY, centerZ, preX, preY, preZ, postX, postY, postZ)

def general_transform(src, dst, mtx):
    """general_transform(vol src, vol dst, vol mtx)"""
    return _pytom_volume.general_transform(src, dst, mtx)

def transformCubic(src, dst, phi, psi, theta, centerX, centerY, centerZ, preX, preY, preZ, postX, postY, postZ):
    """transformCubic(vol src, vol dst, double phi, double psi, double theta, double centerX, double centerY, double centerZ, double preX, double preY, double preZ, double postX, double postY, double postZ)"""
    return _pytom_volume.transformCubic(src, dst, phi, psi, theta, centerX, centerY, centerZ, preX, preY, preZ, postX, postY, postZ)

def transformSpline(src, dst, phi, psi, theta, centerX, centerY, centerZ, preX, preY, preZ, postX, postY, postZ):
    """transformSpline(vol src, vol dst, double phi, double psi, double theta, double centerX, double centerY, double centerZ, double preX, double preY, double preZ, double postX, double postY, double postZ)"""
    return _pytom_volume.transformSpline(src, dst, phi, psi, theta, centerX, centerY, centerZ, preX, preY, preZ, postX, postY, postZ)

def transformFourierSpline(src, phi, psi, theta, shiftX, shiftY, shiftZ):
    """transformFourierSpline(vol_comp src, double phi, double psi, double theta, double shiftX, double shiftY, double shiftZ) -> vol_comp"""
    return _pytom_volume.transformFourierSpline(src, phi, psi, theta, shiftX, shiftY, shiftZ)

def peak(*args):
    """
    peak(vol volume) -> tom::st_idx
    peak(vol volume, vol mask) -> tom::st_idx
    """
    return _pytom_volume.peak(*args)

def conjugate(volume):
    """conjugate(vol_comp volume)"""
    return _pytom_volume.conjugate(volume)

def paste(volume, destination, x, y, z):
    """paste(vol volume, vol destination, long x, long y, long z)"""
    return _pytom_volume.paste(volume, destination, x, y, z)

def pasteCenter(volume, destination):
    """pasteCenter(vol volume, vol destination)"""
    return _pytom_volume.pasteCenter(volume, destination)

def power(*args):
    """
    power(vol volume, float const & exponent)
    power(vol_comp volume, float const & exponent)
    """
    return _pytom_volume.power(*args)

def numberSetVoxels(volume):
    """numberSetVoxels(vol volume) -> std::size_t"""
    return _pytom_volume.numberSetVoxels(volume)

def sum(*args):
    """
    sum(vol volume) -> float
    sum(vol_comp volume) -> std::complex< float >
    """
    return _pytom_volume.sum(*args)

def projectSum(volume, axis):
    """projectSum(vol volume, int axis) -> vol"""
    return _pytom_volume.projectSum(volume, axis)

def mean(volume):
    """mean(vol volume) -> double"""
    return _pytom_volume.mean(volume)

def variance(volume, use_sample_standard_deviation):
    """variance(vol volume, bool use_sample_standard_deviation) -> double"""
    return _pytom_volume.variance(volume, use_sample_standard_deviation)

def min(volume):
    """min(vol volume) -> float"""
    return _pytom_volume.min(volume)

def max(volume):
    """max(vol volume) -> float"""
    return _pytom_volume.max(volume)

def limit(volume, lowerBound, lowerReplacement, upperBound, upperReplacement, doLower, doUpper):
    """limit(vol volume, float lowerBound, float lowerReplacement, float upperBound, float upperReplacement, bool doLower, bool doUpper)"""
    return _pytom_volume.limit(volume, lowerBound, lowerReplacement, upperBound, upperReplacement, doLower, doUpper)

def abs(*args):
    """
    abs(vol volume) -> vol
    abs(vol_comp volume) -> vol_comp
    """
    return _pytom_volume.abs(*args)

def reducedToFull(*args):
    """
    reducedToFull(vol volume) -> vol
    reducedToFull(vol_comp volume) -> vol_comp
    """
    return _pytom_volume.reducedToFull(*args)

def gaussianNoise(volume, mean, sigma):
    """gaussianNoise(vol volume, double mean, double sigma)"""
    return _pytom_volume.gaussianNoise(volume, mean, sigma)

def complexDiv(volume, div):
    """complexDiv(vol_comp volume, vol div) -> vol_comp"""
    return _pytom_volume.complexDiv(volume, div)

def vectorize(source):
    """vectorize(vol source) -> vol"""
    return _pytom_volume.vectorize(source)

def subvolume(*args):
    """
    subvolume(vol source, std::size_t startX, std::size_t startY, std::size_t startZ, std::size_t endX, std::size_t endY, std::size_t endZ) -> vol
    subvolume(vol_comp source, std::size_t startX, std::size_t startY, std::size_t startZ, std::size_t endX, std::size_t endY, std::size_t endZ) -> vol_comp
    """
    return _pytom_volume.subvolume(*args)

def putSubVolume(*args):
    """
    putSubVolume(vol source, vol destination, std::size_t positionX, std::size_t positionY, std::size_t positionZ)
    putSubVolume(vol_comp source, vol_comp destination, std::size_t positionX, std::size_t positionY, std::size_t positionZ)
    """
    return _pytom_volume.putSubVolume(*args)

def writeSubregion(source, filename, positionX, positionY, positionZ):
    """writeSubregion(vol source, std::string filename, std::size_t positionX, std::size_t positionY, std::size_t positionZ)"""
    return _pytom_volume.writeSubregion(source, filename, positionX, positionY, positionZ)

def updateResFromIdx(resV, newV, orientV, orientIdx):
    """updateResFromIdx(vol resV, vol newV, vol orientV, std::size_t const orientIdx)"""
    return _pytom_volume.updateResFromIdx(resV, newV, orientV, orientIdx)

def updateResFromVol(resV, newV, orientV, neworientV):
    """updateResFromVol(vol resV, vol newV, vol orientV, vol neworientV)"""
    return _pytom_volume.updateResFromVol(resV, newV, orientV, neworientV)

def mirrorVolume(src, des):
    """mirrorVolume(vol src, vol des)"""
    return _pytom_volume.mirrorVolume(src, des)

def rescale(source, destination):
    """rescale(vol source, vol destination)"""
    return _pytom_volume.rescale(source, destination)

def rescaleCubic(source, destination):
    """rescaleCubic(vol source, vol destination)"""
    return _pytom_volume.rescaleCubic(source, destination)

def rescaleSpline(source, destination):
    """rescaleSpline(vol source, vol destination)"""
    return _pytom_volume.rescaleSpline(source, destination)

def interpolate(source, x, y, z):
    """interpolate(vol source, float const x, float const y, float const z) -> float"""
    return _pytom_volume.interpolate(source, x, y, z)

def interpolateCubic(source, x, y, z):
    """interpolateCubic(vol source, float const x, float const y, float const z) -> float"""
    return _pytom_volume.interpolateCubic(source, x, y, z)

def interpolateSpline(source, x, y, z):
    """interpolateSpline(vol source, float const x, float const y, float const z) -> float"""
    return _pytom_volume.interpolateSpline(source, x, y, z)

def interpolateFourierSpline(source, x, y, z):
    """interpolateFourierSpline(vol_comp source, float const x, float const y, float const z) -> std::complex< float >"""
    return _pytom_volume.interpolateFourierSpline(source, x, y, z)

def backProject(*args):
    """
    backProject(vol src, vol dst, vol phi, vol theta, vol offsetPaticle, vol offsetProjections)
    backProject(vol src, vol dst, vol phi, vol theta, vol psi, vol offsetPaticle, vol offsetProjections)
    """
    return _pytom_volume.backProject(*args)

def complexRealMult(vol, otherVol):
    """complexRealMult(vol_comp vol, vol otherVol) -> vol_comp"""
    return _pytom_volume.complexRealMult(vol, otherVol)

def real(vol):
    """real(vol_comp vol) -> vol"""
    return _pytom_volume.real(vol)

def imag(vol):
    """imag(vol_comp vol) -> vol"""
    return _pytom_volume.imag(vol)

def mergeRealImag(real, imag):
    """mergeRealImag(vol real, vol imag) -> vol_comp"""
    return _pytom_volume.mergeRealImag(real, imag)

def volTOsf(vol, r, b, m_x, m_y, m_z):
    """volTOsf(vol vol, double const r, int const b, float const m_x, float const m_y, float const m_z) -> vol"""
    return _pytom_volume.volTOsf(vol, r, b, m_x, m_y, m_z)

def fvolTOsf(vol, r, b):
    """fvolTOsf(vol_comp vol, double const r, int const b) -> vol_comp"""
    return _pytom_volume.fvolTOsf(vol, r, b)

def sfFourierShift(vol, r, b, sizex, sizey, sizez, shiftX, shiftY, shiftZ):
    """sfFourierShift(vol_comp vol, double const r, int const b, int sizex, int sizey, int sizez, double shiftX, double shiftY, double shiftZ) -> vol_comp"""
    return _pytom_volume.sfFourierShift(vol, r, b, sizex, sizey, sizez, shiftX, shiftY, shiftZ)

def fullToReduced(*args):
    """
    fullToReduced(vol volume) -> vol
    fullToReduced(vol_comp volume) -> vol_comp
    """
    return _pytom_volume.fullToReduced(*args)
# This file is compatible with both classic and new-style classes.



"""This module defines agnostic correlation functions

Created: unknown date
@author: FF

restructured by @sroet in May 2023
"""
from pytom.gpu.initialize import xp, device
from pytom.gpu.gpuFunctions import argmax
from pytom.agnostic.normalise import (
    meanUnderMask,
    meanVolUnderMask,
    stdUnderMask,
    stdVolUnderMask,
    normaliseUnderMask,
    mean0std1,
)
from pytom.agnostic.tools import create_circle, create_sphere
#Typing imports
from pytom.gpu.initialize import xptyping as xpt
from typing import Optional, Tuple


def flcf(volume: xpt.NDArray, template: xpt.NDArray, mask: Optional[xpt.NDArray]=None, std_v: Optional[float]=None) -> float:
    """Fast local correlation function

    @param volume: target volume
    @param template: template to be searched. It can have smaller size then target volume.
    @param mask: template mask. If not given, a default sphere mask will be used.
    @param std_v: standard deviation of the target volume under mask, which do not need to be
                  calculated again when the mask is identical.

    @return: the local correlation function
    """

    if mask is None:
        radius = volume.shape[0] // 2 - 3
        if len(volume.shape) == 2:
            mask = create_circle(volume.shape, radius, 3)
        elif len(volume.shape) == 3:
            mask = create_sphere(volume.shape, radius, 3)

    p = mask.sum()
    mean_t = meanUnderMask(template, mask, p=p)
    temp = ((template - mean_t) / stdUnderMask(template, mask, mean_t, p=p)) * mask

    if std_v is None:
        mean_v = meanVolUnderMask(volume, mask)
        std_v = stdVolUnderMask(volume, mask, mean_v)

    res = (
        xp.fft.fftshift(
            xp.fft.ifftn(xp.conj(xp.fft.fftn(temp)) * xp.fft.fftn(volume))
        ).real
        / std_v
        / p
    )
    return res


def xcc(volume: xpt.NDArray, template: xpt.NDArray, mask: Optional[xpt.NDArray]=None, volume_is_normalized: bool=False) -> float:
    """
    xcc: Calculates the cross correlation coefficient in real space
    @param volume: A volume
    @type volume:  L{xp.ndarray}
    @param template: A template that is searched in volume. Must be of same size as volume.
    @type template:  L{xp.ndarray}
    @param mask: mask to constrain correlation
    @type mask: L{xp.ndarray}
    @param volume_is_normalized: only used for compatibility with nxcc - not used
    @type volume_is_normalized: L{bool}
    @return: An unscaled value
    @raise exception: Raises a runtime error if volume and template have a different size.
    @author: Thomas Hrabe
    """

    from pytom.agnostic.macros import volumesSameSize

    if not volumesSameSize(volume, template):
        raise RuntimeError("Volume and template must have same size!")

    if mask:  # mask is given
        result = mask * volume * template
    else:
        result = volume * template

    cc = result.sum()

    cc = cc / float(volume.size)

    return cc


def nxcc(volume: xpt.NDArray, template: xpt.NDArray, mask: Optional[xpt.NDArray]=None, volume_is_normalized: bool=False) -> float:
    """
    nxcc: Calculates the normalized cross correlation coefficient in real space
    @param volume: A volume
    @type volume:  L{xp.ndarray}
    @param template: A template that is searched in volume. Must be of same size as volume.
    @type template:  L{xp.ndarray}
    @param mask: mask to constrain correlation
    @type mask: L{xp.ndarray}
    @param volume_is_normalized: speed up if volume is already normalized
    @type volume_is_normalized: L{bool}
    @return: A value between -1 and 1
    @raise exception: Raises a runtime error if volume and template have a different size.
    @author: Thomas Hrabe
    @change: flag for pre-normalized volume, FF
    """

    from pytom.agnostic.macros import volumesSameSize

    if not volumesSameSize(volume, template):
        raise RuntimeError("Volume and template must have same size!")

    if mask is None:
        if not volume_is_normalized:
            v = mean0std1(volume, True)
        else:
            v = volume
        t = mean0std1(template, True)
        p = volume.size
        result = v * t
    else:
        from pytom.agnostic.normalise import normaliseUnderMask

        if not volume_is_normalized:
            (v, p) = normaliseUnderMask(volume, mask)
            (t, p) = normaliseUnderMask(template, mask, p)
            t = t * mask  # multiply with the mask
            result = v * t
        else:
            (t, p) = normaliseUnderMask(template, mask)
            t = t * mask  # multiply with the mask
            result = volume * t

    ncc = result.sum()
    ncc = ncc / float(p)

    return float(ncc)


def xcf(
        volume: xpt.NDArray,
        template: xpt.NDArray,
        mask: Optional[xpt.NDArray]=None,
        std_v: Optional[float]=None,
) -> float:
    """
    XCF: returns the non-normalised cross correlation function. The xcf
    result is scaled only by the square of the number of elements.

    @param volume : The search volume
    @type volume: L{xp.ndarray}
    @param template : The template searched (this one will be used for conjugate complex
                                             multiplication)
    @type template: L{xp.ndarray}
    @param mask: Will be unused, only for compatibility reasons with FLCF
    @param std_v: Will be unused, only for compatibility reasons with FLCF
    @return: XCF volume
    @rtype: L{xp.ndarray}
    @author: Thomas Hrabe
    """

    # if mask:
    #    volume = volume * mask
    # 	 template = template * mask

    # determine fourier transforms of volumes
    if template.dtype in (xp.float64, xp.float32):
        fvolume = xp.fft.fftn(volume)
    else:
        fvolume = volume

    if template.dtype in (xp.float32, xp.float64):
        ftemplate = xp.fft.fftn(template)
    else:
        ftemplate = template

    # perform element wise - conjugate multiplication
    ftemplate = xp.conj(ftemplate)
    fresult = fvolume * ftemplate

    # transform back to real space
    result = abs(xp.fft.ifftn(fresult))

    # xp.fft.iftshift(result)

    return result


def xcf_mult(
        volume: xpt.NDArray,
        template: xpt.NDArray,
        mask: xpt.NDArray,
        std_v: Optional[float]=None,
) -> float:
    """
    XCF: returns the non-normalised cross correlation function. The xcf
    result is scaled only by the square of the number of elements.

    @param volume : The search volume
    @type volume: L{pytom.lib.pytom_volume.vol}
    @param template : The template searched (this one will be used for conjugate complex
                                             multiplication)
    @type template: L{pytom.lib.pytom_volume.vol}
    @param mask: Will be unused, only for compatibility reasons with FLCF
    @param std_v: Will be unused, only for compatibility reasons with FLCF
    @return: XCF volume
    @rtype: L{pytom.lib.pytom_volume.vol}
    @author: Thomas Hrabe
    """

    # if mask:
    #    volume = volume * mask
    # 	 template = template * mask

    # determine fourier transforms of volumes
    fvolume = xp.fft.fftn(volume)
    ftemplate = xp.fft.fftn(template)

    # perform element wise - conjugate multiplication
    ftemplate = xp.conj(ftemplate)
    fresult = fvolume * ftemplate

    # transform back to real space
    result = abs(xp.fft.ifftn(fresult))

    # xp.fft.iftshift(result)

    return result


def norm_xcf(volume: xpt.NDArray, template: xpt.NDArray, mask: Optional[xpt.NDArray]=None, std_v: Optional[float]=None, gpu: bool=False) ->  float:
    """
    nXCF: returns the normalised cross correlation function. Autocorrelation
    of two equal objects would yield a max nxcf peak of 1.

    @param volume: The search volume
    @param template: The template searched (this one will be used for conjugate complex
                                            multiplication)
    @type template: L{xp.ndarray}
    @param mask: template mask. If not given, a default sphere mask will be generated which has the
                 same size with the given template.
    @type mask: L{xp.ndarray}
    @param std_v: Will be unused, only for compatibility reasons with FLCF
    @return: the calculated norm_xcf volume
    @rtype: L{xp.ndarray}
    @author: Thomas Hrabe
    @change: masking of template implemented
    """

    if mask is not None:
        result = xcf_mult(
            normaliseUnderMask(volume=volume, mask=mask, p=None),
            normaliseUnderMask(volume=template, mask=mask, p=None),
            mask,
        )
    else:
        result = xcf_mult(
            mean0std1(volume, True), mean0std1(template, True), mask, True
        )

    return result



def soc(volume: xpt.NDArray, reference: xpt.NDArray, mask: Optional[xpt.NDArray]=None, std_v: Optional[float]=None) -> float:
    """
    soc : Second Order Correlation. Correlation of correlation peaks.
    @param volume: The volume
    @type volume:  L{xp.ndarray}
    @param reference: The reference / template
    @type reference:  L{xp.ndarray}
    @author: Thomas Hrabe
    """

    reference_peak = flcf(reference, reference, mask)
    peaks = flcf(volume, reference, mask)

    return flcf(peaks, reference_peak, mask)


def dev(volume: xpt.NDArray, template: xpt.NDArray, mask: Optional[xpt.NDArray]=None, volume_is_normalized: bool=False) -> float:
    """
    dev: Calculates the squared deviation of volume and template in real space
    @param volume: A volume
    @type volume:  L{xp.ndarray}
    @param template: A template that is searched in volume. Must be of same size as volume.
    @type template:  L{xp.ndarray}
    @param mask: mask to constrain correlation
    @type mask: L{xp.ndarray}
    @param volume_is_normalized: speed up if volume is already normalized
    @type volume_is_normalized: L{bool}
    @return: deviation
    @raise exception: Raises a runtime error if volume and template have a different size.
    @author: FF
    """

    from pytom.agnostic.tools import volumesSameSize

    assert volume.__class__ == xp.ndarray, "dev: volume has to be of type xp.ndarray!"
    assert (
        template.__class__ == xp.ndarray
    ), "dev: template has to be of type xp.ndarray!"
    if not volumesSameSize(volume, template):
        raise RuntimeError("Volume and template must have same size!")

    if not mask:
        p = volume.size
        result = volume - template
    else:
        assert mask.__class__ == xp.ndarray, "dev: mask has to be of type xp.ndarray!"
        p = xp.sum(mask)
        result = (volume - template) * mask

    deviat = sum(result**2)
    deviat = deviat / float(p)

    return deviat


# Band correlation


def band_cc(volume: xpt.NDArray, reference: xpt.NDArray, band: Tuple[float, float], verbose: bool=False, shared: None=None, index: None=None)-> float:
    """
    band_cc: Determines the normalised correlation coefficient within a band
    @param volume: The volume
    @type volume: L{xp.ndarray}
    @param reference: The reference
    @type reference: L{xp.ndarray}
    @param band: [a,b] - specify the lower and upper end of band.
    @return: The correlation coefficient of the two volumes in the specified band.
    @rtype: L{float}
    @author: GS
    """

    if index is not None:
        print(index)

    from pytom.agnostic.filter import bandpass

    # from pytom.agnostic.filter import vol_comp

    if verbose:
        print("lowest freq : ", band[0], " highest freq", band[1])

    vf, m = bandpass(volume, band[0], band[1], returnMask=True, fourierOnly=True)
    rf = bandpass(reference, band[0], band[1], mask=m, fourierOnly=True)  # ,vf[1])

    vf = vf.astype(xp.complex128)
    cc_volume = rf.astype(vf.dtype)

    cc_volume = cc_volume * xp.conj(vf)

    cc = cc_volume.sum()

    cc = cc.real
    v = vf
    r = rf

    abs_v = xp.abs(v)
    abs_r = xp.abs(r)

    sum_v = xp.sum(abs_v**2)
    sum_r = xp.sum(abs_r**2)

    sum_v = xp.abs(sum_v)
    sum_r = xp.abs(sum_r)

    if sum_v == 0:
        sum_v = 1

    if sum_r == 0:
        sum_r = 1

    cc = cc / (xp.sqrt(sum_v * sum_r))

    # numerical errors will be punished with nan
    nan_treshold = 1.1
    if abs(cc) > nan_treshold:
        cc = float("nan")

    return float(cc)


def band_cf(volume: xpt.NDArray, reference: xpt.NDArray, band: Tuple[float, float]=(0, 100)) -> Tuple[xpt.NDArray, xpt.NDArray]:
    """
    band_cf:
    @param volume: The volume
    @param reference: The reference
    @param band: [a,b] - specify the lower and upper end of band. [0,1] if not set.
    @return: First parameter - The correlation of the two volumes in the specified ring.
             Second parameter - The bandpass filter used.
    @rtype: List - [L{xp.ndarray},L{pytom.lib.pytom_freqweight.weight}]
    @author: Thomas Hrabe
    @todo: does not work yet -> test is disabled
    """

    from math import sqrt
    from pytom.agnostic.filter import bandpass
    from pytom.agnostic.transform import fourier_reduced2full

    vf, vfm = bandpass(volume, band[0], band[1], fourierOnly=True, returnMask=True)
    rf = bandpass(reference, band[0], band[1], vf[1], fourierOnly=True)

    v = fourier_reduced2full(vf)
    r = fourier_reduced2full(rf)

    abs_v = xp.abs(v) ** 2
    abs_r = xp.abs(r) ** 2

    sum_v = abs(xp.sum(abs_v))
    sum_r = abs(xp.sum(abs_r))

    if sum_v == 0:
        sum_v = 1

    if sum_r == 0:
        sum_r = 1

    xp.conj(rf[0])

    fresult = vf * rf

    # transform back to real space
    result = xp.fft.ifftshift(xp.fft.ifftn(fresult)).real

    result *= 1 / float(sqrt(sum_v * sum_r))

    return result, vfm


# Fourier Shell Correlation and helper functions


def fsc(volume1, volume2, number_bands=None, mask=None, verbose=False, filename=None):
    """
    FSC - Calculates the Fourier Shell Correlation for two volumes
    @param volume1: volume one
    @type volume1: L{xp.ndarray}
    @param volume2: volume two
    @type volume2: L{xp.ndarray}
    @param number_bands: number of shells for FSC
    @type number_bands: int
    @param mask: mask
    @type mask: L{xp.ndarray}
    @param verbose: flag to activate printing of info
    @type verbose: L{bool}
    @param filename: write FSC to ascii file if specified
    @type filename: string

    @return: Returns a list of cc values
    @author:GS
    @rtype: list[floats]
    """

    from pytom.agnostic.correlation import band_cc
    from pytom.basic.structures import Mask
    from pytom.agnostic.io import read
    from pytom.agnostic.tools import volumesSameSize
    import time

    time.time()

    if not volumesSameSize(volume1, volume2):
        raise RuntimeError("Volumes must have the same size!")

    number_bands = volume1.shape[0] // 2 if number_bands is None else number_bands

    if not mask is None:
        if mask.__class__ == xp.array([0]).__class__:
            volume1 = volume1 * mask
            volume2 = volume2 * mask

        elif mask.__class__ == Mask:
            mask = mask.getVolume()
            volume1 = volume1 * mask
            volume2 = volume2 * mask

        elif mask.__class__ == str:
            mask = read(mask)
            volume1 = volume1 * mask
            volume2 = volume2 * mask

        else:
            raise RuntimeError(
                "FSC: Mask must be a volume OR a Mask object OR a string path to a mask!"
            )

    fsc_result = []
    band = [-1, -1]

    increment = int(volume1.shape[0] / 2 * 1 / number_bands)
    import time

    fvolume1 = xp.fft.fftn(volume1)
    fvolume2 = xp.fft.fftn(volume2)

    for n, i in enumerate(range(0, volume1.shape[0] // 2, increment)):
        band[0] = i
        band[1] = i + increment

        if verbose:
            print("Band : ", band)

        res = band_cc(fvolume1, fvolume2, band, verbose)

        if i == 0 and increment == 1:
            # force a 1 for correlation of the zero frequency
            res = 1

        if verbose:
            print("Correlation ", res)

        fsc_result.append(res)

    if filename:
        f = open(filename, "w")
        for item in fsc_result:
            f.write("%s\n" % item)
        f.close()

    return fsc_result


def fsc_sum(volume, reference, number_of_bands, wedge_angle=-1):
    """
    fsc_sum: Determines the sum of the Fourier Shell Correlation coefficient for a volume and
            reference.
    @param volume: A volume
    @type volume: L{xp.ndarray}
    @param reference: A reference of same size as volume
    @type reference: L{xp.ndarray}
    @param number_of_bands: Number of bands
    @param wedge_angle: A optional wedge angle
    @return: The sum FSC coefficient
    @rtype: float
    @author: Thomas Hrabe
    """

    from pytom.agnostic.correlation import band_cc

    result = 0

    fvolume = xp.fft.fftn(volume)
    freference = xp.fft.fftn(reference)
    numelem = volume.size

    fvolume *= 1 / float(numelem)
    freference *= 1 / float(numelem)

    for i in range(number_of_bands):
        # process bandCorrelation
        band = []
        band[0] = i * volume.shape[0] / number_of_bands
        band[1] = (i + 1) * volume.shape[0] / number_of_bands

        r = band_cc(fvolume, freference, band)
        cc = r[0]
        result = result + cc

    return result * (1 / float(number_of_bands))


def determine_resolution(fsc, resolution_criterion, verbose=False, randomized_fsc=None):
    """
    determine_resolution: Determines frequency and band where correlation drops below the
                         resolution_criterion. Uses linear interpolation between two positions
    @param fsc: The fsc list determined by L{pytom.basic.correlation.FSC}
    @param resolution_criterion: A value between 0 and 1
    @param verbose: Bool that activate writing of info, default=False
    @param randomized_fsc: A value that sets the start of the calculation of randomized FSC. (0-1).
    @return: [resolution,interpolated_band,number_bands]
    @author: Thomas Hrabe
    @todo: Add test!
    """

    fsc = xp.array(fsc)
    number_bands = len(fsc)

    band = number_bands

    if randomized_fsc is None:
        randomized_fsc = xp.ones_like(fsc) * (fsc.min() - 0.1)

    for i in range(number_bands):
        if fsc[i] < resolution_criterion and fsc[i] > randomized_fsc[i]:
            band = i - 1  # select the band that is still larger than criterion
            break

    if verbose:
        print("Band detected at ", band)

    if band == -1:
        raise RuntimeError("Please check your resolution criterion or you FSC!")

    elif band < number_bands:
        fsc1 = fsc[band]
        fsc2 = fsc[band + 1]

        rfsc1 = randomized_fsc[band]
        rfsc2 = randomized_fsc[band + 1]

        try:
            if fsc2 < rfsc2:
                interpolated_band = (fsc1 - rfsc1) / (rfsc2 - rfsc1 + fsc1 - fsc2)
                pass
            else:
                interpolated_band = (resolution_criterion - fsc1) / (fsc2 - fsc1) + band

        except ZeroDivisionError:
            interpolated_band = band

    else:
        interpolated_band = band

    if verbose:
        print("Band interpolated to ", interpolated_band)

    resolution = (interpolated_band + 1) / float(number_bands)

    if resolution < 0:
        resolution = 1
        interpolated_band = number_bands
        print(
            "Warning: PyTom determined a resolution < 0 for your data. "
            'Please check "mass" in data is positive or negative for all cubes.'
        )
        print(f"Warning: Setting resolution to 1 and {interpolated_band}")
        print("")

    return [resolution, interpolated_band, number_bands]


def calc_fsc_true(fsc_t, fsc_n, ring_thickness=1):
    """Calculates the true FSC as defined in Henderson
    @param fsc_t: array with FSC values without randomized phases.
    @type fsc_t: ndarray
    @param fsc_n: array with FSC values without randomized phases.
    @type fsc_n: ndarray

    @return: return the fsc_true array
    @rtype: ndarray
    """

    from numpy import zeros_like

    if len(fsc_t) != len(fsc_n):
        raise Exception("FSC arrays are not of equal length.")
    fsc_true = zeros_like(fsc_t)
    steps = 0
    for i in range(len(fsc_t)):
        if abs(fsc_t[i] - fsc_n[i]) < 1e-1:
            fsc_true[i] = fsc_t[i]
        elif steps < ring_thickness:
            fsc_true[i] = fsc_t[i]
            steps += 1
        else:
            fsc_true[i] = (fsc_t[i] - max(0, fsc_n[i])) / (1 - max(0, fsc_n[i]))

    return fsc_true


def generate_random_phases_3d(shape, reduced_complex=True):
    """This function returns a set of random phases (between -pi and pi), optionally centrosymmetric
    @shape: shape of array
    @type: tuple
    @reduced_complex: is the shape in reduced complex format
    @type: bool, default=True

    @return: returns an array with values between -pi and pi
    @rtype: ndarray
    """

    if len(shape) == 3:
        dx, dy, dz = shape
    else:
        dx, dy = shape
        dz = max(dx, dy)

    rnda = (xp.random.ranf(shape) * xp.pi * 2) - xp.pi

    cc = dx // 2
    ccc = (dx - 1) // 2

    loc = (dz // 2) * (reduced_complex == False)
    centralslice = rnda[:, :, loc]

    centralslice[cc, cc - ccc : cc] = centralslice[cc, -ccc:][::-1] * -1
    centralslice[cc - ccc : cc, cc - ccc :] = (
        xp.rot90(centralslice[-ccc:, cc - ccc :], 2) * -1
    )

    rnda[:, :, loc] = xp.fft.fftshift(centralslice) if reduced_complex else centralslice

    return rnda


def randomize_phase_beyond_freq(volume, frequency):
    """This function randomizes the phases beyond a given frequency,
    while preserving the Friedel symmetry.
    @param volume: target volume
    @type volume: L{xp.ndarray}
    @param frequency: frequency in pixel beyond which phases are randomized.
    @type frequency: int

    @return: 3d volume.
    @rtype: L{xp.ndarray}
    @author: GvdS"""

    dims = len(volume.shape)
    if dims not in {2, 3}:
        raise Exception("Invalid volume dimensions: either supply 3D or 2D ndarray")

    ft = xp.fft.rfftn(volume)
    phase = xp.angle(ft)

    amplitude = xp.abs(ft)
    rnda = generate_random_phases_3d(
        amplitude.shape, reduced_complex=True
    )  # (ranf((dx,dy,dz//2+1)) * pi * 2) - pi
    # TODO: why does 2D have the extra terms "+ dx %2" and "+ dy %2" and casting to int?
    if dims == 2:
        dx, dy = volume.shape
        x, y = xp.meshgrid(
            xp.arange(-dx // 2, dx // 2 + dx % 2), xp.arange(-dy // 2, dy // 2 + dy % 2)
        )
        rf = xp.sqrt(x**2 + y**2).astype(int)
        r = xp.fft.fftshift(rf)[:, : dy // 2 + 1]
    else:
        dx, dy, dz = volume.shape
        x, y, z = xp.meshgrid(
            xp.arange(-dx // 2, dx // 2),
            xp.arange(-dy // 2, dy // 2),
            xp.arange(-dz // 2, dz // 2),
        )
        rf = xp.sqrt(x**2 + y**2 + z**2)  # .astype(int)
        r = xp.fft.fftshift(rf)[:, :, : dz // 2 + 1]
        # centralslice = fftshift(rnda[:,:,0])
        # cc = dx//2
        # ccc= (dx-1)//2
        # centralslice[cc, cc-ccc:cc] = centralslice[cc,-ccc:][::-1]*-1
        # centralslice[cc-ccc:cc,cc-ccc:] = rot90(centralslice[-ccc:,cc-ccc:],2)*-1
        # rnda[:,:,0] = fftshift(centralslice)

    rnda[r <= frequency] = 0
    phase[r > frequency] = 0
    phase += rnda
    image = xp.fft.irfftn((amplitude * xp.exp(1j * phase)), s=volume.shape)

    if xp.abs(image.imag).sum() > 1e-8:
        raise Exception(
            "Imaginary part is non-zero. Failed to centro-summetrize the phases."
        )

    return image.real


# Sub Pixel Peak Methods


def sub_pixel_peak_parabolic(score_volume, coordinates, verbose=False):
    """
    quadratic interpolation of three adjacent samples
    @param score_volume: The score volume
    @param coordinates: (x,y,z) coordinates as tuple (!) where the sub pixel peak will be determined
    @param verbose: be talkative
    @type verbose: bool
    @return: Returns [peak_value,peak_coordinates] with sub pixel accuracy
    @rtype: float, tuple
    """

    if is_border_voxel(score_volume, coordinates):
        if verbose:
            print("sub_pixel_peak_parabolic: peak near borders - no interpolation done")
        return [score_volume[coordinates], coordinates]

    peak_coordinates = coordinates
    l = len(coordinates)
    (x, a1, b1) = qint(
        ym1=score_volume[tuple(xp.array(coordinates) - xp.array([1, 0, 0][:l]))],
        y0=score_volume[coordinates],
        yp1=score_volume[tuple(xp.array(coordinates) + xp.array([1, 0, 0][:l]))],
    )
    (y, a2, b2) = qint(
        ym1=score_volume[tuple(xp.array(coordinates) - xp.array([0, 1, 0][:l]))],
        y0=score_volume[coordinates],
        yp1=score_volume[tuple(xp.array(coordinates) + xp.array([0, 1, 0][:l]))],
    )
    if l > 2:
        (z, a3, b3) = qint(
            ym1=score_volume[tuple(xp.array(coordinates) - xp.array([0, 0, 1]))],
            y0=score_volume[coordinates],
            yp1=score_volume[tuple(xp.array(coordinates) + xp.array([0, 0, 1]))],
        )
        peak_coordinates[0] += x
        peak_coordinates[1] += y
        peak_coordinates[2] += z
        peak_value = (
            score_volume[coordinates]
            + a1 * x**2
            + b1 * x
            + a2 * y**2
            + b2 * y
            + a3 * z**2
            + b3 * z
        )
    else:
        peak_coordinates[0] += x
        peak_coordinates[1] += y
        peak_value = (
            score_volume[coordinates] + a1 * x**2 + b1 * x + a2 * y**2 + b2 * y
        )
    # return coordinates as tuple for indexing
    return [peak_value, tuple(peak_coordinates)]


def sub_pixel_peak(
    score_volume, coordinates, cube_length=8, interpolation="Spline", verbose=False
):
    """
    sub_pixel_peak: Will determine the sub pixel area of peak. Utilizes spline, fourier or
    parabolic interpolation.

    @param verbose: be talkative
    @type verbose: L{str}
    @param score_volume: The score volume
    @param coordinates: [x,y,z] coordinates where the sub pixel peak will be determined
    @param cube_length: length of cube - only used for Spline and Fourier interpolation
    @type cube_length: int (even)
    @param interpolation: interpolation type: 'Spline', 'Quadratic', or 'Fourier'
    @type interpolation: str
    @return: Returns [peak_value,peak_coordinates] with sub pixel accuracy

    change log:
    - 02/07/2013 FF: 2D functionality added
    - 17/05/2023 SR: removed broken 2D functionality
    """
    assert type(interpolation) == str, "sub_pixel_peak: interpolation must be str"
    if (interpolation.lower() == "quadratic") or (interpolation.lower() == "parabolic"):
        (peak_value, peak_coordinates) = sub_pixel_peak_parabolic(
            score_volume=score_volume, coordinates=coordinates, verbose=verbose
        )
        return [peak_value, peak_coordinates]

    from pytom.lib.pytom_volume import vol, subvolume, rescaleSpline, peak # type: ignore
    from pytom.basic.transformations import resize

    cube_start = cube_length // 2
    # This fails for 2D
    size_x, size_y, size_z = score_volume.shape

    if any(
        (
            coordinates[0] - cube_start < 1,
            coordinates[1] - cube_start < 1,
            coordinates[2] - cube_start < 1,
            coordinates[0] - cube_start + cube_length >= size_x,
            coordinates[1] - cube_start + cube_length >= size_y,
            coordinates[2] - cube_start + cube_length >= size_z,
        )
    ):
        if verbose:
            print("SubPixelPeak: position too close to border for sub-pixel")
        return [
            score_volume(coordinates[0], coordinates[1], coordinates[2]),
            coordinates,
        ]

    sub_volume = subvolume(
        score_volume,
        coordinates[0] - cube_start,
        coordinates[1] - cube_start,
        coordinates[2] - cube_start,
        cube_length,
        cube_length,
        cube_length,
    )

    # size of interpolated volume
    scale_size = 10 * cube_length

    # ratio between interpolation area and large volume
    scale_ratio = 1.0 * cube_length / scale_size

    # resize into bigger volume
    if interpolation == "Spline":
        sub_volume_scaled = vol(scale_size, scale_size, scale_size)
        rescaleSpline(sub_volume, sub_volume_scaled)
    else:
        sub_volume_scaled = resize(volume=sub_volume, factor=10)[0]

    peak_coordinates = peak(sub_volume_scaled)

    peak_value = sub_volume_scaled(
        peak_coordinates[0], peak_coordinates[1], peak_coordinates[2]
    )

    # calculate sub pixel coordinates of interpolated peak
    peak_coordinates[0] = (
        peak_coordinates[0] * scale_ratio - cube_start + coordinates[0]
    )
    peak_coordinates[1] = (
        peak_coordinates[1] * scale_ratio - cube_start + coordinates[1]
    )
    peak_coordinates[2] = (
        peak_coordinates[2] * scale_ratio - cube_start + coordinates[2]
    )
    if any(
        (
            peak_coordinates[0] > score_volume.sizeX(),
            peak_coordinates[1] > score_volume.sizeY(),
            peak_coordinates[2] > score_volume.sizeZ(),
        )
    ):
        if verbose:
            print("SubPixelPeak: peak position too large :( return input value")
        # something went awfully wrong here. return regular value
        return [
            score_volume(coordinates[0], coordinates[1], coordinates[2]),
            coordinates,
        ]
    return [peak_value, peak_coordinates]


def sub_pixel_max_3d(
    volume,
    k=0.01,
    ignore_border=50,
    interpolation="filt_bspline",
    plan=None,
    profile=True,
    num_threads=1024,
    zoomed=None,
    fast_sum=None,
    max_id=None,
):
    """
    Function to find the highest point in a 3D array, with subpixel accuracy using cubic spline
    interpolation.

    @param inp: A 3D numpy/cupy array containing the data points.
    @type inp: numpy/cupy array 3D
    @param k: The interpolation factor used in the spline interpolation,
              k < 1 is zoomed in, k>1 zoom out.
    @type k: float
    @return: A list of maximal value in the interpolated volume and a list of x position,
             the y position and the value.
    @returntype: list
    """

    from pytom.voltools import transform

    ox, oy, oz = volume.shape
    ib = ignore_border
    cropped_volume = volume[ib : ox - ib, ib : oy - ib, ib : oz - ib].astype(xp.float32)

    if profile:
        stream = xp.cuda.Stream.null
        t_start = stream.record()

    # x,y,z = xp.array(max_index(cropped_volume)) + ignore_border
    x, y, z = (
        xp.array(xp.unravel_index(cropped_volume.argmax(), cropped_volume.shape))
        + ignore_border
    )

    dx, dy, dz = volume.shape
    translation = [dx // 2 - x, dy // 2 - y, dz // 2 - z]

    if profile:
        t_end = stream.record()
        t_end.synchronize()

        time_took = xp.cuda.get_elapsed_time(t_start, t_end)
        print(f"initial find max time: \t{time_took:.3f}ms")
        t_start = stream.record()

    b = max(0, int(volume.shape[0] // 2 - 4 / k))
    zx, zy, zz = volume.shape
    out = volume[b : zx - b, b : zy - b, b : zz - b]

    transform(
        out,
        output=zoomed,
        scale=(k, k, k),
        device=device,
        translation=translation,
        interpolation=interpolation,
    )

    if profile:
        t_end = stream.record()
        t_end.synchronize()

        time_took = xp.cuda.get_elapsed_time(t_start, t_end)
        print(f"transform finished in \t{time_took:.3f}ms")
        t_start = stream.record()

    nblocks = int(xp.ceil(zoomed.size / num_threads / 2))
    argmax(
        (
            nblocks,
            1,
        ),
        (num_threads, 1, 1),
        (zoomed, fast_sum, max_id, zoomed.size),
        shared_mem=8 * num_threads,
    )
    x2, y2, z2 = xp.unravel_index(max_id[fast_sum.argmax()], zoomed.shape)

    peak_value = zoomed[x2][y2][z2]
    peak_shift = [
        x + (x2 - zoomed.shape[0] // 2) * k - volume.shape[0] // 2,
        y + (y2 - zoomed.shape[1] // 2) * k - volume.shape[1] // 2,
        z + (z2 - zoomed.shape[2] // 2) * k - volume.shape[2] // 2,
    ]

    if profile:
        t_end = stream.record()
        t_end.synchronize()

        time_took = xp.cuda.get_elapsed_time(t_start, t_end)
        print(f"argmax finished in \t{time_took:.3f}ms")
        t_start = stream.record()
        print()

    return [peak_value, peak_shift]


def max_index(volume, num_threads=1024):
    nblocks = int(xp.ceil(volume.size / num_threads / 2))
    fast_sum = -1000000 * xp.ones((nblocks), dtype=xp.float32)
    max_id = xp.zeros((nblocks), dtype=xp.int32)
    argmax(
        (
            nblocks,
            1,
        ),
        (num_threads, 1, 1),
        (volume, fast_sum, max_id, volume.size),
        shared_mem=16 * num_threads,
    )
    mm = min(max_id[fast_sum.argmax()], volume.size - 1)
    indices = xp.unravel_index(mm, volume.shape)
    return indices


def is_border_voxel(volume, coordinates):
    shape = volume.shape
    for n, pos in enumerate(coordinates):
        if pos < 1 or pos >= shape[n] - 1:
            return True
    return False


def qint(ym1, y0, yp1):
    """
    quadratic interpolation of three adjacent samples

    [p,y,a] = qint(ym1,y0,yp1)

    returns the extremum location p, height y, and half-curvature a
    of a parabolic fit through three points.
    Parabola is given by y(x) = a*(x-p)^2+b,
    where y(-1)=ym1, y(0)=y0, y(1)=yp1.
    @param ym1: y(-1)
    @type ym1: float
    @param y0: y(0)
    @type y0: float
    @param yp1: y(+1)
    @type yp1: float
    @return: peak-ordinate, peak-value
    """
    p = (yp1 - ym1) / (2 * (2 * y0 - yp1 - ym1))
    y = y0 - 0.25 * (ym1 - yp1) * p
    a = 0.5 * (ym1 - 2 * y0 + yp1)
    # b = y0 - a*(p**2)
    return p, y, a

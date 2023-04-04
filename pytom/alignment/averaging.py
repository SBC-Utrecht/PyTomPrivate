from pytom.gpu.initialize import xp
import numpy as np


analytWedge = False


def splitParticleList(particleList, setParticleNodesRatio=3, numberOfNodes=10):
    """
    @param particleList: The particle list
    @param setParticleNodesRatio: minimum number of particles per node
    @type setParticleNodesRatio: L{int}
    @return: list of particle lists, splitFactor (number of processors or smaller for few particles)
    @rtype: list, L{int}
    @author: FF
    """

    particleNodesRatio = float(len(particleList)) / float(numberOfNodes)
    splitFactor = numberOfNodes
    # make sure each node gets at least setParticleNodesRatio particles.
    if particleNodesRatio < setParticleNodesRatio:
        splitFactor = len(particleList) / int(setParticleNodesRatio)
    splitLists = particleList.splitNSublists(splitFactor)  # somehow ...
    return splitLists


def average(particleList, averageName, showProgressBar=False, verbose=False,
            createInfoVolumes=False, weighting=False, norm=False, previousReference=None, wiener_filter=False,
            gpuId=None):
    """
    average : Creates new average from a particleList
    @param particleList: The particles
    @param averageName: Filename of new average
    @param verbose: Prints particle information. Disabled by default.
    @param createInfoVolumes: Create info data (wedge sum, inverted density) too? False by default.
    @param weighting: apply weighting to each average according to its correlation score
    @param norm: apply normalization for each particle
    @return: A new Reference object
    @rtype: L{pytom.basic.structures.Reference}
    @author: Thomas Hrabe
    @change: limit for wedgeSum set to 1% or particles to avoid division by small numbers - FF
    """
    from pytom.lib.pytom_volume import read, vol, reducedToFull, fullToReduced, complexRealMult, real, abs
    from pytom.basic.filter import lowpassFilter, rotateWeighting, profile2FourierVol
    from pytom.lib.pytom_volume import transformSpline as transform
    from pytom.basic.fourier import convolute, fft, powerspectrum, radialaverage, iftshift
    from pytom.basic.structures import Reference
    from pytom.basic.normalise import mean0std1
    from pytom.tools.ProgressBar import FixedProgBar
    from pytom.alignment.alignmentFunctions import invert_WedgeSum
    from math import exp
    import os
    if len(particleList) == 0:
        raise RuntimeError('The particle list is empty. Aborting!')

    # initialize to None to ensure variable exists
    numberAlignedParticles, progressBar = None, None
    if showProgressBar:
        progressBar = FixedProgBar(0, len(particleList), 'Particles averaged ')
        progressBar.update(0)
        numberAlignedParticles = 0

    # pre-check that scores != 0
    if weighting:
        wsum = 0.
        for particleObject in particleList:
            wsum += particleObject.getScore().getValue()
        if wsum < 0.00001:
            weighting = False
            print("Warning: all scores have been zero - weighting not applied")

    # read box dims and interpolation center from the first particle
    particleObject = particleList[0]
    wedgeInfo = particleObject.getWedge()
    particle = read(particleObject.getFilename())
    sizeX, sizeY, sizeZ = particle.sizeX(), particle.sizeY(), particle.sizeZ()
    centerX, centerY, centerZ = sizeX // 2, sizeY // 2, sizeZ // 2

    # allocate boxes for averaging
    newParticle = vol(sizeX, sizeY, sizeZ)
    result = vol(sizeX, sizeY, sizeZ)
    result.setAll(0.0)
    wedgeSum = wedgeInfo.returnWedgeVolume(sizeX, sizeY, sizeZ)
    wedgeSum.setAll(0)

    # bug check: size dim should be sizeZ / 2 + 1 as pytom works with fourier reduced volumes
    assert wedgeSum.sizeX() == sizeX and wedgeSum.sizeY() == sizeY and wedgeSum.sizeZ() == sizeZ / 2 + 1, \
        "wedge initialization result in wrong dims :("

    # initialize to None to ensure variable exists
    noise_spectrum, previous_average, previous_average_rotated, ssnr_filter = None, None, None, None
    if wiener_filter:
        noise_spectrum = vol(wedgeSum.sizeX(), wedgeSum.sizeY(), wedgeSum.sizeZ())
        noise_spectrum.setAll(0.0)
        previous_average = previousReference.getVolume()
        previous_average_rotated = vol(sizeX, sizeY, sizeZ)

    for particleObject in particleList:
        if 0 and verbose:
            print(particleObject)

        if not os.path.exists(particleObject.getFilename()):
            continue
        particle = read(particleObject.getFilename())
        if norm:  # normalize the particle
            mean0std1(particle)  # happen inplace

        # get data about particle
        wedgeInfo = particleObject.getWedge()
        shiftV = particleObject.getShift()
        rotation = particleObject.getRotation()
        rotinvert = rotation.invert()

        # load wedges and for 3d ctf calculate noise spectrum
        if wiener_filter:
            wedge = wedgeInfo.returnWedgeVolume(sizeX, sizeY, sizeZ, False)

            # =====> first calculate the noise spectrum with unrotated wedge
            previous_average_rotated.setAll(0.0)
            transform(previous_average, previous_average_rotated, rotation[0], rotation[1], rotation[2],
                      centerX, centerY, centerZ, shiftV[0], shiftV[1], shiftV[2], 0, 0, 0)

            p_noise = abs(
                fft(particle, scaling='sqrtN') - complexRealMult(fft(previous_average_rotated,
                                                                     scaling='sqrtN'), wedge))
            # power spectrum is abs square of fourier transform
            noise_spectrum = noise_spectrum + real(p_noise * p_noise)

            # =====> calculate the CTF^2 for the rotated particle
            wedge = rotateWeighting(weighting=wedge,
                                    z1=rotinvert[0], z2=rotinvert[1], x=rotinvert[2], mask=None,
                                    isReducedComplex=True, returnReducedComplex=True)
            wedge = wedge * wedge
        elif analytWedge:
            # TODO McHaillet: analytWedge should first rotate the fourier grid and then calculate the wedge
            # TODO instead of calculate the wedge and then interpolating, doubt if this is ever getting a fix
            # > analytical buggy version
            wedge = wedgeInfo.returnWedgeVolume(sizeX, sizeY, sizeZ, False, rotinvert)
        else:
            # > FF: interpol bugfix
            wedge = rotateWeighting(weighting=wedgeInfo.returnWedgeVolume(sizeX, sizeY, sizeZ, False),
                                    z1=rotinvert[0], z2=rotinvert[1], x=rotinvert[2], mask=None,
                                    isReducedComplex=True, returnReducedComplex=True)
            # wedge = wedgeInfo.returnWedgeVolume(sizeX, sizeY, sizeZ, False, rotation=rotinvert)

            # < FF
            # > TH bugfix
            # wedgeVolume = wedgeInfo.returnWedgeVolume(wedgeSizeX=sizeX, wedgeSizeY=sizeY, wedgeSizeZ=sizeZ,
            #                                    humanUnderstandable=True, rotation=rotinvert)
            # wedge = rotate(volume=wedgeVolume, rotation=rotinvert, imethod='linear')
            # < TH

        # convolute the particle with the unrotated wedge
        particle = wedgeInfo.apply(particle)

        # shift and rotate particle
        newParticle.setAll(0)
        # pass rotinvert here, why otherwise calculate rotinvert before...
        # previously: -rotation[1], -rotation[0], -rotation[2]
        transform(particle, newParticle, rotinvert[0], rotinvert[1], rotinvert[2],
                  centerX, centerY, centerZ, -shiftV[0], -shiftV[1], -shiftV[2], 0, 0, 0)

        if weighting:
            weight = 1. - particleObject.getScore().getValue()
            weight = exp(-1. * weight)
            result = result + newParticle * weight
            wedgeSum = wedgeSum + wedge * weight
        else:
            result = result + newParticle
            wedgeSum = wedgeSum + wedge

        if showProgressBar:
            numberAlignedParticles = numberAlignedParticles + 1
            progressBar.update(numberAlignedParticles)

    # get file name format
    root, ext = os.path.splitext(averageName)

    # calculate the signal power spectrum from the previous average
    if wiener_filter:
        noise_spectrum = noise_spectrum / len(particleList)

        # write estimated signal and noise 1d vectors to disk
        signal_vector = radialaverage(powerspectrum(previous_average), isreduced=False)
        noise_vector = radialaverage(noise_spectrum, isreduced=True)

        # write as a 1d array to text file is most efficient
        np.savetxt(f'{root}-signal-spectrum-1d.txt', signal_vector)
        np.savetxt(f'{root}-noise-spectrum-1d.txt', noise_vector)

        # calculate the inverse ssnr for the wiener filter
        ssnr_vector = [s / n for s, n in zip(signal_vector, noise_vector)]
        ssnr_filter = fullToReduced(iftshift(profile2FourierVol(ssnr_vector, dim=sizeX, reduced=False)))
        invert_WedgeSum(invol=ssnr_filter, r_max=sizeX / 2 - 2, lowlimit=1e-6,
                        lowval=1e-6)

    # remove unknown fourier areas
    result = lowpassFilter(result, sizeX / 2 - 1, 0.)[0]

    # if createInfoVolumes:
    result.write(f'{root}-PreWedge{ext}')

    # wedgeSum = wedgeSum*0+len(particleList)
    wedgeSum.write(f'{root}-WedgeSumUnscaled{ext}')

    if wiener_filter:
        wedgeSum = wedgeSum + ssnr_filter
        invert_WedgeSum(invol=wedgeSum, r_max=sizeX / 2 - 2, lowlimit=1e-6,
                        lowval=1e-6)
        wedgeSum.write(f'{root}-wiener-filter{ext}')
    else:
        invert_WedgeSum(invol=wedgeSum, r_max=sizeX / 2 - 2, lowlimit=.05 * len(particleList),
                        lowval=.05 * len(particleList))

    if createInfoVolumes:
        w1 = reducedToFull(wedgeSum)
        w1.write(f'{root}-WedgeSumInverted{ext}')

    result = convolute(v=result, k=wedgeSum, kernel_in_fourier=True)
    result.write(averageName)

    if createInfoVolumes:
        resultINV = result * -1
        # write sign inverted result to disk (good for chimera viewing ... )
        resultINV.write(f'{root}-INV{ext}')

    newReference = Reference(averageName, particleList)

    return newReference


def average_wiener(previous_average=None):
    """
    If no previous average is provided will use Grigorieff method to estimate the signal.
    :param previous_average:
    :return:
    """
    pass


def averageGPU(particleList, averageName, showProgressBar=False, verbose=False,
               createInfoVolumes=False, weighting=False, norm=False, previousReference=None,
               wiener_filter=False, gpuId=None, profile=False):
    """
    average : Creates new average from a particleList
    @param particleList: The particles
    @param averageName: Filename of new average
    @param verbose: Prints particle information. Disabled by default.
    @param createInfoVolumes: Create info data (wedge sum, inverted density) too? False by default.
    @param weighting: apply weighting to each average according to its correlation score
    @param norm: apply normalization for each particle
    @return: A new Reference object
    @rtype: L{pytom.basic.structures.Reference}
    @author: Thomas Hrabe
    @change: limit for wedgeSum set to 1% or particles to avoid division by small numbers - FF
    """
    from pytom.agnostic.io import read, write, read_size
    from pytom.agnostic.filter import bandpass as lowpassFilter, applyFourierFilterFull, profile2FourierVol
    from pytom.voltools import transform
    from pytom.basic.structures import Reference
    from pytom.agnostic.normalise import mean0std1
    from pytom.agnostic.tools import invert_WedgeSum, create_sphere
    from pytom.agnostic.transform import fourier_full2reduced
    from pytom.simulation.microscope import radial_average
    from pytom.tools.ProgressBar import FixedProgBar
    import cupy as xp
    import os

    try:
        from cupyx.scipy.fftpack.fft import fftn as fftnP
        from cupyx.scipy.fftpack.fft import ifftn as ifftnP
        from cupyx.scipy.fftpack.fft import get_fft_plan
    except ModuleNotFoundError:  # TODO go to 'from cupyx.scipy.fftpack import ...' once fully moved to cupy > 8.3
        from cupyx.scipy.fftpack import fftn as fftnP
        from cupyx.scipy.fftpack import ifftn as ifftnP
        from cupyx.scipy.fftpack import get_fft_plan

    if not gpuId is None:
        device = f'gpu:{gpuId}'
        xp.cuda.Device(gpuId).use()
    else:
        print(gpuId)
        raise Exception('Running gpu code on non-gpu device')

    cstream = xp.cuda.Stream()
    if profile:
        stream = xp.cuda.Stream.null
        t_start = stream.record()

    if len(particleList) == 0:
        raise RuntimeError('The particle list is empty. Aborting!')

    # pre-check that scores != 0
    if weighting:
        wsum = 0.
        for particleObject in particleList:
            wsum += particleObject.getScore().getValue()
        if wsum < 0.00001:
            weighting = False
            print("Warning: all scores have been zero - weighting not applied")

    # extract root path and extension for writing files
    root, ext = os.path.splitext(averageName)

    shape = tuple(read_size(particleList[0].getFilename()))
    center = shape[0] // 2, shape[1] // 2, shape[2] // 2

    # TODO ifftshift to shift centered spectrum back to corner!
    wedge = xp.zeros(shape, dtype=xp.float32)
    wedgeSum = xp.zeros(shape, dtype=xp.float32)
    newParticle = xp.zeros(shape, dtype=xp.float32)
    result = xp.zeros(shape, dtype=xp.float32)
    fftplan = get_fft_plan(wedge.astype(xp.complex64))

    ortho_scaling, low_pass, noise_variance_sum, signal_variance_sum, signal_raw, signal = (None, None, None,
                                                                                            None, None, None)
    if wiener_filter:
        # TODO mask is necessary for wiener filter spectra estimation
        # prev_average_rotated = xp.zeros(shape)
        ortho_scaling = xp.sqrt(result.size, dtype=xp.float32)
        low_pass = xp.fft.ifftshift(create_sphere(shape, radius=shape[0] // 2 - 1))
        noise_variance_sum = xp.zeros(shape)
        signal_variance_sum = xp.zeros(shape)
        signal_raw = xp.zeros(shape, dtype=xp.complex64)
        signal = fftnP(read(previousReference.getFilename()), plan=fftplan) * (1 / ortho_scaling)

        # create signal spectrum
        # prev_average = read(previousReference.getFilename())
        # signal = xp.abs(xp.fft.fftshift(xp.fft.fftn(prev_average, norm='ortho'))) ** 2 / 2
        # _, r_signal = radial_average(signal)
        # signal_spectrum = profile2FourierVol(r_signal * t_fudge, dim=shape) * low_pass

    # for profiling
    total = len(particleList)
    if profile:
        t_end = stream.record()
        t_end.synchronize()

        time_took = xp.cuda.get_elapsed_time(t_start, t_end)
        print(f'startup time: \t{time_took:.3f}ms')
        t_start = stream.record()

    if showProgressBar:
        progressBar = FixedProgBar(0, len(particleList), 'Particles averaged ')
        progressBar.update(0)

    for n, particleObject in enumerate(particleList):
        # get the transformation parameters from the particle object
        rotation = particleObject.getRotation()
        rotinvert = rotation.invert()
        zxz_backward = (rotinvert[0], rotinvert[2], rotinvert[1])
        shift_forward = tuple(particleObject.getShift().toVector())
        shift_backward = (-shift_forward[0], -shift_forward[1], -shift_forward[2])

        particle = read(particleObject.getFilename(), deviceID=device)

        if norm:  # normalize the particle
            mean0std1(particle)  # happen inplace

        # get the wedge per particle because the wedge can differ
        wedgeInfo = particleObject.getWedge().convert2numpy()

        if wiener_filter:
            # align the particle and wedge to the reference orientation
            transform(particle, rotation=zxz_backward, rotation_order='rzxz', translation=shift_backward,
                      center=center, interpolation='filt_bspline', output=newParticle, device=device)  # * mask
            particle_ft = fftnP(newParticle, plan=fftplan) * (1 / ortho_scaling)
            ctf = xp.fft.ifftshift(
                transform(wedgeInfo.returnWedgeVolume(humanUnderstandable=True), rotation=zxz_backward,
                          rotation_order='rzxz', center=center, device=device))
            ctf_squared = ctf ** 2

            # calculate signal, ctf sum, noise and signal variance
            signal_raw += (ctf * particle_ft)
            wedgeSum += ctf_squared
            signal_variance_sum += ctf_squared * xp.abs(particle_ft) ** 2
            noise_variance_sum += xp.abs(particle_ft - ctf * signal) ** 2

        else:
            wedgeZero = xp.fft.ifftshift(wedgeInfo.returnWedgeVolume(*shape, True))

            # apply wedge to particle
            particle = (ifftnP(fftnP(particle, plan=fftplan) * wedgeZero, plan=fftplan)).real

            # create spectral wedge weighting
            wedge *= 0
            transform(xp.fft.fftshift(wedgeZero), rotation=zxz_backward,
                      rotation_order='rzxz', center=center, output=wedge, device=device,
                      interpolation='linear')

            # shift and rotate particle
            newParticle *= 0
            transform(particle, output=newParticle, rotation=zxz_backward,
                      center=center, translation=shift_backward,
                      device=device, interpolation='filt_bspline', rotation_order='rzxz')

            # add to average and wedgeweighting
            result += newParticle
            wedgeSum += xp.fft.ifftshift(wedge)

        # for profiling
        if n % total == 0:
            if profile:
                t_end = stream.record()
                t_end.synchronize()

                time_took = xp.cuda.get_elapsed_time(t_start, t_end)
                print(f'total time {n:5d}: \t{time_took:.3f}ms')
                t_start = stream.record()

        if showProgressBar:
            progressBar.update(n)

        cstream.synchronize()

    if wiener_filter:
        # write out for combining parallel procs
        write(f'{root}-PreWedge{ext}', ifftnP(signal_raw * ortho_scaling, plan=fftplan).real)
        write(f'{root}-WedgeSumUnscaled{ext}', wedgeSum)
        write(f'{root}-signal-spectrum{ext}', signal_variance_sum)
        write(f'{root}-noise-spectrum{ext}', noise_variance_sum)

        # calculate radial profiles
        ctf_radial = radial_average(xp.fft.fftshift(wedgeSum))[1]
        signal_variance = radial_average(xp.fft.fftshift(signal_variance_sum))[1] / ctf_radial
        noise_variance = (radial_average(xp.fft.fftshift(noise_variance_sum))[1] /
                          (len(particleList) - 1))
        ssnr = (signal_variance / noise_variance) - (1 / ctf_radial)

        # determine wiener filter and filter result
        wiener_filter = xp.ones(shape)
        denom = (wedgeSum + profile2FourierVol(1 / ssnr, dim=shape) * low_pass)
        wiener_filter[denom != 0] = (wiener_filter[denom != 0] / denom[denom != 0])
        wiener_filter *= low_pass
        result = ifftnP((signal_raw * wiener_filter * ortho_scaling).astype(np.complex64), plan=fftplan).real
        # * mask

        write(f'{root}-wiener-filter{ext}', wiener_filter)
        write(averageName, result)

    else:
        # apply spectral weighting to sum
        result = lowpassFilter(result, high=shape[0] / 2 - 1, sigma=0)
        # write the unfiltered average, needed for parallel pooling of results
        write(f'{root}-PreWedge{ext}', result)

        write(f'{root}-WedgeSumUnscaled{ext}', fourier_full2reduced(wedgeSum))
        # prev: wedgeSumINV =
        invert_WedgeSum(wedgeSum, r_max=shape[0] // 2 - 2., lowlimit=.05 * len(particleList), lowval=.05 * len(
            particleList))

        if createInfoVolumes:
            # write(f'{root}-WedgeSumInverted{ext}', xp.fft.fftshift(wedgeSumINV))
            write(f'{root}-WedgeSumInverted{ext}', xp.fft.fftshift(wedgeSum))

        # the wedge sum should already be centered in the corner ?
        result = applyFourierFilterFull(result, wedgeSum)  # prev wedgeSumINV

        # do a low pass filter
        result = lowpassFilter(result, shape[0] / 2 - 2, (shape[0] / 2 - 1) / 10.)[0]
        write(averageName, result)

        if createInfoVolumes:
            resultINV = result * -1
            # write sign inverted result to disk (good for chimera viewing ... )
            write(f'{root}-INV{ext}', resultINV)

    newReference = Reference(averageName, particleList)

    return newReference


def average_wiener_gpu(previous_average=None, gpu_id=0):
    """
    If no previous average is provided will use Grigorieff method to estimate the signal.
    :param previous_average:
    :return:
    """
    pass


def average_wiener_parallel(mpi_class=None):
    pass


def average_parallel(mpi_class=None):
    pass


def averageParallel(particleList, averageName, showProgressBar=False, verbose=False,
                    createInfoVolumes=False, weighting=None, norm=False,
                    setParticleNodesRatio=3, cores=6, gpuID=None):
    """
    compute average using parfor
    @param particleList: The particles
    @param averageName: Filename of new average
    @param verbose: Prints particle information. Disabled by default.
    @param createInfoVolumes: Create info data (wedge sum, inverted density) too? False by default.
    @param weighting: weight particles by exp CC in average
    @type weighting: bool
    @param setParticleNodesRatio: minimum number of particles per node
    @type setParticleNodesRatio: L{int}
    @return: A new Reference object
    @rtype: L{pytom.basic.structures.Reference}
    @author: FF

    """
    from pytom.lib.pytom_volume import read, complexRealMult
    from pytom.basic.fourier import fft, ifft
    from pytom.basic.filter import lowpassFilter
    from pytom.basic.structures import Reference
    from pytom.alignment.alignmentFunctions import invert_WedgeSum
    from pytom.basic.files import em2mrc
    from multiprocessing import Process
    import os
    import time

    splitLists = splitParticleList(particleList, setParticleNodesRatio=setParticleNodesRatio, numberOfNodes=cores)
    splitFactor = len(splitLists)

    avgNameList = []
    preList = []
    wedgeList = []
    for ii in range(splitFactor):
        avgName = averageName + '_dist' + str(ii) + '.em'
        avgNameList.append(avgName)
        preList.append(averageName + '_dist' + str(ii) + '-PreWedge.em')
        wedgeList.append(averageName + '_dist' + str(ii) + '-WedgeSumUnscaled.em')

    procs = []
    for i in range(splitFactor):
        proc = Process(target=average, args=(
        splitLists[i], avgNameList[i], showProgressBar, verbose, createInfoVolumes, weighting, norm))
        procs.append(proc)
        proc.start()

    while procs:
        procs = [proc for proc in procs if proc.is_alive()]
        time.sleep(.1)

    # averageList = mpi.parfor( average, list(zip(splitLists, avgNameList, [showProgressBar]*splitFactor,
    #                                       [verbose]*splitFactor, [createInfoVolumes]*splitFactor,
    #                                            [weighting]*splitFactor, [norm]*splitFactor)), verbose=True)

    # collect results from files
    unweiAv = read(preList[0])
    wedgeSum = read(wedgeList[0])
    os.system('rm ' + wedgeList[0])
    os.system('rm ' + avgNameList[0])
    os.system('rm ' + preList[0])
    for ii in range(1, splitFactor):
        av = read(preList[ii])
        unweiAv += av
        os.system('rm ' + preList[ii])
        w = read(wedgeList[ii])
        wedgeSum += w
        os.system('rm ' + wedgeList[ii])
        os.system('rm ' + avgNameList[ii])

    if createInfoVolumes:
        root, ext = os.path.splitext(averageName)

        unweiAv.write(f'{root}-PreWedge{ext}')
        wedgeSum.write(f'{root}-WedgeSumUnscaled{ext}')

    # convolute unweighted average with inverse of wedge sum
    invert_WedgeSum(invol=wedgeSum, r_max=unweiAv.sizeX() / 2 - 2., lowlimit=.05 * len(particleList),
                    lowval=.05 * len(particleList))
    invert_WedgeSum(invol=wedgeSum, r_max=unweiAv.sizeX() / 2 - 2., lowlimit=.05 * len(particleList),
                    lowval=.05 * len(particleList))

    if createInfoVolumes:
        wedgeSum.write(averageName[:len(averageName) - 3] + '-WedgeSumINV.em')

    fResult = fft(unweiAv)
    r = complexRealMult(fResult, wedgeSum)
    unweiAv = ifft(r)
    unweiAv.shiftscale(0.0, 1 / float(unweiAv.sizeX() * unweiAv.sizeY() * unweiAv.sizeZ()))
    # low pass filter to remove artifacts at fringes
    unweiAv = lowpassFilter(volume=unweiAv, band=unweiAv.sizeX() / 2 - 2, smooth=(unweiAv.sizeX() / 2 - 1) / 10.)[0]

    if averageName.endswith("mrc"):
        averageNameEM = averageName[:-3] + 'em'
        unweiAv.write(averageNameEM)
        em2mrc(averageNameEM, './' if not os.path.dirname(averageName) else os.path.dirname(averageName))
        os.remove(averageNameEM)

    else:
        unweiAv.write(averageName)

    return Reference(averageName, particleList)


def averageParallelGPU(particleList, averageName, showProgressBar=False, verbose=False,
                       createInfoVolumes=False, weighting=None, norm=False,
                       setParticleNodesRatio=3, cores=6, gpuID=None):
    """
    compute average using parfor
    @param particleList: The particles
    @param averageName: Filename of new average
    @param verbose: Prints particle information. Disabled by default.
    @param createInfoVolumes: Create info data (wedge sum, inverted density) too? False by default.
    @param weighting: weight particles by exp CC in average
    @type weighting: bool
    @param setParticleNodesRatio: minimum number of particles per node
    @type setParticleNodesRatio: L{int}
    @return: A new Reference object
    @rtype: L{pytom.basic.structures.Reference}
    @author: FF

    """
    from pytom.agnostic.tools import invert_WedgeSum
    from pytom.agnostic.io import write, read
    import os

    splitLists = splitParticleList(particleList, setParticleNodesRatio=setParticleNodesRatio, numberOfNodes=cores)
    splitFactor = len(splitLists)

    avgNameList = []
    preList = []
    wedgeList = []
    for ii in range(splitFactor):
        avgName = averageName + '_dist' + str(ii) + '.em'
        avgNameList.append(avgName)
        preList.append(averageName + '_dist' + str(ii) + '-PreWedge.em')
        wedgeList.append(averageName + '_dist' + str(ii) + '-WedgeSumUnscaled.em')

    #####
    averageGPU(splitLists[0], avgNameList[0], showProgressBar, verbose, createInfoVolumes, weighting, norm, gpuID)
    # averageList = mpi.parfor( average, list(zip(splitLists, avgNameList, [showProgressBar]*splitFactor,
    #                                       [verbose]*splitFactor, [createInfoVolumes]*splitFactor,
    #                                            [weighting]*splitFactor, [norm]*splitFactor)), verbose=True)

    unweiAv = read(preList[0])
    wedgeSum = read(wedgeList[0])
    os.system('rm ' + wedgeList[0])
    os.system('rm ' + avgNameList[0])
    os.system('rm ' + preList[0])
    for ii in range(1, splitFactor):
        # print(preList[ii], wedgeList[ii], avgNameList[ii])
        av = read(preList[ii])
        unweiAv += av
        os.system('rm ' + preList[ii])
        w = read(wedgeList[ii])
        wedgeSum += w
        # os.system('rm ' + wedgeList[ii])
        # os.system('rm ' + avgNameList[ii])

    if createInfoVolumes:
        root, ext = os.path.splitext(averageName)

        write(f'{root}-PreWedge{ext}', unweiAv)
        write(f'{root}-WedgeSumUnscaled{ext}', wedgeSum)

    # convolute unweighted average with inverse of wedge sum
    wedgeINV = invert_WedgeSum((wedgeSum), r_max=unweiAv.shape[0] / 2 - 2., lowlimit=.05 * len(particleList),
                               lowval=.05 * len(particleList))

    if createInfoVolumes:
        write(averageName[:len(averageName) - 3] + '-WedgeSumINV.em', wedgeINV)

    r = xp.fft.rfftn(unweiAv) * wedgeINV
    unweiAv = (xp.fft.irfftn(r)).real
    # unweiAv.shiftscale(0.0,1/float(unweiAv.sizeX()*unweiAv.sizeY()*unweiAv.sizeZ()))
    # low pass filter to remove artifacts at fringes
    # unweiAv = lowpassFilter(volume=unweiAv, band=unweiAv.sizeX()/2-2, smooth=(unweiAv.sizeX()/2-1)/10.)[0]

    write(averageName, unweiAv)
    return 1
from pytom.gpu.initialize import xp


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
            createInfoVolumes=False, weighting=False, norm=False, gpuId=None):
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
    from pytom.lib.pytom_volume import read, vol, reducedToFull
    from pytom.basic.filter import lowpassFilter, rotateWeighting
    from pytom.lib.pytom_volume import transformSpline as transform
    from pytom.basic.fourier import convolute
    from pytom.basic.structures import Reference
    from pytom.basic.normalise import mean0std1
    from pytom.tools.ProgressBar import FixedProgBar
    from pytom.alignment.alignmentFunctions import invert_WedgeSum
    from math import exp
    import os
    if len(particleList) == 0:
        raise RuntimeError('The particle list is empty. Aborting!')

    if showProgressBar:
        progressBar = FixedProgBar(0, len(particleList), 'Particles averaged ')
        progressBar.update(0)
        numberAlignedParticles = 0

    result = []
    wedgeSum = []

    newParticle = None
    # pre-check that scores != 0
    if weighting:
        wsum = 0.
        for particleObject in particleList:
            wsum += particleObject.getScore().getValue()
        if wsum < 0.00001:
            weighting = False
            print("Warning: all scores have been zero - weighting not applied")

    n = 0

    for particleObject in particleList:
        if 0 and verbose:
            print(particleObject)

        if not os.path.exists(particleObject.getFilename()):
            continue
        particle = read(particleObject.getFilename())
        if norm:  # normalize the particle
            mean0std1(particle)  # happen inplace

        wedgeInfo = particleObject.getWedge()
        # apply its wedge to itself

        rotation = particleObject.getRotation()
        rotinvert = rotation.invert()

        if not result:
            sizeX = particle.sizeX()
            sizeY = particle.sizeY()
            sizeZ = particle.sizeZ()
            newParticle = vol(sizeX, sizeY, sizeZ)

            centerX = sizeX // 2
            centerY = sizeY // 2
            centerZ = sizeZ // 2

            result = vol(sizeX, sizeY, sizeZ)
            result.setAll(0.0)

            if analytWedge:
                wedgeSum = wedgeInfo.returnWedgeVolume(wedgeSizeX=sizeX, wedgeSizeY=sizeY, wedgeSizeZ=sizeZ)
            else:
                # > FF bugfix
                wedgeSum = wedgeInfo.returnWedgeVolume(sizeX, sizeY, sizeZ)
                # < FF
                # > TH bugfix
                # wedgeSum = vol(sizeX,sizeY,sizeZ)
                # < TH
                # wedgeSum.setAll(0)
            assert wedgeSum.sizeX() == sizeX and wedgeSum.sizeY() == sizeY and wedgeSum.sizeZ() == sizeZ / 2 + 1, \
                "wedge initialization result in wrong dims :("
            wedgeSum.setAll(0)
            # wedgeFilter = wedgeInfo.returnWedgeFilter(particle.sizeX(), particle.sizeY(), particle.sizeZ())

        if wedgeInfo._type in ['SingleTiltWedge', 'DoubleTiltWedge']:
            particle = wedgeInfo.apply(particle)  # dont for 3d ctf, because particle is premult with ctf

        ### create spectral wedge weighting
        if analytWedge:
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

        if wedgeInfo._type == 'Wedge3dCTF':  # square the 3d ctf for the weighting
            wedge = wedge * wedge

        ### shift and rotate particle
        shiftV = particleObject.getShift()
        newParticle.setAll(0)

        transform(particle, newParticle, -rotation[1], -rotation[0], -rotation[2],
                  centerX, centerY, centerZ, -shiftV[0], -shiftV[1], -shiftV[2], 0, 0, 0)

        if weighting:
            weight = 1. - particleObject.getScore().getValue()
            # weight = weight**2
            weight = exp(-1. * weight)
            result = result + newParticle * weight
            wedgeSum = wedgeSum + wedge * weight
        else:
            result = result + newParticle
            wedgeSum = wedgeSum + wedge

        if showProgressBar:
            numberAlignedParticles = numberAlignedParticles + 1
            progressBar.update(numberAlignedParticles)

        n += 1
    ###apply spectral weighting to sum
    result = lowpassFilter(result, sizeX / 2 - 1, 0.)[0]

    root, ext = os.path.splitext(averageName)

    # if createInfoVolumes:
    result.write(f'{root}-PreWedge{ext}')

    # wedgeSum = wedgeSum*0+len(particleList)
    wedgeSum.write(f'{root}-WedgeSumUnscaled{ext}')
    invert_WedgeSum(invol=wedgeSum, r_max=sizeX / 2 - 2., lowlimit=.05 * len(particleList),
                    lowval=.05 * len(particleList))

    if createInfoVolumes:
        w1 = reducedToFull(wedgeSum)
        w1.write(f'{root}-WedgeSumInverted{ext}')

    result = convolute(v=result, k=wedgeSum, kernel_in_fourier=True)

    # do a low pass filter
    # result = lowpassFilter(result, sizeX/2-2, (sizeX/2-1)/10.)[0]
    result.write(averageName)

    if createInfoVolumes:
        resultINV = result * -1
        # write sign inverted result to disk (good for chimera viewing ... )
        resultINV.write('{root}-INV{ext}')
    newReference = Reference(averageName, particleList)

    return newReference


def averageGPU(particleList, averageName, showProgressBar=False, verbose=False,
               createInfoVolumes=False, weighting=False, norm=False, gpuId=None, profile=False):
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
    from pytom.agnostic.filter import bandpass as lowpassFilter, applyFourierFilterFull
    from pytom.voltools import transform
    from pytom.basic.structures import Reference
    from pytom.agnostic.normalise import mean0std1
    from pytom.agnostic.tools import invert_WedgeSum
    from pytom.agnostic.transform import fourier_full2reduced
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

    sx, sy, sz = read_size(particleList[0].getFilename())
    wedgeInfo = particleList[0].getWedge().convert2numpy()

    # TODO ifftshift to shift centered spectrum back to corner!
    wedgeZero = xp.fft.ifftshift(wedgeInfo.returnWedgeVolume(sx, sy, sz, True))
    wedge = xp.zeros_like(wedgeZero, dtype=xp.float32)
    wedgeSum = xp.zeros_like(wedge, dtype=xp.float32)

    newParticle = xp.zeros((sx, sy, sz), dtype=xp.float32)

    centerX = sx // 2
    centerY = sy // 2
    centerZ = sz // 2

    result = xp.zeros((sx, sy, sz), dtype=xp.float32)

    fftplan = get_fft_plan(wedge.astype(xp.complex64))

    n = 0

    total = len(particleList)

    if profile:
        t_end = stream.record()
        t_end.synchronize()

        time_took = xp.cuda.get_elapsed_time(t_start, t_end)
        print(f'startup time {n:5d}: \t{time_took:.3f}ms')
        t_start = stream.record()

    for particleObject in particleList:

        rotation = particleObject.getRotation()
        rotinvert = rotation.invert()
        shiftV = particleObject.getShift()

        particle = read(particleObject.getFilename(), deviceID=device)

        if norm:  # normalize the particle
            mean0std1(particle)  # happen inplace

        # get the wedge per particle because the wedge can differ
        wedgeInfo = particleObject.getWedge().convert2numpy()
        wedgeZero = xp.fft.ifftshift(wedgeInfo.returnWedgeVolume(sx, sy, sz, True))

        # apply wedge to particle
        particle = (ifftnP(fftnP(particle, plan=fftplan) * wedgeZero, plan=fftplan)).real

        ### create spectral wedge weighting
        wedge *= 0
        transform(xp.fft.fftshift(wedgeZero), rotation=(rotinvert[0], rotinvert[2], rotinvert[1]),
                  rotation_order='rzxz', center=(centerX, centerY, centerZ), output=wedge, device=device,
                  interpolation='linear')

        ### shift and rotate particle
        newParticle *= 0
        transform(particle, output=newParticle, rotation=(-rotation[1], -rotation[2], -rotation[0]),
                  center=(centerX, centerY, centerZ), translation=(-shiftV[0], -shiftV[1], -shiftV[2]),
                  device=device, interpolation='filt_bspline', rotation_order='rzxz')

        # add to average and wedgeweighting
        result += newParticle
        wedgeSum += xp.fft.ifftshift(wedge)

        if n % total == 0:
            if profile:
                t_end = stream.record()
                t_end.synchronize()

                time_took = xp.cuda.get_elapsed_time(t_start, t_end)
                print(f'total time {n:5d}: \t{time_took:.3f}ms')
                t_start = stream.record()
        cstream.synchronize()
        n += 1

    ###apply spectral weighting to sum

    root, ext = os.path.splitext(averageName)

    result = lowpassFilter(result, high=sx / 2 - 1, sigma=0)
    # if createInfoVolumes:
    write(f'{root}-PreWedge{ext}', result)
    write(f'{root}-WedgeSumUnscaled{ext}', fourier_full2reduced(wedgeSum))

    # prev: wedgeSumINV =
    invert_WedgeSum(wedgeSum, r_max=sx // 2 - 2., lowlimit=.05 * len(particleList), lowval=.05 * len(particleList))

    if createInfoVolumes:
        # write(f'{root}-WedgeSumInverted{ext}', xp.fft.fftshift(wedgeSumINV))
        write(f'{root}-WedgeSumInverted{ext}', xp.fft.fftshift(wedgeSum))

    # the wedge sum should already be centered in the corner ?
    result = applyFourierFilterFull(result, wedgeSum)  # prev wedgeSumINV

    # do a low pass filter
    result = lowpassFilter(result, sx / 2 - 2, (sx / 2 - 1) / 10.)[0]
    write(averageName, result)

    if createInfoVolumes:
        resultINV = result * -1
        # write sign inverted result to disk (good for chimera viewing ... )
        write(f'{root}-INV{ext}', resultINV)

    newReference = Reference(averageName, particleList)

    return newReference


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
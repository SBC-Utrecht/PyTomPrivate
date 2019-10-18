#!/usr/bin/env pytom
'''
Created on Nov 4, 2014

@author: FF
@lastchange: added spherical mask option
'''

if __name__ == '__main__':
    # parse command line arguments

    import sys
    from pytom.tools.script_helper import ScriptHelper, ScriptOption
    from pytom.tools.parse_script_options import parse_script_options
    from pytom.tools.files import checkFileExists,checkDirExists
    from pytom.alignment.GLocalSampling import GLocalSamplingJob, mainAlignmentLoop
    from pytom.basic.structures import ParticleList, Reference, Mask, SampleInformation, PointSymmetry
    from pytom.score.score import FLCFScore
    # from pytom.frontend.serverpages.createAlignmentJob import createRunscripts
    from pytom.angles.localSampling import LocalSampling
    from pytom.alignment.preprocessing import Preprocessing
    
    helper = ScriptHelper(sys.argv[0].split('/')[-1], 
                          description='Create an GLocalSampling job. Documentation is available at\n\
                          http://www.pytom.org/doc/pytom/alignment.html',
                          authors='Friedrich Foerster',
                          options=[ScriptOption(['-p','--particleList'],
                                                'Particle list : xml file storing information to all subvolumes',
                                                'string', 'required'),
                                   ScriptOption(['-r','--reference'],
                                                'Reference : the initial reference - if none provided average of particle list',
                                                'string', 'optional'),
                                   ScriptOption(['-m','--mask'], 'Mask : a path to a mask ', 'string', 'required'),
                                   ScriptOption(['--SphericalMask'], 'Mask is spherical / speed up!', 'no arguments',
                                                'optional', False),
                                   ScriptOption(['--angleShells'], '# angle shells used for angular refinement.',
                                                'uint', 'optional', 3),
                                   ScriptOption(['--angleIncrement'], 'Angular increment for refinement.',
                                                'uint', 'optional', 3),
                                   ScriptOption(['--symmetry'], 'PointSymmetry : specify n-fold symmetry (n)',
                                                'int', 'optional'),
                                   ScriptOption(['--symmetryAngleZ'], 'PointSymmetry axis tilt around Z axis',
                                                'float', 'optional'),
                                   ScriptOption(['--symmetryAngleX'], 'PointSymmetry axis tilt around X axis',
                                                'float', 'optional'),
                                   ScriptOption(['-d', '--destination'], 'Destination : destination directory',
                                                'string', 'required'),
                                   ScriptOption(['-n', '--numberIterations'], 'Number of iterations',
                                                'string', 'required'),
                                   ScriptOption(['-b','--binning'],
                                                'Perform binning (downscale) of subvolumes by factor.',
                                                'uint', 'optional', 1),
                                   ScriptOption(['--pixelSize'], 'Pixelsize in Angstrom', 'float', 'required'),
                                   ScriptOption(['--particleDiameter'], 'Particle diameter in Angstrom', 'float',
                                                'required', -1),
                                   ScriptOption(['-w', '--weighting'], 'Weight particles by exp of CC', 'no arguments',
                                                'optional', False),
                                   ScriptOption(['-c', '--compound'], 'Use compound weighting in Fourier space',
                                                'no arguments', 'optional', False),
                                   ScriptOption(['-j','--jobName'], 'Specify job.xml output filename', 'string',
                                                'required'),
                                   ScriptOption(['--noShift'], 'Remove all shifts from the particlelist, useful for particle polishing', 'no arguments', 'optional', False)])

    particleList, reference, mask, isSphere, angShells, angleInc, symmetryN, symmetryAxisZ, symmetryAxisX,\
    destination, numberIterations, binning,\
    pixelSize, diameter, weighting, compound, jobName, noshift = parse_script_options(sys.argv[1:], helper)

    #particleList
    if not checkFileExists(particleList):
        raise RuntimeError('ParticleList file ' + particleList + ' does not exist!')
    pl = ParticleList()
    pl.fromXMLFile(particleList)

    # Remove all shifts, useful for particle polishing
    if noshift:
        for particle in pl:
            particle.getShift().scale(0)

    if reference:
        if not checkFileExists(reference):
            raise RuntimeError('Reference file ' + reference + ' does not exist!')
        ref = Reference(referenceFile=reference)
    else:
        ref = Reference()
    
    if not checkFileExists(mask):
        raise RuntimeError('Mask file ' + mask + ' does not exist!')

    m = Mask(filename=mask, isSphere=isSphere)
    
    if not checkDirExists(destination):
        raise RuntimeError('Destination directory ' + destination + ' does not exist!')

    rot = LocalSampling(angShells,angleInc)

    sampleI = SampleInformation(pixelSize=pixelSize, particleDiameter=diameter)

    if symmetryN is None or symmetryAxisZ is None or symmetryAxisX is None:
        sym = None
    else:
        sym = PointSymmetry(nfold=symmetryN,z2=symmetryAxisZ,x=symmetryAxisX)

    ################# fixed parameters #################
    score   = FLCFScore()
    score.setRemoveAutocorrelation(flag=False)
    #pre     = Preprocessing(lowestFrequency = float(lowestFrequency),highestFrequency = float(highestFrequency))
    pre     = Preprocessing()
    adaptive_res = 0.1
    fsc_criterion = 0.143

    locJob = GLocalSamplingJob(pl=pl, ref=ref, mask=m,
                               sample_info=sampleI, rotations=rot,
                               preprocessing=pre,
                               dest=destination, max_iter=int(numberIterations), score=score, binning=int(binning),
                               weighting=weighting, compoundWedge=compound,
                               symmetries=None, adaptive_res=adaptive_res, fsc_criterion=fsc_criterion)
    locJob.toXMLFile(jobName)
    # run script
    mainAlignmentLoop( alignmentJob=locJob, verbose=False)
    

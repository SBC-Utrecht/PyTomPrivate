#!/usr/bin/env pytom
from pytom.basic.structures import ParticleList
from pytom.alignment.averaging import averageParallel, averageParallelGPU
import os


if __name__=='__main__':
    import sys
    from pytom.tools.script_helper import ScriptHelper, ScriptOption
    from pytom.tools.parse_script_options import parse_script_options

    options = [ScriptOption(['-p','--particleList'], 'Particle List', True, False),
               ScriptOption(['-a','--averageName'], 'Filename of output average.', True, False),
               ScriptOption(['-c','--numberOfCores'],'Number of Cores used for average calculation', True, False),
               ScriptOption(['-w','--weighting'],'Weight particles by exp CC in average. False by default.', False, True),
               ScriptOption(['-v','--verbose'],'Print particle information. False by default.', False, True),
               ScriptOption(['-s','--showProgressBar'],'Show progress bar. False by default.', False, True),
               ScriptOption(['-i','--createInfoVolumes'],'Create Info data (wedge sum, inverted density) too? False by default.', False, True),
               ScriptOption(['-n','--normalize'],'Normalize average. False by default.', False, True),
               ScriptOption(['-g', '--gpuID'], 'Provide a gpu if you want to use one (one only)', True, True),
               ScriptOption(['-h', '--help'], 'Help.', False, True)]


    helper = ScriptHelper(sys.argv[0].split('/')[-1], # script name
                          description='Convert coordinate list to particle list.',
                          authors='Friedrich Foerster',
                          options=options)
    if len(sys.argv) == 1:
        print(helper)
        sys.exit()
    try:
        plName, outname, cores, weighting, verbose, showProgressBar, createInfoVol, norm, gpuID, help = dd = parse_script_options(sys.argv[1:], helper)
    except Exception as e:
        print(e)
        sys.exit()

    if help is True:
        print(helper)
        sys.exit()

    try: cores = int(cores)
    except: cores = 1

    if not os.path.exists(plName):
        print('Please provide an existing particle list')
        sys.exit()

    try:

        gpuId = None if gpuID is None else int(gpuID)
        pnr = 3 if gpuID is None else 1
    except Exception as e:
        print(e)
        if ',' in gpuID:
            print('\n\nPlease provide only one gpu')
        sys.exit()

    even = ParticleList()
    even.fromXMLFile(plName)

    if gpuID is None:
        averageParallel(particleList=even,
                           averageName=outname,
                           showProgressBar=showProgressBar, verbose=verbose, createInfoVolumes=createInfoVol,
                           weighting=weighting, norm=norm,
                           setParticleNodesRatio=pnr, cores=cores)

    else:
        averageParallelGPU(particleList=even,
                           averageName=outname,
                           showProgressBar=showProgressBar, verbose=verbose, createInfoVolumes=createInfoVol,
                           weighting=weighting, norm=norm,
                           setParticleNodesRatio=pnr, cores=1, gpuID=gpuID)
    

#!/usr/bin/env pytom
'''
Created on Jul 19, 2011

@author: hrabe
'''

    

if __name__ == '__main__':
    
    #this script is significantly linked to 
    #pytom/frontend/serverpages/createReconstructionJob.py
    #any changes like parameter names must be changed in the other script, too!

    import sys,getopt
    from pytom.tools.script_helper import ScriptHelper, ScriptOption
    from pytom.tools.parse_script_options import parse_script_options
    from pytom.basic.structures import Particle, ParticleList
    from pytom.reconstruction.reconstructionStructures import ProjectionList
    from pytom.tools.files import checkFileExists,checkDirExists
    from pytom.basic.files import read_em_header
    
    helper = ScriptHelper(sys.argv[0].split('/')[-1],
                      description='Reconstruct particles in a particle list. Documentation is available at\n\
                          http://www.pytom.org/doc/pytom/resonstructTomograms.html or\n\
                          http://www.pytom.org/doc/pytom/resonstructSubtomograms.html',
                      authors='Thomas Hrabe, FF',

                      options= [ScriptOption(['-t','--tomogram'], 'Reconstruct a tomogram. Specify name of tomogam here. You do not need a particle list for that!', 'has arguments', 'optional'),
                                ScriptOption(['-p','--particleList'], 'XML particle list.', 'has arguments', 'required'),
                                ScriptOption(['--projectionList'], 'XML projection list.', 'has arguments', 'optional'),
                                ScriptOption(['--projectionDirectory'], 'Directory containing the projections.', 'has arguments', 'required'),
                                ScriptOption(['-w','--applyWeighting'], 'If projections are not weighted, apply weighting before. If omited, no weighting.', 'no arguments', 'optional'),
                                ScriptOption(['-s','--size'], 'Size of particle cube / tomogram.', 'has arguments', 'required'),
                                ScriptOption(['-b','--coordinateBinning'], 'Binning factor of coordinates. If particle coordinates are determined in binned volume (with respect to projections) this binning factor needs to be specified.', 'has arguments', 'optional'),
                                ScriptOption(['-o','--recOffset'], 'Cropping offset of the binned tomogram.', 'has arguments', 'required'),
                                ScriptOption(['--projBinning'], 'Bin projections BEFORE reconstruction. 1 is no binning, 2 will merge two voxels to one, 3 -> 1, 4 ->1 ...', 'has arguments', 'optional'),
                                ScriptOption(['-a', '--alignResultFile'], 'Use an alignmentResultFile to generate the aligned files in memory.', 'has arguments', 'optional'),
                                ScriptOption(['-r', '--particlePolishFile'], 'Use an particlePolishFile to generate the polished cutouts.', 'has arguments', 'optional'),
                                ScriptOption(['-n', '--numProcesses'], 'The number of processes to use. Default: 10', 'has arguments', 'optional'),
                                ScriptOption(['--numReadProcesses'], 'The maximum number of reading processes. BEWARE!'
                                              ' a single read thread easily consumes up to 4 Gb (for a 4k picture) so '
                                              'keep that in mind to not overload the nodes. If not specified this will '
                                              'be the same as numProcesses.', 'has arguments', 'optional')])
        
    particleList = None
    projectionDirectory = None
    aw = False
    
    try:
        tomogram, particleListXMLPath, projectionList, projectionDirectory, aw, size, coordinateBinning, recOffset, \
        projBinning, alignmentResultFile, particlePolishFile, numProcesses, numReadProcesses = parse_script_options(sys.argv[1:], helper)
    
    except Exception as e:
        print e
        sys.exit()

    print(projectionDirectory, projectionList)
   
    size = [int(i) for i in size.split(',')]
    if len(size) == 1:
        tmp = size[0]
        size.append(tmp)
        size.append(tmp)

    if projBinning:
        projBinning = int(projBinning)
    else:
        projBinning = 1

    if numProcesses is None:
        numProcesses = 10
    else:
        numProcesses = int(numProcesses)

    if numReadProcesses is None:
        numReadProcesses = numProcesses
    else:
        numReadProcesses = int(numReadProcesses)

    if coordinateBinning:
        coordinateBinning = float(coordinateBinning)
    else:
        coordinateBinning = 1
    if not aw:
        aw = False

    if recOffset:
        recOffset = [int(i) for i in recOffset.split(',')]
    else:
        recOffset = [0.,0.,0.]
    
        
    projections = ProjectionList()
    if checkFileExists(projectionList):
        projections.fromXMLFile(projectionList)
    elif checkDirExists(projectionDirectory):
        projections.loadDirectory(projectionDirectory)
    else:
        raise RuntimeError('Neither projectionList existed nor the projectionDirectory you specified! Abort')


    if checkFileExists(particlePolishFile):
        from pytom.basic.datatypes import LOCAL_ALIGNMENT_RESULTS
        import numpy

        polishedCoordinates = numpy.loadtxt(particlePolishFile,dtype=LOCAL_ALIGNMENT_RESULTS)


    if alignmentResultFile and not checkFileExists(alignmentResultFile):
        raise Exception('alignmentResultFile does not exists. Please provide an existing file or omit this flag.')

    if tomogram:
        vol = projections.reconstructVolume( dims=size, reconstructionPosition=recOffset,
            binning=projBinning, applyWeighting=aw)
        vol.write(tomogram)
        
    else:
        # transform the cropping offset
        tmp = projections[0]
        sx = tmp.getXSize() # here should be the size of original projection!
        sy = tmp.getYSize()
        recOffset[0] = -sx/2 + recOffset[0]*coordinateBinning
        recOffset[1] = -sy/2 + recOffset[1]*coordinateBinning
        recOffset[2] = -sx/2 + recOffset[2]*coordinateBinning

        # set particle list in order to reconstruct subtomograms
        particleList = ParticleList()

        try:
            particleList.fromXMLFile(particleListXMLPath)
        except RuntimeError:
            print 'Error reading particleList XML file! Abort'
            sys.exit()

        if not len(polishedCoordinates['AlignmentTransX']) == len(particleList) * len(projections):
            raise Exception("The length of the polished alignment list does not correspond to the theoretical length, "
                            "are you sure every parameter is correct? polished list len {:d}, particle  list len {:d}"
                            " and projection len {:d}".format(len(polishedCoordinates['AlignmentTransX']), len(particleList), len(projections)))

        from pytom.basic.structures import PickPosition
        for n, particle in enumerate(particleList):
            pickPosition = particle.getPickPosition()
            x = (pickPosition.getX()*coordinateBinning+ recOffset[0])/projBinning
            y = (pickPosition.getY()*coordinateBinning+ recOffset[1])/projBinning
            if particlePolishFile:
                x += polishedCoordinates['AlignmentTransX'][n] / float(projBinning)
                y += polishedCoordinates['AlignmentTransY'][n] / float(projBinning)
                #particle.setFilename(particle.getFilename()[:-3]+'_polished.em')
            z = (pickPosition.getZ()*coordinateBinning+ recOffset[2])/projBinning
            particle.setPickPosition( PickPosition(x=x, y=y, z=z))

        projections.reconstructVolumes(particles=particleList, cubeSize=int(size[0]), \
                                       binning=projBinning, applyWeighting=aw, \
                                       showProgressBar=True, verbose=False, \
                                       preScale=projBinning, postScale=1, alignResultFile=alignmentResultFile,
                                       num_procs=numProcesses, num_procs_read=numReadProcesses)

            




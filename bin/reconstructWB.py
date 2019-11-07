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

                      options= [ScriptOption(['-t','--tomogram'], 'Reconstruct a tomogram. Specify name of tomogam here. You do not need a particle list for that!', 'string', 'optional'),
                                ScriptOption(['-p','--particleList'], 'XML particle list.', 'string', 'required'),
                                ScriptOption(['--projectionList'], 'XML projection list.', 'string', 'optional'),
                                ScriptOption(['--projectionDirectory'], 'Directory containing the projections.', 'string', 'required'),
                                ScriptOption(['-w','--applyWeighting'], 'If projections are not weighted, apply weighting before. If omited, no weighting.', 'no arguments', 'optional', False),
                                ScriptOption(['-s','--size'], 'Size of particle cube / tomogram.', 'has arguments', 'required'),
                                ScriptOption(['-b','--coordinateBinning'], 'Binning factor of coordinates. If particle coordinates are determined in binned volume (with respect to projections) this binning factor needs to be specified.', 'float', 'optional', 1),
                                ScriptOption(['-o','--recOffset'], 'Cropping offset of the binned tomogram.', 'int,int,int', 'required'),
                                ScriptOption(['--projBinning'], 'Bin projections BEFORE reconstruction. 1 is no binning, 2 will merge two voxels to one, 3 -> 1, 4 ->1 ...', 'uint', 'optional', 1),
                                ScriptOption(['-a', '--alignResultFile'], 'Use an alignmentResultFile to generate the aligned files in memory.', 'string', 'optional'),
                                ScriptOption(['-r', '--particlePolishFile'], 'Use an particlePolishFile to generate the polished cutouts.', 'string', 'optional'),
                                ScriptOption(['-n', '--numProcesses'], 'The number of processes to use.', 'uint', 'optional', 10),
                                ScriptOption(['--numReadProcesses'], 'The maximum number of reading processes. BEWARE!'
                                              ' a single read thread easily consumes up to 4 Gb (for a 4k picture) so '
                                              'keep that in mind to not overload the nodes. If not specified this will '
                                              'be the same as numProcesses.', 'uint', 'optional'),
                                ScriptOption(['--dimZ'], 'The dimension on the Z axis, default is the same as dimension X', 'uint', 'optional'),
                                ScriptOption(['--notPolished'], "Do not apply the shifts in the particlelist, for testing purposes", 'no arguments', 'optional', False)])
        
    particleList = None

    tomogram, particleListXMLPath, projectionList, projectionDirectory, aw, size, coordinateBinning, recOffset, \
    projBinning, alignmentResultFile, particlePolishFile, numProcesses, numReadProcesses, dimz, notpolished = parse_script_options(sys.argv[1:], helper)

    print(projectionDirectory, projectionList)
   
    size = [int(i) for i in size.split(',')]
    if len(size) == 1:
        tmp = size[0]
        size.append(tmp)
        size.append(tmp)
    elif len(size) == 2 or len(size) > 3: raise RuntimeError("You should specify 1 or 3 sizes for the 3 dimensional cube.")

    if numReadProcesses is None:
        numReadProcesses = numProcesses
    else:
        numReadProcesses = numReadProcesses
        
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

        polishedCoordinates = numpy.loadtxt(particlePolishFile, dtype=LOCAL_ALIGNMENT_RESULTS)
    else:
        polishedCoordinates = None

    if alignmentResultFile and not checkFileExists(alignmentResultFile):
        raise Exception('alignmentResultFile does not exists. Please provide an existing file or omit this flag.')

    if tomogram:
        vol = projections.reconstructVolume( dims=size, reconstructionPosition=recOffset,
            binning=projBinning, applyWeighting=aw)
        vol.write(tomogram)
        
    else:
        # transform the cropping offset
        tmp = projections[0]
        sx = tmp.getXSize()  # here should be the size of original projection!
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

        if polishedCoordinates is not None:
            if not len(polishedCoordinates['AlignmentTransX']) == len(particleList) * len(projections):
                raise Exception("The length of the polished alignment list does not correspond to the theoretical length, "
                                "are you sure every parameter is correct? polished list len {:d}, particle  list len {:d}"
                                " and projection len {:d}".format(len(polishedCoordinates['AlignmentTransX']), len(particleList), len(projections)))

        from pytom.basic.structures import PickPosition
        for n, particle in enumerate(particleList):
            pickPosition = particle.getPickPosition()
            x = (pickPosition.getX() * coordinateBinning + recOffset[0]) / projBinning
            y = (pickPosition.getY() * coordinateBinning + recOffset[1]) / projBinning
            z = (pickPosition.getZ() * coordinateBinning + recOffset[2]) / projBinning
            particle.setPickPosition(PickPosition(x=x, y=y, z=z))

            # Shifting each particle (by the shift of another particle) does not really help with aligning to every
            # specific projection, I should build support for the particlepolishfile into reconstructVolumes
            # to shift every single projection in a unique way for every single particle

            # if particlePolishFile:
            #    x -= polishedCoordinates['AlignmentTransX'][n] / float(projBinning)
            #    y -= polishedCoordinates['AlignmentTransY'][n] / float(projBinning)

        projections.reconstructVolumes(particles=particleList, cubeSize=int(size[0]),
                                       binning=projBinning, applyWeighting=aw,
                                       showProgressBar=True, verbose=False,
                                       preScale=projBinning, postScale=1, alignResultFile=alignmentResultFile,
                                       num_procs=numProcesses, num_procs_read=numReadProcesses, particle_polish_file=polishedCoordinates, dimz=dimz, notpolished=notpolished, coordbinning=coordinateBinning, offset=recOffset)

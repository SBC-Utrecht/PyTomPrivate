#!/usr/bin/env pytom
"""
Created on Jul 20, 2013

@author: FF
"""

if __name__ == '__main__':
    import sys
    #from pytom.reconstruction.TiltAlignmentStructures import TiltAlignmentParameters, TiltSeries, TiltAlignment
    from pytom.tools.script_helper import ScriptHelper, ScriptOption
    from pytom.tools.parse_script_options import parse_script_options
    from pytom.reconstruction.reconstructionFunctions import alignWeightReconstruct

    options=[ScriptOption(['--tiltSeriesName'], 'Name tilt series - either prefix of sequential tilt series files \
             expected as "tiltSeriesName_index.em/mrc" or full name of stack "tiltSeriesName.st"',
                          'string', 'required'),
             ScriptOption(['--tiltSeriesFormat'], 'Format of tilt series (series of "em" or "mrc" images or "st" stack).',
                          'string', 'optional', 'em'),
             ScriptOption(['--firstIndex'], 'Index of first projection.', 'uint', 'optional', 1),
             ScriptOption(['--lastIndex'], 'Index of last projection.', 'uint', 'required'),
             ScriptOption(['--tltFile'], 'tltFile containing tilt angles.', 'string', 'optional'),
             ScriptOption(['--prexgFile'], 'prexgFile containing pre-shifts from IMOD.', 'string', 'optional'),
             ScriptOption(['--preBin'], 'pre-Binning in IMOD prior to marker determination.', 'int', 'optional'),
             ScriptOption(['--referenceIndex'], 'Index of reference projection used for alignment.', 'int',
                          'required'),
             ScriptOption(['--markerFile'], 'Name of EM markerfile or IMOD wimp File containing marker coordinates.',
                          'string', 'required'),
             ScriptOption(['--referenceMarkerIndex'], 'Index of reference marker to set up coordinate system.',
                          'int', 'required', 1),
             ScriptOption(['--handflip'], 'Is your tilt series outside of 0-180deg (Specify if yes).', 'no arguments',
                          'optional', False),
             ScriptOption(['--projectionTargets'],
                          'Relative or absolute path to the aligned projections that will be generated + file prefix.',
                          'string', 'optional', 'align/myTilt'),
             ScriptOption(['--fineAlignFile'],
                          'Relative or absolute path to the file with fineAlign parameters (type should be *.dat).',
                          'string', 'optional'),
             ScriptOption(['--projectionBinning'], 'Binning of projections during read.', 'uint',
                          'optional', 1),
             ScriptOption(['--lowpassFilter'], 'Lowpass filter in Nyquist after binning.', 'float', 'optional', 1.0),
             ScriptOption(['--tomogramFile'],
                          'Relative or absolute path to final tomogram (no tomogram written if not specified).',
                          'string', 'optional'),
             ScriptOption(['--fileType'], 'File type (can be em or mrc - no tomogram written if not specified).',
                          'string', 'optional'),
             ScriptOption(['--tomogramSize'], 'Size of tomogram (no tomogram written if not specified).',
                          'int,int,int', 'optional', [0, 0, 0]),
             ScriptOption(['--reconstructionCenter'],
                          'Center where tomogram will be reconstructed (no tomogram written if not specified).',
                          'int,int,int', 'optional', [0, 0, 0]),
             ScriptOption(['--weightingType'], 'Type of weighting (-1 default r-weighting, 0 no weighting)', 'int',
                          'optional', -1),
             ScriptOption(['--verbose'], 'Enable verbose mode', 'no arguments', 'optional')]
    
    helper = ScriptHelper(sys.argv[0].split('/')[-1],
                          description='Align and weight projections, save them and reconstruct tomogram (optional). \n\
                                      See http://pytom.org/doc/pytom/reconstructTomograms.html for documentation.',
                          authors='Friedrich Foerster',
                          options = options)

    tiltSeriesName, tiltSeriesFormat, firstProj, lastProj, \
    tltFile, prexgFile, preBin, referenceIndex, markerFileName, referenceMarkerIndex, handflip, \
    alignedTiltSeriesName, fineAlignFile, projBinning, lowpassFilter, \
    volumeName, filetype, \
    voldims, \
    reconstructionPosition, \
    weightingType, verbose = parse_script_options(sys.argv[1:], helper)

    # input parameters
    #tiltSeriesName = tiltSeriesPath + tiltSeriesPrefix  # "../projections/tomo01_sorted" # ending is supposed to be tiltSeriesName_index.em (or mrc)

    # only write projections and do NOT reconstruct tomogram (following parameters would be obsolete)
    onlyWeightedProjections = False
    if not volumeName:
        onlyWeightedProjections = True
         
    if filetype is None:
        if volumeName:
            filetype = volumeName.split('.')[-1]
    else: #TODO bug? shouldn't this be one identation level in? otherwise it overwrites userinput
        filetype = 'em'

    if voldims == [0, 0, 0]:
        onlyWeightedProjections = True

    outMarkerFileName = 'MyMarkerFile.em'
    if verbose:
        print "Tilt Series: "+str(tiltSeriesName)+", "+str(firstProj)+"-"+str(lastProj)
        print "Index of Reference Projection: "+str(referenceIndex)
        print "Marker Filename: "+str(markerFileName)
        print "TltFile: "+str(tltFile)
        print "prexgFile: "+str(prexgFile)
        print "Index of Reference Marker: "+str(referenceMarkerIndex)
        print "Handflip: "+str(handflip)
        print "Projection Targets: "+str(alignedTiltSeriesName)
        print "FineAlignmentFile: "+str(fineAlignFile)
        print "Binning Factor of Projections: "+str(projBinning)+", lowpass filter (in Ny): "+str(lowpassFilter)
        print "Name of Reconstruction Volume: "+str(volumeName)+" of Filetype: "+str(filetype)
        print "Reconstruction size: "+str(voldims)
        print "Reconstruction center: "+str(reconstructionPosition)
        print "write only aligned projections out: "+str(onlyWeightedProjections)

    alignWeightReconstruct(tiltSeriesName=tiltSeriesName, markerFileName=markerFileName, lastProj=lastProj,
                           tltfile=tltFile, prexgfile=prexgFile, preBin=preBin,
                           volumeName=volumeName, volumeFileType=filetype,
                           voldims=voldims, recCent=reconstructionPosition,
                           tiltSeriesFormat=tiltSeriesFormat, firstProj=firstProj, irefmark=referenceMarkerIndex, ireftilt=referenceIndex,
                           handflip=handflip,
                           alignedTiltSeriesName=alignedTiltSeriesName,
                           weightingType=weightingType,
                           lowpassFilter=lowpassFilter, projBinning=projBinning,
                           outMarkerFileName=outMarkerFileName, verbose=True)

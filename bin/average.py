#!/usr/bin/env pytom

"""
Created on Aug 24, 2011

@author: hrabe
"""



if __name__ == '__main__':
# parse command line arguments
    import sys
    from pytom.tools.script_helper import ScriptHelper, ScriptOption
    from pytom.tools.parse_script_options import parse_script_options
    from pytom.basic.structures import ParticleList
    from pytom.alignment.alignmentFunctions import distributeAverage
    from pytom_mpi import isInitialised
    helper = ScriptHelper(sys.argv[0].split('/')[-1],
                          description='Averge a particle. Documentation is available at\n\
                          http://www.pytom.org/doc/pytom/average.html',
                          authors='Thomas Hrabe',
                          options= [ScriptOption(['-p','--particleList'], 'XML particle list to be averaged.', 'string', 'required'),
                                    ScriptOption(['-a','--average'], 'Resulting average filename.', 'string', 'required'),
                                    ScriptOption(['--startIndex'], 'Average a subset of particles: start particle index.', 'int', 'optional'),
                                    ScriptOption(['--endIndex'], 'Average a subset of particles: last particle index.', 'int', 'optional'),
                                    ScriptOption(['--minimumScore'], 'Average a subset of particles: minimum score of particle. Does not work when start & end index are specified!', 'float', 'optional'),
                                    ScriptOption(['--maximumScore'], 'Average a subset of particles: maximum score of particle. Does not work when start & end index are specified!', 'float', 'optional'),
                                    ScriptOption(['--infoVolumes'], 'Generate info volumes like wedge volume.', 'no arguments', 'optional', False),
                                    ScriptOption(['--progressbarOff'], 'Display a progressbar. On by default', 'no arguments', 'optional', False),
                                    ScriptOption(['--fromAlignmentList'], 'Average from alignment list XML instead from particleListXML. Optional, off by default.', 'string', 'optional'),
                                    ScriptOption(['--subregion'], 'Average of particles whose coordinates are in [startX,startY,startZ,endX,endY,endZ].', '[int,int,int,int,int,int]', 'optional')])

    particleListName, averageName, startIndex, endIndex, minimum, maximum, infoVolumes, progressbarOff, alignmentFileName ,subRegion = parse_script_options(sys.argv[1:], helper)
        
    pl = None
    
    if particleListName:
        pl = ParticleList('/')
        pl.fromXMLFile(particleListName)
    
    if startIndex is not None and endIndex is not None:
        pl = pl[startIndex:endIndex]

    if subRegion:
        newPl = ParticleList()

        for aParticle in pl:
            position = aParticle.getPickPosition()

            xOK = subRegion[0] <= position.getX() <= subRegion[3]
            yOK = subRegion[1] <= position.getY() <= subRegion[4]
            zOK = subRegion[2] <= position.getZ() <= subRegion[5]

            if xOK and yOK and zOK:
                newPl.append(aParticle)

        if len(newPl) > 0:
            pl = newPl
        
    if startIndex is None and endIndex is None and minimum is not None and maximum is not None:
        newPl = ParticleList() 
    
        for aParticle in pl:
            if minimum <= aParticle.getScore().getValue() <= maximum:
                newPl.append(aParticle)
                
        if len(newPl) > 0:
            pl = newPl

    if not particleListName and alignmentFileName:
        pl = ParticleList('/')
        pl.fromAlignmentList(alignmentFileName)
        
    distributeAverage(pl, averageName, progressbarOff, False, infoVolumes, sendEndMessage=True)

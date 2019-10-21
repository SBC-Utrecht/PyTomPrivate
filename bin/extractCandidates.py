#!/usr/bin/env pytom

'''
Created on Jun 9, 2010

@author: chen
'''


def extractCandidates(jobFilename='', resultFilename='', orientFilename='', sizeParticle=None, maxNumParticle=0, minScore=-1, write2disk=0, margin=None):
    # construct the original job from the xml file
    if jobFilename=='':
        jobFilename='JobInfo.xml'
        
    from pytom.localization.peak_job import PeakJob
    job = PeakJob()
    job.fromXMLFile(jobFilename)
    
    from pytom.localization.extraction_peak_result import ExPeakResult
    res = ExPeakResult()
    res.volFilename = job.volume.getFilename()
    
    if resultFilename=='':
        res.resultFilename='scores.em'
    else:
        res.resultFilename = resultFilename
    if orientFilename=='':
        res.orientFilename='angles.em'
    else:
        res.orientFilename = orientFilename
        
    res.angleList = job.rotations[:]
    res.score = job.score.__class__
    
    if maxNumParticle <= 0:
        return None
    
    res.readAll()
    
    if sizeParticle==None:
        ref = job.reference.getVolume()
        sizeParticle = [ref.sizeX(),ref.sizeY(),ref.sizeZ()]
    
    particleList = res.findParticles(sizeParticle,maxNumParticle,minScore,write2disk,margin)
    
    return particleList

if __name__ == '__main__':
    import sys, os
    from pytom.tools.script_helper import ScriptHelper, ScriptOption
    from pytom.tools.parse_script_options import parse_script_options
     
    helper = ScriptHelper(sys.argv[0].split('/')[-1],
                          description='Extract candidate molecules from localization result.',
                          authors='Yuxiang Chen',
                          options= [ScriptOption(['-j', '--jobFile'], 'Localization job XML file.', 'has arguments', 'required'),
                                    ScriptOption(['-r', '--result'], 'File with score coefficients (score.em).', 'has arguments', 'required'),
                                    ScriptOption(['-o', '--orientation'], 'File with orientation indices (angles.em).', 'has arguments', 'required'),
                                    ScriptOption(['-n', '--numberCandidates'], 'Number of candidates to extract.', 'has arguments', 'required'),
                                    ScriptOption(['-s', '--size'], 'Radius around potential candidate that will be ignored during further processing.', 'has arguments', 'required'),
                                    ScriptOption(['-p', '--particleList'], 'Name of particle list XML file.', 'has arguments', 'required'),
                                    ScriptOption(['-t', '--particlePath'], 'Path prepended to each particle.', 'has arguments', 'optional'),
                                    ScriptOption(['-v', '--minimalScoreValue'], 'Minimal score value to which to extract.', 'float', 'optional', -1),
                                    ScriptOption(['-m', '--motlList'], 'Write a MOTL file with candidates. The name of the file will be an extension of the particle list with .em.', 'has arguments','optional'),
                                    ScriptOption(['-g', '--margin'], 'Size of outer margin that will be ignored for potential candidates.', 'has arguments','optional'),
                                    ScriptOption(['-w', '--sizeCubes'], 'If specified, it will cut out candidates from the original tomogram with the specified size.', 'has arguments','optional'),
                                    ScriptOption(['--scale'], 'Scale coordinates by a factor. Set > 1 to adjust to larger volumes. Use 2 if the localization tomo was 1x binned.', 'has arguments', 'optional')])

    jobFilename, resultFilename, orientFilename, maxNumParticle, sizeParticle, plFilename, particlePath, minScore, motlFilename, margin, write2disk, scale = parse_script_options(sys.argv[1:], helper)

    if write2disk == None:
        write2disk = 0
    if margin.__class__ == str:
        margin = int(margin)
    if particlePath is None:
        particlePath = './'
    if not scale is None:
        scale = float(scale)
    else:
        scale = 1.0

    res=extractCandidates(jobFilename,resultFilename,orientFilename,int(sizeParticle),int(maxNumParticle),minScore,int(write2disk),margin)
    if not plFilename and not motlFilename:
        raise RuntimeError('You must specify at least a particle list or a motl file as result of this script!')
    
    if plFilename:
        from pytom.basic.structures import ParticleList
        pl = ParticleList()
        
        from pytom.localization.peak_job import PeakJob
        job = PeakJob()
        job.fromXMLFile(jobFilename)

        wedge = job.wedge

        if particlePath[-1] != os.sep:
            particlePath += os.sep
        
        for particle in res:
            newParticle = particle.toParticle()
            newParticle.setWedge(wedge)
            newParticle.setFilename(particlePath + newParticle.getFilename())

            if scale != 1.0:
                pi = newParticle.getPickPosition()
                pi.setX(pi.getX() * scale)
                pi.setY(pi.getY() * scale)
                pi.setZ(pi.getZ() * scale)
                newParticle.setPickPosition(pi)

            pl.append(newParticle)
        
        pl.toXMLFile(plFilename)
        
    
    if motlFilename:
        from pytom.basic.structures import ParticleList
        
        pl = ParticleList()
        for newParticle in res:

            if scale != 1.0:
                pi = newParticle.getPickPosition()
                pi.setX(pi.getX() * scale)
                pi.setY(pi.getY() * scale)
                pi.setZ(pi.getZ() * scale)
                newParticle.setPickPosition(pi)

            pl.append(newParticle.toParticle())
            
        pl.toMOTL(motlFilename)

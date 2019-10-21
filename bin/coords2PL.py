#!/usr/bin/env pytom
"""
    little script to convert coordinate list to particle list xml file
    FF Jan 2013
"""
from pytom.basic.structures import ParticleList

def convertCoords2PL(coordinate_file, particleList_file, subtomoPrefix=None, 
        wedgeAngle=None):
    pl = ParticleList()
    pl.loadCoordinateFile( filename=coordinate_file, name_prefix=subtomoPrefix, 
        wedgeAngle=wedgeAngle)
    pl.toXMLFile(particleList_file)


if __name__ == '__main__':
    import sys
    from pytom.tools.script_helper import ScriptHelper, ScriptOption
    from pytom.tools.parse_script_options import parse_script_options
    
    helper = ScriptHelper(
        sys.argv[0].split('/')[-1], # script name
        description='Convert coordinate list to particle list.',
        authors='Friedrich Foerster',
        options=[ScriptOption(['-p', '--particleList'], 'Particle List', 'string', 'required'),
                 ScriptOption(['-c', '--coords'], 'Coordinate List (ascii file from EMAN2)', 'has arguments', 'required'),
                 ScriptOption(['-s', '--subtomoPrefix'], 'path and filename for subtomogram files'
                                                         ' (e.g., MyPath/particle_)', 'string', 'optional'),
                 ScriptOption(['-w', '--wedgeAngles'], 'missing wedge angle(s) [counter-clock, clock] or single angle',
                              'has arguments', 'optional')])

    plName, coordName, subtomoPrefix, w = parse_script_options(sys.argv[1:], helper)

    if w:
        if len(w.split(',')) > 1:
            wedgeAngle = []
	    for kk in w.split(','):
	        wedgeAngle.append(float(kk))
        else:
	    wedgeAngle = float(w)
    else:
        wedgeAngle = None

    convertCoords2PL(coordinate_file=coordName, particleList_file=plName, 
        subtomoPrefix=subtomoPrefix, wedgeAngle=wedgeAngle)


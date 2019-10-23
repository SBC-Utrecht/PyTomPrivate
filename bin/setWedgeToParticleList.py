#!/usr/bin/env pytom
'''
Created on Sep 25, 2012

@author: yuxiangchen
'''

if __name__ == '__main__':
    # parse command line arguments
    import sys
    from pytom.tools.script_helper import ScriptHelper, ScriptOption
    from pytom.tools.parse_script_options import parse_script_options
    from pytom.basic.structures import ParticleList, SingleTiltWedge
    
    helper = ScriptHelper(sys.argv[0].split('/')[-1], # script name
                          description='Set the same wedge to all the particles in the particle list.',
                          authors='Yuxiang Chen',
                          options=[ScriptOption(['-p'], 'Particle list', 'string', 'required'),
                                   ScriptOption(['-w'], 'Wedge Angle (degree)', 'int', 'required'),
                                   ScriptOption(['-o'], 'Output particle list', 'string', 'optional')])

    pl_name, wedge_angle, output = parse_script_options(sys.argv[1:], helper)
    
    pl = ParticleList()
    pl.fromXMLFile(pl_name)
    w = SingleTiltWedge(wedge_angle)
    pl.setWedgeAllParticles(w)
    
    if output is None:
        pl.toXMLFile(pl_name)
    else:
        pl.toXMLFile(output)
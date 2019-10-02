#!/usr/bin/env pytom

'''
Created on Feb 6, 2014

@author: yuxiangchen
'''

if __name__ == '__main__':
    # parse command line arguments
    import sys
    from pytom.tools.script_helper import ScriptHelper, ScriptOption
    from pytom.tools.parse_script_options import parse_script_options
    from pytom.basic.structures import Particle, ParticleList, SingleTiltWedge
    import os
    
    helper = ScriptHelper(sys.argv[0].split('/')[-1], # script name
                          description='Extract certain classes from the particle list.',
                          authors='Yuxiang Chen',
                          options=[ScriptOption(['-p'], 'Particle name prefix', 'has arguments', 'required'),
                                   ScriptOption(['-c'], 'Class labels to be extracted. If multiple, please separate with comma.', 'has arguments', 'required'),
                                   ScriptOption(['-o'], 'Output particle list', 'has arguments', 'required')])

    try:
        pl_name, class_names, output = parse_script_options(sys.argv[1:], helper)
    except:
        sys.exit()
    
    from pytom.basic.structures import ParticleList
    pl = ParticleList()
    pl.fromXMLFile(pl_name)
    
    pls = pl.splitByClass()
    class_labels = class_names.split(',')
    
    res = ParticleList()
    for pp in pls:
        if pp[0].getClass() in class_labels:
            res += pp
    
    res.toXMLFile(output)
    
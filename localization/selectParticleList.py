#!/usr/bin/env python

'''
Created on Mar 20, 2012

@author: yuxiangchen
'''

if __name__ == '__main__':
    # parse command line arguments
    import sys
    from pytom.tools.script_helper import ScriptHelper, ScriptOption
    from pytom.tools.parse_script_options import parse_script_options
    helper = ScriptHelper(sys.argv[0].split('/')[-1], # script name
                          description='Select the particle list for alignment according to the class list.',
                          authors='Yuxiang Chen',
                          options=[ScriptOption('-p', 'Particle list.', 'string', 'required'),
                                   ScriptOption('-c', 'Class label file.', 'string', 'required'),
                                   ScriptOption('-o', 'Output particle list.', 'string', 'required'),
                                   ScriptOption('-t', 'True positive class.', 'int', 'optional', 0)])

    pl_filename, class_label_filename, output, tp_label = parse_script_options(sys.argv[1:], helper)
    
    from pytom.basic.structures import ParticleList
    pl = ParticleList(".")
    pl.fromXMLFile(pl_filename)
    
    pl.setClassFromLocalizationClassList(class_label_filename)
    new_pl = pl.particlesFromClass(tp_label)
    new_pl.toXMLFile(output)
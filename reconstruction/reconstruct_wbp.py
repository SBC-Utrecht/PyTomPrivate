#!/usr/bin/env python

'''
Created on Feb 28, 2013

@author: yuxiangchen
'''

if __name__ == '__main__':
    # parse command line arguments
    import sys
    from pytom.tools.script_helper import ScriptHelper, ScriptOption
    from pytom.tools.parse_script_options import parse_script_options
    helper = ScriptHelper(sys.argv[0].split('/')[-1], # script name
                          description='Reconstruct a volume from its projections using weighted back projection.',
                          authors='Yuxiang Chen',
                          options=[ScriptOption('-p', 'Projection prefix.', 'has arguments', 'required'),
                                   ScriptOption('-a', 'Projection suffix.', 'has arguments', 'required'),
                                   ScriptOption('-w', 'Apply weighting to the projections or not', 'no arguments', 'optional'),
                                   ScriptOption('-f', 'Projection start index (int).', 'has arguments', 'required'),
                                   ScriptOption('-t', 'Projection end index (int).', 'has arguments', 'required'),
                                   ScriptOption('-s', 'Volume size.', 'has arguments', 'required'),
                                   ScriptOption('-o', 'Output filename.', 'has arguments', 'required')])

    try:
        prefix, suffix, weighting, start_idx, end_idx, vol_size, output = parse_script_options(sys.argv[1:], helper)
    except Exception as e:
        print e
        sys.exit()
    
    # parse the argument
    start_idx = int(start_idx)
    end_idx = int(end_idx)
    vol_size = [int(i) for i in vol_size.split(",")]
    
    from pytom.reconstruction.reconstructionStructures import Projection, ProjectionList
    projections = ProjectionList()
    for i in xrange(start_idx, end_idx+1):
        p = Projection(prefix+str(i)+suffix)
        projections.append(p)
    
    vol = projections.reconstructVolume( dims=vol_size, reconstructionPosition=[0,0,0], binning=1, applyWeighting=weighting)
    vol.write(output)
    
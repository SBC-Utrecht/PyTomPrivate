#!/usr/bin/env pytom

'''
Created on Apr 4, 2014

@author: yuxiangchen
'''

if __name__ == '__main__':
    # parse command line arguments
    import sys
    from pytom.tools.script_helper import ScriptHelper, ScriptOption
    from pytom.tools.parse_script_options import parse_script_options
    helper = ScriptHelper(sys.argv[0].split('/')[-1], # script name
                          description='Mirror the volume.',
                          authors='Yuxiang Chen',
                          options=[ScriptOption(['-v'], 'Input volume.', 'string', 'required'),
                                   ScriptOption(['-o'], 'Output volume, defaults to the input name appended with "_mirror.em".', 'string', 'optional')])

    input_filename, output_filename = parse_script_options(sys.argv[1:], helper)
    
    # process the arguments
    if output_filename is None:
        output_filename = input_filename.split('.')[0] + '_mirror.em'

    from pytom_volume import read, vol, mirrorVolume
    v = read(input_filename)
    res = vol(v)
    mirrorVolume(v, res)

    res.write(output_filename)

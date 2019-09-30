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
                          options=[ScriptOption(['-v'], 'Input volume.', 'has arguments', 'required'),
                                   ScriptOption(['-o'], 'Output volume.', 'has arguments', 'optional'),
                                   ScriptOption(['-h', '--help'], 'Help.', 'no arguments', 'optional')])
    if len(sys.argv) == 1:
        print helper
        sys.exit()
    try:
        input_filename, output_filename, help = parse_script_options(sys.argv[1:], helper)
    except Exception as e:
        print e
        sys.exit()
    if help is True:
        print helper
        sys.exit()
    
    # process the arguments
    if output_filename is None:
        output_filename = input_filename.split('.')[0] + '_mirror.em'

    from pytom_volume import read, vol, mirrorVolume
    v = read(input_filename)
    res = vol(v)
    mirrorVolume(v, res)

    res.write(output_filename)

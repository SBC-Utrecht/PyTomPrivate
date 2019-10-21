#!/usr/bin/env pytom
'''
Created on Jul 21, 2011

@author: hrabe
'''





if __name__ == '__main__':
    # parse command line arguments
    import sys
    from pytom.tools.script_helper import ScriptHelper, ScriptOption
    from pytom.tools.parse_script_options import parse_script_options
    from pytom.basic.filter import bandpassFilter
    from pytom_volume import read
    
    helper = ScriptHelper(sys.argv[0].split('/')[-1], # script name
                          description='Will filter (band/low/high pass) a file. Choose values between 0 and 1.',
                          authors='Thomas Hrabe',
                          options=[ScriptOption(['-f', '--file'], 'Filename', 'string', 'required'),
                                   ScriptOption(['-t', '--target'], 'Name of filtered file', 'string', 'required'),
                                   ScriptOption(['-l', '--lowestFrequency'], 'The lowest frequency. 0 if not set (in pixels)', 'uint', 'optional', 0),
                                   ScriptOption(['-h', '--highestFrequency'], 'The highest frequency. Volume size / 2 if not set (in pixels)', 'uint', 'optional'),
                                   ScriptOption(['-s', '--smooth'], 'Smoothing of bandpass (in voxels).', 'uint', 'optional', 0)])

    filename, target, lowestFrequency, highestFrequency, smooth = parse_script_options(sys.argv[1:], helper)
    
    if not filename or not target:
        print helper
        sys.exit()

    v = read(filename)
    
    if highestFrequency:
        highestFrequency = int(highestFrequency)
    else:
        highestFrequency = v.sizeX()/2

    r = bandpassFilter(v, str(lowestFrequency) + ':' + str(highestFrequency)+';', None, smooth)
    
    r[0].write(target)

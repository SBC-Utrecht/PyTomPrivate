#!/usr/bin/env pytom
'''
Created on Nov 8, 2013

@author: thrabe
'''

if __name__ == '__main__':
 # parse command line arguments
    import sys
    from pytom.tools.script_helper import ScriptHelper, ScriptOption
    from pytom.tools.parse_script_options import parse_script_options 
    from pytom.basic.files import mmCIF2em
    
    helper = ScriptHelper(sys.argv[0].split('/')[-1], # script name
                          description='Compile a electron density from PDB file\n\
                          http://pytom.org/doc/pytom/files.html',
                          authors='Thomas Hrabe',
                          options=[ScriptOption(['-m', '--mmCIFFile'], 'A mmCIF file', 'string', 'required'),
                                   ScriptOption(['-c', '--chain'], 'A Chain', 'has arguments', 'optional'),
                                   ScriptOption(['-s', '--pixelSize'], 'Pixel size of output volume (in Angstrom)', 'float', 'required'),
                                   ScriptOption(['-v', '--volumeSize'], 'Volume length (size) in all dimensions', 'uint', 'required'),
                                   ScriptOption(['-o', '--outputVolumePath'], 'Path to output volume ', 'string', 'required'),
                                   ScriptOption(['-i', '--invertDensity'],'Set if density should be negative', 'no arguments', 'optional', False)])

    mmCIF, chain, pixelSize, cubeSize, volumePath ,densityNegative = parse_script_options(sys.argv[1:], helper)
    
    volume = mmCIF2em(mmCIF, pixelSize, cubeSize, chain = chain,densityNegative = densityNegative)
    
    volume.write(volumePath)
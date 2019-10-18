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
    
    helper = ScriptHelper(sys.argv[0].split('/')[-1], # script name
                          description='Convert ccp4 file to em.',
                          authors='Thomas Hrabe',
                          options=[ScriptOption(['-f', '--file'], 'Filename', 'string', 'optional'),
                                   ScriptOption(['-d', '--directory'], 'A directory of files.', 'string', 'optional'),
                                   ScriptOption(['-t', '--targetPath'], 'Path to new file.', 'string', 'required')])

    filename, directory, target = parse_script_options(sys.argv[1:], helper)
    
    if filename:
        #convert only one file
        ccp42em(filename,target)
    elif directory:
        import os
        
        fileList = os.listdir(directory)
        for file in fileList:
            if file[len(file)-3:len(file)] == '.ccp4':
                print directory + os.sep + file , target
                ccp42em(directory + os.sep + file,target)


def ccp42em(filename, target):
    from pytom_volume import read
    from pytom.tools.files import checkFileExists, checkDirExists
    import os

    if not checkFileExists(filename):
        raise RuntimeError('CCP4 file not found! ', filename)

    if not checkDirExists(target):
        raise RuntimeError('Destination directory not found! ', target)

    emfile = read(filename)

    splitName = filename.split(os.sep)
    filename = splitName[len(splitName) - 1]

    newFilename = target + os.sep + filename[0:len(filename) - 3] + '.em'

    emfile.write(newFilename, 'em')

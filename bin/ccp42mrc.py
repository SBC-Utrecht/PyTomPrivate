#!/usr/bin/env pytom

'''
Created on Jul 21, 2011

@author: hrabe
'''


from pytom.basic.files import cpp42mrc


if __name__ == '__main__':
    # parse command line arguments
    import sys
    from pytom.tools.script_helper import ScriptHelper, ScriptOption
    from pytom.tools.parse_script_options import parse_script_options
    
    helper = ScriptHelper(sys.argv[0].split('/')[-1], # script name
                          description='Convert ccp4 file to mrc.',
                          authors='Thomas Hrabe',
                          options=[ScriptOption(['-f', '--file'], 'Filename', 'string', 'optional'),
                                   ScriptOption(['-d', '--directory'], 'A directory of files.', 'string', 'optional'),
                                   ScriptOption(['-t', '--targetPath'], 'Path to new file.', 'string', 'required')])

    filename, directory, target = parse_script_options(sys.argv[1:], helper)
    
    if filename:
        #convert only one file
        ccp42mrc(filename, target)
    elif directory:
        import os
        
        fileList = os.listdir(directory)
        for file in fileList:
            if file[len(file)-3:len(file)] == '.ccp4':
                print directory + os.sep + file , target
                ccp42mrc(directory + os.sep + file,target)

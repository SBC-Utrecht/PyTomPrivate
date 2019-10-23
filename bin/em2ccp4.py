#!/usr/bin/env pytom

'''
Created on Jul 21, 2011

@author: hrabe
'''


from pytom.basic.files import em2cpp4


if __name__ == '__main__':
    # parse command line arguments
    import sys
    from pytom.tools.script_helper import ScriptHelper, ScriptOption
    from pytom.tools.parse_script_options import parse_script_options
    
    helper = ScriptHelper(sys.argv[0].split('/')[-1], # script name
                          description='Convert em file to ccp4.',
                          authors='Thomas Hrabe',
                          options=[ScriptOption(['-f', '--file'], 'Filename', 'string', 'optional'),
                                   ScriptOption(['-d', '--directory'], 'A directory of files.', 'string', 'optional'),
                                   ScriptOption(['-t', '--targetPath'], 'Path to new file.', 'string', 'required')])

    filename, directory, target, help = parse_script_options(sys.argv[1:], helper)
    
    if filename:
        #convert only one file
        em2ccp4(filename,target)
    elif directory:
        import os
        
        fileList = os.listdir(directory)
        for file in fileList:
            if file[len(file)-3:len(file)] == '.em':
                print directory + os.sep + file, target
                em2ccp4(directory + os.sep + file, target)
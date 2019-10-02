'''
Created on Jun 29, 2011

@author: yuxiangchen
'''

def parse_script_options(args, helper):
    """
    parse script options
    @param args: array of input arguments (typically from command line)
    @type args: L{list}
    @param helper: script helper
    @type helper: L{pytom.tools.script_helper.ScriptHelper}
    @return: result of parsing
    @rtype: L{list}
    """
    import getopt

    exception = "" # To buffer all exceptions to present them in one batch to the user so all errors can be resolved in one go
    
    res = [] # the result of parsing
    opt_str = ""
    long_opt = []
    for opt in helper.options:
        res.append(None)
        for name in opt.option_str:
            if name[:2] == "--":
                if opt.arg:
                    long_opt.append(name[2:]+'=')
                else:
                    long_opt.append(name[2:])
            elif name[:1] == "-":
                if opt.arg:
                    opt_str += name[1]+":"
                else:
                    opt_str += name[1]
            else:
                exception += "Option format invalid: {:s}\n".format(name)
    
    try:
        # getopt raises exceptions when an option is passed that needs an arguments but has none
        opts, args = getopt.getopt(args, opt_str, long_opt)

        if len(opts) == 0:
            import sys
            print helper
            sys.exit()

        for o, a in opts:

            # Shortcut because otherwise the code will complain when no required options are passed
            if o == "--help":
                import sys
                print helper
                sys.exit()

            for i in range(len(helper.options)):
                if o in helper.options[i].option_str:
                    if helper.options[i].arg:
                        if a is None:
                            exception += "This option ({:s}) needs an argument, none given\n".format(o)
                        res[i] = a
                    else:
                        res[i] = True
                    break

        for i, _ in enumerate(res):
            if res[i] is None and helper.options[i].required:
                exception += "Required argument not passed, use any of the following: " + " ".join(
                    helper.options[i].option_str) + "\n"

    except Exception as e:
        exception += str(e) + "\n"

    if exception != "":
        raise Exception("Some exception(s) while parsing the command line arguments, see below for details:\n" + exception + "\nUse '--help' to get information on how this script can be called.")

    return res[:-1] # Because the last option is by definition the help option
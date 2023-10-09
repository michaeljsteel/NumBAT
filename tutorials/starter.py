# startup arg parser for all literature examples

def read_args(enum, argv, sub='', refine=5):

    prefix = 'tut_{0:02d}'.format(enum)
    prefix +=sub 

    if len(argv)>1 and argv[1]=='fast=1':  # choose between faster or more accurate calculation
        refine_fac=1 
        #prefix = 'f'+prefix
        s_fmode = ' - fast mode'
    else: 
        s_fmode = ''
        refine_fac=refine 

    print('\n\nCommencing NumBAT tutorial %d%s%s'%(enum,sub, s_fmode))
    return prefix, refine_fac


# startup arg parser for all literature examples

def read_args(enum, argv, sub='', dflt_refine=5):

    prefix_str = f'josab_{enum:02d}'
    prefix_str+=sub

    if len(argv)>1 and argv[1]=='fast=1':  # choose between faster or more accurate calculation
        refine_fac=1
        #prefix_str = 'f'+prefix_str
        s_fmode = ' - fast mode'
    else:
        s_fmode = ''
        refine_fac=dflt_refine

    print('\n\nCommencing NumBAT JOSA B example %d%s%s'%(enum,sub, s_fmode))
    return prefix_str, refine_fac


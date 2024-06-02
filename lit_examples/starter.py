# startup arg parser for all literature examples

def read_args(enum, argv, sub=''):

    prefix_str = f'lit_{enum:02d}'
    prefix_str+=sub 

    if len(argv)>1 and argv[1]=='fast=1':  # choose between faster or more accurate calculation
        #prefix_str = 'f'+prefix_str
        s_fmode = ' - fast mode'
        refine_fac=1 
    else: 
        s_fmode = ''
        refine_fac=5 

    print('\n\nCommencing NumBAT literature example %d%s%s'%(enum,sub, s_fmode))
    return prefix_str, refine_fac


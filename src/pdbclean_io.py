import os
import datetime
#
def check_project(projdir=None, level='top', verbose=True):
    """
    check_project
    """
    if projdir is None:
        print("Please provide a project directory path")
    else:
        if(level=='top'):
            create_dir(projdir, verbose=verbose)
        else:
            create_dir(projdir+'/'+level, verbose=verbose)

def create_dir(dirpath, verbose=True):
    """
    """
    if not os.path.exists(dirpath):
        os.mkdir(dirpath)
        if verbose:
            now=datetime.datetime.now()
            f = open(dirpath+'/info.txt', 'w')
            f.write('directory created on {0}'.format(now))
            f.close()
    else:
        if verbose:
            print('{0} already exists, with content:'.format(dirpath))
            print(os.listdir(dirpath))

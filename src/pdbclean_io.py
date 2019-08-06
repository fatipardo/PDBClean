import sys, os, glob, shutil, datetime
#

def check_project(projdir=None, level='top', action='create', verbose=True):
    """
    check_project
    """
    if projdir is None:
        print("Please provide a project directory path")
    else:
        dirname=projdir
        if(level!='top'):
            dirname=dirname+'/'+level
        if(action=='create'):
            create_dir(dirname, verbose=verbose)
        elif(action=='clean'):
            clean_dir(dirname, verbose=verbose)
        elif(action=='delete'):
            delete_dir(dirname, verbose=verbose)

#################
# FILE HANDLING #
#################

def new_filepath(input_file, output_dir, force_new=True):
    """
    """
    filename = os.path.basename(input_file)
    output_file = output_dir+'/'+filename
    if force_new:
        if os.path.exists(output_file):
            os.remove(output_file)
    return output_file

def list_files_in_dir(path, ext):
    """
    list_files_in_dir
    """
    filelist = ()
    if os.path.exists(path):
        filelist = glob.glob(path+'/*'+ext)
    return filelist

######################
# DIRECTORY HANDLING #
######################

def define_dirs(project_dir=None, source=None, target=None):
    """
    define_dirs

    Input  : project_dir, source and target
    Output : source_dir = project_dir+/+source
             target_dir = project_dir+/+target
    """
    if project_dir is None:
        source_dir = None
        target_dir = None
    else:
        if source is None:
            source_dir = None
        else:
            source_dir = project_dir+'/'+source
        if target is None:
            target_dir = None
        else:
            target_dir = project_dir+'/'+target
    return source_dir, target_dir

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
            
def clean_dir(dirpath, verbose=True):
    """
    """
    if os.path.exists(dirpath):
        listfile = (file for file in os.listdir(dirpath) if os.path.isfile(os.path.join(dirpath, file)))
        if verbose:
            print('Cleaning {0}...'.format(dirpath))
        for f in listfile:
            os.remove(dirpath+'/'+f)

def delete_dir(dirpath, verbose=True):
    """
    """
    if os.path.exists(dirpath):
        shutil.rmtree(dirpath)
        if verbose:
            print('Deleting {0}...'.format(dirpath))

##########################
# COLABORATORY INTERFACE #
##########################

def define_rundir(local_relative_path, colab=False):
    """
    """
    if not colab:
        RUNDIR=local_relative_path
    else:
        RUNDIR='./'
    return RUNDIR

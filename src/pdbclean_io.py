import sys
import os
import glob
import shutil
import datetime
import pickle
import numpy as np
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

def new_filepath(input_file, output_dir, subdir='.', force_new=True):
    """
    """
    filename = os.path.basename(input_file)
    output_file = output_dir+'/'+subdir+'/'+filename
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
        log(logdir=dirpath, action='create')
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

def log(logdir='.', action='create', field=None, fieldname=None):
    """
    :param logdir:
    :param action:
    :param field:
    :return:
    """
    logfile=logdir+'/log.pkl'
    if(action=='create'):
        now = datetime.datetime.now()
        dictionary = {}
        dictionary['directory_location'] = logdir
        dictionary['time_created'] = now
        with open(logfile, 'wb') as fp:
            pickle.dump(dictionary, fp)
    elif(action=='add'):
        if field is None or fieldname is None:
            print('LOG WARNING: please provide a field and fieldname to add to the dictionary')
        else:
            with open(logfile, 'rb') as fp:
                dictionary = pickle.load(fp)
            dictionary[fieldname] = field
            with open(logfile, 'wb') as fp:
                pickle.dump(dictionary, fp)
    elif(action=='inspect'):
        print('The following fields have been filled in {0}: '.format(logfile))
        with open(logfile, 'rb') as fp:
            dictionary = pickle.load(fp)
        print(dictionary.keys())

def load_field_from_dictionary(logdir='.', fieldname=None):
    """
    :param logdir:
    :param fieldname:
    :return:
    """
    logfile = logdir + '/log.pkl'
    if fieldname is None:
        print('WARNING: please provide a field name to load from {0}'.format(logfile))
    else:
        with open(logfile, 'rb') as fp:
            dictionary = pickle.load(fp)
        field = dictionary[fieldname]
        stripped_field = np.copy(field)
        for i in np.arange(field.shape[0]):
            stripped_field[i] = field[i].strip('\n')
        return stripped_field

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

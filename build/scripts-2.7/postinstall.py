#!/usr/bin/python
# -*- coding: utf-8 -*-
from __future__ import absolute_import
import os
import sys
import shutil
import offspec_databrowser as mypackage

DESKTOP_FOLDER = get_special_folder_path("CSIDL_DESKTOPDIRECTORY")
NAME = 'offspec-databrowser.lnk'

if sys.argv[1] == '-install':
    create_shortcut(
        os.path.join(sys.prefix, 'pythonw.exe'), # program
        mypackage.__file__, # description
        os.path.join(DESKTOP_FOLDER, NAME),
        #os.path.join(DESKTOP_FOLDER, NAME), # filename
        '-moffspec_databrowser.databrowser', #parameters
        #mypackage.__file__, # parameters
        '', # workdir
        os.path.join(os.path.dirname(mypackage.__file__), 'favicon.ico'), # iconpath
    )
    ## move shortcut from current directory to DESKTOP_FOLDER
    #shutil.move(os.path.join(DESKTOP_FOLDER, NAME),
    #            os.path.join(DESKTOP_FOLDER, NAME))
    # tell windows installer that we created another 
    # file which should be deleted on uninstallation
    file_created(os.path.join(DESKTOP_FOLDER, NAME))

if sys.argv[1] == '-remove':
    pass
    # This will be run on uninstallation. Nothing to do.

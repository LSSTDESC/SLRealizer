#!/usr/bin/python2
# -*- coding: utf-8 -*-

import dropbox, sys, os
import requests

def dropbox_upload(path_to_files, save_name):
    """
    This method saves output to dropbox.
    path_to_files: relative path to the output file
    save_name: saving name in dropbox
    """

    access_token = 'mwiZY0MbbvAAAAAAAAAAB3qk-8Gk9zGkaypCKGNzOeYXn0iQ1S-JHHGgee3dUru6'
    dbx = dropbox.Dropbox(access_token)
    print(type(path_to_files))
    print(path_to_files)
    f = open(path_to_files)
    save_full_name = '/'+save_name # it's just dropbox syntax
    dbx.files_upload(f.read(), save_full_name, mode=dropbox.files.WriteMode.overwrite)

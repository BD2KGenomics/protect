#!/usr/bin/env python2.7
# Copyright (C) 2016 UCSC Computational Genomics Lab
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""
Author : Arjun Arkal Rao
Affiliation : UCSC BME, UCSC Genomics Institute
File : encrypt_files_in_dir_to_s3.py

SOURCE: https://github.com/jvivian/one_off_scripts/blob/master/
        encrypt_files_in_dir_to_s3.py
ORIGINAL AUTHOR: John Vivian

Move files in a directory, or entire directory structures to S3 with (or without) encryption.
"""

import argparse
import base64
import hashlib
import os
import subprocess
import sys
import re

from boto.s3.connection import S3Connection


class InputParameterError(Exception):
    """
    This Error Class will be raised  in the case of a bad parameter provided.
    """
    __module__ = Exception.__module__


def generate_unique_key(master_key, url):
    """
    This module will take a master key and a url, and then make a new key specific to the url, based
    off the master.
    :param str master_key: Path to the master key used for encryption.
    :param str url: Full URL to the potential file location in S3.
    :returns new_key: The new key that is obtained by hashing the master key:url combination.
    """
    with open(master_key, 'r') as keyfile:
        master_key = keyfile.read()
    assert len(master_key) == 32, 'Invalid Key! Must be 32 characters.  Key: %s' % master_key + \
                                  ', Length: %s' % len(master_key)
    new_key = hashlib.sha256(master_key + url).digest()
    assert len(new_key) == 32, 'New key is invalid and is not ' + \
        '32 characters: {}'.format(new_key)
    return new_key


class BucketInfo(object):
    """
    This class contains all the functions relevant to this script for working with a given S3
    bucket.
    """
    def __init__(self, bucket_name):
        self.bucket_name = bucket_name
        # Set up the https url base. A typical s3 endpoint url base would look like
        #                       https://s3-ENDPOINT.amazonaws.com/
        # However the https endpoints for us-east-1 are slightly different and use
        #                       https://s3.amazonaws.com/
        # REF: http://docs.aws.amazon.com/general/latest/gr/rande.html#s3_region
        endpoint = self._get_bucket_endpoint()
        endpoint = '-' + endpoint if endpoint else ''
        self._https_url_base = 'https://s3' + endpoint + '.amazonaws.com/' + self.bucket_name
        # Set up the s3 url base
        self._s3_url_base = 'S3://' + self.bucket_name

    def _get_bucket_endpoint(self):
        """
        Queries S3 to identify the region hosting the provided bucket.
        """
        conn = S3Connection()
        bucket = conn.lookup(self.bucket_name)
        if not bucket:
            # TODO: Make the bucket here?
            raise InputParameterError('The provided bucket %s doesn\'t exist' % self.bucket_name)
        endpoint = str(bucket.get_location())
        return endpoint

    def object_https_url(self, key):
        """
        Returns the full https url for key given this bucket.
        :param key: the remote filename
        :return: Full https url to the file in S3
        """
        return os.path.join(self._https_url_base, key)

    def object_s3_url(self, key):
        """
        Returns the full s3 url for key given this bucket.
        :param key: the remote filename
        :return: Full https url to the file in S3
        """
        return os.path.join(self._s3_url_base, key)


def write_to_s3(datum, master_key, bucket_name, remote_dir=''):
    """
    This module will take in some datum (a file, or a folder) and write it to
    S3.  It requires a master key to encrypt the datum with, and a bucket to
    drop the results into.  If remote dir is set, the datum is dropped into the
    provided directory.
    :param str datum: File or folder that needs to be transferred to S3
    :param str master_key: Path to the master key used for encryption.
    :param str bucket_name: AWS bucket to store the remote data
    :param str remote_dir: An optional parameter describing a remote pseudo directory in the bucket
                           where the data will be stored.
    """
    # Instantiate the bucket info class to set up the https and s3 url bases.
    bucket_info = BucketInfo(bucket_name)

    # Retain the base dir separately from the file name / folder structure of DATUM.  This way it
    # can be easily joined into an AWS filename
    folder_base_dir = os.path.split(datum)[0]
    # Ensure files are either "regular files" or folders
    if os.path.isfile(datum):
        files = [os.path.basename(datum)]
    elif os.path.isdir(datum):
        files = ['/'.join([re.sub(folder_base_dir, '', folder), filename]).lstrip('/')
                 for folder, _, files in os.walk(datum) for filename in files]
    else:
        raise RuntimeError(datum + 'was neither regular file nor folder.')

    # Write each file to S3
    for file_path in files:
        # key, here, refers to the key or token used to access the file once it's in S3.
        # THIS IS NOT RELATED TO THE ENCRYPTION KEY.
        key = os.path.join(remote_dir, file_path)
        # base command call
        command = ['s3am', 'upload']
        if master_key:
            new_key = generate_unique_key(master_key, bucket_info.object_https_url(key))
            #  Add base64 encoded key
            command.extend(['--sse-key-base64', base64.b64encode(new_key)])
        # Add source path info to the call
        command.extend(['file://' + os.path.join(folder_base_dir, file_path)])
        # Add destination to the call
        command.append(bucket_info.object_s3_url(key))
        subprocess.call(command)
    return None


def main():
    """
    This is the main module for the script.  The script will accept a file, or a directory, and then
    encrypt it with a provided key before pushing it to S3 into a specified bucket.
    """
    parser = argparse.ArgumentParser(description=main.__doc__, add_help=True)
    parser.add_argument('-M', '--master_key', dest='master_key', help='Path  to the master key ' +
                        'used for the encryption.  Data is transferred without encryption if this' +
                        'is not provided.', type=str, required=False, default=None)
    parser.add_argument('-B', '--bucket', dest='bucket', help='S3 bucket.', type=str, required=True)
    parser.add_argument('-R', '--remote_dir', dest='remote_dir', help='Pseudo directory within ' +
                        'the bucket to store the file(s).  NOTE: Folder structure below ' +
                        'REMOTE_DIR will be retained.', type=str, required=False, default='')
    parser.add_argument('data', help='File(s) or folder(s) to transfer to S3.', type=str, nargs='+')
    params = parser.parse_args()
    # Input handling
    if params.master_key and not os.path.exists(params.master_key):
        raise InputParameterError('The master key was not found at ' +
                                  params.master_key)
    #  If the user doesn't have ~/.boto , it doesn't even make sense to go ahead
    if not os.path.exists(os.path.expanduser('~/.boto')):
        raise RuntimeError('~/.boto not found')
    # Ensure that the remote directory doesn't start with a /
    if params.remote_dir.startswith('/'):
        raise InputParameterError('The remote dir cannot start with a \'/\'')

    # Process each of the input arguments.
    for datum in params.data:
        datum = os.path.abspath(datum)
        if not os.path.exists(datum):
            print('ERROR: %s could not be found.' % datum, file=sys.stderr)
            continue
        write_to_s3(datum, params.master_key, params.bucket, params.remote_dir)
    return None


if __name__ == '__main__':
    main()

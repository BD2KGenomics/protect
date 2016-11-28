#!/usr/bin/env python2.7
# Copyright 2016 Arjun Arkal Rao
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
"""
Author : Arjun Arkal Rao
Affiliation : UCSC BME, UCSC Genomics Institute
File : protect/ProTECT.py

Program info can be found in the docstring of the main function.
Details can also be obtained by running the script with -h .
"""
from __future__ import print_function
from urlparse import urlparse

import errno
import gzip
import os
import re
import shutil
import subprocess
import sys
import tarfile
import time


def get_files_from_filestore(job, files, work_dir, docker=False):
    """
    This is adapted from John Vivian's return_input_paths from the RNA-Seq pipeline.

    Returns the paths of files from the FileStore if they are not present.
    If docker=True, return the docker path for the file.
    If the file extension is tar.gz, then tar -zxvf it.

    files is a dict with:
        keys = the name of the file to be returned in toil space
        value = the input value for the file (can be toil temp file)
    work_dir is the location where the file should be stored
    cache indiciates whether caching should be used
    """
    for name in files.keys():
        outfile = job.fileStore.readGlobalFile(files[name], '/'.join([work_dir, name]))
        # If the files will be sent to docker, we will mount work_dir to the container as /data and
        # we want the /data prefixed path to the file
        if docker:
            files[name] = docker_path(outfile)
        else:
            files[name] = outfile
    return files


def docker_path(filepath):
    """
    Given a path, returns that files path inside the docker mount directory
    (/data).
    """
    return os.path.join('/data', os.path.basename(filepath))


def docker_call(tool, tool_parameters, work_dir, java_opts=None, outfile=None,
                dockerhub='aarjunrao', interactive=False):
    """
    Makes subprocess call of a command to a docker container. work_dir MUST BE AN ABSOLUTE PATH or
    the call will fail.  outfile is an open file descriptor to a writeable file.
    """
    # If an outifle has been provided, then ensure that it is of type file, it is writeable, and
    # that it is open.
    if outfile:
        assert isinstance(outfile, file), 'outfile was not passsed a file'
        assert outfile.mode in ['w', 'a', 'wb', 'ab'], 'outfile not writeable'
        assert not outfile.closed, 'outfile is closed'
    # If the call is interactive, set intereactive to -i
    if interactive:
        interactive = '-i'
    else:
        interactive = ''
    # If a tag is passed along with the image, use it.
    if ':' in tool:
        docker_tool = '/'.join([dockerhub, tool])
    # Else use 'latest'
    else:
        docker_tool = ''.join([dockerhub, '/', tool, ':latest'])
    # Get the docker image on the worker if needed
    call = ['docker', 'images']
    dimg_rv = subprocess.check_output(call)
    existing_images = [':'.join(x.split()[0:2]) for x in dimg_rv.splitlines()
                       if x.startswith(dockerhub)]
    if docker_tool not in existing_images:
        try:
            call = ' '.join(['docker', 'pull', docker_tool]).split()
            subprocess.check_call(call)
        except subprocess.CalledProcessError as err:
            raise RuntimeError('docker command returned a non-zero exit status ' +
                               '(%s)' % err.returncode + 'for command \"%s\"' % ' '.join(call),)
        except OSError:
            raise RuntimeError('docker not found on system. Install on all' +
                               ' nodes.')
    # If java options have been provided, it needs to be in the docker call
    if java_opts:
        base_docker_call = ' docker run -e JAVA_OPTS=-Xmx{} '.format(java_opts) + '--rm=true ' + \
            '-v {}:/data --log-driver=none '.format(work_dir) + interactive
    else:
        base_docker_call = ' docker run --rm=true -v {}:/data '.format(work_dir) + \
            '--log-driver=none ' + interactive
    call = base_docker_call.split() + [docker_tool] + tool_parameters
    try:
        subprocess.check_call(call, stdout=outfile)
    except subprocess.CalledProcessError as err:
        raise RuntimeError('docker command returned a non-zero exit status (%s)' % err.returncode +
                           'for command \"%s\"' % ' '.join(call),)
    except OSError:
        raise RuntimeError('docker not found on system. Install on all nodes.')


def untargz(input_targz_file, untar_to_dir):
    """
    This module accepts a tar.gz archive and untars it.

    RETURN VALUE: path to the untar-ed directory/file

    NOTE: this module expects the multiple files to be in a directory before
          being tar-ed.
    """
    assert tarfile.is_tarfile(input_targz_file), 'Not a tar file.'
    tarball = tarfile.open(input_targz_file)
    return_value = os.path.join(untar_to_dir, tarball.getmembers()[0].name)
    tarball.extractall(path=untar_to_dir)
    tarball.close()
    return return_value


def gunzip(input_gzip_file, block_size=1024):
    """
    Gunzips the input file to the same directory
    :param input_gzip_file: File to be gunzipped
    :return: path to the gunzipped file
    """
    assert os.path.splitext(input_gzip_file)[1] == '.gz'
    assert is_gzipfile(input_gzip_file)
    with gzip.open(input_gzip_file) as infile:
        with open(os.path.splitext(input_gzip_file)[0], 'w') as outfile:
            while True:
                block = infile.read(block_size)
                if block == '':
                    break
                else:
                    outfile.write(block)
    return outfile.name


def is_gzipfile(filename):
    """
    This function attempts to ascertain the gzip status of a file based on the "magic signatures" of
    the file. This was taken from the stack overflow
    http://stackoverflow.com/questions/13044562/python-mechanism-to-identify-compressed-file-type\
        -and-uncompress
    """
    assert os.path.exists(filename), 'Input {} does not '.format(filename) + \
        'point to a file.'
    with open(filename, 'rb') as in_f:
        start_of_file = in_f.read(3)
        if start_of_file == '\x1f\x8b\x08':
            return True
        else:
            return False


def get_file_from_s3(job, s3_url, encryption_key=None, per_file_encryption=True,
                     write_to_jobstore=True):
    """
    Downloads a supplied URL that points to an unencrypted, unprotected file on Amazon S3. The file
    is downloaded and a subsequently written to the jobstore and the return value is a the path to
    the file in the jobstore.

    :param str s3_url: URL for the file (can be s3 or https)
    :param str encryption_key: Path to the master key
    :param bool per_file_encryption: If encrypted, was the file encrypted using the per-file method?
    :param bool write_to_jobstore: Should the file be written to the job store?
    """
    work_dir = job.fileStore.getLocalTempDir()

    parsed_url = urlparse(s3_url)
    if parsed_url.scheme == 'https':
        download_url = 'S3:/' + parsed_url.path  # path contains the second /
    elif parsed_url.scheme == 's3':
        download_url = s3_url
    else:
        raise RuntimeError('Unexpected url scheme: %s' % s3_url)

    filename = '/'.join([work_dir, os.path.basename(s3_url)])
    # This is common to encrypted and unencrypted downloads
    download_call = ['s3am', 'download', '--download-exists', 'resume']
    # If an encryption key was provided, use it.
    if encryption_key:
        download_call.extend(['--sse-key-file', encryption_key])
        if per_file_encryption:
            download_call.append('--sse-key-is-master')
    # This is also common to both types of downloads
    download_call.extend([download_url, filename])
    attempt = 0
    while True:
        try:
            with open(work_dir + '/stderr', 'w') as stderr_file:
                subprocess.check_call(download_call, stderr=stderr_file)
        except subprocess.CalledProcessError:
            # The last line of the stderr will have the error
            exception = ''
            with open(stderr_file.name) as stderr_file:
                for line in stderr_file:
                    line = line.strip()
                    if line:
                        exception = line
            if exception.startswith('boto'):
                exception = exception.split(': ')
                if exception[-1].startswith('403'):
                    raise RuntimeError('s3am failed with a "403 Forbidden" error  while obtaining '
                                       '(%s). Did you use the correct credentials?' % s3_url)
                elif exception[-1].startswith('400'):
                    raise RuntimeError('s3am failed with a "400 Bad Request" error while obtaining '
                                       '(%s). Are you trying to download an encrypted file without '
                                       'a key, or an unencrypted file with one?' % s3_url)
                else:
                    raise RuntimeError('s3am failed with (%s) while downloading (%s)' %
                                       (': '.join(exception), s3_url))
            elif exception.startswith('AttributeError'):
                exception = exception.split(': ')
                if exception[-1].startswith("'NoneType'"):
                    raise RuntimeError('Does (%s) exist on s3?' % s3_url)
                else:
                    raise RuntimeError('s3am failed with (%s) while downloading (%s)' %
                                       (': '.join(exception), s3_url))
            else:
                if attempt < 3:
                    attempt += 1
                    continue
                else:
                    raise RuntimeError('Could not diagnose the error while downloading (%s)' %
                                       s3_url)
        except OSError:
            raise RuntimeError('Failed to find "s3am". Install via "apt-get install --pre s3am"')
        else:
            break
        finally:
            os.remove(stderr_file.name)
    assert os.path.exists(filename)
    if write_to_jobstore:
        filename = job.fileStore.writeGlobalFile(filename)
    return filename


def get_file_from_cghub(job, cghub_xml, cghub_key, univ_options):
    """
    This function will download the file from cghub using the xml specified by cghub_xml

    ARGUMENTS
    1. cghub_xml: Path to an xml file for cghub.
    2. cghub_key: Credentials for a cghub download operation.
    3. write_to_jobstore: Flag indicating whether the final product should be written to jobStore.

    RETURN VALUES
    1. A path to the prefix for the fastqs that is compatible with the pipeline.
    """
    work_dir = job.fileStore.getLocalTempDir()
    # Get from S3 if required
    if cghub_xml.startswith('http'):
        assert cghub_xml.startswith('https://s3'), 'Not an S3 link'
        cghub_xml = get_file_from_s3(job, cghub_xml, encryption_key=univ_options['sse_key'],
                                     write_to_jobstore=False)
    else:
        assert os.path.exists(cghub_xml), 'Could not find file: %s' % cghub_xml
    shutil.copy(cghub_xml, os.path.join(work_dir, 'cghub.xml'))
    cghub_xml = os.path.join(work_dir, 'cghub.xml')
    assert os.path.exists(cghub_key), 'Could not find file: %s' % cghub_key
    shutil.copy(cghub_key, os.path.join(work_dir, 'cghub.key'))
    cghub_key = os.path.join(work_dir, 'cghub.key')
    temp_fastqdir = os.path.join(work_dir, 'temp_fastqdir')
    os.mkdir(temp_fastqdir)
    base_parameters = ['-d', docker_path(cghub_xml),
                       '-c', docker_path(cghub_key),
                       '-p', docker_path(temp_fastqdir)]
    attempt_number = 0
    while True:
        # timeout increases by 10 mins per try
        parameters = base_parameters + ['-k', str((attempt_number + 1) * 10)]
        try:
            docker_call('genetorrent', tool_parameters=parameters, work_dir=work_dir,
                        dockerhub=univ_options['dockerhub'])
        except RuntimeError as err:
            time.sleep(600)
            job.fileStore.logToMaster(err.message)
            attempt_number += 1
            if attempt_number == 3:
                raise
            else:
                continue
        else:
            break
    analysis_id = [x for x in os.listdir(temp_fastqdir)
                   if not (x.startswith('.') or x.endswith('.gto'))][0]
    files = [x for x in os.listdir(os.path.join(temp_fastqdir, analysis_id))
             if not x.startswith('.')]
    if len(files) == 2:
        prefixes = [os.path.splitext(x)[1] for x in files]
        if {'.bam', '.bai'} - set(prefixes):
            raise RuntimeError('This is probably not a TCGA archive for WXS or RSQ. If you are ' +
                               'sure it is, email aarao@ucsc.edu with details.')
        else:
            bamfile = os.path.join(temp_fastqdir, analysis_id,
                                   [x for x in files if x.endswith('.bam')][0])
            return bam2fastq(bamfile, univ_options)
    elif len(files) == 1:
        if not files[0].endswith('.tar.gz'):
            raise RuntimeError('This is probably not a TCGA archive for WXS or RSQ. If you are ' +
                               'sure it is, email aarao@ucsc.edu with details.')
        else:
            out_fastq_dir = os.path.join(work_dir, 'fastqs')
            os.mkdir(out_fastq_dir)
            fastq_file = untargz(os.path.join(temp_fastqdir, analysis_id, files[0]), out_fastq_dir)
            if fastq_file.endswith(('.fastq', '.fastq.gz')):
                return re.sub('_2.fastq', '_1.fastq', fastq_file)
            else:
                raise RuntimeError('This is probably not a TCGA archive for WXS or RSQ. If you ' +
                                   'are sure it is, email aarao@ucsc.edu with details.')
    else:
        raise RuntimeError('This is probably not a TCGA archive for WXS or RSQ. If you are sure ' +
                           'it is, email aarao@ucsc.edu with details.')


def bam2fastq(bamfile, univ_options):
    """
    split an input bam to paired fastqs.

    ARGUMENTS
    1. bamfile: Path to a bam file
    2. univ_options: Dict of universal arguments used by almost all tools
         univ_options
                |- 'dockerhub': <dockerhub to use>
                +- 'java_Xmx': value for max heap passed to java
    """
    work_dir = os.path.split(bamfile)[0]
    base_name = os.path.split(os.path.splitext(bamfile)[0])[1]
    parameters = ['SamToFastq',
                  ''.join(['I=', docker_path(bamfile)]),
                  ''.join(['F=/data/', base_name, '_1.fastq']),
                  ''.join(['F2=/data/', base_name, '_2.fastq']),
                  ''.join(['FU=/data/', base_name, '_UP.fastq'])]
    docker_call(tool='picard', tool_parameters=parameters, work_dir=work_dir,
                dockerhub=univ_options['dockerhub'], java_opts=univ_options['java_Xmx'])
    first_fastq = ''.join([work_dir, '/', base_name, '_1.fastq'])
    assert os.path.exists(first_fastq)
    return first_fastq


def export_results(job, fsid, file_path, univ_options, subfolder=None):
    """
    Write out a file to a given location. The location can be either a directory on the local
    machine, or a folder with a bucket on AWS.
    TODO: Azure support
    :param str fsid: The file store id for the file to be exported
    :param str file_path: The path to the file that neeeds to be exported
    :param dict univ_options: A dict of the universal options passed to this script. The important
           dict entries are ouput_folder and storage_location.
                * storage_location: 'Local' or an 'aws:<bucket_name>'.
                * output_folder: The folder to store the file. This must exist on the local machine
                                 if storage_location is 'Local'. If the storage_location is an aws
                                 bucket,  this string represents the path to the file in the bucket.
                                 To keep it in the base directory, specify 'NA'.
    :param str subfolder: A sub folder within the main folder where this data should go
    :return: None
    """
    job.fileStore.logToMaster('Exporting %s to output location' % fsid)
    try:
        assert univ_options['output_folder'], 'Need a path to a folder to write out files'
        assert univ_options['storage_location'], 'Need to know where the files need to go. ' + \
                                                 'Local or AWS/Azure, etc.'
    except AssertionError as err:
        # This isn't a game killer.  Continue the pipeline without erroring out but do inform the
        # user about it.
        print('ERROR:', err.message, file=sys.stderr)
        return
    if univ_options['output_folder'] == 'NA':
        output_folder = ''
    else:
        output_folder = univ_options['output_folder']
    output_folder = os.path.join(output_folder, univ_options['patient'])
    output_folder = os.path.join(output_folder, subfolder) if subfolder else output_folder
    if univ_options['storage_location'] == 'local':
        # Handle Local
        try:
            # Create the directory if required
            os.makedirs(output_folder, 0755)
        except OSError as err:
            if err.errno != errno.EEXIST:
                raise
        output_url = 'file://' + os.path.join(output_folder, os.path.basename(file_path))
    elif univ_options['storage_location'].startswith('aws'):
        # Handle AWS
        bucket_name = univ_options['storage_location'].split(':')[-1]

        file_name = os.path.basename(file_path)
        output_url = os.path.join('S3://', bucket_name, output_folder.strip('/'), file_name)
    # Can't do Azure or google yet.
    else:
        print("Currently doesn't support anything but Local and aws.")
        return
    job.fileStore.exportFile(fsid, output_url)


def write_to_s3(file_path, key_path, bucket_name, output_folder, overwrite=True):
    """
    Write the file to S3.

    :param file_path: The file to be written
    :param key_path: Path to the encryption Key
    :param bucket_name: The bucket where the data will be written
    :param output_folder: The location in the bucket for the output data
    :param overwrite: Should the data be overwritten if it exists in the bucket?
    """
    assert overwrite in (True, False)
    overwrite = 'overwrite' if overwrite else 'skip'
    file_name = os.path.basename(file_path)
    output_file = os.path.join('S3://', bucket_name, output_folder.strip('/'), file_name)
    subprocess.check_call(['s3am', 'upload', '--exists=' + overwrite, '--sse-key-file', key_path,
                           file_path, output_file])


def file_xext(filepath):
    """
    Get the file extension wrt compression from the filename (is it tar or targz)
    :param str filepath: Path to the file
    :return str ext: Compression extension name
    """
    ext = os.path.splitext(filepath)[1]
    if ext == '.gz':
        xext = os.path.splitext(os.path.splitext(filepath)[0])[1]
        if xext == '.tar':
            ext = xext + ext
    elif ext == '.tar':
        pass  # ext is already .tar
    else:
        ext = ''
    return ext


def strip_xext(filepath):
    """
    Strips the compression extension from the filename
    :param filepath: Path to compressed file.
    :return str filepath: Path to the file with the compression extension stripped off.
    """
    ext_size = len(file_xext(filepath).split('.')) - 1
    for i in xrange(0, ext_size):
        filepath = os.path.splitext(filepath)[0]
    return filepath


def delete_fastqs(job, fastqs):
    """
    This module will delete the fastqs from the job Store once their purpose has been achieved (i.e.
    after all mapping steps)

    ARGUMENTS
    :param dict fastqs: Dict of list of input fastqs
    """
    for key in fastqs.keys():
        if key == 'patient_id':
            continue
        job.fileStore.logToMaster('Deleting fastq files for "%s".' % key)
        for i in fastqs[key]:
            job.fileStore.deleteGlobalFile(i)
    return None


# Exception for bad parameters provided
class ParameterError(Exception):
    """
    This Error Class will be raised  in the case of a bad parameter provided.
    """
    pass

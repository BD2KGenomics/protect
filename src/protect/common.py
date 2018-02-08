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

from collections import defaultdict
from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText
from urlparse import urlparse

import errno
import gzip
import logging
import os
import re
import smtplib
import socket
import subprocess
import sys
import tarfile
import urllib2
import uuid


def get_files_from_filestore(job, files, work_dir, docker=False):
    """
    Download a dict of files to the given directory and modify the path to a docker-friendly one if
    requested.

    :param dict files: A dictionary of filenames: fsIDs
    :param str work_dir: The destination directory
    :param bool docker: Should the file path be converted to our standard docker '/data/filename'?
    :return: Dict of files: (optionallly docker-friendly) fileepaths
    :rtype: dict
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


def docker_path(filepath, work_dir=None):
    """
    Given a path, return that files path inside the docker mount directory (/data).

    :param str filepath: The path to a file
    :param str work_dir: The part of the path to replace with /data
    :return: The docker-friendly path for `filepath`
    :rtype: str
    """
    if work_dir:
        return re.sub(work_dir, '/data', filepath)

    else:
        return os.path.join('/data', os.path.basename(filepath))


def docker_call(tool, tool_parameters, work_dir, java_xmx=None, outfile=None,
                dockerhub='aarjunrao', interactive=False, tool_version='latest'):
    """
    Make a subprocess call of a command to a docker container.

    :param str tool: The tool to run
    :param list tool_parameters: Parameters passed to `tool`
    :param str work_dir: The absolute path to the working directory to be mounted into the container
    :param str java_xmx: The heap space in human readable format to provide java ('20G' will pass
           -Xmx20G to java)
    :param file outfile: The file object to dump stdout
    :param str dockerhub: The dockerhub from where the tool will be pulled
    :param bool interactive: Should the docker container be run in interactive mode?
    :param str tool_version: What dockerised tool version should be used?
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
    # Set the tool version
    docker_tool = ''.join([dockerhub, '/', tool, ':', tool_version])
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
    if java_xmx:
        base_docker_call = ' docker run -e JAVA_OPTS=-Xmx{} '.format(java_xmx) + '--rm=true ' + \
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
    Accept a tar.gz archive and untar it to the given location.  The archive can have either one
    file, or many files in a single directory.

    :param str input_targz_file: Path to a tar.gz archive
    :param str untar_to_dir: The directory where untared files will be dumped

    :return: path to the untar-ed directory/file
    :rtype: str
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
    :rtype: str
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
    Attempt to ascertain the gzip status of a file based on the "magic signatures" of the file.

    This was taken from the stack overflow post
    http://stackoverflow.com/questions/13044562/python-mechanism-to-identify-compressed-file-type\
        -and-uncompress

    :param str filename: A path to a file
    :return: True if the file appears to be gzipped else false
    :rtype: bool
    """
    assert os.path.exists(filename), 'Input {} does not '.format(filename) + \
        'point to a file.'
    with open(filename, 'rb') as in_f:
        start_of_file = in_f.read(3)
        if start_of_file == '\x1f\x8b\x08':
            return True
        else:
            return False


def get_file_from_gdc(job, gdc_url, gdc_download_token, write_to_jobstore=True):
    """
    Download a supplied "URL" that points to a file in the NCBI GDC database.  The path to the gdc
    download token must be provided.  The file is downloaded and written to the jobstore if
    requested.

    :param str gdc_url: URL for the file (in the form of gdc://<UUID>)
    :param str gdc_download_token: Path to the gdc download token
    :param bool write_to_jobstore: Should the file be written to the job store?
    :return: Path to the downloaded file or fsID (if write_to_jobstore was True)
    :rtype: list(str|toil.fileStore.FileID)
    """
    work_dir = job.fileStore.getLocalTempDir()

    parsed_url = urlparse(gdc_url)
    assert parsed_url.scheme == 'gdc', 'Unexpected url scheme: %s' % gdc_url

    file_dir = '/'.join([work_dir, parsed_url.netloc])

    # This is common to encrypted and unencrypted downloads
    currwd = os.getcwd()
    os.chdir(work_dir)
    try:
        download_call = ['gdc-client', 'download', '-t', gdc_download_token, parsed_url.netloc]
        subprocess.check_call(download_call)
    finally:
        os.chdir(currwd)

    assert os.path.exists(file_dir)
    output_files = [os.path.join(file_dir, x) for x in os.listdir(file_dir)
                    if not x.endswith('logs')]
    # NOTE: We only handle vcf and bam+bai
    if len(output_files) == 1:
        assert output_files[0].endswith('vcf')
    else:
        if not {os.path.splitext(x)[1] for x in output_files} >= {'.bam', '.bai'}:
            raise ParameterError('Can currently only handle pre-indexed GDC bams.')
        # Always [bam, bai]
        output_files = [x for x in output_files if x.endswith(('bam', 'bai'))]
        output_files = sorted(output_files, key=lambda x: os.path.splitext(x)[1], reverse=True)
    if write_to_jobstore:
        output_files = [job.fileStore.writeGlobalFile(f) for f in output_files]
    return output_files


def get_file_from_s3(job, s3_url, encryption_key=None, per_file_encryption=True,
                     write_to_jobstore=True):
    """
    Download a supplied URL that points to a file on Amazon S3.  If the file is encrypted using
    sse-c (with the user-provided key or with a hash of the usesr provided master key) then the
    encryption keys will be used when downloading.  The file is downloaded and written to the
    jobstore if requested.

    :param str s3_url: URL for the file (can be s3, S3 or https)
    :param str encryption_key: Path to the master key
    :param bool per_file_encryption: If encrypted, was the file encrypted using the per-file method?
    :param bool write_to_jobstore: Should the file be written to the job store?
    :return: Path to the downloaded file or fsID (if write_to_jobstore was True)
    :rtype: str|toil.fileStore.FileID
    """
    work_dir = job.fileStore.getLocalTempDir()

    parsed_url = urlparse(s3_url)
    if parsed_url.scheme == 'https':
        download_url = 'S3:/' + parsed_url.path  # path contains the second /
    elif parsed_url.scheme in ('s3', 'S3'):
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
    exception = ''
    while True:
        try:
            with open(work_dir + '/stderr', 'w') as stderr_file:
                subprocess.check_call(download_call, stderr=stderr_file)
        except subprocess.CalledProcessError:
            # The last line of the stderr will have the error
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


def get_file_from_url(job, any_url, encryption_key=None, per_file_encryption=True,
                      write_to_jobstore=True):
    """
    Download a supplied URL that points to a file on an http, https or ftp server.  If the file is
    found to be an https s3 link then the file is downloaded using `get_file_from_s3`. The file is
    downloaded and written to the jobstore if requested.
    Encryption arguments are for passing to `get_file_from_s3` if required.

    :param str any_url: URL for the file
    :param str encryption_key: Path to the master key
    :param bool per_file_encryption: If encrypted, was the file encrypted using the per-file method?
    :param bool write_to_jobstore: Should the file be written to the job store?
    :return: Path to the downloaded file or fsID (if write_to_jobstore was True)
    :rtype: str|toil.fileStore.FileID
    """
    work_dir = job.fileStore.getLocalTempDir()

    filename = '/'.join([work_dir, str(uuid.uuid4())])
    url = any_url
    parsed_url = urlparse(any_url)
    try:
        response = urllib2.urlopen(url)
    except urllib2.HTTPError:
        if parsed_url.netloc.startswith(('s3', 'S3')):
            job.fileStore.logToMaster("Detected https link is for an encrypted s3 file.")
            return get_file_from_s3(job, any_url, encryption_key=encryption_key,
                                    per_file_encryption=per_file_encryption,
                                    write_to_jobstore=write_to_jobstore)
        else:
            raise
    else:
        with open(filename, 'w') as f:
            f.write(response.read())

    if write_to_jobstore:
        filename = job.fileStore.writeGlobalFile(filename)
    return filename


def bam2fastq(bamfile, univ_options, picard_options):
    """
    Split an input bam to paired fastqs.

    :param str bamfile: Path to a bam file
    :param dict univ_options: Dict of universal options used by almost all tools
    :param dict picard_options: Dict of options specific to Picard
    :return: Path to the _1.fastq file
    :rtype: str
    """
    work_dir = os.path.split(bamfile)[0]
    base_name = os.path.split(os.path.splitext(bamfile)[0])[1]
    parameters = ['SamToFastq',
                  ''.join(['I=', docker_path(bamfile)]),
                  ''.join(['F=/data/', base_name, '_1.fastq']),
                  ''.join(['F2=/data/', base_name, '_2.fastq']),
                  ''.join(['FU=/data/', base_name, '_UP.fastq'])]
    docker_call(tool='picard', tool_parameters=parameters, work_dir=work_dir,
                dockerhub=univ_options['dockerhub'], java_xmx=univ_options['java_Xmx'],
                tool_version=picard_options['version'])
    first_fastq = ''.join([work_dir, '/', base_name, '_1.fastq'])
    assert os.path.exists(first_fastq)
    return first_fastq


def export_results(job, fsid, file_name, univ_options, subfolder=None):
    """
    Write out a file to a given location. The location can be either a directory on the local
    machine, or a folder with a bucket on AWS.

    :param str fsid: The file store id for the file to be exported
    :param str file_name: The name of the file that neeeds to be exported (path to file is also
           acceptable)
    :param dict univ_options: Dict of universal options used by almost all tools
    :param str subfolder: A sub folder within the main folder where this data should go
    :return: None
    """
    job.fileStore.logToMaster('Exporting %s to output location' % fsid)
    file_name = os.path.basename(file_name)
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
        output_url = 'file://' + os.path.join(output_folder, file_name)
    elif univ_options['storage_location'].startswith('aws'):
        # Handle AWS
        bucket_name = univ_options['storage_location'].split(':')[-1]
        output_url = os.path.join('S3://', bucket_name, output_folder.strip('/'), file_name)
    # Can't do Azure or google yet.
    else:
        # TODO: Azure support
        print("Currently doesn't support anything but Local and aws.")
        return
    job.fileStore.exportFile(fsid, output_url)


def delete_fastqs(job, patient_dict):
    """
    Delete the fastqs from the job Store once their purpose has been achieved (i.e. after all
    mapping steps)

    :param dict patient_dict: Dict of list of input fastqs
    """
    for key in patient_dict.keys():
        if 'fastq' not in key:
            continue
        job.fileStore.logToMaster('Deleting "%s:%s" ' % (patient_dict['patient_id'], key) +
                                  'from the filestore.')
        job.fileStore.deleteGlobalFile(patient_dict[key])
    return None


def delete_bams(job, bams, patient_id):
    """
    Delete the bams from the job Store once their purpose has been achieved (i.e. after all
    mutation calling steps). Will also delete the chimeric junction file from Star.

    :param dict bams: Dict of bam and bai files
    :param str patient_id: The ID of the patient for logging purposes.
    """
    bams = {b: v for b, v in bams.items() if b.endswith('.bam') or b.endswith('.bai')}
    if bams:
        for key, val in bams.items():
            job.fileStore.logToMaster('Deleting "%s" for patient "%s".' % (key, patient_id))
            job.fileStore.deleteGlobalFile(val)
    elif 'rna_genome' in bams:
        delete_bams(job, bams['rna_genome'], patient_id)
        job.fileStore.logToMaster('Deleting "rna_transcriptome.bam" for patient "%s".' % patient_id)
        job.fileStore.deleteGlobalFile(bams['rna_transcriptome.bam'])

    elif 'rnaChimeric.out.junction' in bams:
        job.fileStore.logToMaster('Deleting "rnaChimeric.out.junction" for patient "%s".' % patient_id)
        job.fileStore.deleteGlobalFile(bams['rnaChimeric.out.junction'])

    else:
        assert False


# Exception for bad parameters provided
class ParameterError(Exception):
    """
    This Error Class will be raised  in the case of a bad parameter provided.
    """
    pass


def read_peptide_file(in_peptfile):
    """
    Reads an input peptide fasta file into memory in the form of a dict of fasta record: sequence

    :param str in_peptfile: Path to a peptide fasta
    :return: Dict of fasta record: sequence
    :rtype: dict
    """
    peptides = defaultdict()
    pept = None
    with open(in_peptfile, 'r') as peptfile:
        for line in peptfile:
            if line.startswith('>'):
                pept = line.strip().lstrip('>')
                peptides[pept] = ''
            else:
                peptides[pept] = line.strip()
    return peptides


def parse_chromosome_string(job, chromosome_string):
    """
    Parse a chromosome string into a list.

    :param chromosome_string: Input chromosome string
    :return: list of chromosomes to handle
    :rtype: list
    """
    if chromosome_string is None:
        return []
    else:
        assert isinstance(chromosome_string, str)
        chroms = [c.strip() for c in chromosome_string.split(',')]
        if 'canonical' in chroms:
            assert 'canonical_chr' not in chroms, 'Cannot have canonical and canonical_chr'
            chr_prefix = False
            chroms.remove('canonical')
            out_chroms = [str(c) for c in range(1, 23)] + ['X', 'Y']
        elif 'canonical_chr' in chroms:
            assert 'canonical' not in chroms, 'Cannot have canonical and canonical_chr'
            chr_prefix = True
            chroms.remove('canonical_chr')
            out_chroms = ['chr' + str(c) for c in range(1, 23)] + ['chrX', 'chrY']
        else:
            chr_prefix = None
            out_chroms = []
        for chrom in chroms:
            if chr_prefix is not None and chrom.startswith('chr') is not chr_prefix:
                job.fileStore.logToMaster('chromosome %s does not match the rest that %s begin '
                                          'with "chr".' % (chrom,
                                                           'all' if chr_prefix else 'don\'t'),
                                          level=logging.WARNING)
            out_chroms.append(chrom)
        return chrom_sorted(out_chroms)


def chrom_sorted(in_chroms):
    """
    Sort a list of chromosomes in the order 1..22, X, Y, M, <others in alphabetical order>.

    :param list in_chroms: Input chromosomes
    :return: Sorted chromosomes
    :rtype: list[str]
    """
    in_chroms.sort()
    canonicals = [str(c) for c in range(1, 23)] + ['X', 'Y', 'M', 'MT']
    canonical_chr = ['chr' + c for c in canonicals]
    out_chroms_dict = {
        'can': [c for c in in_chroms if c in canonicals],
        'can_chr': [c for c in in_chroms if c in canonical_chr],
        'others': [c for c in in_chroms if c not in canonicals + canonical_chr]}

    assert not (out_chroms_dict['can'] and out_chroms_dict['can_chr'])
    assert not ('M' in out_chroms_dict['can']and 'MT' in out_chroms_dict['can'])
    assert not ('chrM' in out_chroms_dict['can_chr'] and 'chrMT' in out_chroms_dict['can_chr'])

    out_chroms_dict['can'] = canonical_chrom_sorted(out_chroms_dict['can'])
    out_chroms_dict['can_chr'] = canonical_chrom_sorted(out_chroms_dict['can_chr'])

    out_chroms = out_chroms_dict['can'] or out_chroms_dict['can_chr']
    out_chroms.extend(out_chroms_dict['others'])
    return out_chroms


def canonical_chrom_sorted(in_chroms):
    """
    Sort a list of chromosomes in the order 1..22, X, Y, M/MT

    :param list in_chroms: Input chromosomes
    :return: Sorted chromosomes
    :rtype: list[str]
    """
    if len(in_chroms) == 0:
        return []
    chr_prefix = False
    mt = False
    if in_chroms[0].startswith('chr'):
        in_chroms = [x.lstrip('chr') for x in in_chroms]
        chr_prefix = True
    if 'MT' in in_chroms:
        in_chroms[in_chroms.index('MT')] = 'M'
        mt = True
    in_chroms = sorted(in_chroms, key=lambda c: int(c) if c not in ('X', 'Y', 'M') else c)
    try:
        m_index = in_chroms.index('M')
    except ValueError:
        pass
    else:
        in_chroms.pop(m_index)
        in_chroms.append('M')
    # At this point it should be nicely sorted
    if mt:
        in_chroms[in_chroms.index('M')] = 'MT'
    if chr_prefix:
        in_chroms = [''.join(['chr', x]) for x in in_chroms]
    return in_chroms


def email_report(job, univ_options):
    """
    Send an email to the user when the run finishes.

    :param dict univ_options: Dict of universal options used by almost all tools
    """
    fromadd = "results@protect.cgl.genomics.ucsc.edu"
    msg = MIMEMultipart()
    msg['From'] = fromadd
    if  univ_options['mail_to'] is None:
        return
    else:
        msg['To'] = univ_options['mail_to']
    msg['Subject'] = "Protect run for sample %s completed successfully." % univ_options['patient']
    body = "Protect run for sample %s completed successfully." % univ_options['patient']
    msg.attach(MIMEText(body, 'plain'))
    text = msg.as_string()

    try:
        server = smtplib.SMTP('localhost')
    except socket.error as e:
        if e.errno == 111:
            print('No mail utils on this maachine')
        else:
            print('Unexpected error while attempting to send an email.')
        print('Could not send email report')
    except:
        print('Could not send email report')
    else:
        server.sendmail(fromadd, msg['To'], text)
        server.quit()

from __future__ import print_function

import tarfile
import argparse
import json
import logging
import os
import shutil
import subprocess
import sys
import textwrap

log = logging.getLogger(__name__)
 #TODO: add in comments about default arguments

class PipelineWrapperBuilder(object):
    """
    This class can be used to define wrapper scripts to run specific Toil pipelines in docker
    containers. The purpose of this class is to provide a convenient way to define a command line
    interface for a wrapper script and the logic to run the pipeline with this interface.
    """
    def __init__(self, name, desc, config, add_no_clean=True, add_resume=True):
        """
        :param str name: The name of the command to start the workflow.
        :param str desc: The description of the workflow.
        :param str config: A format string where each key exactly matches the name of an argument
            defined in the arg_builder context manager. Note that dashes in argument names should be
            changed to underscores in this string e.g. 'no-clean' should be 'no_clean'.
        """
        self._name = name
        self._desc = desc
        self._config = config
        self._no_clean = add_no_clean
        self._resume = add_resume

    def run(self, args, pipeline_command):
        """
        Invokes the pipeline with the defined command. Command line arguments, and the command need
        to be set with arg_builder, and command_builder respectively before this method can be
        invoked.
        """

        # output that must be moved but not renamed
        consistentNaming = ['alignments/normal_dna_fix_pg_sorted.bam',
                           'alignments/normal_dna_fix_pg_sorted.bam.bai',
                           'alignments/rna_genome_sorted.bam',
                           'alignments/rna_genome_sorted.bam.bai',
                           'alignments/rna_transcriptome.bam',
                           'alignments/tumor_dna_fix_pg_sorted.bam',
                           'alignments/tumor_dna_fix_pg_sorted.bam.bai',
                           'mutations/merged/all_merged.vcf',
                           'rankboost/mhcii_rankboost_concise_results.tsv',
                           'rankboost/mhci_rankboost_concise_results.tsv',
                           ]

        # output that must be renamed as well as moved
        # map of the original name to the final name
        renamingNeeded = {'binding_predictions': 'binding_predictions.tar',
                         'expression': 'expression.tar',
                         'haplotyping': 'haplotyping.tar',
                         'peptides': 'peptides.tar',
                         'rankboost': 'rankboost.tar',
                         'reports': 'reports.tar',
                         'mutations/snpeffed/mutations.vcf': 'all_snpeffed.vcf',
                         'mutations/transgened/mutations.vcf': 'all_transgened.vcf',
                         'mutations/merged': 'merged_perchrom.tar',
                         'mutations/muse': 'muse_perchrom.tar',
                         'mutations/mutect': 'mutect_perchrom.tar',
                         'mutations/radia': 'radia_perchrom.tar',
                         'mutations/somaticsniper': 'somaticsniper_perchrom.tar',
                         'mutations/strelka/snv': 'strelka_snv_perchrom.tar',
                         'mutations/strelka/indel': 'strelka_indel_perchrom.tar'}

        def make_output(output_dir, source_dir):
            """
            :param output_dir: dir to write the output to
            :param source_dir: dir containing the directory structure to be parsed
            :return:
            """
            def make_tar(dir, tar):
                with tarfile.open(tar, "w:gz") as tar:
                    tar.add(dir)
            # the output dir is where the real output directories are written
            protect_outputs = os.listdir(source_dir)
            for protectOut in protect_outputs:
                def getName(fileName):
                    return os.path.join(os.path.join(source_dir, protectOut), fileName)
                # move individual files out
                for fileName in consistentNaming:
                    shutil.copyfile(getName(fileName), os.path.join(output_dir, os.path.basename(fileName)))
                for src, dst in renamingNeeded.iteritems():
                    if dst.endswith('.tar'):
                        make_tar(getName(src), os.path.join(output_dir, dst))
                    else:
                        shutil.copyfile(getName(src), os.path.join(output_dir, dst))
            shutil.rmtree(source_dir)


        # prepare workdir
        mount = self._prepare_mount(args)
        self._workdir = os.path.join(mount, 'Toil-' + self._name)
        # insure the pairs are in the same directory, as protect expects
        # This is made more complicated by the fact CWLTool mounts inputs into random, read-only dirs
        # to get around this we copy all inputs into their own directories that we own
        tumor_dna_dir = os.path.expanduser('~/tumorDNA')
        tumor_rna_dir = os.path.expanduser('~/tumorRNA')
        normal_dna_dir = os.path.expanduser('~/normalDNA')
        os.mkdir(tumor_dna_dir)
        os.mkdir(tumor_rna_dir)
        os.mkdir(normal_dna_dir)
        shutil.copy(args.tumor_dna, tumor_dna_dir)
        shutil.copy(args.tumor_rna, tumor_rna_dir)
        shutil.copy(args.normal_dna, normal_dna_dir)
        shutil.copy(args.tumor_dna2, tumor_dna_dir)
        shutil.copy(args.tumor_rna2, tumor_rna_dir)
        shutil.copy(args.normal_dna2, normal_dna_dir)
        args.tumor_dna = os.path.join(tumor_dna_dir, os.path.basename(args.tumor_dna))
        args.tumor_dna2 = os.path.join(tumor_dna_dir, os.path.basename(args.tumor_dna2))
        args.tumor_rna = os.path.join(tumor_rna_dir, os.path.basename(args.tumor_rna))
        args.tumor_rna2 = os.path.join(tumor_rna_dir, os.path.basename(args.tumor_rna2))
        args.normal_dna = os.path.join(normal_dna_dir, os.path.basename(args.normal_dna))
        args.normal_dna2 = os.path.join(normal_dna_dir, os.path.basename(args.normal_dna2))

        # prepare config
        args_dict = vars(args)
        args_dict['output_dir'] = mount
        self._config = textwrap.dedent(self._config.format(**args_dict))
        self._sample_name = args_dict["sample_name"]
        config_path = os.path.join(self._workdir, 'config')
        command = self._make_prefix(os.path.join(self._workdir, 'jobStore'),
                                          config_path,
                                          self._workdir) + pipeline_command
        if self._resume and args.resume:
            command.append('--restart')
        self._create_workdir(args)
        with open(config_path, 'w') as f:
            f.write(self._config)

        try:
            subprocess.check_call(command)
        except subprocess.CalledProcessError as e:
            print(e, file=sys.stderr)
        finally:
            log.info('Pipeline terminated, changing ownership of output files from root to user.')
            stat = os.stat(self._mount)
            subprocess.check_call(['chown', '-R', '{}:{}'.format(stat.st_uid, stat.st_gid),
                                   self._mount])
            make_output(self._mount, os.path.join(self._mount, 'output'))
            if self._no_clean and args.no_clean:
                log.info('Flag "--no-clean" was used, therefore %s was not deleted.', self._workdir)
            else:
                log.info('Cleaning up temporary directory: %s', self._workdir)
                shutil.rmtree(self._workdir)

    def get_args(self):
        """
        Use this context manager to add arguments to an argparse object with the add_argument
        method. Arguments must be defined before the command is defined. Note that
        no-clean and resume are added upon exit and should not be added in the context manager. For
        more info about these default arguments see below.
        """
        parser = argparse.ArgumentParser(description=self._desc,
                                         formatter_class=MyUniversalHelpFormatter)
        # default args
        if  self._no_clean:
            parser.add_argument('--no-clean', action='store_true',
                                help='If this flag is used, temporary work directory is not '
                                     'cleaned.')
        if self._resume:
            parser.add_argument('--resume', action='store_true',
                                help='If this flag is used, a previously uncleaned workflow in the'
                                     ' same directory will be resumed')
        return parser

    def _prepare_mount(self, args):
        assert args is not None
        # Get name of most recent running container. If socket is mounted, should be this one.
        name_command = ['docker', 'ps', '--format', '{{.Names}}']
        try:
            name = subprocess.check_output(name_command).split('\n')[0]
        except subprocess.CalledProcessError as e:
            raise RuntimeError('No container detected, ensure Docker is being run with: '
                               '"-v /var/run/docker.sock:/var/run/docker.sock" as an argument.'
                               ' \n\n{}'.format(e.message))
        # Get name of mounted volume
        blob = json.loads(subprocess.check_output(['docker', 'inspect', name]))
        mounts = blob[0]['Mounts']
        # Ensure docker.sock is mounted correctly
        sock_mnt = [x['Source'] == x['Destination'] for x in mounts if 'docker.sock' in x['Source']]
        require(len(sock_mnt) == 1, 'Missing socket mount. Requires the following: '
                                      'docker run -v /var/run/docker.sock:/var/run/docker.sock')
        self._mount = args.work_mount
        return self._mount

    def _make_prefix(self, jobstore_path, config_path, workdir_path):
        return [self._name, jobstore_path,
                '--config', config_path,
                '--workDir', workdir_path,
                '--retryCount', '1']

    def _create_workdir(self, args):
        if os.path.exists(self._workdir):
            if self._resume and args.resume:
                 log.info('Reusing temporary directory: %s', self._workdir)
            else:
                raise UserError('Temporary directory {} already exists. Run with --resume option or'
                                ' remove directory.'.format(self._workdir))
        else:
            os.makedirs(self._workdir)
            log.info('Temporary directory created: %s', self._workdir)

class UserError(Exception):
    pass


class DefinitionError(Exception):
    pass


def require(expression, message):
    if not expression:
        raise UserError('\n\n' + message + '\n\n')

def check_for_input(tool_input, name):
    require(tool_input, 'Cannot find {0} input, please use --{0} and provide full path to file.'
            .format(name))
    return tool_input[0]

class MyUniversalHelpFormatter(argparse.HelpFormatter):
    '''
    This formatter formats both the description and argument defaults formatting
    of the argparse help string.
    '''
    def _fill_text(self, text, width, indent):
        '''
        This module was taken from the ArgumentDefaultsHelpFormatter class
        within argparse.  It deals with the formatting of arguments in that it
        appends the default value to teh description of each argument.
        '''
        return ''.join(indent + line for line in text.splitlines(True))
    def _get_help_string(self, action):
        '''
        This module was taken from the RawDescriptionHelpFormatter class
        within argparse.  It deals with the formatting of the description string
        and allows properly formatted descriptions ot be printed without line
        wrapping.
        '''
        help = action.help
        if '%(default)' not in action.help:
            if action.default is not argparse.SUPPRESS:
                defaulting_nargs = [argparse.OPTIONAL, argparse.ZERO_OR_MORE]
                if action.option_strings or action.nargs in defaulting_nargs:
                    help += ' (default: %(default)s)'
        return help

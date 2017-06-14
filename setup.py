from __future__ import print_function
from pkg_resources import parse_version, SetuptoolsLegacyVersion
from setuptools import find_packages, setup
from setuptools.command.test import test as TestCommand
from version import version

import errno
import subprocess
import sys


toil_version = '3.8.0'
s3am_version = '2.0.1'
gdc_version = 'v1.1.0'


def check_tool_version(tool, required_version, binary=False):
    """
    This will ensure that the required_version of `tool` is at least `required_version`.
    :param str tool: The tool under review
    :param str required_version: The version of the tool required by ProTECT
    :param bool binary: Is the tool a binary
    :return: None
    """
    if binary:
        try:
            installed_version = subprocess.check_output([tool, '--version'],
                                                        stderr=subprocess.STDOUT)
        except OSError as err:
            if err.errno == errno.ENOENT:
                raise RuntimeError('Is %s installed as a binary and present on your $PATH?' % tool)
            else:
                raise
        installed_version = installed_version.rstrip()
    else:
        try:
            module = __import__(tool + '.version')
        except ImportError:
            raise RuntimeError('Is %s installed as a library, and is it accessible in the same '
                               'environment as ProTECT?' % tool)
        try:
            installed_version = getattr(module, 'version').version
        except AttributeError:
            raise RuntimeError('Does %s have a version.py?' % tool)

    if type(parse_version(installed_version)) == SetuptoolsLegacyVersion:
        print('Detecting that the installed version of "%s"(%s) is probably based off a git commit '
              'and assuming this build is for testing purposes.  If this is not the case, please '
              'try again with a valid version of "%s".' % (tool, installed_version, tool))
    elif parse_version(installed_version) < parse_version(required_version):
        raise RuntimeError('%s was detected to be version (%s) but ProTECT requires (%s)' %
                           (tool, installed_version, required_version))


# Check Toil version
check_tool_version('toil', toil_version, binary=False)
# Check S3am version
check_tool_version('s3am', s3am_version, binary=True)
# Check gdc-client version
check_tool_version('gdc-client', gdc_version, binary=True)


# Set up a test class
# noinspection PyAttributeOutsideInit
class PyTest(TestCommand):
    user_options = [('pytest-args=', 'a', "Arguments to pass to py.test")]

    def initialize_options(self):
        TestCommand.initialize_options(self)
        self.pytest_args = []

    def finalize_options(self):
        TestCommand.finalize_options(self)
        self.test_args = []
        self.test_suite = True

    def run_tests(self):
        import pytest
        # Sanitize command line arguments to avoid confusing Toil code attempting to parse them
        sys.argv[1:] = []
        err_number = pytest.main(self.pytest_args)
        sys.exit(err_number)


setup(name='protect',
      version=version,
      description='Prediction of T-Cell Epitopes for Cancer Therapy',
      url='http://github.com/BD2KGenomics/protect',
      author='Arjun Arkal Rao',
      author_email='aarao@ucsc.edu',
      license='Apache',
      install_requires=[
          'PyYAML',
          'pandas==0.19.2'
      ],
      tests_require=[
          'pytest==2.8.3'],
      test_suite='protect',
      entry_points={
          'console_scripts': [
              'ProTECT = protect.pipeline.ProTECT:main']},
      cmdclass={'test': PyTest},
      package_dir={'': 'src'},
      packages=find_packages('src', exclude=['*.test']),
      package_data={'':['src/protect/pipeline/input_parameters.yaml']},
      include_package_data=True,
      zip_safe=False)

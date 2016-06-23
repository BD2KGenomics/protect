from __future__ import print_function
from distutils.version import LooseVersion

import sys

from setuptools import find_packages, setup
from setuptools.command.test import test as TestCommand
import errno
import subprocess


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

    if LooseVersion(installed_version) < LooseVersion(required_version):
        raise RuntimeError('%s was detected to be version (%s) but ProTECT requires (%s)' %
                           (tool, installed_version, required_version))


# Check Toil version
check_tool_version('toil', '3.2.0', binary=False)
# Check S3am version
check_tool_version('s3am', '2.0', binary=True)


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
      version='2.0',
      description='Prediction of T-Cell Epitopes for Cancer Therapy',
      url='http://github.com/BD2KGenomics/protect',
      author='Arjun Arkal Rao',
      author_email='aarao@ucsc.edu',
      license='Apache',
      install_requires=[
          'PyYAML'
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
      zip_safe=False)

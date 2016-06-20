from __future__ import print_function
from distutils.version import LooseVersion
from setuptools import find_packages, setup
from setuptools.command.install import install as _install
from setuptools.command.develop import develop as _develop
import errno
import subprocess


class install(_install):
    """
    This is a custom install class that first checks that the Installed versions of Toil (Hopefully
    in the same venv as ProTECT) and s3am (Hopefull in its own venv) are compatible with ProTECT.
    This way, even if we are using s3am and/or Toil from --system-site-packages we won't attempt to
    overwrite it if the versions are incompatible.
    """
    def run(self):
        # Check Toil version
        check_tool_version('toil', '3.2.0', binary=False)
        # Check S3am version
        check_tool_version('s3am', '2.0', binary=True)
        _install.run(self)


class develop(_develop):
    """
    This is a custom develop class that first checks that the Installed versions of Toil (Hopefully
    in the same venv as ProTECT) and s3am (Hopefull in its own venv) are compatible with ProTECT.
    This way, even if we are using s3am and/or Toil from --system-site-packages we won't attempt to
    overwrite it if the versions are incompatible.
    """
    def run(self):
        # Check Toil version
        check_tool_version('toil', '3.2.0', binary=False)
        # Check S3am version
        check_tool_version('s3am', '2.0', binary=True)
        _develop.run(self)


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
            raise RuntimeError('Is %s installed as a library in the same environment as ProTECT?' %
                               tool)
        try:
            installed_version = getattr(module, 'version').version
        except AttributeError:
            raise RuntimeError('Does %s have a version.py?' % tool)

    if LooseVersion(installed_version) < LooseVersion(required_version):
        raise RuntimeError('%s was detected to be version (%s) but ProTECT requires (%s)' %
                           (tool, installed_version, required_version))


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
      cmdclass={'install': install,
                'develop': develop},
      package_dir={'': 'src'},
      packages=find_packages('src', exclude=['*.test']),
      zip_safe=False)

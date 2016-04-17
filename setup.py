from setuptools import setup

setup(name='protect',
      version='1.0',
      description='Prediction of T-Cell Epitopes for Cancer Therapy',
      url='http://github.com/BD2KGenomics/protect',
      author='Arjun Arkal Rao',
      author_email='aarao@ucsc.edu',
      license='Apache',
      install_requires=[
          'toil',
          'pysam'
      ],
      zip_safe=False)
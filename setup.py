from setuptools import find_packages, setup

setup(name='protect',
      version='2.0',
      description='Prediction of T-Cell Epitopes for Cancer Therapy',
      url='http://github.com/BD2KGenomics/protect',
      author='Arjun Arkal Rao',
      author_email='aarao@ucsc.edu',
      license='Apache',
      install_requires=[
          'toil>=3.2.0a2.dev184',
          'PyYAML'
      ],
      package_dir={'': 'src'},
      packages=find_packages('src', exclude=['*.test']),
      zip_safe=False)

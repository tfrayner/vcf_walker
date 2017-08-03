import os
import re

import multiprocessing # workaround for Python bug. See http://bugs.python.org/issue15881#msg170215

from setuptools import setup, find_packages

README  = open(os.path.join(os.path.dirname(__file__), 'README.rst')).read()
SCRIPTS = [ os.path.join('bin', x) for x in os.listdir('bin') if re.search('\.py$', x) ]

# Allow setup.py to be run from any path
os.chdir(os.path.normpath(os.path.join(os.path.abspath(__file__), os.pardir)))

from vcf_walker import __package_version__

setup(
  name='vcf_walker',
  version=__package_version__,
  packages=find_packages(),
  include_package_data=True,
  license='GPLv3 License',
  description='VCF annotation and file format conversion tools.',
  long_description=README,
  url='https://xp-dev.com/svn/tfrayner-nursery',
  author='Tim Rayner',
  author_email='tim.rayner@cruk.cam.ac.uk',
  classifiers=[
    'Environment :: Command Line',
    'Intended Audience :: Developers',
    'License :: OSI Approved :: GPLv3 License',
    'Operating System :: OS Independent',
    'Programming Language :: Python',
    'Programming Language :: Python :: 2',
    'Programming Language :: Python :: 2.7',
  ],
  test_suite='nose.collector',
  scripts=SCRIPTS,
  install_requires=[
    'bx-python>=0.7.3',
    'biopython>=1.66',
    'PyVCF>=0.6.7',
    'bcbio-gff>=0.6.2',
  ],
)

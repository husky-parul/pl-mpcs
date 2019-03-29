import sys
import os
# Make sure we are running python3.5+
if 10 * sys.version_info[0]  + sys.version_info[1] < 35:
    sys.exit("Sorry, only Python 3.5+ is supported.")

from setuptools import setup


def readme():
    print("Current dir = %s" % os.getcwd())
    print(os.listdir())
    with open('README.rst') as f:
        return f.read()

setup(
      name             =   'mpcs',
      version          =   '1.0.3',
      description      =   'This app simulates an MPC compute call and creates a z-score file.', 
      long_description =   readme(),
      author           =   'Rudolph Pienaar',
      author_email     =   'rudolph.pienaar@gmail.com',
      url              =   'https://github.com/FNNDSC/pl-mpcs',
      packages         =   ['mpcs'],
      install_requires =   ['chrisapp', 'pudb', 'numpy'],
      test_suite       =   'nose.collector',
      tests_require    =   ['nose'],
      scripts          =   ['mpcs/mpcs.py'],
      license          =   'MIT',
      zip_safe         =   False
     )

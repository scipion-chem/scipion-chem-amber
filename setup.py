"""A setuptools based setup module.

See:
https://packaging.python.org/en/latest/distributing.html
https://github.com/pypa/sampleproject
"""

# Always prefer setuptools over distutils
from setuptools import setup, find_packages
from codecs import open
from os import path

here = path.abspath(path.dirname(__file__))

# Get the long description from the README file
with open(path.join(here, 'README.rst'), encoding='utf-8') as f:
    long_description = f.read()

with open('requirements.txt') as f:
    requirements = f.read().splitlines()

# Arguments marked as "Required" below must be included for upload to PyPI.
# Fields marked as "Optional" may be commented out.

setup(
    name='scipion-chem-amber',  # Required
    version='0.1',  # Required
    description='Scipion amber plugin.',  # Required
    long_description= 'Scipion AMBER plugin for Scipion. Aida Pinacho Master Thesis',  # Optional
    url='https://github.com/scipion-chem/scipion-chem-amber',  # Optional
    author='Aida Pinacho',  # Optional
    #author_email='you@yourinstitution.email',  # Optional
    #keywords='',  # Optional
    packages=find_packages(),
    install_requires=[requirements],
    entry_points={'pyworkflow.plugin': 'amber=amber'},
    package_data={  # Optional
       'amber': ['icon.png', 'protocols.conf'],
    }
)

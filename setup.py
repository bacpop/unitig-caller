from setuptools import setup, find_packages
from codecs import open
from os import path
import os
import re
import io


def read(*names, **kwargs):
    with io.open(
        os.path.join(os.path.dirname(__file__), *names),
        encoding=kwargs.get("encoding", "utf8")
    ) as fp:
        return fp.read()


def find_version(*file_paths):
    version_file = read(*file_paths)
    version_match = re.search(r"^__version__ = ['\"]([^'\"]*)['\"]",
                              version_file, re.M)
    if version_match:
        return version_match.group(1)
    raise RuntimeError("Unable to find version string.")


here = path.abspath(path.dirname(__file__))


with open(path.join(here, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setup(
    name='unitig_caller',
    version=find_version("pyseer/__init__.py"),
    description='unitig-caller: wrapper around mantis to detect presence of sequence elements',
    long_description=long_description,
    url='https://github.com/johnlees/unitig-caller',
    author='John Lees',
    author_email='john@johnlees.me',
    license='Apache Software License',
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'License :: OSI Approved :: Apache Software License',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7'
    ],
    keywords='gwas bacteria k-mer unitig',
    packages=['unitig_caller'],
    entry_points={
        "console_scripts": [
            'unitig_caller = unitig_caller.__main__:main'
            ]
    }
)

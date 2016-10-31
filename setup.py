#!/usr/bin/env python

try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

config = {
    'name': 'vcfped',
    'version': "1.1.0",
    'description': 'Identification of trios and other close relationships in multisample VCF files',
    'long_description': open('README.md').read(),
    'author': 'Magnus Dehli Vigeland',
    'author_email': 'magnusdv@medisin.uio.no',
    'license': 'GPL-2',
    'url': 'https://github.com/magnusdv/vcfped',
    'packages': [
        'vcfped'
    ],
    'classifiers': [
        'Development Status :: 4 - Beta',
        'Programming Language :: Python :: 2.7',
        'License :: OSI Approved :: GNU General Public License v2 or later (GPLv2+)'
    ],
    'keywords': [
        'relatedness',
        'variant files'
        'trio'
        'exome',
        'vcf',
        'family-based sequencing'
    ],
    'entry_points': {
        'console_scripts': [
            'vcfped = vcfped.vcfped:main'
        ]
    },
    'package_data': {
        'vcfped': [
            'testfiles/*'
        ]
    }
}

setup(**config)

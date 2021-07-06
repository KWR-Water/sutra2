import os
from setuptools import setup, find_packages

def read(fname):
    with open(os.path.join(os.path.dirname(__file__), fname)) as f:
        return f.read()

setup(
    name='greta',
    version='0.2.3',
    packages=find_packages(exclude=['tests*']),
    license='MIT',
    description='A python package for correction, validation and analysis of ground water quality samples',
    long_description=read('README.rst'),
    long_description_content_type="text/x-rst",
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'Intended Audience :: Other Audience',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Topic :: Scientific/Engineering :: Hydrology',
    ],
    python_requires='>=3.6',
    project_urls={
    'Source': 'https://github.com/KWR-Water/greta',
    'Documentation': 'http://greta.readthedocs.io/en/latest/',
    'Tracker': 'https://github.com/KWR-Water/greta/issues',
    'Help': 'https://github.com/KWR-Water/greta/issues',
    # 'Help': 'https://stackoverflow.com/questions/tagged/greta'
    },
    install_requires=[
        'pandas>=0.23',
        'molmass',

        ],
    include_package_data=True,
    url='https://github.com/KWR-Water/greta',
    author='KWR Water Research Institute',
    author_email='martin.korevaar@kwrwater.nl, martin.van.der.schans@kwrwater.nl, alex.hocking@kwrwater.nl, steven.ros@kwrwater.nl',
)

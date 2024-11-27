from setuptools import setup, find_packages

setup(
    name='seqret',
    version='1.0.0',
    packages=find_packages(),
    install_requires=[
        'biopython==1.78',
        'dash==2.12.1',
        'dash-bootstrap-components==1.4.2',
        'plotly==5.16.1',
        'dash-bio==1.0.2',
        'dash-core-components==2.0.0',
        'dash-html-components==2.0.0',
        'dash-table==5.0.0',
        'waitress==2.0.0'
    ],
    entry_points={
        'console_scripts': [
            #'seqret=seqret.app:start_app'
            'seqret=seqret.cli:main'
        ],
    },
    package_data={
        # Specify package and the files to include
        'seqret': ['assets/*.ico'],
    },
    author='Nicholas Freitas',
    author_email='nicholas.freitas@ucsf.edu',
    description='GUI for optimizing DNA sequences for improved protein expression.',
    long_description='https://github.com/pinneylab/SeqRET',
    url='',
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
)
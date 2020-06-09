from setuptools import setup, find_packages

data_files_to_include = [('', ['README.md', 'LICENSE'])]

setup(name='immprint',
      version='0.0.1',
      description='Guess if T-cell samples come from the same individual.',
      long_description='',
      author='Thomas Dupic',
      author_email='dupic.thomas@gmail.com',
      license='GPLv3',
      classifiers=[
            'Development Status :: 3 - Alpha',
            'Intended Audience :: Healthcare Industry',
            'Intended Audience :: Science/Research',
            'Topic :: Scientific/Engineering :: Bio-Informatics',
            'Topic :: Scientific/Engineering :: Physics',
            'Topic :: Scientific/Engineering :: Medical Science Apps.',
            'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
            'Natural Language :: English',
            'Programming Language :: Python :: 3.6',
            ],
      setup_requires=['wheel'],
      entry_points={'console_scripts': 'immprint = immprint.main:parse_arguments'},
      python_requires='>=3.6',
      packages=find_packages())

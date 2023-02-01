from setuptools import setup,find_packages

def readme():
    with open('README.md') as f:
        return f.read()

setup(name='salty',
      version='1.0.0',
      description='In silico lineage typing of Staphylococcus aureus.',
      long_description=readme(),
      classifiers=[
          'License :: OSI Approved :: GPLv3',
          'Programming Language :: Python :: 3.7',
          'Topic :: Scientific/Engineering :: Bio-Informatics',
          'Topic :: Scientific/Engineering :: Medical Science Apps.',
          'Intended Audience :: Science/Research',
      ],
      keywords='microbial genomics typing',
      url='https://github.com/LanLab/SaLTy',
      author='Liam Cheney',
      author_email='liam.cheney1@gmail.com',
      license='GPLv3',
      packages=find_packages(),
      include_package_data=True,
      entry_points={
          'console_scripts': ['salty=salty.salty:main'],
      },
      zip_safe=False)

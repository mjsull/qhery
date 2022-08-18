from setuptools import setup, find_packages

setup(name='qhery',
      version='0.1.0',
      description='SARS-CoV-2 resistance analysis tool.',
      author='Mitchell Sullivan',
      author_email='mjsull@gmail.com',
      url='https://github.com/mjsull/qhery',
      python_requires='>=3.6, <4',
      packages=find_packages(include=['qhery', 'qhery.*']),
      entry_points={
            "console_scripts": [
                  "qhery = qhery.__main__:main"
                  ]
      },
      install_requires=[
            'pysam>=0.19.1'
      ],
      setup_requires=['pytest-runner', 'flake8'],
      tests_require=['pytest'],
      package_data={
            'qhery': ['data/nCoV-2019.reference.fasta',
                      'data/nsp_coords.tsv',
                      'data/proteins.faa',
                      'data/resistance_table.tsv',
                      'data/Sars_cov_2.ASM985889v3.101.gff3'],
      })
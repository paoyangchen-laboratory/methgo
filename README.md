# MethGo
**MethGo** is a comprehensive tool for analyzing whole genome bisulfite sequencing data.

## Dependencies
It needs Python 2.7, samtools, bedtools, and the packages listed below:
- [NumPy](http://www.numpy.org/)
- [matplotlib](http://matplotlib.org/)
- [pandas](http://pandas.pydata.org/)
- [pysam](http://pysam.readthedocs.org/)
- [biopython](http://biopython.org/)
- [pyfasta](https://pypi.python.org/pypi/pyfasta/)
- [Cython](http://cython.org/)
- [pybedtools](https://pythonhosted.org/pybedtools/)

## Setting up environment
1. Installing virtualenv

    - Installing globally with pip (if you have pip 1.3 or greater installed globally)

            $ [sudo] pip install virtualenv
            $ virtualenv --no-site-packages --python=python2.7 .venv

    - Installing locally from source

            $ curl -O https://pypi.python.org/packages/source/v/virtualenv/virtualenv-12.1.1.tar.gz
            $ tar xvfz virtualenv-12.1.1.tar.gz
            $ python virtualenv-12.1.1/virtualenv.py --no-site-packages --python=python2.7 .venv

2. Activating the virtual environment

        $ source .venv/bin/activate

3. Installing packages

        $ pip install -r methgo/requirements.txt

## License
MIT

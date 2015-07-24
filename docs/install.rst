Installation
============

1. Obtain Python 2.7 and virturalenv.

  .. note::
    MethGo depends on SAMtools and BEDtools, so please make sure you already
    have them on your server.

2. Create a virtual environment somewhere on your disk, and then activate it.

  ::

  $ virturalenv --no-site-packages --python=python2.7 methgo_env
  $ cd methgo_env
  $ source bin/activate


3. Download the source code and install the requirements.

  ::

  $ git clone https://github.com/wwliao/methgo.git
  $ pip install -r methgo/requirements/base.txt
  $ pip install -r methgo/requirements/addition.txt

  .. note::
    If you're using Mac and the OS version is larger than 10.8, run the
    following line before you install the requirements

    ::

    $ export CFLAGS=-Qunused-arguments

  pip will install the following packages:

  * NumPy
  * SciPy
  * matplotlib
  * pandas
  * PySAM == 0.8.0
  * Biopython
  * pyfasta
  * Cython
  * pybedtools

4. Add <your methgo path> to the PATH environment variable.

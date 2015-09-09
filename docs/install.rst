Installation
============

1. Obtain Python 2.7 and virturalenv.

  .. note::
    MethGo depends on `SAMtools <http://www.htslib.org/>`_ and
    `BEDtools <http://bedtools.readthedocs.org/>`_, so please make sure you
    already have them on your server.

2. Create a virtual environment somewhere on your disk, and then activate it.

  ::

  $ virtualenv --no-site-packages --python=python2.7 methgo_env
  $ source methgo_env/bin/activate


3. Download the source code and install the requirements.

  ::

  $ git clone https://github.com/wwliao/methgo.git
  $ pip install -r methgo/requirements/base.txt
  $ pip install -r methgo/requirements/addition.txt

  .. note::
    If you're using Mac and the OS version is larger than 10.8, run the
    following line before you install the requirements:

    ::

    $ export CFLAGS=-Qunused-arguments

  pip will install the following packages:

  * `NumPy <http://www.numpy.org/>`_
  * `SciPy <http://www.scipy.org/>`_
  * `matplotlib <http://matplotlib.org/>`_
  * `pandas <http://matplotlib.org/>`_
  * `PySAM (0.8.0) <http://matplotlib.org/>`_
  * `Biopython <http://biopython.org/>`_
  * `pyfasta <https://pypi.python.org/pypi/pyfasta/>`_
  * `Cython <http://cython.org/>`_
  * `pybedtools <https://pythonhosted.org/pybedtools/>`_

4. Add your MethGo path to the PATH environment variable.

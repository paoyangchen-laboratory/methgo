Installation
============

1. Obtain Python 2.7 and virturalenv.

  .. note::
    MethGo depends on SAMtools and BEDtools, so please make sure you already
    have them on your server.

2. Create a virtual environment somewhere on your disk, and then activate it.

  ::

  $ virturalenv methgo_env
  $ cd methgo_env
  $ source bin/activate


3. Download the source code and install the requirements.

  ::

  $ git clone https://github.com/wwliao/methgo.git
  $ cd methgo
  $ chmod +x methgo
  $ echo 'export PATH=$PATH:~/.local/methgo' >> ~/.bashrc
  $ source ~/.bashrc
  $ pip install -r requirements/base.txt
  $ pip install -r requirements/addition.txt


  pip will install the following packages:

  * NumPy
  * matplotlib
  * pandas
  * PySAM == 0.8.0
  * Biopython
  * pyfasta
  * Cython
  * pybedtools


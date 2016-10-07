Installation Guide
==================

You need a working python (both python2 and 3) installation together with some standard modules.
You can use for example the `anaconda3 installer <https://www.continuum.io/downloads/>`_.

The advantage of the anaconda3 installer is that you get a lot of additional modules and programs,
that make it really easy to work with python. 
For example `Ipython <http://ipython.org/>`_ and the `jupyter notebooks <http://jupyter.org/>`_.
I highly recommend to use those.

Unix
++++

Being in the root directory of this package, where ``setup.py`` lies, just type in your terminal::

    python setup.py install --user

This should also resolve all dependencies automatically.


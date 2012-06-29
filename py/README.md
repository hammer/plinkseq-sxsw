To build:

    $ cd pyplinkseqint
    $ make
    $ cp pyplinkseqint.h /usr/local/include
    $ cp libpyplinkseqint.dylib /usr/local/lib
    $ cd ../pyplinkseq
    $ python setup.py build_ext --inplace
    $ cp pyplinkseq.so /usr/local/lib/pythonX.Y/site-packages


To use:

    >>> import pyplinkseq
    >>> pyplinkseq.gstore_version()
    >>> pyplinkseq.set_project("/path/to/plinkseq/project")
    >>> pyplinkseq.summary()




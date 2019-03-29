pl-mpcs
================================

.. image:: https://badge.fury.io/py/mpcs.svg
    :target: https://badge.fury.io/py/mpcs

.. image:: https://travis-ci.org/FNNDSC/mpcs.svg?branch=master
    :target: https://travis-ci.org/FNNDSC/mpcs

.. image:: https://img.shields.io/badge/python-3.5%2B-blue.svg
    :target: https://badge.fury.io/py/pl-mpcs

.. contents:: Table of Contents


Abstract
--------

This app simulates  a call to a remote Multi-Party Compute (MPC) in the context of a FreeSurfer workflow.

This particular application simply returns a z-score file to be consumed by a downstream plugin, typciall ``pl-z2labelmap``.

NOTE: The <inputDir> is largely ignored by this plugin.


Synopsis
--------

.. code:: bash

    python mpcs.py                                                  \
        [--random] [--seed <seed>]                                  \
        [-p <f_posRange>] [--posRange <f_posRange>]                 \
        [-n <f_negRange>] [--negRange <f_negRange>]                 \
        [-z <zFile>] [--zFile <zFile>]                              \
        [-v <level>] [--verbosity <level>]                          \
        [--version]                                                 \
        [--man]                                                     \
        [--meta]                                                    \
        <inputDir>
        <outputDir> 


Run
----

This ``plugin`` can be run in two modes: natively as a python package or as a containerized docker image.

Using PyPI
~~~~~~~~~~

To run from PyPI, simply do a 

.. code:: bash

    pip install mpcs

and run with

.. code:: bash

    mpcs.py --man /tmp /tmp

to get inline help. The app should also understand being called with only two positional arguments

.. code:: bash

    mpcs.py /some/input/directory /destination/directory


Using ``docker run``
~~~~~~~~~~~~~~~~~~~~

To run using ``docker``, be sure to assign an "input" directory to ``/incoming`` and an output directory to ``/outgoing``. *Make sure that the* ``$(pwd)/out`` *directory is world writable!*

Now, prefix all calls with 

.. code:: bash

    docker run --rm -v $(pwd)/out:/outgoing                             \
            fnndsc/pl-mpcs mpcs.py                                      \

Thus, getting inline help is:

.. code:: bash

    mkdir in out && chmod 777 out
    docker run --rm -v $(pwd)/in:/incoming -v $(pwd)/out:/outgoing      \
            fnndsc/pl-mpcs mpcs.py                                      \
            --man                                                       \
            /incoming /outgoing

Arguments
---------

.. code::

    [--random] [--seed <seed>]
    If specified, generate a z-score file based on <posRange> and 
    <negRange>. In addition, if a further optional <seed> is passed,
    then initialize the random generator with that seed, otherwise
    system time is used.

    [-p <f_posRange>] [--posRange <f_posRange>]
    Positive range for random max deviation generation.

    [-n <f_negRange>] [--negRange <f_negRange>]
    Negative range for random max deviation generation.

    [-z <zFile>] [--zFile <zFile>]
    z-score file to save in output directory. Defaults to 'zfile.csv'.

    [-v <level>] [--verbosity <level>]
    Verbosity level for app. Not used currently.

    [--version]
    If specified, print version number. 
    
    [--man]
    If specified, print (this) man page.

    [--meta]
    If specified, print plugin meta data.


Examples
--------

Create a z-file with values between -3.0 and +3.0
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: bash

    mkdir in out && chmod 777 out
    docker run --rm -v $(pwd)/in:/incoming -v $(pwd)/out:/outgoing      \
            fnndsc/pl-mpcs mpcs.py                                      \
            -random --seed 1                                            \
            --posRange 3.0 --negRange -3.0                              \
            in out







cromosim
========

CROMOSIM is a Python Library for Crowd Motion Simulation.

The aim of this open source project is to make it possible to users to run
simulations based from several models (cellular automata, microscopic
simulations or using compartments), to test new configurations, and even
to investigate the possibility to program their own model in complex
geometry: do-it yourself !

This package proposes Python implementations of the numerical methods detailed
in the book “Crowds in equations: an introduction to the microscopic modeling
of crowds” by B. Maury (ENS Ulm & Univ. Paris-Saclay) and S. Faure (CNRS), World
Scientific 2018, Advanced textbooks in mathematics.

Github repository
=================

Download the source code and the examples here:

https://github.com/sylvain-faure/cromosim/

How to use cromosim ?
=====================

First you have to install cromosim, either by using pip

.. code-block:: javascript

  ~$ pip install cromosim

or by manually installing the package

.. code-block:: python

  cromosim$ python setup.py install

Once cromosim is installed, you can verify that it is possible to import it
into Python:

.. code-block:: python

  ~$ python
  Python 3.7.7 (default, Mar 10 2020, 15:43:33)
  [Clang 11.0.0 (clang-1100.0.33.17)] on darwin
  Type "help", "copyright", "credits" or "license" for more information.
  >>> import cromosim
  >>> print(cromosim.__version__)
  2.0.0
  >>>

Now to make a first simulation, you can download one of the examples found in:

https://github.com/sylvain-faure/cromosim/tree/master/examples

or retrieve all the examples available using the following command (``svn``
is the Subversion command):

.. code-block:: python

  ~$ svn export https://github.com/sylvain\-faure/cromosim/trunk/examples my-cromosim

and then run a first example:

.. code-block:: python

  ~$ cd my-cromosim/micro/granular
  granular$ python micro_granular.py --json input_room.json

These examples will allow you to start your own calculations.


Documentation:
==============

`<http://www.cromosim.fr>`_

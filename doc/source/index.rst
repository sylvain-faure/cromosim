.. cromosim documentation master file, created by
   sphinx-quickstart on Wed Sep 13 07:38:52 2017.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to cromosim's documentation!
====================================

A Python Library for microscopic Crowd Motion Simulation.

The aim of this open source project is to allow users to run
simulations based on several models (cellular automata, microscopic simulations
or using compartments), to test new configurations, and even to investigate the
possibility to program its own model in complex geometry : do-it yourself !

This package proposes Python implementations of the numerical methods detailed
in the following book:

.. [MF2018] "*Crowds in equations: an introduction to the microscopic modeling of
    crowds*" by  B. Maury (ENS Ulm & Univ. Paris-Saclay) and S. Faure (CNRS), World
    Scientific 2018, Advanced textbooks in mathematics.

This documentation concerns the **Release 2.0** of cromosim. The old documentation
corresponding to the **Release 1.0** is available there: http://www.cromosim.fr/1.0


Main new features of the **Release 2.0** of cromosim (mainly concerning
microscopic crowd modeling):

  * addition of a Destination class allowing to use a more elaborate color
    code: a door can be represented by a red line, another by a yellow line, etc ...
  * possibility of having several elements of the Domain class. It is for
    example possible to define a domain per floor.
  * possibility of having several possible destinations for people and
    chaining them
  * possibility of defining stairs, thanks to the Destination class
  * during the initialization of groups of people, their radius and speed of
    movement can be defined with probability laws (normal or uniform). They are
    also assigned an initial destination.

Unfortunately these big changes break the backward compatibility with
**Release 1.0** scripts, at least for microscopic crowd modeling.

This documentation is hosted on a server of the Mathematics Department of Orsay
(Univ. Paris-Saclay, France).

.. raw:: html

   <div align=center>
   <img src="_static/logo_lmo.png" alt="Laboratoire de Mathématiques d'Orsay" width="25%" height="auto"/>
   <img src="_static/logo_Paris-Saclay.png" alt="Université Paris-Saclay" width="25%" height="auto"/>
   </div>

.. currentmodule:: cromosim

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

.. code-block:: bash

  cromosim$ python setup.py install

Once cromosim is installed, you can verify that it is possible to import it
into Python:

.. code-block:: bash

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

.. code-block:: bash

  ~$ svn export https://github.com/sylvain\-faure/cromosim/trunk/examples my-cromosim

and then run a first example:

.. code-block:: bash

  ~$ cd my-cromosim/micro/granular
  granular$ python micro_granular.py --json input_room.json

These examples will allow you to start your own calculations.


Examples for users
==================

.. toctree::
   :maxdepth: 2

   Build a domain <example_dom>
   Understand one-dimensional microscopic models : the Follow-the-Leader model <example_ftl>
   Try a cellular automata <example_ca>
   Make a microscopic simulation <example_micro>
   Discover the compartment model <example_comp>

Documentation of the code
=========================

.. toctree::
   :maxdepth: 2
   :caption: Contents :

   Classes used to build the computational domain <domain>
   Module for the Follow the Leader model <ftl>
   Module for cellular automata <ca>
   Module for microscopic simulations with the granular or the social force models <micro>
   Module for the compartment model <comp>

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

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

This documentation concerns the **Release 3.X** of cromosim. 

The old documentations corresponding to the **Release 1.0 or 2.0** are available there: https://www.cromosim.fr/1.0
or https://www.cromosim.fr/2.0

Main new features of the **Release 3.X** of cromosim (mainly concerning
microscopic crowd modeling):

  * The addition of a visualization module to provide more efficient 
    Matplotlib visualization as well as visualization using PyVista for 
    rendering.
  * The option to follow up one destination with a new destination chosen 
    at random from among those listed (the selection may be weighted).
  * For the velocity projection step of the granular microscopic model, 
    the Nesterov method has been added in addition to the Uzawa method. 
    Two implementations are available: one in Python and the other in C++.
  * The option to add a multiplicative coefficient—beta—in front of the 
    desired speeds of certain individuals in microscopic models. This 
    coefficient is calculated by taking into account these individuals’ 
    frustration; they become frustrated if other individuals obstruct 
    their movement. A frustrated individual will have a small value 
    for beta, which reduces their desired speed; they will not “push” 
    the individuals in front of them.
  * The ability to set the positions of a group of characters on colored 
    pixels. For example, to place characters on chairs, simply add a 
    pixel of a specific color to each chair and specify that color in 
    the lines corresponding to the group within the JSON file.

Main new features of the **Release 2.X** of cromosim (mainly concerning
microscopic crowd modeling):

  * addition of a Destination class allowing to use a more elaborate color
    code: a door can be represented by a red line, another by a yellow line, etc ...
  * possibility of having several elements of the Domain class. It is for
    example possible to define a domain per floor.
  * possibility of having several possible destinations for people and
    chaining them.
  * possibility of defining stairs, thanks to the Destination class.
  * during the initialization of groups of people, their radius and speed of
    movement can be defined with probability laws (normal or uniform). They are
    also assigned an initial destination.

Unfortunately these big changes break the backward compatibility with
**Release 1.0 or 2.0** scripts, at least for microscopic crowd modeling.

This documentation is hosted on a server of the Mathematics Department of Orsay
(Univ. Paris-Saclay, France).

.. raw:: html

   <div align=center>
   <img src="_static/logo_lmo.png" alt="Laboratoire de Mathématiques d'Orsay" width="25%" height="auto"/>
   <img src="_static/logo_Paris-Saclay.png" alt="Université Paris-Saclay" width="25%" height="auto"/>
   <img src="_static/cnrs.png" alt="CNRS" width="9%" height="auto"/>
   </div>

.. currentmodule:: cromosim

Github repository
=================

Download the source code and the examples here:

https://github.com/sylvain-faure/cromosim/

How to use cromosim ?
=====================

You can install the main modules with pip:

.. code-block:: bash

  ~$ pip install cromosim

Then, you can either run one of the examples in the “examples” folder 
from the command line, or start with a notebook from the “notebooks” folder.

If you'd rather install Cromosim manually, you can download the entire 
repository from GitHub and use pip:

.. code-block:: bash

  ~$ git clone git@github.com:sylvain-faure/cromosim.git
  ~$ cd cromosim
  cromosim$ pip install .

To try one of the scripts in the 'examples' folder:

.. code-block:: bash

  cromosim$ cd examples/micro/social
  cromosim/examples/micro/social$ python micro_social.py --json input_room.json

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
   Module for visualizing individuals in microscopic models using Matplotlib or PyVista rendering <visu>

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

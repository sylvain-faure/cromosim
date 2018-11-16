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
    crowds*" by  B. Maury (ENS Ulm & Univ. Paris-Sud) and S. Faure (CNRS), World
    Scientific 2018, Advanced textbooks in mathematics.


This documentation is hosted on a server of the Mathematics Department of Orsay
(Univ. Paris-Sud, France).

.. raw:: html

   <div align=center>
   <img src="_static/logo_lmo.png" alt="Laboratoire de Mathématiques d'Orsay" width="25%" height="auto"/>
   <img src="_static/logo_psud.png" alt="Université Paris-Sud" width="15%" height="auto"/>
   </div>

.. currentmodule:: cromosim

Github repository
=================

Download the source code and the examples here:

https://github.com/sylvain-faure/cromosim/

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

   Domain class used to build the computational domain <domain>
   Module for the Follow the Leader model <ftl>
   Module for cellular automata <ca>
   Module for microscopic simulations with the granular or the social force models <micro>
   Module for the compartment model <comp>

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

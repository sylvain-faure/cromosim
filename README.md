# cromosim

CROMOSIM is a Python Library for Crowd Motion Simulation.

- **Website:** http://www.cromosim.fr
- **Source code:** https://github.com/sylvain-faure/cromosim

The aim of this open source project is to make it possible to users to run simulations based from several models (cellular automata, microscopic simulations or using compartments), to test new configurations, and even to investigate the possibility to program their own model in complex geometry : do-it yourself !

This package proposes Python implementations of the numerical methods detailed in the book “Crowds in equations: an introduction to the microscopic modeling of crowds” by B. Maury (ENS Ulm & Univ. Paris-Sud) and S. Faure (CNRS), World Scientific 2018, Advanced textbooks in mathematics.



How to use cromosim ?
---------------------

- You can either install the main modules with pip:

    ~$ pip install cromosim

and then download one of the examples:

    ~$ mkdir my_first_example
    ~$ cd my_first_example
    my_first_example$ wget https://github.com/sylvain-faure/cromosim/blob/master/examples/domain/domain_room.py
    my_first_example$ wget https://github.com/sylvain-faure/cromosim/blob/master/examples/domain/room.png
    my_first_example$ python domain_room.py

- Or clone the entire folder with git and install it manually:

    ~$ git clone git@github.com:sylvain-faure/cromosim.git
    ~$ cd cromosim
    cromosim$ python setup.py install

and then try one of the scripts in the 'examples' folder :
    
    cromosim$ cd examples/domain
    cromosim/examples/domain$ python domain_room.py

For more information
---------------------

https://www.cromosim.fr

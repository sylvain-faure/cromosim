# cromosim

CROMOSIM is a Python Library for Crowd Motion Simulation.

- **Website:** https://www.cromosim.fr
- **Source code:** https://github.com/sylvain-faure/cromosim

The aim of this open source project is to make it possible to users to run simulations based from several models (cellular automata, microscopic simulations or using compartments), to test new configurations, and even to investigate the possibility to program their own model in complex geometry : do-it yourself !

This package proposes Python implementations of the numerical methods detailed in the book “Crowds in equations: an introduction to the microscopic modeling of crowds” by Bertrand Maury (ENS Ulm & Univ. Paris-Saclay) and Sylvain Faure (CNRS), World Scientific 2018, Advanced textbooks in mathematics.



How to use cromosim ?
---------------------

You can either install the main modules with pip

    ~$ pip install cromosim

and then manually download the files for one of the examples you're interested in, such as

    https://github.com/sylvain-faure/cromosim/blob/master/examples/micro/social/micro_social.py
    https://github.com/sylvain-faure/cromosim/blob/master/examples/micro/social/input_room.json
    https://github.com/sylvain-faure/cromosim/blob/master/examples/micro/social/room.png

and finally, once you've retrieved the files, you can run the script:

    ~$ python micro_social.py --json input_room.json

Or clone the entire folder with git and install it manually:

    ~$ git clone git@github.com:sylvain-faure/cromosim.git
    ~$ cd cromosim
    cromosim$ python setup.py install

and then try one of the scripts in the 'examples' folder:
    
    cromosim$ cd examples/micro/social
    cromosim/examples/micro/social$ python micro_social.py --json input_room.json

For more information
---------------------

https://www.cromosim.fr

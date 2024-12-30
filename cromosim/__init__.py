# Authors:
#     Sylvain Faure <sylvain.faure@universite-paris-saclay.fr>
#     Bertrand Maury <bertrand.maury@universite-paris-saclay.fr>
#
# License: GPL

from .version import version as __version__

try:
    from .domain import Domain
    # if Domain:
    #     print("cromosim.domain ==> 'Domain' class is ok")
    # else:
    #     print("cromosim.domain ==> there is a problem with 'Domain' class")
    from .domain import Destination
    # if Destination:
    #     print("cromosim.domain ==> 'Destination' class is ok")
    # else:
    #     print("cromosim.domain ==> there is a problem with 'Destination' class")
except Exception as error:
    print("cromosim.domain ==> an exception occurred on the 'domain' module:", error)

try:
    from . import ca
    # print("cromosim.ca ==> 'cellular automata' module is ok")
except Exception as error:
    print("cromosim.ca ==> an exception occurred on the 'cellular automata' module:", error)

try:
    from . import ftl
    # print("cromosim.ftl ==> 'follow the leader' module is ok")
except Exception as error:
    print("cromosim.ftl ==> an exception occurred on the 'follow the leader':", error)

try:
    from . import comp
    # print("cromosim.comp ==> 'compartment model' module is ok")
except Exception as error:
    print("cromosim.comp ==> an exception occurred on the 'compartment model' module:", error)

try:
    from . import micro
    # print("cromosim.micro ==> 'micro' module is ok")
except Exception as error:
    print("cromosim.micro ==> an exception occurred on the 'micro' module:", error)

name = "cromosim"

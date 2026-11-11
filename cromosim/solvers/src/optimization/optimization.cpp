#include "optimization/minimizer.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

PYBIND11_MODULE(_optimization, module)
{
  module.doc() = "Portable multithreaded optimization solver";
  pybind11::class_<Minimizer>(module, "Minimizer")
    .def(pybind11::init<long long, double, double, double, double, double, bool>(),
         pybind11::arg("max_iterations"), pybind11::arg("rho"),
         pybind11::arg("tolerance"), pybind11::arg("time_step"),
         pybind11::arg("minimum_people_distance"),
         pybind11::arg("minimum_wall_distance"), pybind11::arg("nesterov") = false)
    .def("run", &Minimizer::run, pybind11::arg("xyrv"), pybind11::arg("contacts"),
         pybind11::arg("desired_velocity"), pybind11::arg("mass"),
         pybind11::call_guard<pybind11::gil_scoped_release>())
    .def("set_thread_count", &Minimizer::set_thread_count, pybind11::arg("thread_count"))
    .def_property_readonly("thread_count", &Minimizer::thread_count);
}

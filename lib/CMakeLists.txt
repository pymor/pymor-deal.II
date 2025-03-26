cmake_minimum_required(VERSION 2.8.12)

INCLUDE(dealii_macro)
DEAL_II_INITIALIZE_CACHED_VARIABLES()

# Create the binding library
pybind11_add_module(pymor_dealii_bindings bindings.cc)
# Don't add a 'lib' prefix to the shared library
# set_target_properties(pymor_dealii_bindings PROPERTIES PREFIX "")

PYMOR_DEAL_II_SETUP_TARGET(pymor_dealii_bindings)
install(TARGETS pymor_dealii_bindings LIBRARY DESTINATION .)

# Create the example library
pybind11_add_module(dealii_elasticity py_elasticity.cc elasticity.cc)
# set_target_properties(dealii_elasticity PROPERTIES PREFIX "")
PYMOR_DEAL_II_SETUP_TARGET(dealii_elasticity)
install(TARGETS dealii_elasticity LIBRARY DESTINATION .)

add_executable(eoc eoc.cc elasticity.cc)
target_compile_options(eoc PUBLIC ${PYBIND11_CPP_STANDARD})
PYMOR_DEAL_II_SETUP_TARGET(eoc)

add_subdirectory(pybind11)

FIND_PACKAGE(deal.II 8.5 QUIET CONFIG REQUIRED)
DEAL_II_INITIALIZE_CACHED_VARIABLES()

# Create the binding library
pybind11_add_module(pymor_dealii_bindings bindings.cc)
# Don't add a 'lib' prefix to the shared library
# set_target_properties(pymor_dealii_bindings PROPERTIES PREFIX "")
DEAL_II_SETUP_TARGET(pymor_dealii_bindings)
# target_link_libraries(pymor_dealii_bindings ${PYTHON_LIBRARY})

# Create the example library
pybind11_add_module(dealii_elasticity py_elasticity.cc elasticity.cc)
# set_target_properties(dealii_elasticity PROPERTIES PREFIX "")
DEAL_II_SETUP_TARGET(dealii_elasticity)

add_executable(eoc eoc.cc elasticity.cc)
target_compile_options(eoc PUBLIC ${PYBIND11_CPP_STANDARD})
DEAL_II_SETUP_TARGET(eoc)

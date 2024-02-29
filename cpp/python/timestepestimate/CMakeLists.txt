# SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
# SPDX-License-Identifier: LGPL-3.0-or-later


add_python_targets(timestepestimate __init__ )

dune_add_pybind11_module(NAME _timestepestimate)
set_property(
  TARGET _timestepestimate
  PROPERTY LINK_LIBRARIES Dune::Common ikarus APPEND
)
if(SKBUILD)
  install(TARGETS _timestepestimate LIBRARY DESTINATION python/timestepestimate)
endif()

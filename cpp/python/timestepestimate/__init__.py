# SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
# SPDX-License-Identifier: LGPL-3.0-or-later
from dune.common.hashit import hashIt
from dune.common.checkconfiguration import assertHave, ConfigurationError

try:
    from dune.packagemetadata import registerExternalModule
    import pathlib

    # register timestepestimate to be recognized by dune-py (code generation module)
    # as a module of the dune universe
    registerExternalModule(
        moduleName="timestepestimate",
        modulePath=str(pathlib.Path(__file__).parent.resolve()),
    )

except ImportError:
    pass

from ._timestepestimate import *

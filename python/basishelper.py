# SPDX-FileCopyrightText: 2022 The Ikarus Developers mueller@ibb.uni-stuttgart.de
# SPDX-License-Identifier: LGPL-3.0-or-later

from dune.common.hashit import hashIt
from dune.iga.basis import preBasisTypeName,defaultGlobalBasis
from dune.generator.generator import SimpleGenerator


def globalBasis(gv, tree):
    generator = SimpleGenerator("Basis", "Ikarus::Python")

    pbfName = preBasisTypeName(tree, gv.cppTypeName)
    element_type = f"Ikarus::Basis<{pbfName}>"
    includes = []
    includes += list(gv.cppIncludes)
    includes += ["dune/iga/nurbsbasis.hh"]
    includes += ["ikarus/python/basis/basis.hh"]
    
    moduleName = "Basis_" + hashIt(element_type)
    module = generator.load(
        includes=includes, typeName=element_type, moduleName=moduleName
    )
    basis = defaultGlobalBasis(gv, tree)
    return module.Basis(basis)

# def sparseFlatAssembler(fes, dirichletValues):
#     element_type = f"Ikarus::SparseFlatAssembler<std::vector<{fes[0].cppTypeName}>,{dirichletValues.cppTypeName}>"
#     generator = SimpleGenerator("SparseFlatAssembler", "Ikarus::Python")

#     includes = []
#     includes += ["ikarus/assembler/simpleAssemblers.hh"]
#     includes += fes[0]._includes  # include header of finite element
#     includes += ["ikarus/python/assembler/flatAssembler.hh"]
#     includes += ["dune/iga/nurbsbasis.hh"]
#     moduleName = "SparseFlatAssembler_" + hashIt(element_type)
#     module = generator.load(
#         includes=includes, typeName=element_type, moduleName=moduleName
#     )
#     return module.SparseFlatAssembler(fes, dirichletValues)
# ---------------------------------------------------------------------------- #
# emDNA / src / DNA
#
# Nicolas Clauvelin (n.clauvelin@gmail.com)
#
# ---------------------------------------------------------------------------- #


set(DNA_SOURCES
    ${CMAKE_CURRENT_LIST_DIR}/AnchoredBpCollection.cpp
    ${CMAKE_CURRENT_LIST_DIR}/BpCollection.cpp
    ${CMAKE_CURRENT_LIST_DIR}/BpCollectionElasticEnergy.cpp
    ${CMAKE_CURRENT_LIST_DIR}/BpCollectionElectrostaticEnergy.cpp
    ${CMAKE_CURRENT_LIST_DIR}/BpCollectionFactory.cpp
    ${CMAKE_CURRENT_LIST_DIR}/BpCollectionGeometry.cpp
    ${CMAKE_CURRENT_LIST_DIR}/BpCollectionGradient.cpp
    ${CMAKE_CURRENT_LIST_DIR}/BpCollectionHessian.cpp
    ${CMAKE_CURRENT_LIST_DIR}/BpCollection_Interface.cpp
    ${CMAKE_CURRENT_LIST_DIR}/BpCollection_x3DNA.cpp
    ${CMAKE_CURRENT_LIST_DIR}/BpGeometryFunctions.cpp
    ${CMAKE_CURRENT_LIST_DIR}/BpStepDofs.cpp
    ${CMAKE_CURRENT_LIST_DIR}/BpStepParams.cpp
    ${CMAKE_CURRENT_LIST_DIR}/ClampedBpCollection.cpp
    ${CMAKE_CURRENT_LIST_DIR}/FreeBpCollection.cpp
    ${CMAKE_CURRENT_LIST_DIR}/HessianFunctions.cpp
    ${CMAKE_CURRENT_LIST_DIR}/IndexManager.cpp
    ${CMAKE_CURRENT_LIST_DIR}/IndexPair.cpp
    ${CMAKE_CURRENT_LIST_DIR}/PullingBpCollection.cpp
    ${CMAKE_CURRENT_LIST_DIR}/StepBlock.cpp
    ${CMAKE_CURRENT_LIST_DIR}/StepGradient.cpp
    PARENT_SCOPE)

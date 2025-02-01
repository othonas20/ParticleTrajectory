# Include header files
include_directories(${CMAKE_CURRENT_LIST_DIR}/FlowField)
include_directories(${CMAKE_CURRENT_LIST_DIR}/ODE)
include_directories(${CMAKE_CURRENT_LIST_DIR}/SolidWall)

list(
    APPEND LIB_SRC
    ${CMAKE_CURRENT_LIST_DIR}/FlowField/FlowField.cpp
    ${CMAKE_CURRENT_LIST_DIR}/FlowField/VortexFlowField.cpp
    ${CMAKE_CURRENT_LIST_DIR}/FlowField/FlowCylinder.cpp
    ${CMAKE_CURRENT_LIST_DIR}/FlowField/FlowPipeLaminar.cpp
    ${CMAKE_CURRENT_LIST_DIR}/FlowField/FlowWall.cpp
    ${CMAKE_CURRENT_LIST_DIR}/ODE/ODE.cpp
    ${CMAKE_CURRENT_LIST_DIR}/SolidWall/SolidWall.cpp
    ${CMAKE_CURRENT_LIST_DIR}/SolidWall/Cylinder.cpp
    ${CMAKE_CURRENT_LIST_DIR}/SolidWall/WallStraight.cpp
    ${CMAKE_CURRENT_LIST_DIR}/SolidWall/Pipe.cpp
)

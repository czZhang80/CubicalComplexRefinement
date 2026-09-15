# CMakeList.txt: CubicalComplexRefinement 的 CMake 项目，在此处包括源代码并定义
# 项目特定的逻辑。
#
cmake_minimum_required (VERSION 3.20)

# 如果支持，请为 MSVC 编译器启用热重载。
if (POLICY CMP0141)
  cmake_policy(SET CMP0141 NEW)
  set(CMAKE_MSVC_DEBUG_INFORMATION_FORMAT "$<IF:$<AND:$<C_COMPILER_ID:MSVC>,$<CXX_COMPILER_ID:MSVC>>,$<$<CONFIG:Debug,RelWithDebInfo>:EditAndContinue>,$<$<CONFIG:Debug,RelWithDebInfo>:ProgramDatabase>>")
endif()

project ("CubicalComplexRefinement")

find_package(Eigen3 3.3 REQUIRED NO_MODULE)

# 将源代码添加到此项目的可执行文件。
add_executable (CubicalComplexRefinement "CubicalComplexRefinement.cpp" "CubicalComplexRefinement.h" "main.cpp" "test.h"
  test.cpp)

target_link_libraries(CubicalComplexRefinement PRIVATE Eigen3::Eigen)

if (CMAKE_VERSION VERSION_GREATER 3.12)
  set_property(TARGET CubicalComplexRefinement PROPERTY CXX_STANDARD 20)
endif()

# TODO: 如有需要，请添加测试并安装目标。

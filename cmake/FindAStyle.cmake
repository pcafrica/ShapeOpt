# - Try to find AStyle
# Once done this will define
# ASTYLE_FOUND - System has AStyle libraries.
# ASTYLE_EXECUTABLE - The AStyle executable.
# ASTYLE_INCLUDE_DIRS - The AStyle include directories.
# ASTYLE_LIBRARIES - The libraries needed to use AStyle.

# find_path(ASTYLE_INCLUDE_DIR astyle.h)

# find_library(ASTYLE_LIBRARY NAMES astyle)

find_program(ASTYLE_EXECUTABLE NAMES astyle)

# set(ASTYLE_LIBRARIES ${ASTYLE_LIBRARY})
# set(ASTYLE_INCLUDE_DIRS ${ASTYLE_INCLUDE_DIR})

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set AStyle_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(AStyle DEFAULT_MSG
                                  # ASTYLE_LIBRARY ASTYLE_INCLUDE_DIR
                                  ASTYLE_EXECUTABLE)

mark_as_advanced( # ASTYLE_LIBRARY ASTYLE_INCLUDE_DIR
    ASTYLE_EXECUTABLE)

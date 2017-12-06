# Copyright 2017 Manuel Fasching <manuel.fasching@tum.de>
# Distributed under the MIT License
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is furnished
# to do so, subject to the following conditions:
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
# INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
# PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE
# FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
# ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.

if(NOT EXISTS "/home/hpc/pr63so/di29zux/ExaHyPE-Engine/3d/easi/libs/ImpalaJIT_build/install_manifest.txt")
    message(FATAL_ERROR "Cannot find install manifest: /home/hpc/pr63so/di29zux/ExaHyPE-Engine/3d/easi/libs/ImpalaJIT_build/install_manifest.txt")
endif(NOT EXISTS "/home/hpc/pr63so/di29zux/ExaHyPE-Engine/3d/easi/libs/ImpalaJIT_build/install_manifest.txt")

file(READ "/home/hpc/pr63so/di29zux/ExaHyPE-Engine/3d/easi/libs/ImpalaJIT_build/install_manifest.txt" files)
string(REGEX REPLACE "\n" ";" files "${files}")
foreach(file ${files})
    message(STATUS "Uninstalling $ENV{DESTDIR}${file}")
    if(IS_SYMLINK "$ENV{DESTDIR}${file}" OR EXISTS "$ENV{DESTDIR}${file}")
        exec_program(
                "/lrz/mnt/sys.x86_64/tools/cmake/3.4.0/bin/cmake" ARGS "-E remove \"$ENV{DESTDIR}${file}\""
                OUTPUT_VARIABLE rm_out
                RETURN_VALUE rm_retval
        )
        if(NOT "${rm_retval}" STREQUAL 0)
            message(FATAL_ERROR "Problem when removing $ENV{DESTDIR}${file}")
        endif(NOT "${rm_retval}" STREQUAL 0)
    else(IS_SYMLINK "$ENV{DESTDIR}${file}" OR EXISTS "$ENV{DESTDIR}${file}")
        message(STATUS "File $ENV{DESTDIR}${file} does not exist.")
    endif(IS_SYMLINK "$ENV{DESTDIR}${file}" OR EXISTS "$ENV{DESTDIR}${file}")
endforeach(file)

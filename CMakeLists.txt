cmake_minimum_required(VERSION 3.23)
project(Omega_K)

set(CMAKE_CXX_STANDARD 14)


set(SOURCE_LIB extra_tools.cpp compression_algorithm.cpp tester.h tester.cpp CubicSlpine.h)			# Тоже самое, но для библиотеки

set(SOURCE_EXE main.cpp)		# Установка переменной со списком исходников для исполняемого файла



add_library(extra STATIC ${SOURCE_LIB})	# Создание статической библиотеки с именем extra

add_executable(Omega_K ${SOURCE_EXE})	# Создает исполняемый файл с именем untitled

#target_link_libraries(untitled fftw3 extra)		# Линковка программы с библиотекой

find_library(FFTW_LIBRARY
        NAMES fftw3 fftw)
set(FFTW_LIBRARIES "${FFTW_LIBRARY}")
if(UNIX AND NOT WIN32)
    find_library(FFTW_libm_LIBRARY
            NAMES m)
    list(APPEND FFTW_LIBRARIES "${FFTW_libm_LIBRARY}")
endif()


target_link_libraries(extra ${FFTW_LIBRARIES})
target_link_libraries(Omega_K extra)
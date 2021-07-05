set(CMAKE_CXX_COMPILER "/vol6/software/mpi/mpi-intel2013/bin/mpicxx")
set(CMAKE_CXX_COMPILER_ARG1 "")
set(CMAKE_CXX_COMPILER_ID "Intel")
set(CMAKE_CXX_COMPILER_VERSION "13.0.0.20120731")
set(CMAKE_CXX_COMPILER_VERSION_INTERNAL "")
set(CMAKE_CXX_COMPILER_WRAPPER "")
set(CMAKE_CXX_STANDARD_COMPUTED_DEFAULT "98")
set(CMAKE_CXX_COMPILE_FEATURES "cxx_std_98;cxx_binary_literals;cxx_std_11;cxx_alias_templates;cxx_attributes;cxx_auto_type;cxx_decltype;cxx_default_function_template_args;cxx_defaulted_functions;cxx_deleted_functions;cxx_explicit_conversions;cxx_extended_friend_declarations;cxx_extern_templates;cxx_func_identifier;cxx_lambdas;cxx_local_type_template_args;cxx_long_long_type;cxx_nullptr;cxx_range_for;cxx_right_angle_brackets;cxx_rvalue_references;cxx_static_assert;cxx_template_template_parameters;cxx_trailing_return_types;cxx_uniform_initialization;cxx_variadic_macros;cxx_variadic_templates")
set(CMAKE_CXX98_COMPILE_FEATURES "cxx_std_98;cxx_binary_literals")
set(CMAKE_CXX11_COMPILE_FEATURES "cxx_std_11;cxx_alias_templates;cxx_attributes;cxx_auto_type;cxx_decltype;cxx_default_function_template_args;cxx_defaulted_functions;cxx_deleted_functions;cxx_explicit_conversions;cxx_extended_friend_declarations;cxx_extern_templates;cxx_func_identifier;cxx_lambdas;cxx_local_type_template_args;cxx_long_long_type;cxx_nullptr;cxx_range_for;cxx_right_angle_brackets;cxx_rvalue_references;cxx_static_assert;cxx_template_template_parameters;cxx_trailing_return_types;cxx_uniform_initialization;cxx_variadic_macros;cxx_variadic_templates")
set(CMAKE_CXX14_COMPILE_FEATURES "")
set(CMAKE_CXX17_COMPILE_FEATURES "")
set(CMAKE_CXX20_COMPILE_FEATURES "")

set(CMAKE_CXX_PLATFORM_ID "Linux")
set(CMAKE_CXX_SIMULATE_ID "")
set(CMAKE_CXX_SIMULATE_VERSION "")



set(CMAKE_AR "/usr/bin/ar")
set(CMAKE_CXX_COMPILER_AR "")
set(CMAKE_RANLIB "/usr/bin/ranlib")
set(CMAKE_CXX_COMPILER_RANLIB "")
set(CMAKE_LINKER "/usr/bin/ld")
set(CMAKE_MT "")
set(CMAKE_COMPILER_IS_GNUCXX )
set(CMAKE_CXX_COMPILER_LOADED 1)
set(CMAKE_CXX_COMPILER_WORKS TRUE)
set(CMAKE_CXX_ABI_COMPILED TRUE)
set(CMAKE_COMPILER_IS_MINGW )
set(CMAKE_COMPILER_IS_CYGWIN )
if(CMAKE_COMPILER_IS_CYGWIN)
  set(CYGWIN 1)
  set(UNIX 1)
endif()

set(CMAKE_CXX_COMPILER_ENV_VAR "CXX")

if(CMAKE_COMPILER_IS_MINGW)
  set(MINGW 1)
endif()
set(CMAKE_CXX_COMPILER_ID_RUN 1)
set(CMAKE_CXX_IGNORE_EXTENSIONS inl;h;hpp;HPP;H;o;O;obj;OBJ;def;DEF;rc;RC)
set(CMAKE_CXX_SOURCE_FILE_EXTENSIONS C;M;c++;cc;cpp;cxx;mm;CPP)
set(CMAKE_CXX_LINKER_PREFERENCE 30)
set(CMAKE_CXX_LINKER_PREFERENCE_PROPAGATES 1)

# Save compiler ABI information.
set(CMAKE_CXX_SIZEOF_DATA_PTR "8")
set(CMAKE_CXX_COMPILER_ABI "ELF")
set(CMAKE_CXX_LIBRARY_ARCHITECTURE "")

if(CMAKE_CXX_SIZEOF_DATA_PTR)
  set(CMAKE_SIZEOF_VOID_P "${CMAKE_CXX_SIZEOF_DATA_PTR}")
endif()

if(CMAKE_CXX_COMPILER_ABI)
  set(CMAKE_INTERNAL_PLATFORM_ABI "${CMAKE_CXX_COMPILER_ABI}")
endif()

if(CMAKE_CXX_LIBRARY_ARCHITECTURE)
  set(CMAKE_LIBRARY_ARCHITECTURE "")
endif()

set(CMAKE_CXX_CL_SHOWINCLUDES_PREFIX "")
if(CMAKE_CXX_CL_SHOWINCLUDES_PREFIX)
  set(CMAKE_CL_SHOWINCLUDES_PREFIX "${CMAKE_CXX_CL_SHOWINCLUDES_PREFIX}")
endif()





set(CMAKE_CXX_IMPLICIT_INCLUDE_DIRECTORIES "/usr/local/mpi-intel2013/include;/opt/intel/composer_xe_2013.0.079/ipp/include;/vol6/software/valgrind-3.10.0/include;/vol6/software/io_tools/netcdf/mpi/4.1.2/include;/vol6/software/io_tools/hdf5/mpi/1.8.11/include;/vol6/software/libraries/fftw/mpi/3.3.3/include;/opt/intel/composer_xe_2013.0.079/mkl/include;/opt/intel/composer_xe_2013.0.079/tbb/include;/opt/intel/composer_xe_2013.0.079/compiler/include/intel64;/opt/intel/composer_xe_2013.0.079/compiler/include;/usr/include/c++/4.4.7;/usr/include/c++/4.4.7/x86_64-kylin-linux;/usr/include/c++/4.4.7/backward;/usr/local/include;/usr/lib/gcc/x86_64-kylin-linux/4.4.7/include;/usr/include")
set(CMAKE_CXX_IMPLICIT_LINK_LIBRARIES "mpichcxx;mpich;pmi;opa;mpl;rt;glex;pmi;pthread;imf;svml;irng;m;ipgo;decimal;cilkrts;stdc++;gcc;gcc_s;irc;svml;c;gcc;gcc_s;irc_s;dl;c")
set(CMAKE_CXX_IMPLICIT_LINK_DIRECTORIES "/usr/local/glex-2.9b-cn-lib/lib64;/usr/local/glex-2.9b-cn-lib/lib;/usr/local/mpi-intel2013/lib;/vol6/software/lapack-3.4.2;/vol6/software/mpi/mpi-intel2013/lib;/vol6/software/BLAS;/vol6/software/io_tools/netcdf/mpi/4.1.2/lib;/vol6/software/io_tools/hdf5/mpi/1.8.11/lib;/vol6/software/libraries/fftw/mpi/3.3.3/lib;/opt/intel/composer_xe_2013.0.079/compiler/lib/intel64;/opt/intel/composer_xe_2013.0.079/mkl/lib/intel64;/opt/intel/composer_xe_2013.0.079/ipp/lib/intel64;/opt/intel/composer_xe_2013.0.079/tbb/lib/intel64;/usr/lib/gcc/x86_64-kylin-linux/4.4.7;/usr/lib64;/lib64;/usr/lib;/lib")
set(CMAKE_CXX_IMPLICIT_LINK_FRAMEWORK_DIRECTORIES "")

cmake_minimum_required (VERSION 2.8)
project (DEPENDENCIES)

include(ExternalProject)

ExternalProject_Add(
        EP_zlib
        URL http://zlib.net/zlib-1.2.5.tar.gz
        CONFIGURE_COMMAND <SOURCE_DIR>/configure --prefix=<INSTALL_DIR>
        BUILD_IN_SOURCE 1
	PREFIX ${DEPENDENCIES_SOURCE_DIR}
)

ExternalProject_Add(
        EP_protobuf
        URL http://protobuf.googlecode.com/files/protobuf-2.3.0.tar.gz
        CONFIGURE_COMMAND <SOURCE_DIR>/configure
        --prefix=<INSTALL_DIR>
        PREFIX ${DEPENDENCIES_SOURCE_DIR}
)

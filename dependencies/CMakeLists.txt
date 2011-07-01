# todo: will build command none work?

cmake_minimum_required (VERSION 2.8) 
project (DEPENDENCIES)

include(ExternalProject)

ExternalProject_Add(
	EP_sqlite3
	URL http://www.sqlite.org/sqlite-autoconf-3070603.tar.gz
	CONFIGURE_COMMAND <SOURCE_DIR>/configure --prefix=<INSTALL_DIR>
	PREFIX ${DEPENDENCIES_SOURCE_DIR}
)

ExternalProject_Add(
	EP_protobuf
	URL http://protobuf.googlecode.com/files/protobuf-2.3.0.tar.gz
	CONFIGURE_COMMAND <SOURCE_DIR>/configure --prefix=<INSTALL_DIR>
	PREFIX ${DEPENDENCIES_SOURCE_DIR}
)

ExternalProject_Add(
	EP_readline
	URL ftp://ftp.cwru.edu/pub/bash/readline-6.2.tar.gz
	CONFIGURE_COMMAND <SOURCE_DIR>/configure --prefix=<INSTALL_DIR>
	PREFIX ${DEPENDENCIES_SOURCE_DIR}
)


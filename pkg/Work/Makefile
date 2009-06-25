RUNIT_DIR=$(shell pwd)
PKG_PATH=$(shell dirname ${RUNIT_DIR})

SOURCES=$(wildcard ${PKG_PATH}/src/*.c)
OBJECTS=$(patsubst %.c,%.o,$(SOURCES))

UNAME=$(shell uname)
ifeq (${UNAME},Linux)
	MKVARS=${PKG_PATH}/src/Makevars
	TARGET=Rcplex.so
else
	MKVARS=${PKG_PATH}/src/Makevars.win
	TARGET=Rcplex.dll
endif

CONF=${PKG_PATH}/configure

# This has to change according to architecture
CONFIG_ARGS=-with-cplex-dir=/s/cplex/cplex90
#CONFIG_ARGS=-with-cplex-dir=/home/theussl/lib/CPLEX/cplex112

all: build test

build: ${TARGET}

${CONF}: ${PKG_PATH}/configure.ac
	cd ${PKG_PATH}; \
	autoconf
	cd ${RUNIT_DIR}

${MKVARS}: ${CONF}
	cd ${PKG_PATH}; \
	./configure ${CONFIG_ARGS}
	cd ${RUNIT_DIR}

${TARGET}: ${MKVARS} ${SOURCES}
	cd $(shell dirname $<); \
	R CMD SHLIB -o ${TARGET} $(foreach SRC,${SOURCES},$(shell basename ${SRC}))
	cd ${RUNIT_DIR}

#%.o: %.c

cleanso:
	rm -f ${PKG_PATH}/src/*.so ${PKG_PATH}/src/*.o

cleanconf:
	rm -rf ../config.log ../config.status ../autom4te.cache ../src/Makevars

clean: cleanso cleanconf

test:
	R --slave < runtests.R

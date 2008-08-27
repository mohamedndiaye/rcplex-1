#!/bin/bash

TEST_FILE=$1

if test $# -ne 1 ; then
    echo "Usage: runonefile.sh SOMETEST.R"
    exit 1
fi

R_CODE="library(RUnit)"
R_CODE="${R_CODE};source(\"setuptests.R\")"
R_CODE="${R_CODE};res <- runTestFile('${TEST_FILE}',rngKind='default',rngNormalKind='default')"
R_CODE="${R_CODE};printTextProtocol(res)"

echo ${R_CODE} | R --slave

library(RUnit)

# load R files and shared object
source("setuptests.R")

suite <- defineTestSuite("Rcplex Unit Tests",
                        dirs=c("."),
                        testFileRegexp="^runit.+\\.R$",
                        testFuncRegexp="^test.+",
                        rngKind="default",
                        rngNormalKind="default")
res <- runTestSuite(suite)
printTextProtocol(res,file="rcplex.test.out")
warnings()
Rcplex.close()
q()

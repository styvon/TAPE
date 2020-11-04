library(TAPE)
bfile = strsplit(system.file("extdata", "toy.bim", package = "TAPE"), split=".bim")[[1]]
print(bfile)
test_that("linkplink:setgeno", {
  expect_equal(setgeno(bfile, 1:10000, memoryChunk=1,isDiagofKinSetAsOne=T),10000)
})
closegeno()


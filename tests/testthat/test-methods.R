test_that("Testing getPConsMtx function", {
    input1 <-list(c(1L,1L))
    output <- matrix(c(0,1,1,0), nrow = 2, ncol = 2)
    expect_equal(getPConsMtx(input1, .5), output)
    
    input2 <-list(c(1L,1L), c(1L,1L))
    expect_equal(getPConsMtx(input2, .5), output)
})
test_that("Testing findNDim function", {
    expect_equal(findNDim(100, NULL, NULL, NULL), seq(4,7))
    expect_equal(findNDim(100, NULL, 0.01, 0.02), seq(1,2))
})
test_that("Testing findDropouts function", {
    expect_equal(findDropouts(matrix(1, 4, 4) - diag(4), matrix(1, 4, 4)), diag(4)==TRUE)
})
# Test

test_that("Test motif empty", {
  seq = "AAATTTCCC"
  code = codes.code(c("GGG"))
  ml = ccmotif.lengths(seq, code)
  expect_equal(length(ml), 2)
  expect_equal(length(ml$incode), 0)
  expect_equal(length(ml$outcode), 1)
  expect_equal(ml$outcode[1], 3)
})

test_that("Test motif 1", {
  seq = "AAATTTCCC"
  code = codes.code(c("AAA"))
  ml = ccmotif.lengths(seq, code)
  expect_equal(length(ml), 2)
  expect_equal(length(ml$incode), 1)
  expect_equal(length(ml$outcode), 1)
  expect_equal(ml$incode[1], 1)
  expect_equal(ml$outcode[1], 2)
})

test_that("Test motif 2", {
  seq = "AAATTTAAACCC"
  code = codes.code(c("AAA"))
  ml = ccmotif.lengths(seq, code)
  expect_equal(length(ml), 2)
  expect_equal(length(ml$incode), 2)
  expect_equal(length(ml$outcode), 2)
  expect_equal(ml$incode[1], 1)
  expect_equal(ml$incode[2], 1)
  expect_equal(ml$outcode[1], 1)
  expect_equal(ml$outcode[2], 1)
})

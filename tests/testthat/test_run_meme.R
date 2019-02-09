context("run_meme()")

test_that("run_meme() gives correct error without bin", {

  skip_on_os("windows")

  expect_error(run_meme(bin = "__not__meme__"),
               "could not find the MEME binary")

})

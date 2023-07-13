# Just checking the code run
expect_silent(score.split(
  n1 = 5,
  n2 = 5,
  W1area = 5,
  W2area = 5,
  score = "lcv"
))

expect_silent(score.split(
  n1 = 5,
  n2 = 5,
  W1area = 5,
  W2area = 5,
  score = "lcv2"
))

expect_silent(score.pp(
  X = spatstat.random::rpoispp(10)
))

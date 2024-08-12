data.GaussKronrod <- function(data.id, a, b, k = 15){

  wk <- gaussKronrod()$wk
  sk <- gaussKronrod()$sk
  K <- length(sk)
  P <- (b-a)/2
  st <- outer(P, sk + 1)+a
  id.GK <- rep(seq_along(data.id$id), each = K)
  data.id2 <- data.id[id.GK, ]

  list(K = K, P = P, st = st, wk = wk, data.id2 = data.id2,
       id.GK = id.GK)

}

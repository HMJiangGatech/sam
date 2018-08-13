# SAM

library(SAM)

total_t = 0
total_l = 0
nlamb = 30
for (i in 1:t) {
  t0 = proc.time()
  out.trn = samQL(X, y, nlambda=nlamb)
  total_t = total_t + proc.time() - t0
  out.tst = predict(out.trn, Xt)
  total_l = total_l + mean((out.tst$values[,nlamb]-yt)^2/2)
}
print("sam lin-reg:")
print(total_t / t)
print(total_l / t)


total_t = 0
total_l = 0
nlamb = 30
for (i in 1:t) {
  t0 = proc.time()
  out.trn = samQL(X, y, nlambda=nlamb, regfunc="MCP")
  total_t = total_t + proc.time() - t0
  out.tst = predict(out.trn, Xt)
  total_l = total_l + mean((out.tst$values[,nlamb]-yt)^2/2)
}
print("sam lin-reg with MCP:")
print(total_t / t)
print(total_l / t)
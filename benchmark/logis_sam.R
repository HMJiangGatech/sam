# SAM

library(SAM)

total_t = 0
total_l = 0
nlamb = 20
for (i in 1:t) {
  t0 = proc.time()
  out.trn = samLL(X, y, nlambda=nlamb)
  total_t = total_t + proc.time() - t0
  out.tst = predict(out.trn, Xt)
  total_l = total_l + mean(out.tst$labels[,nlamb]==yt)
}
print("sam log-reg:")
print(total_t / t - genZ_t)
print(total_l / t)

total_t = 0
total_l = 0
nlamb = 20
for (i in 1:t) {
  t0 = proc.time()
  out.trn = samLL(X, y, nlambda=nlamb, regfunc="MCP")
  total_t = total_t + proc.time() - t0
  out.tst = predict(out.trn, Xt)
  total_l = total_l + mean(out.tst$labels[,nlamb]==yt)
}
print("sam log-reg with MCP:")
print(total_t / t - genZ_t)
print(total_l / t)
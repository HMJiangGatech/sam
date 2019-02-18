# grpreg

library(grpreg)

total_t = 0
total_l = 0
nlamb = 120
for (i in 1:t) {
  t0 = proc.time()
  out.trn = grpreg(Z[,2:ncol(Z)], y, group=index[2:length(index)], nlambda=nlamb, eps=1e-7)
  total_t = total_t + proc.time() - t0
  out.tst = predict(out.trn, Zt[,2:ncol(Zt)])
  total_l = total_l + mean((out.tst[,nlamb]-yt)^2/2)
}
print("grpreg lin-reg:")
print(total_t / t)
print(total_l / t)


total_t = 0
total_l = 0
nlamb = 120
for (i in 1:t) {
  t0 = proc.time()
  out.trn = grpreg(Z[,2:ncol(Z)], y, group=index[2:length(index)], penalty="grMCP", nlambda=nlamb, eps=1e-7)
  total_t = total_t + proc.time() - t0
  out.tst = predict(out.trn, Zt[,2:ncol(Zt)])
  total_l = total_l + mean((out.tst[,nlamb]-yt)^2/2)
}
print("grpreg lin-reg with MCP:")
print(total_t / t)
print(total_l / t)
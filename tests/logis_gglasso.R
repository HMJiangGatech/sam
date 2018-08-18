# gglasso

library(gglasso)


total_t = 0
total_l = 0
nlamb = 100
for (i in 1:t) {
  t0 = proc.time()
  out.trn = gglasso(Z[,2:ncol(Z)], y*2-1, group=rep(1:d,each=p), loss="logit", nlambda=nlamb)
  total_t = total_t + proc.time() - t0
  out.tst = predict(out.trn, Zt[,2:ncol(Zt)])
  total_l = total_l + mean(out.tst[,nlamb]==(yt*2-1))
}
print("gglasso log-reg:")
print(total_t / t)
print(total_l / t)
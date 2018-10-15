# grplasso
nlab = 10

library(grplasso)

total_t = 0
total_l = 0
for (i in 1:t) {
  t0 = proc.time()
  lambda <- lambdamax(Z, y = y, index = index, penscale = sqrt,
                      model = LogReg()) * 0.5^(0:5)
  out.trn <- grplasso(Z, y = y, index = index, lambda = lambda, model = LogReg(),
                      penscale = sqrt,
                      control = grpl.control(update.hess = "lambda", trace = 0))
  total_t = total_t + proc.time() - t0
  out.tst = predict(out.trn, Zt, type="response")
  total_l = total_l + mean((out.tst[,6]>0.5)==yt)
}



print("grplasso log-reg:")
print(total_t / t)
print(total_l / t)
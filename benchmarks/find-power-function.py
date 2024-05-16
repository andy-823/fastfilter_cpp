import numpy as np
import math as m

file = open('input.txt')
file_data = file.readlines()
file.close()

xvalues = []
yvalues = []
for row in file_data:
  a, b = row.split()
  xvalues.append(a)
  yvalues.append(b)

xvalues = np.array(xvalues, dtype=float)
yvalues = np.array(yvalues, dtype=float)

best_l = []
best_pow = []
sqr_errors = []

# print(xvalues, yvalues)

step = 0.01
cur_pow = 0
while cur_pow <= 1:
  l = 0
  r = 100
  while r - l > 0.001:
    mid = (l + r) / 2
    is_good = True
    for i in range(0, len(xvalues)):
      n = xvalues[i]
      got = mid * m.pow(n, cur_pow)
      expected = yvalues[i] * 0.85
      is_good = is_good and expected >= got
    if is_good:
      l = mid
    else:
      r = mid
  
  sqr_error = 0
  for i in range(0, len(xvalues)):
    n = xvalues[i]
    got = mid * m.pow(n, cur_pow)
    expected = yvalues[i] * 0.85
    sqr_error += abs(expected - got)**2
  best_l.append(l)
  best_pow.append(cur_pow)
  sqr_errors.append(sqr_error)
  cur_pow += step

min_error = 1e9
min_pow = -1
min_l = -1
for i in range(0, len(best_l)):
  if sqr_errors[i] < min_error:
    min_error = sqr_errors[i]
    min_pow = best_pow[i]
    min_l = best_l[i]

print(min_pow, min_l)
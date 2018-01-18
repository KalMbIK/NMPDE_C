import matplotlib.pyplot as plt
import csv

Y = []
with open('cmake-build-debug/test.csv', 'rb') as f:
    reader = csv.reader(f)
    for row in reader:
        Y.append(map(float,row))
plt.plot(Y[0], Y[1], label='density')
plt.plot(Y[0], Y[2], label='velocity')
plt.plot(Y[0], Y[3], label='pressure')
plt.legend(loc='best')
plt.show()

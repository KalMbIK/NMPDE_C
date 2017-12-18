import matplotlib.pyplot as plt
import csv

Y = []
with open('cmake-build-debug/test.csv', 'rb') as f:
    reader = csv.reader(f)
    for row in reader:
        Y.append(map(float,row))
plt.plot(Y[0],Y[1])
plt.plot(Y[0],Y[2])
plt.show()

import matplotlib.pyplot as plt
import csv

Y = []
with open('cmake-build-debug/test.csv', 'rb') as f:
    reader = csv.reader(f)
    for row in reader:
        Y.append(map(float,row))
plt.plot(Y)
plt.show()

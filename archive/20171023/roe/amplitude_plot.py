import matplotlib.pyplot as plt
import numpy as np

data = np.loadtxt("sine_max.txt")

x,amp = data[0:230,0], data[0:230,1]

plt.plot(x,amp)
plt.xlabel("Time [s]")
plt.ylabel("Amplitude [kg/m^3]")
plt.show()
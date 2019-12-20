import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib import style
import numpy as np

style.use('fivethirtyeight')

fig = plt.figure()
ax1 = fig.add_subplot(1,1,1)

def animate(i):
    filename = '../data/int1/quantumdot/energy/RBMPJ/2D/2P/1.000000w/ADAM_MC524288.dat'
    data = np.loadtxt(filename)
    ax1.plot(data, label="Simulated")
    
ani = animation.FuncAnimation(fig, animate, interval=1000)

plt.axhline(3, label="Exact", color='r' , linestyle='--')
plt.xlabel("Iteration")
plt.ylabel("Energy")
plt.legend(loc='best')
plt.show()

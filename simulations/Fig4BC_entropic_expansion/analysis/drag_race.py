import matplotlib.pyplot as plt
import numpy as np

L = np.linspace(1,10)

r= 0.0125
drag_MT = 3 * L / (np.log(0.5 * L / r) + 0.312) * np.pi
for d0 in [0.04]:
    term = d0-r
    print term
    hydro_drag = 2*np.pi*r/term
    plt.axhline(y=hydro_drag)


plt.plot(L,drag_MT/L)
plt.ylim(ymin=0)
plt.show()
import matplotlib.pylab as plt
import numpy as np
plt.close('all')

N = 25
y = np.random.randn(N)
x = np.arange(N)

y2 = np.random.randn(25)

# serie A
p1a, = plt.plot(x, y,       "ro", ms=10, mfc="r", mew=2, mec="r")
p1b, = plt.plot(x[:5], y[:5] ,  "w+", ms=10, mec="w", mew=2)
p1c, = plt.plot(x[5:10], y[5:10], "w*", ms=10, mec="w", mew=2)

# serie B
p2a, = plt.plot(x, y2,       "bo", ms=10, mfc="b", mew=2, mec="b")
p2b, = plt.plot(x[15:20], y2[15:20] ,  "w+", ms=10, mec="w", mew=2)
p2c, = plt.plot(x[10:15], y2[10:15], "w*", ms=10, mec="w", mew=2)

line_columns = [
                p1a, p2a,
                (p1a, p1b), (p2a, p2b),
                (p1a, p1c), (p2a, p2c)
                ]


plt.legend(line_columns, ['a']*2 + ['']*(len(line_columns)-2),
             title=' No Prop    Prop +    Prop *',
             ncol=3, numpoints=1, handletextpad=1.5,markerfirst=False)


# plt.gca().add_artist(leg2)
# plt.gca().add_artist(leg3)

plt.show()
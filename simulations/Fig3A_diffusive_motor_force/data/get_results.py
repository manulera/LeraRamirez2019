import numpy as np

with open('single_force.txt','r') as ins:

    points = [[], []]
    for line in ins:
        ls = line.split()
        if len(ls) and ls[0]!="%":
            ind = int(ls[1]) - 1
            points[ind].append(np.array(ls[3], dtype=float))

    dist = np.array(points[0])-np.array(points[1])
    out = str()
    for i in dist:
        out += str(i) + " "

    print out[:-1]
import numpy as np

with open('fiber.txt') as ins:

    points = [[], []]
    for line in ins:
        ls = line.split()
        if len(ls) and ls[0]!="%":
            ind = int(ls[1]) - 1
            points[ind].append(np.array(ls[4], dtype=float))

    if not len(points[0])!=len(points[1]):
        dist = np.abs(np.array(points[0])-np.array(points[1]))
    else:
        dist = np.zeros(101)

    out = str()
    for i in dist:
        out += str(i) + " "

    print out[:-1]
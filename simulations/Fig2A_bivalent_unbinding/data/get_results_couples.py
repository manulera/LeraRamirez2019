
from numpy import mean
with open('couples.txt') as ins:

    ase = list()
    mot = list()
    for line in ins:
        ls = line.split()
        if len(ls) and ls[0]!="%":
            if ls[0]=="cross_couple":
                ase.append(float(ls[-1]))
            if ls[0]=="motor_couple":
                mot.append(float(ls[-1]))
    if len(ase):
        print mean(ase), mean(mot)
    else:
        print "nan", "nan"
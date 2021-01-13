with open('couples.txt') as ins:

    count=0
    line = ins.readline()
    while count<3:
        if "frame" in line:
            count+=1
        line = ins.readline()

    out = list()
    for line in ins:
        ls = line.split()
        if len(ls)>3 and ls[0]!="%":
            out.append(ls[6])

    print " ".join(out)
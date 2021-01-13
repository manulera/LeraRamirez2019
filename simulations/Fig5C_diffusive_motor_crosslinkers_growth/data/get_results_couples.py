

with open('couples.txt') as ins:

    out = list()
    for line in ins:
        if len(line)>3 and line[0]!="%":
            ls = line.split()
            out.append(ls[-1])

    print " ".join(out)
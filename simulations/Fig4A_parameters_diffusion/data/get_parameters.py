#!/usr/bin/env python2
# Meant to read parameters from the simulation.


from leraramirez2019 import read_config

def get_vals(pile,extract):
    out = list()
    for i in range(0, len(extract), 3):
        com = extract[i]
        what = extract[i + 1]
        obj = read_config.get_command(pile, com)

        if len(what)==1:
            out.append(obj.value(*what))
        elif len(what)>1:
            out.append(obj[what[1]].value(what[0]))
        else:
            out.append(obj.cnt)
    return out

pile = read_config.parse("properties.cmo")
par_names = []
extract =[
    ['set', 'couple', 'couple1'],  ["stiffness"],          "stifness",
    ['set', 'hand', 'motor'],  ["diffusion"],          "diffusion",
]
par_names+=extract[2::3]
out_list = list()
out_list += get_vals(pile,extract)

pile = read_config.parse("config.cym")
extract =[['new', 'fiber', 'fiber1'],  ["length",1],          "length",]
par_names+=extract[2::3]
out_list += get_vals(pile,extract)
print "|".join(par_names)
print "|".join(map(str,out_list))
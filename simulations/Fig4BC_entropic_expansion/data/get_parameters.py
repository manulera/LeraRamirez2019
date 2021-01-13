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
    ['set', 'hand', 'cross_hand'],  ["diffusion"],          "cl_diff",
    ['set', 'couple', 'cross_couple'],  ["stiffness"],          "stiffness",
]
par_names+=extract[2::3]
out_list = list()
out_list += get_vals(pile,extract)

pile = read_config.parse("config.cym")

extract = [
    ['set', 'simul', '*'],          ["viscosity"],          "viscosity",
    ['new','couple', 'cross_couple'],   [],    "cl_number",
    ['change', 'simul', '*'], ["time_step"], "time_step",
    ['run', 'simul', '*'], ["nb_steps", 1], "nb_steps",
    ['run', 'simul', '*'], ["nb_frames", 1], "nb_frames",
]
par_names+=extract[2::3]
out_list += get_vals(pile,extract)


with open("config.cym") as ins:
    for line in ins:
        if "occup" in line:
            occup = line.split()[3]
            break
        if "overlap_length" in line:
            overlap_length = line.split()[3]
par_names.append("occup")
out_list.append(occup)
par_names.append("overlap_length")
out_list.append(overlap_length)
print "|".join(par_names)
print "|".join(map(str,out_list))
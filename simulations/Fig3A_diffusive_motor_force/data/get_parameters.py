#!/usr/bin/env python2
# Meant to read parameters from the simulation.


from leraramirez2019 import read_config

def get_vals(pile,extract):
    out = list()
    for i in range(0, len(extract), 3):
        com = extract[i]
        what = extract[i + 1]
        obj = read_config.get_command(pile, com)
        if len(what):
            out.append(obj.value(*what))
        else:
            out.append(obj.cnt)
    return out

pile = read_config.parse("properties.cmo")
par_names = []
extract =[
    ['set', 'hand', 'motor_hand'],  ["stall_force"],        "stall_force",
    ['set', 'hand', 'motor_hand'],  ["unloaded_speed"],     "mt_speed",
    ['set', 'hand', 'cross_hand'],  ["diffusion"],          "cl_diff",
    ['set', 'single', 'anchor_single'],  ["stiffness"],          "anchor_stiffness",
    ['set', 'couple', 'motor_couple'],  ["stiffness"],          "couple_stiffness",
]
par_names+=extract[2::3]
out_list = list()
out_list += get_vals(pile,extract)

pile = read_config.parse("config.cym")

extract = [
    ['set', 'simul', '*'],          ["viscosity"],          "viscosity",
    ['new','couple', 'motor_couple'],   [],    "mt_number",
]
par_names+=extract[2::3]
out_list += get_vals(pile,extract)

print "|".join(par_names)
print "|".join(map(str,out_list))
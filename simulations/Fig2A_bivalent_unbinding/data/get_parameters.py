
# Meant to read parameters from the simulation.


from leraramirez2019 import read_config



pile = read_config.parse("properties.cmo")
par_names = []
extract =[
    ['set', 'hand', 'motor_hand'],  ["stall_force"],        "stall_force",
    ['set', 'hand', 'motor_hand'],  ["unloaded_speed"],     "mt_speed",
    ['set', 'hand', 'motor_hand'],  ["binding"],          "mt_binding",
    ['set', 'hand', 'motor_hand'],  ["unbinding"],        "mt_unbinding",
    ['set', 'hand', 'cross_hand'],  ["binding"],         "cl_binding",
    ['set', 'hand', 'cross_hand'],  ["unbinding"],       "cl_unbinding",
    ['set', 'couple', 'cross_couple'],  ["stiffness"],         "cl_stiffness",
    ['set','simul', '*'],   ["time_step"],    "time_step",
    ['set', 'simul', '*'],          ["viscosity"],          "viscosity",
]
par_names+=extract[2::3]
out_list = list()
out_list += read_config.get_vals(pile,extract)

pile = read_config.parse("config.cym")

extract = [
    ['new','couple', 'cross_couple'],   [],    "cl_number",
    ['new','couple', 'motor_couple'],   [],    "mt_number",
    ['run','simul', '*'],   ["nb_steps",1],    "nb_steps",
    ['run','simul', '*'],   ["nb_frames",1],    "nb_frames",
]
par_names+=extract[2::3]
out_list += read_config.get_vals(pile,extract)

print "|".join(par_names)
print "|".join(map(str,out_list))
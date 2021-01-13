% Cytosim config file
% Manuel Lera Ramirez - 01 December, 2017

[[diffusion = [0.011 + (0.085- 0.011)/49.*i for i in range(50) ]]]
[[stiffness = [10. + (300.- 10.)/49.*i for i in range(50) ]]]


set simul ase1
{
time_step = 0.00001
viscosity = 0.05
}

set space cell
{
geometry = ( sphere 1 )
}

new space cell

set fiber fiber1
{
rigidity = 30
segmentation = 1
lattice = 1, 0.008
}

set hand motor
{
binding_rate = inf
binding_range = 0.1
unbinding_rate = 0
unbinding_force = inf

is_gillespie = 0

activity = ase_walk
use_lattice = 0
step_size = 0.008
dangling_chance = 1
diffusion = [[diffusion]];
}


set couple couple1
{
hand1 = motor;
hand2 = motor;
stiffness = [[stiffness]];
}

new 1 fiber fiber1
{
length = 10
orientation = ( 1 0 0 )
position = ( 0 0 0 )
}
new 1 fiber fiber1
{
length = 10
orientation = ( 1 0 0 )
position = ( 0 0 0 )
}

new 1000 single couple1
{
	attach1 = fiber1, 0, center
	attach2 = fiber2, 0, center
}

run simul *
{
nb_steps  = 100000
nb_frames = 2
solve = 0
}


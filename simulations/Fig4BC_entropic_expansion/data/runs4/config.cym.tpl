% Cytosim config file
% Manuel Lera Ramirez - 01 December, 2017

%% Parameters

[[cross_diff = 0.043]]
[[stiffness = 250.]]
[[occup = [0.1]]]


[[viscosity = [0.01]]]
[[mt_length = 20.]]
[[overlap_length = random.random()*14.5+0.5]]

[[nb_cross = int(occup*overlap_length/0.008)]]

set simul mysim
{
    time_step = 0.001
    viscosity = [[viscosity]]
}

set space bind
{
	geometry = (circle [[overlap_length/2.]])
}

set space move
{
	geometry = (circle [[mt_length]])
}

new space move
new space bind

set fiber mt
{
    rigidity = 30
    segmentation = 1
    lattice = 1, 0.008
	display = (explode=1,1)
}
set fiber mt2
{
    rigidity = 30
    segmentation = 1
    lattice = 1, 0.008
	display = (explode=1,1)
	viscosity=10000
}

%% cross 1 -----------------------------------------------------------------
set hand cross_hand
{
	binding_rate = inf
	unbinding_rate = 0
    binding_range = 0.05 % Has to be bigger than the length
    unbinding_force = inf

	is_gillespie = 0
    hold_growing_end = 0
    activity = ase_walk
    use_lattice = 1
    step_size = 0.008
    dangling_chance = 1
	diffusion = [[cross_diff]];

    display = ( color=purple; )

}

set couple cross_couple
{
	activity=diffuse
	hand1 = cross_hand;
	hand2 = cross_hand;
	stiffness = [[stiffness]];
	diffusion = 100;
	confine = inside, , bind
}


%% Creating fibers

new 1 fiber mt
{
    position_ends = [[-overlap_length/2.]],[[mt_length-overlap_length/2.]]
}
new 1 fiber mt2
{
    position_ends = [[+overlap_length/2.]],[[-mt_length+overlap_length/2.]]
}



new [[nb_cross]] couple cross_couple

% Equilibration step
run 10000 simul
{
    solve = 0
}

change couple cross_couple
{
    confine = inside, , move
}
change simul mysim
{
    time_step = 0.00001
}


run simul *
{
    nb_steps  = 1500000
    nb_frames = 100
    solve = 1
}



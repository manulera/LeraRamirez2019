% Cytosim config file
% Manuel Lera Ramirez - 01 December, 2017

%% Parameters

[[nb_motor = 35]]
[[cross_diff = 10**-(random.random()*7)]]
[[motor_speed = [0.2]]]
[[mt_length = 6.]]

set simul mysim 
{
    time_step = 0.00001
    viscosity = 0.1
}

set space grow
{
	geometry = (circle [[mt_length/2]])
}


new space grow

set fiber mt
{
    rigidity = 30
    segmentation = 1
    lattice = 1, 0.008
	activity =grow
    growing_speed    = [[motor_speed*2]] % We force the speed to keep the overlap length constant this value is irrelevant
    growing_force    = inf
    min_length       = 0.008
	growth_space     = grow
	display = (explode=1,1)
}


%% motor -----------------------------------------------------------------
set hand motor_hand
{	
	% To reach equilibrium fast
    binding = 1, 0.05
    unbinding = 0.05, inf
    
	stall_force = 6
    unloaded_speed  = [[motor_speed]]
    activity = move
    display = ( color=green; )
}

%% cross 1 -----------------------------------------------------------------
set hand cross_hand
{
	binding_rate = 1
	unbinding_rate = 0.05
    binding_range = 0.05 % Has to be bigger than the length
    unbinding_force = inf
	
	is_gillespie = 0
    hold_growing_end = 0
    activity = ase_walk
    use_lattice = 0
    step_size = 0.008
    dangling_chance = 0
	diffusion = [[cross_diff]];
    
    display = ( color=purple; )

}
% cross 2
set hand anchor_hand
{
	binding = inf, 0.1
	unbinding = 0, inf
}

set single anchor_single
{
	hand = anchor_hand
	activity = fixed
	% Stiffness should scale with the force, otherwise the microtubule slides too much and the 
	% setup is not comparable
	stiffness = [[40./(cross_diff/0.0042/motor_speed+1./6.)/(0.008/motor_speed)]]
}

set couple motor_couple
{
	activity=diffuse
	hand1 = motor_hand;
	hand2 = cross_hand;
	stiffness = 100;
	diffusion = 100;
}

%% Creating fibers

new 1 fiber mt
{
    length = [[mt_length]]
    orientation = ( 1 0 0 )
    position = (0 0 0)
}
new 1 fiber mt
{
    length = [[mt_length]]
    orientation = ( -1 0 0 )
    position = (0 0 0)
}

new [[nb_motor]] couple motor_couple


new single anchor_single
{
	position = 0
	attach = fiber1, 0, center 
}
new single anchor_single
{
	position = 0
	attach = fiber2, 0, center 
}

run simul *
{
    nb_steps  = 4000000
    nb_frames = 100
    solve = 1
}


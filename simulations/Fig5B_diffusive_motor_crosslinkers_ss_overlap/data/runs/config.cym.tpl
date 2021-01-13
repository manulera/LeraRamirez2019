% Cytosim config file
% Manuel Lera Ramirez - 01 December, 2017

%% Parameters

[[nb_motor = [100,200,300]]]
[[nb_cross = random.randint(1,200)]]
[[cross_diff = [0.1]]]
[[mot_diff = [0.1]]]
[[viscosity = [0.1]]]
[[motor_speed = [0.2]]]
[[mt_length = 3.]]

set simul mysim
{
    time_step = 0.001
    viscosity = [[viscosity]]
}

set space grow
{
	geometry = (circle [[mt_length/2.]])
}

set space contain
{
	geometry = (rectangle [[mt_length]])
}

new space grow
new space contain

set fiber mt
{
    rigidity = 30
    segmentation = 1
    lattice = 1, 0.008
	display = (explode=1,1)
	confine = inside, 100, contain
}


%% motor -----------------------------------------------------------------
set hand motor_hand
{
	% To reach equilibrium fast
    binding = 1, 0.05
    unbinding = 1, inf

	stall_force = 6
    unloaded_speed  = [[motor_speed]]
    activity = move
    display = ( color=green; )
}
set hand motor_hand_diff
{
	binding_rate = 1
	unbinding_rate = 1
    binding_range = 0.05 % Has to be bigger than the length
    unbinding_force = inf

	is_gillespie = 0
    hold_growing_end = 0
    activity = ase_walk
    use_lattice = 0
    step_size = 0.008
    dangling_chance = 0
	diffusion = [[mot_diff]];

    display = ( color=yellow; )

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

set couple motor_couple
{
	activity=diffuse
	hand1 = motor_hand;
	hand2 = motor_hand_diff;
	stiffness = 100;
	diffusion = 100;
	confine = inside, ,grow
}
set couple cross_couple
{
	activity=diffuse
	hand1 = cross_hand;
	hand2 = cross_hand;
	stiffness = 100;
	diffusion = 100;
	confine = inside, ,grow
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
new [[nb_cross]] couple cross_couple



% Equilibration step for the crosslinkers
run 1000 simul
{
    solve = 0
}

% Equilibration step for the motors

change couple motor_couple
{
	confine = inside, ,contain
}
change couple cross_couple
{
	confine = inside, ,contain
}
run 1000 simul
{
    solve = 0
}

change simul mysim
{
    time_step = 0.00001
}

% We run the simulation
run simul *
{
    nb_steps  = 20000000
    nb_frames = 100
    solve = 1
}


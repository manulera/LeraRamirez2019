% Cytosim config file
% Manuel Lera Ramirez - 01 December, 2017

%% Parameters

[[nb_cross = 150]]
[[nb_motor = [100]]]
[[cross_diff = [0.1]]]
[[motor_speed = [0.05]]]
[[growth_speed = random.random()*motor_speed*0.8]]
[[mt_length = 3.]]
set simul mysim 
{
    time_step = 0.001
    viscosity = 1
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
    activity = grow
    growing_speed    = 0
    growing_force    = inf
    min_length       = 0.008
	display = (explode=1,1)
}

%% motor -----------------------------------------------------------------
set hand motor_hand
{	
	% To reach equilibrium fast
    binding = 1, 0.045
    unbinding = 0.1, inf
    
	stall_force = 6
    unloaded_speed  = 0
    activity = move
    display = ( color=green; )
}

set couple motor_couple
{
	activity=diffuse
    hand1 = motor_hand
    hand2 = motor_hand
    diffusion = 100
    stiffness = 100
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
	% No way to reach high occupancy otherwise
	diffusion = 0.01;
    
    display = ( color=purple; )

}

set couple cross_couple
{
	activity=diffuse
	hand1 = cross_hand;
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

new [[nb_cross]] couple cross_couple

% Equilibration step
run 10000 simul
{
    solve = 0
}


delete couple cross_couple { state1 = 0; }

delete couple cross_couple { state2 = 0; }


change hand motor_hand
{
	unloaded_speed = [[motor_speed]]
}
change hand cross_hand
{
    diffusion = [[cross_diff]]
}
change simul mysim 
{
    time_step = 0.00001
}
change fiber mt
{
    growing_speed = [[growth_speed]]
}

run simul *
{
    nb_steps  = 10000000
    nb_frames = 100
    solve = 1
}


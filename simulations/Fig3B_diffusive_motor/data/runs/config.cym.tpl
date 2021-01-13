% Cytosim config file
% Manuel Lera Ramirez - 01 December, 2017

%% Parameters

[[nb_motor = random.randint(1,300)]]
[[cross_diff = [0.1]]]
[[viscosity = [1,0.1,0.01]]]
[[motor_speed = [0.2]]]
[[mt_length = 3.]]

set simul mysim 
{
    time_step = 0.0001
    viscosity = [[viscosity]]
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
	display = (explode=1,1)
	activity = grow
	growing_speed    = [[motor_speed*2]] % We force the speed to keep the overlap length constant this value is irrelevant
    growing_force    = inf
    min_length       = 0.008
	growth_space     = grow
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

%% cross 1 -----------------------------------------------------------------
set hand cross_hand
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
	diffusion = [[cross_diff]];
    
    display = ( color=purple; )

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



% Equilibration step
run 1000 simul
{
    solve = 0
}


change simul mysim
{
    time_step = 0.00001
}

run simul *
{
    nb_steps  = 4000000
    nb_frames = 100
    solve = 1
}


% Cytosim config file
% Manuel Lera Ramirez - 01 December, 2017

%% Parameters

[[nb_cross = random.randint(1,375)]]
[[nb_motor = [30,60,90]]]
[[motor_speed = [0.05]]]
[[mt_length = 3.]]

set simul mysim
{
    time_step = 0.01
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
    activity = grow
    growing_speed    = [[motor_speed*2]] % We force the speed to keep the overlap length constant this value is irrelevant
    growing_force    = inf
    min_length       = 0.008
	growth_space     = grow
	display = (explode=1,1)
}

%% motor -----------------------------------------------------------------
set hand motor_hand
{
    binding = 0.01, 0.01
    unbinding = 0.01, inf
    
    activity = move
    max_speed = 0
    stall_force = 6
    
    display = ( size=2; color=blue; )
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
    binding = 1, 0.01
    unbinding = 2.38, inf
    display = (size=9; color=green )
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

change simul mysim
{
    time_step = 0.001
    viscosity = 1
}

% Equilibration step
run 100000 simul
{
    solve = 0
}

change hand motor_hand
{
	unloaded_speed = [[motor_speed]]
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


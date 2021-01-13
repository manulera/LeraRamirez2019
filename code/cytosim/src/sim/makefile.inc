# Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.


OBJ_HANDS := hand.o hand_prop.o\
             motor.o motor_prop.o\
             slider.o slider_prop.o\
             actor.o actor_prop.o\
             nucleator.o nucleator_prop.o\
             tracker.o tracker_prop.o\
             cutter.o cutter_prop.o\
             chewer.o chewer_prop.o\
             rescuer.o rescuer_prop.o\
             regulator.o regulator_prop.o\
             mighty.o mighty_prop.o

OBJ_DIGITS := digit.o digit_prop.o\
              walker.o walker_prop.o\
              ase_walker.o ase_walker_prop.o\
              myosin.o myosin_prop.o\
              dynein.o dynein_prop.o\
              kinesin.o kinesin_prop.o

OBJ_SINGLES := single.o single_prop.o single_set.o\
               picket.o picket_long.o wrist.o wrist_long.o

OBJ_COUPLES := couple.o couple_prop.o couple_set.o couple_long.o\
               crosslink.o crosslink_long.o crosslink_prop.o\
               bridge.o bridge_prop.o\
               duo.o duo_long.o duo_prop.o\
               shackle.o shackle_long.o shackle_prop.o\
               flipper.o flipper_prop.o\
               trapper.o trapper_long.o trapper_prop.o\
               swapper.o swapper_prop.o\
               fork.o fork_prop.o
         
OBJ_FIBERS := filament.o fiber.o rigid_fiber.o fiber_prop.o fiber_set.o\
              fiber_locus.o fiber_binder.o lattice.o\
              dynamic_fiber.o dynamic_fiber_prop.o\
              classic_fiber.o classic_fiber_prop.o\
              treadmilling_fiber.o treadmilling_fiber_prop.o\
              growing_fiber.o growing_fiber_prop.o\
              tubule.o tubule_prop.o

OBJ_ORGANIZERS := organizer.o organizer_set.o \
                  aster.o aster_prop.o nucleus.o nucleus_prop.o\
                  fake.o fake_prop.o bundle.o bundle_prop.o

OBJ_SIM := movable.o mecable.o point_set.o object.o object_set.o\
           sphere_prop.o sphere.o sphere_set.o\
           bead_prop.o bead.o bead_set.o meca.o\
           solid_prop.o solid.o solid_set.o\
           field.o field_prop.o field_set.o\
           event.o event_prop.o event_set.o\
           point_exact.o point_interpolated.o\
           fiber_grid.o point_grid.o space_set.o\
           simul_prop.o simul.o interface.o parser.o


OBJ_CYTOSIM:=$(OBJ_SIM) $(OBJ_HANDS) $(OBJ_DIGITS) $(OBJ_FIBERS) $(OBJ_SINGLES) $(OBJ_COUPLES) $(OBJ_ORGANIZERS)

OBJ_CYTOSIM1D:=$(OBJ_CYTOSIM:.o=1D.o)
OBJ_CYTOSIM2D:=$(OBJ_CYTOSIM:.o=2D.o)
OBJ_CYTOSIM3D:=$(OBJ_CYTOSIM:.o=3D.o)

#----------------------------targets--------------------------------------------

INC_SIM=$(addprefix -Isrc/, SFMT math base sim sim/fibers sim/hands sim/singles sim/couples sim/organizers sim/spaces)

# need to fix the paths of libraries:
SIM_OBJS=$< $(addprefix build/, $(notdir $(filter-out %.cc, $^)))

#--- normal build:


sim: sim.cc libcytosim.a libcytospace.a libcytomath.a libcytobase.a | bin
	$(CXXFLAGS) $(INC_SIM) $(CODE_VERSION) $(SIM_OBJS) $(LINK) -o bin/$@
	$(DONE)
vpath sim  bin


libcytosim.a: $(OBJ_CYTOSIM) splash.o
	$(MAKELIB) build/$@ $(addprefix build/, $(notdir $^))
	$(DONE)

#--- build for debuging:

simd: sim.cc libcytosim.a libcytospace.a libcytomath.a libcytobase.a | bin
	$(CXX) $(FlagsD) $(INC_SIM) $(CODE_VERSION) $(SIM_OBJS) $(LINK) -o bin/$@
	$(DONE)
vpath simd  bin


#--- builds with specific dimension:


bin1/sim: sim.cc libcytosim1D.a libcytospace1D.a libcytomath1D.a libcytobase.a
	$(CXXFLAGS) $(INC_SIM) $(CODE_VERSION) -DDIM=1 $(SIM_OBJS) $(LINK) -o bin1/sim
	$(DONE)

libcytosim1D.a: $(OBJ_CYTOSIM1D) splash.o
	$(MAKELIB) build/$@ $(addprefix build/, $(notdir $^))
	$(DONE)


bin2/sim: sim.cc libcytosim2D.a libcytospace2D.a libcytomath2D.a libcytobase.a
	$(CXXFLAGS) $(INC_SIM) $(CODE_VERSION) -DDIM=2 $(SIM_OBJS) $(LINK) -o bin2/sim
	$(DONE)

libcytosim2D.a: $(OBJ_CYTOSIM2D) splash.o
	$(MAKELIB) build/$@ $(addprefix build/, $(notdir $^))
	$(DONE)


bin3/sim: sim.cc libcytosim3D.a libcytospace3D.a libcytomath3D.a libcytobase.a
	$(CXXFLAGS) $(INC_SIM) $(CODE_VERSION) -DDIM=3 $(SIM_OBJS) $(LINK) -o bin3/sim
	$(DONE)

libcytosim3D.a: $(OBJ_CYTOSIM3D) splash.o
	$(MAKELIB) build/$@ $(addprefix build/, $(notdir $^))
	$(DONE)


#----------------------------rules----------------------------------------------

$(OBJ_CYTOSIM): %.o: %.cc %.h | build
	$(CXXFLAGS) -c $(INC_SIM) $< -o build/$@

$(OBJ_CYTOSIM1D): %1D.o: %.cc %.h | build
	$(CXXFLAGS) -c -DDIM=1 $(INC_SIM) $< -o build/$@

$(OBJ_CYTOSIM2D): %2D.o: %.cc %.h | build
	$(CXXFLAGS) -c -DDIM=2 $(INC_SIM) $< -o build/$@

$(OBJ_CYTOSIM3D): %3D.o: %.cc %.h | build
	$(CXXFLAGS) -c -DDIM=3 $(INC_SIM) $< -o build/$@

splash.o: %.o: %.cc %.h | build
	$(CXXFLAGS) -c $(CODE_VERSION) $(COMPILER_VERSION) -Isrc/base -Isrc/math $< -o build/$@

mecaGL.o: %.o: %.cc %.h | build
	$(CXXFLAGS) -c $(INC_SIM) -Isrc/disp $< -o build/$@


// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
/**
 @file
 @brief Global Compile Switches
*/

#ifndef SIM_H
#define SIM_H

/**
 Enable code to be able to read old trajectory files
 Option normally ON
 */
#define BACKWARD_COMPATIBILITY


/**
 Enables special feature to simulate fluid flow 
 Option normally OFF
 */
//#define NEW_CYTOPLASMIC_FLOW


/**
 Use alternative binding scheme in which Hands have their own counters
 */
//#define TRICKY_HAND_ATTACHMENT


/**
 If the keyword below is defined, the viscous drag of the fibers
 will be different in the transverse and parallel directions, such that
 it will be 2x easier to move a fiber along it longitudinal direction.
 
 This is unpublished development, and the keyword should NOT be enabled
 */
//#define NEW_ANISOTROPIC_FIBER_DRAG


/**
 Enables Myosin, Kinesin and Dynein
 Option normally OFF
 */
//#define NEW_HANDS


/** 
 Enables advanced Space
 */
#define NEW_SPACES


/// a macro to print some text only once
#define PRINT_ONCE(a) { static bool virgin=true; if (virgin) { virgin=false; MSG(a); } }


/// Lattice with values instead of 0 1
#define MANU_LATTICE
/// A temporary solution to have multiple lattices in the same MT, representing the different sides
//#define MULTI_LATTICE


// TRAP_SINGLES
#define TRAP_SINGLES 2

#define GILLESPIE_DIGITS


#endif


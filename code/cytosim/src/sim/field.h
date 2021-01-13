// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include "field_base.h"
#include "field_values.h"


/// type of Field used in Cytosim
typedef FieldBase<FieldScalar> Field;


/// initialize diffusion matrix (only for FieldScalar)
template < > 
void FieldBase<FieldScalar>::prepare();

template < >
void FieldBase<FieldScalar>::prepareDiffusion(real);

template < >
void FieldBase<FieldScalar>::prepareDiffusion(real, unsigned char *);

/// diffusion step
template < > 
void FieldBase<FieldScalar>::step(FiberSet&);

/// calculate second derivative of scalar field
template < >
void FieldBase<FieldScalar>::laplacian(const real*, real*) const;

/// perform diffusion in X-direction
template < >
void FieldBase<FieldScalar>::diffuseX(real*, real);

/// set values on X edges of field
template < >
void FieldBase<FieldScalar>::setEdgesX(real*, real);

/// set values on Y edges of field
template < >
void FieldBase<FieldScalar>::setEdgesY(real*, real);

/// set values on Z edges of field
template < >
void FieldBase<FieldScalar>::setEdgesZ(real*, real);

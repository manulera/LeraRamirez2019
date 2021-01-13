// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#ifndef FIELD_H
#define FIELD_H

#include "dim.h"
#include "real.h"
#include "grid.h"
#include "space.h"
#include "object.h"
#include "iowrapper.h"
#include "messages.h"
#include "exceptions.h"
#include "matsparsesym1.h"
#include "matsparsesym2.h"
#include "field_prop.h"
class FiberSet;

#ifdef DISPLAY
   #include "gle.h"
   #include "grid_display.h"
#endif


/// value of type VAL defined as a function of position over the simulation Space
/**
 A field represents a value that is varying with position over space.
 It does so by storing a different value for each cell in a regular grid.
 
 Each cell holds the amount of molecules.
 The local concentration can be obtained by dividing by the cell volume:
 @code
 concentration = FieldGrid::cell(position) / FieldGrid::cellVolume();
 @endcode

 Note that the field is build on a Grid with square cells, because diffusion/reaction are then
 easier to implement. The Grid does not necessarily match the edges of the Space exactly,
 but instead extends outside, such as to cover the 'inside' region entirely.
 
 */
template < typename VAL >
class FieldBase : public Grid<VAL, DIM, unsigned int>, public Object
{
public:
    
    /// the type of Grid from which the Field is derived
    typedef Grid<VAL, DIM, unsigned int> FieldGrid;
    
    /// forward the index type:
    typedef typename FieldGrid::index_type index_type;
    
    /// type of value
    typedef VAL value_type;
    
    /// property
    FieldProp const* prop;
    
    ///
    static real display_amp;
    
    ///
    static Space const* display_spc;
    
private:
    
    /// disabled default constructor
    FieldBase();

    /// duplicate field
    real*    fiTMP;
    
    /// allocated size of fiTMP
    unsigned fiTMPSize;
    
    /// matrix for diffusion
    MatrixSparseSymmetric1 fiDiffusionMatrix;
    
    /// initialize to cover the given Space with squares of size 'step'
    void setGrid(Vector inf, Vector sup, real step, bool tight)
    {
        assert_true( step > REAL_EPSILON );
        // we add a safety border (in micro-meters)
        const real extra = tight ? 0 : 1;
        
        int size[3] = { 0, 0, 0 };
        // we use square cells:
        for ( int d = 0; d < DIM; ++d )
        {
            size[d] = (int)ceil( (sup[d]-inf[d]+extra) / step );
            real mid = 0.5 * ( inf[d] + sup[d] );
            inf[d] = mid - 0.5 * step * size[d];
            sup[d] = mid + 0.5 * step * size[d];
        }
        
        FieldGrid::setDimensions(inf, sup, size);
        
        //verify the cell size:
        for ( int d = 0; d < DIM; ++d )
        {
            real dif = FieldGrid::delta(d) * step - 1.0;
            if ( fabs(dif) > 1e-6 )
            {
                MSG.warning("Field:step[%i] is not as expected:\n", d);
                MSG.warning("  field: %f  prop: %f\n", FieldGrid::cellWidth(d), step);
            }
        }            
    }
    
    /// allocate memory for the scalar field (setGrid() must be called before)
    void createGrid()
    {
        assert_true( FieldGrid::hasDimensions() );
        // delete preexisting grid if necessary:
        FieldGrid::destroy();
        // create the grid using the calculated dimensions:
        FieldGrid::createCells();
        // set all values to zero (already done in the constructor of FieldValue)
        // FieldGrid::clear();
        // report dimensions:
        //FieldGrid::printSummary(std::clog, "Field");
    }

public:
    #pragma mark -
    
    /// constructor
    FieldBase(FieldProp const* p)
    {
        prop=p;
        fiTMP=0;
        fiTMPSize=0;
    }
    
    /// destructor
    ~FieldBase()
    {
        if ( fiTMP )
            delete[] fiTMP;
    }
    
    /// initialize with squares of size 'step'
    void setField()
    {
        assert_true( prop );

        if ( ! FieldGrid::hasCells() )
        {            
            if ( prop->confine_space_ptr == 0 )
                throw InvalidParameter("A space must be defined to set a field");
            
            Vector inf, sup;
            prop->confine_space_ptr->boundaries(inf, sup);
            
            if ( prop->periodic )
            {
                FieldGrid::periodic(true);
                setGrid(inf, sup, prop->step, true);
            }
            else
            {
                setGrid(inf, sup, prop->step, false);
            }
            createGrid();
            //std::clog << "setField() step="<< prop->step<< " nCells="<< FieldGrid::nbCells()<<std::endl;
            MSG(4, "Field %lx set with %i cells of size %.3f um\n", this, FieldGrid::nbCells(), prop->step);
        }
    }
    
        
    /// true if field was set
    unsigned hasField()  const { return FieldGrid::hasCells(); }
    
    /// size of cell
    real cellWidth() const { return FieldGrid::cellWidth(0); }
    
    //------------------------------ simulation --------------------------------
    #pragma mark -

    /// set all cells to value = volume * conc
    void setConcentration(value_type conc)
    {
        FieldGrid::setValues( conc * FieldGrid::cellVolume() );
    }
 
    
    /// set cells that are inside `spc` to value = volume * conc
    void setConcentration(Space const* spc, value_type in, value_type ou)
    {
        real i = in * FieldGrid::cellVolume();
        real o = ou * FieldGrid::cellVolume();
        
        for ( index_type c = 0; c < FieldGrid::nbCells(); ++c )
        {
            Vector w;
            FieldGrid::setPositionFromIndex(w, c, 0.5);
            if ( spc->inside(w) )
                FieldGrid::icell(c) = i;
            else
                FieldGrid::icell(c) = o;
        }
    }
    
    /// simulation step 
    void step(FiberSet&) {}

    /// calculate second derivative of field
    void laplacian(const real*, real*) const {}
    
    /// calculate second derivative of field
    void diffuseX(real*, real) {}
    
    /// set values of field on its edges
    void setEdgesX(real*, real) {}
    
    /// set values of field on its edges
    void setEdgesY(real*, real) {}

    /// set values of field on its edges
    void setEdgesZ(real*, real) {}

    /// initialize diffusion matrix (only for FieldScalar)
    void prepareDiffusion(real) {}

    /// initialize diffusion matrix (only for FieldScalar)
    void prepareDiffusion(real, unsigned char *) {}

    /// initialize Field
    void prepare() {}

    //------------------------------- object -----------------------------------
    #pragma mark -

    /// a unique character identifying the class
    static const Tag TAG = 'i';
    
    /// return unique character identifying the class
    Tag    tag() const { return TAG; }
    
    /// return index of 'prop' in corresponding PropertyList
    Property const* property() const { return prop; }
    
    /// a static_cast<> of Node::next()
    FieldBase<VAL>* next()  const  { return static_cast<FieldBase<VAL>*>(nNext); }
    
    /// a static_cast<> of Node::prev()
    FieldBase<VAL>* prev()  const  { return static_cast<FieldBase<VAL>*>(nPrev); }

    //------------------------------ read/write --------------------------------
    #pragma mark -
    
    /// print total, minimum and maximum value
    void   writeInfo(std::ostream& out) const
    {
        real vol = FieldGrid::cellVolume();
        value_type sum, mn, mx;
        FieldGrid::infoValues(sum, mn, mx);
        out << prop->name() << " sum " << sum << " min " << mn/vol << " max " << mx/vol << std::endl;
    }

    /// write Field to file using VAL::write()
    /** Some of this should be moved to Grid */
    void   write(Outputter& out) const
    {
        if ( FieldGrid::hasCells() && prop->save )
        {        
            out.writeUInt16(DIM);
            for ( int d = 0; d < DIM; ++d )
            {
                out.writeSoftSpace();
                out.writeUInt32(FieldGrid::dim(d));
                out.writeFloat(FieldGrid::inf(d));
                out.writeFloat(FieldGrid::sup(d));
            }
            out.writeSoftSpace();
            out.writeUInt32(FieldGrid::nbCells());
            for ( index_type c = 0; c < FieldGrid::nbCells(); ++c )
                FieldGrid::icell(c).write(out);
            out.writeSoftNewline();
        }
        
        if ( prop->positive )
        {
            if ( FieldGrid::hasNegativeValue() )
                throw Exception("Aborting because Field has negative values");
        }
    }
    
    
    /// read Field from file using VAL::read()
    void   read_data(Inputter& in, Simul&)
    {
        int  size[DIM] = { 0 };
        real minB[DIM] = { 0 }, maxB[DIM] = { 0 };
        
        try {
            unsigned int dim = in.readUInt16();
            if ( dim != DIM )
                throw InvalidIO("field::dimensionality mismatch");
            
            for ( unsigned int d = 0; d < dim; ++d )
            {
                size[d] = in.readUInt32();
                minB[d] = in.readFloat();
                maxB[d] = in.readFloat();
            }
            
            FieldGrid::setDimensions(minB, maxB, size);
            createGrid();
            
            index_type nbc = in.readUInt32();
            if ( nbc != FieldGrid::nbCells() )
            {
                printf("file: %u field:%u\n", nbc, FieldGrid::nbCells());
                throw InvalidIO("mismatch in Field::size");
            }
            //std::clog << "Field::read() nb_cells=" << nbc << std::endl;

            for ( index_type c = 0; c < nbc; ++c )
                FieldGrid::icell(c).read(in);
        }
        catch( Exception & e ) {
            e << ", in Field::read()";
            throw;
        }
    }
    
    /// read Field and checks that the Grid::step has not changed
    void   read(Inputter& in, Simul& sim, Tag)
    {
        read_data(in, sim);
        
        if ( prop )
        {
            for ( unsigned int d = 0; d < DIM; ++d )
            {
                real dif = FieldGrid::delta(d)*prop->step - 1.0;
                if ( fabs(dif) > 1e-6 )
                {
                    MSG.warning("Field:step[%i] has changed:\n", d);
                    MSG.warning("  file: %f  prop: %f\n", FieldGrid::cellWidth(d), prop->step);
                }
            }
            
            /*
             we should extrapolate the data that were read to a grid with the
             resolution specified by prop->step
             */
            
        }
    }
    
    //------------------------------- display ----------------------------------

#ifdef DISPLAY
    
    static bool field_set_color(VAL const& val, Vector const& pos)
    {
        if ( display_spc && ! display_spc->inside(pos) )
            return false;
        val.setColor(display_amp);
        return true;
    }
    
    /// openGL display function
    void display() const
    {
        display_amp = 1.0 / ( prop->display_scale * FieldGrid::cellVolume() );
        
        glPushAttrib(GL_ENABLE_BIT);
        glDisable(GL_LIGHTING);
        glDisable(GL_DEPTH_TEST);
        glDisable(GL_CULL_FACE);
        drawValues(*this, field_set_color);
        if ( 0 )
        {
            glColor4f(1, 0, 1, 1);
            glLineWidth(0.5);
            drawEdges(*this);
        }
        glPopAttrib();
    }

    
    /// openGL display function
    /**
     display all cells that are inside field:confine_space
     */
    void display(bool all, Vector3 const& dir, const real pos) const
    {
        display_amp = 1.0 / ( prop->display_scale * FieldGrid::cellVolume() );
        if ( all )
            display_spc = 0;
        else
            display_spc = prop->confine_space_ptr;
        
        //glPushAttrib(GL_ENABLE_BIT|GL_POLYGON_BIT);
        glPushAttrib(GL_ENABLE_BIT);
        glDisable(GL_LIGHTING);
        glDisable(GL_DEPTH_TEST);
        glDisable(GL_CULL_FACE);
        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
        //glLineWidth(1);
        //glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
#if ( DIM == 3 )
        drawValues(*this, field_set_color, dir, pos);
#else
        drawValues(*this, field_set_color);
#endif
        glPopAttrib();
    }
    
#endif
};


#endif

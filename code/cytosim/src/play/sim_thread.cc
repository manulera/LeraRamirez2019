// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#include <sys/time.h>
#include <time.h>
#include "sim_thread.h"
#include "exceptions.h"
#include "picket.h"


//------------------------------------------------------------------------------

/**
 This uses a Parser that cannot write to disc.
 The function callback is called when Parser::hold() is reached.
 */
SimThread::SimThread(void (*callback)(void))
: Parser(simul, 1, 1, 1, 1, 0), hold_callback(callback)
{
    mFlag   = 0;
    mState  = 0;
    mHold   = 0;
    mPeriod = 1;
    pthread_mutex_init(&mMutex, 0);
    pthread_cond_init(&mCondition, 0);
}

/**
 Possible issue:
 When quitting the application, this destructor might be called after the
 destructor of Simul(), in which case it will access non-existent data,
 most likely causing a crash().
 */
SimThread::~SimThread()
{
    stop();
    pthread_cond_destroy(&mCondition);
    pthread_mutex_destroy(&mMutex);
    pthread_detach(mThread);
}

//------------------------------------------------------------------------------
#pragma mark -

void SimThread::hold()
{
    //assert_true( pthread_equal(pthread_self(), mThread) );

    if ( ++mHold >= mPeriod )
    {
        mHold = 0;
        hold_callback();
        if ( mFlag )
        {
            // gently terminate:
            //std::clog << pthread_self() << " exit\n";
            mFlag = 0;
            mState = 0;
            unlock();
            pthread_exit(0);
        }
        wait();
        //std::clog << pthread_self() << " continue\n";
    }
}


void SimThread::self_start()
{
    //assert_true( pthread_equal(pthread_self(), mThread) );
    mState = 1;
    
    try {
        if ( Parser::readConfig() )
            std::cerr << "You must specify a config file\n";
    }
    catch( Exception & e ) {
        std::cerr << std::endl << "Error: " << e.what() << std::endl;
        //flashText("Error: simulation died");
    }
    // the thread terminates normally:
    hold_callback();
    mState = 0;
    unlock();
}


void SimThread::self_extend()
{
    //assert_true( pthread_equal(pthread_self(), mThread) );
    mState = 1;
    
    try {
        simul.prepare();
        // enter an infinite loop. Thread will need to be killed
        while ( 1 )
        {
            simul.step();
            simul.solve();
            hold();
        }
    }
    catch( Exception & e ) {
        std::cerr << std::endl << "Error: " << e.what() << std::endl;
        simul.relax();
        mState = 0;
        unlock();
        //flashText("Error: %s", e.what());
    }
}



void* start_helper(void * arg)
{
    SimThread * lt = static_cast<SimThread*>(arg);
    lt->self_start();
    return 0;
}


void* extend_helper(void * arg)
{
    SimThread * lt = static_cast<SimThread*>(arg);
    lt->self_extend();
    return 0;
}

//------------------------------------------------------------------------------
#pragma mark -

/**
 This attempts to start the live simulation
 */
int SimThread::start()
{
    if ( !mState )
    {
        if ( 0 == trylock() )
        {
            if ( pthread_create(&mThread, 0, start_helper, this) )
                return 2;
            return 0;
        }
    }
    return 1;
}


int SimThread::extend()
{
    if ( !mState )
    {
        if ( 0 == trylock() )
        {
            if ( pthread_create(&mThread, 0, extend_helper, this) )
                return 2;
            return 0;
        }
    }
    return 1;
}


/**
 ask the live-thread to exit at the next spontaneous halt
*/ 
void SimThread::stop()
{
    if ( mState )
    {
        // request clean termination:
        mFlag = 1;
        signal();
        // wait for termination:
        pthread_join(mThread, 0);
        mState = 0;
    }
}

/**
 kill the live-thread immediately
 */
void SimThread::cancel()
{
    if ( mState )
    {
        // request hard termination:
        if ( 0 == pthread_cancel(mThread) )
        {
            // wait for termination:
            pthread_join(mThread, 0);
            mState = 0;
            unlock();
        }
    }
}

//------------------------------------------------------------------------------
#pragma mark -


SingleProp * SimThread::getHandleProperty() const
{
    return simul.findProperty<SingleProp*>("single", "mouse_single");
}


SingleProp * SimThread::makeHandleProperty(real range)
{
    // Create a Hand that attaches fast and never detach:
    HandProp * hap = new HandProp("mouse_hand");
    hap->binding_range   = range;
    hap->binding_rate    = 1000;
    hap->unbinding_rate  = 0;
    hap->unbinding_force = INFINITY;
    hap->complete(&simul);
    simul.properties.deposit(hap);

    SingleProp * sip = new SingleProp("mouse_single");
    sip->hand = "mouse_hand";
    sip->stiffness = 256;
    sip->complete(&simul);
    simul.properties.deposit(sip);
    
    return sip;
}


Single * SimThread::createHandle(const Vector & pos, real range)
{
    SingleProp * sip = getHandleProperty();
    if ( sip == 0 )
        sip = makeHandleProperty(range);
    Single * res = new Picket(sip, pos);
    simul.singles.add(res);
    mHandle = res;
    return res;
}


ObjectList SimThread::allHandles(SingleProp const* sip) const
{
    return simul.singles.collect(match_property, sip);
}


bool SimThread::selectClosestHandle(Vector const& pos, real range)
{
    SingleProp * sip = getHandleProperty();
    
    if ( sip )
    {
        ObjectList objs = allHandles(sip);
    
        real dsm = 0;
        Single * res = 0;
        for ( ObjectList::iterator oi = objs.begin(); oi < objs.end(); ++oi )
        {
            Single * s = static_cast<Single*>(*oi);
            real d = ( s->posFoot() - pos ).normSqr();
            if ( res == 0  ||  d < dsm )
            {
                res = s;
                dsm = d;
            }
        }
        if ( res && dsm < range )
        {
            mHandle = res;
            return 1;
        }
    }
    return 0;
}


Single * SimThread::handle()
{
    SingleProp * sip = getHandleProperty();
    if ( sip && mHandle )
    {
        ObjectList objs = allHandles(sip);
        for ( ObjectList::iterator oi = objs.begin(); oi < objs.end(); ++oi )
            if ( *oi == mHandle )
                return mHandle;
    }
    mHandle = 0;
    return 0;
}


void SimThread::detachHandle()
{
    if ( mHandle )
    {
        if ( mHandle->attached() )
            mHandle->detach();
    }
}

void SimThread::moveHandle(const Vector & pos)
{
    if ( mHandle )
    {
        mHandle->setPosition(pos);
    }
}


void SimThread::moveHandles(const Vector & vec)
{
    SingleProp * sip = getHandleProperty();
    if ( sip )
        ObjectSet::translateObjects(allHandles(sip), vec);
}


void SimThread::deleteHandles()
{
    lock();
    SingleProp * sip = getHandleProperty();
    if ( sip )
        simul.erase(allHandles(sip));
    mHandle = 0;
    unlock();
}

void SimThread::clear()
{
    simul.erase();
    mHandle = 0;
}

//------------------------------------------------------------------------------
#pragma mark -

/**
 Read config file from the start, allowing parameters to be changed, 
 while simulation objects remain as they are.
 
 If the simulation is running live, this will pause it,
 read the config file, and allow it to proceed.
 */
void SimThread::reloadParameters()
{
    lock();
    // the parser can only change properties:
    try {
        if ( Parser(simul, 0, 1, 0, 0, 0).readConfig() )
            std::cerr << "Error : File not found";
    }
    catch( Exception & e ) {
        std::cerr << "Error : " << e.what();
    }
    // do not allow display parameters to be changed:
    simul.prop->display_fresh = false;
    unlock();
}


/**
 This will execute the code specified in `iss`.
 
 If the simulation is running live, the SimThread is paused,
 and the code is executed with another Parser.
 When that is completed, the SimThread is released.
 
 The parser has full rights during the execution
 */
int  SimThread::execute(std::istream& iss)
{
    lock();
    try {
        Parser(simul, 1, 1, 1, 1, 1).parse(iss, ", while executing magic code");
    }
    catch( Exception & e ) {
        std::cerr << "Error : " << e.what();
    }
    unlock();
    return 0;
}

/**
 Save current state in two files
 */
void SimThread::exportObjects(bool binary)
{
    lock();
    try {
 
        char str[64] = { '\0' };
        
        snprintf(str, sizeof(str), "properties%04i.cmo", reader.currFrame());
        simul.writeProperties(str, true);
        
        snprintf(str, sizeof(str), "objects%04i.cmo", reader.currFrame());
        Outputter out(str, false, binary);
        simul.writeObjects(out);
    }
    catch( Exception & e ) {
        std::cerr << "Error in Simul::exportObjects(): " << e.what();
    }
    unlock();
}


void SimThread::writeProperties(std::ostream& os, bool prune)
{
    lock();
    try {
        simul.writeProperties(os, prune);
    }
    catch( Exception & e ) {
        std::cerr << "Error in Simul::writeProperties(): " << e.what();
    }
    unlock();
}



// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#ifndef LIVE_THREAD_H
#define LIVE_THREAD_H

#include <pthread.h>
#include "parser.h"
#include "simul.h"
#include "frame_reader.h"


/// SimThread is used to run a simulation in a dedicated thread
class SimThread : private Parser
{
    /// disabled default constructor
    SimThread();
    
private:
    
    /// Simulation object
    Simul           simul;

    /// reader used to access frames in trajectory file
    FrameReader     reader;

    
    /// callback invoked when the thread is halted
    void           (*hold_callback)(void);
    
    /// slave thread
    pthread_t       mThread;
    
    /// a flag reflecting if mThread is running or not
    bool            mState;
    
    /// a flag to indicate that mThread should terminate
    bool            mFlag;
    
    /// mutex to access simulation state
    pthread_mutex_t mMutex;

    /// condition variable used to control the thread execution
    pthread_cond_t  mCondition;
    
    /// counter for hold()
    unsigned int    mHold;
    
    /// period for hold()
    unsigned int    mPeriod;

    
    /// the current Single being controlled with the mouse
    Single *        mHandle;
    
    /// return the SingleProp used for the handles
    SingleProp *  getHandleProperty() const;

    /// make a new SingleProp for the handles with given attachment range
    SingleProp *  makeHandleProperty(real range);
    
    /// return list of Handles
    ObjectList    allHandles(SingleProp const*) const;

public:
    
    /// run the simulation live
    void          self_start();
    
    /// continue to run a simulation beyond its normal termination
    void          self_extend();

    /// redefines Interface::hold(), will be called between commands
    void          hold();
    
    /// create a SimThread with given holding function callback
    SimThread(void (*callback)(void));
    
    /// destructor
    ~SimThread();

#if ( 1 )

    /// lock access to the Simulation data
    void       lock()    { pthread_mutex_lock(&mMutex); }
    
    /// unlock access to the Simulation data
    void       unlock()  { pthread_mutex_unlock(&mMutex);}
    
    /// try to lock access to the Simulation data
    int        trylock() { return pthread_mutex_trylock(&mMutex); }

    /// wait for the condition
    int        wait()    { return pthread_cond_wait(&mCondition, &mMutex); }
    
    /// send signal to other threads
    void       signal()  { if ( mState ) pthread_cond_signal(&mCondition); }

#else
    
    /// lock access to the Simulation data
    void       lock()    { pthread_mutex_lock(&mMutex); std::clog << pthread_self() << "   lock\n"; }
    
    /// unlock access to the Simulation data
    void       unlock()  { pthread_mutex_unlock(&mMutex); std::clog << pthread_self() << " unlock\n"; }
    
    /// try to lock access to the Simulation data
    int        trylock() { int r=pthread_mutex_trylock(&mMutex); std::clog << pthread_self() << " trylock " << r << "\n"; return r; }
    
    /// wait for the condition to be
    int        wait()    { std::clog << pthread_self() << " wait\n"; return pthread_cond_wait(&mCondition, &mMutex); }
    
    /// signal other thread to continue
    void       signal()  { std::clog << pthread_self() << " signal\n"; pthread_cond_signal(&mCondition); }
    
#endif
    
    /// Simul reference
    Simul&     sim() { return simul; }
    
    /// set how many 'hold()' are necessary to halt the thread
    void       period(unsigned int c) { mPeriod = c; }

    
    /// true if child thread is running
    bool       alive() { return mState; }
    
    /// start the thread that will run a simulation
    int        start();
    
    /// continue to run the simulation after its normal termination
    int        extend();
    
    /// gently stop the simulation
    void       stop();

    /// stop the simulation
    void       cancel();

    /// clear the simulation world
    void       clear();
    
    /// halt the live simulation, read the config file and change the object parameters
    void       reloadParameters();
    
    /// execute given code
    int        execute(std::istream&);
    
    /// export simulation Propertes and Objects to file
    void       exportObjects(bool binary);
    
    /// export properties to file
    void       writeProperties(std::ostream&, bool prune);

    
    /// open file for input
    void       openFile()       { reader.openFile(simul.prop->trajectory_file); }
    
    /// true if ready to read from file
    bool       goodFile()       { return reader.good(); }
    
    /// index of current frame
    int        eof()            { return reader.eof(); }
    
    /// index of current frame
    int        currFrame()      { return reader.currFrame(); }
    
    /// attempt to load specified frame
    int        loadFrame(int f) { lock(); int r=reader.loadFrame(simul, f); unlock(); return r; }

    /// load next frame in file
    int        loadNextFrame()  { lock(); int r=reader.loadNextFrame(simul); unlock(); return r; }

    
    /// return current handle
    Single *   handle();

    /// make a new Single that can be controlled by the user
    Single *   createHandle(const Vector&, real range);
    
    /// switch current handle
    bool       selectClosestHandle(const Vector&, real range);
    
    /// detach current handle
    void       detachHandle();
    
    /// move the current handle
    void       moveHandle(const Vector&);
    
    /// move all handles
    void       moveHandles(const Vector&);
    
    /// delete all handles
    void       deleteHandles();
    
    /// detach current handle from mouse control
    void       releaseHandle() { mHandle = 0; }
    
};


#endif


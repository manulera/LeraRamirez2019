// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#include "hand_prop.h"

/**
 Callback for mouse clicks
 */
void Player::processMouseClick(int, int, const Vector3& pos3, int)
{
    // distance in pixels where mouse-Hand binds:
    const int pixrad = 6;
    
    const real range = pixrad * glApp::currentView().pixelSize();
    Vector pos(pos3.XX, pos3.YY, pos3.ZZ);

    simThread.lock();
    if ( simThread.selectClosestHandle(pos, range) )
        simThread.moveHandle(pos);
    else
    {
        if ( simThread.handle() )
        {
            simThread.detachHandle();
            simThread.moveHandle(pos);
        }
        else
        {
            Single * s = simThread.createHandle(pos, range);
            PointDisp *& pd = s->prop->hand_prop->disp;
            if ( pd == 0 )
            {
                pd         = new PointDisp("hand:display", "mouse");
                pd->size   = 2 * pixrad;
                pd->color  = glApp::currentView().front_color;
                pd->color2 = pd->color;
                dproperties.deposit(pd);
            }
        }
    }
    simThread.unlock();
}


//------------------------------------------------------------------------------
/**
 Processes mouse motion
 */
void Player::processMouseDrag(int, int, Vector3& ori3, const Vector3& pos3, int mode)
{
    Vector pos(pos3.XX, pos3.YY, pos3.ZZ);
    Vector ori(ori3.XX, ori3.YY, ori3.ZZ);
    
    simThread.lock();
    if ( mode )
    {
        simThread.moveHandles(pos-ori);
        ori3 = pos3;
    }
    else
        simThread.moveHandle(pos);
    simThread.unlock();
}



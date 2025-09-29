#File: plotfield.py
#Purpose: helper function for plotting vector and scalar field
#Author: Abhikalp Shekhar

import numpy as np 
from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
from matplotlib import cm
plt.style.use('_mpl-gallery')
from vector import *
from electrostatics import *
from Utils import *
import math

def getBoudingBox(charges):
    verify(isinstance(charges,(list,tuple)) and len(charges)>0, "charges should be a tuple")
    [verify(isinstance(x,(Charge,Point)),"bad input") for x in charges]
    if isinstance(charges[0], Charge):
        points = [x.p for x in charges]
    else:
        points = charges
    xmin = min([p.x for p in points])
    ymin = min([p.y for p in points])
    xmax = max([p.x for p in points])
    ymax = max([p.y for p in points])
    if xmax-xmin < 1:
        d = 1.0
        xmax=xmax+d/2.0
        xmin=xmin-d/2.0
    if ymax-ymin < 1:
        d=1.0
        ymax=ymax+d/2.0
        ymin=ymin-d/2.0
    return BoundingBox(xmin,xmax,ymin,ymax)

def showBox(node):
    def drawbox(box,x1=[],y1=[],x2=[],y2=[]):
        x1=[box.xmin,box.xmax]
        y1=[box.ymin,box.ymin]
        plt.plot(x1,y1,'b-')
        x1=[box.xmin,box.xmin]
        y1=[box.ymin,box.ymax]
        plt.plot(x1,y1,'b-')
        x1=[box.xmin,box.xmax]
        y1=[box.ymax,box.ymax]
        plt.plot(x1,y1,'b-')
        x1=[box.xmax,box.xmax]
        y1=[box.ymin,box.ymax]
        plt.plot(x1,y1,'b-')

    drawbox(node.box)
    if node.ne: showBox(node.ne)
    if node.nw: showBox(node.nw)
    if node.se: showBox(node.se)
    if node.sw: showBox(node.sw)

def showTree(node,wts,pts):
    showBox(node)
    sz=len(pts)
    for i in range(sz):
        if wts[i]<0:
            plt.scatter(pts[i].x,pts[i].y,c='red')
        else:
            plt.scatter(pts[i].x,pts[i].y,c='green')

def showCharges(charges,box,showvals,fig):
    xx=[c.p.x for c in charges]
    yy=[c.p.y for c in charges]
    mx=max([abs(x.q) for x in charges ])
    sz=[30+5*int(abs(x.q)/mx) for x in charges]
    vp=np.ma.masked_where(np.asarray([c.q for c in charges]) <0.0,sz)
    vn=np.ma.masked_where(np.asarray([c.q for c in charges]) >0.0,sz)
    axs = fig.add_subplot(221)
    axs.set_title('Charge Distribution')
    axs.scatter(xx,yy,s=vp,marker="+",c="red")
    axs.scatter(xx,yy,s=vn,marker="o",c="green")

def showPotentialsAndFields(charges,box,options,fig,plotfield=False):
    nxpoints = options.get('nxpoints',20)
    nypoints = options.get('nypoints',20)
    xx = np.linspace(box.xmin,box.xmax,nxpoints)
    yy = np.linspace(box.ymin,box.ymax,nypoints)
    XX,YY=np.meshgrid(xx,yy)
    potential = np.zeros((nxpoints,nypoints))
    def potentialfn(p):
        return Charge.getPotentialAtPoint(charges,p,options)
    
    if plotfield:
        ex=np.zeros((nxpoints,nypoints))
        ey=np.zeros((nxpoints,nypoints))
        scalarfld = ScalarField(potentialfn)

    for i in range(XX.shape[0]):
        for j in range(XX.shape[1]):
            p = Point(XX[i,j],YY[i,j])
            potential[i,j] = potentialfn(p)
            if plotfield:
                e = scalarfld.gradient(p).scale(-1.0)
                ex[i,j],ey[i,j]=(e.x,e.y)

    axs = fig.add_subplot(222, projection='3d')
    axs.plot_surface(XX, YY, potential, cmap='cool', alpha=0.8)
    #axs.set_zlim3d(-1, 1)
    axs.set_title('Potential', fontsize=14)
    axs.set_xlabel('x', fontsize=12)
    axs.set_ylabel('y', fontsize=12)
    axs.set_zlabel('z', fontsize=12)

    if plotfield:
        color = np.log(np.hypot(ex, ey))
        axs2 = fig.add_subplot(223)
        axs2.quiver(XX,YY,ex,ey,color='b', linewidth=0.5, cmap=plt.get_cmap('gist_earth'))
        axs2.set_aspect('equal')
        axs3 = fig.add_subplot(224)
        axs3.streamplot(XX,YY,ex,ey,color=color,linewidth=0.5, cmap=plt.cm.inferno, density = 2, arrowstyle='->', arrowsize=1)
        axs2.set_title('Electric Field', fontsize=14)
        axs3.set_title('Electric Field', fontsize=14)


def plot2DChargeAndFields(charges,box=None,options=None):
    if options is None:
        options = {}
    verify(isinstance(charges,(list,tuple)), "charges should be a tuple")
    [verify(isinstance(x,Charge),"bad input") for x in charges]
    if box is None:
        box = getBoudingBox(charges)
    fig = plt.figure()
    
    showCharges(charges,box,options.get("showvalues",True),fig)
    showPotentialsAndFields(charges,box,options,fig,True)
    plt.show()


if __name__ == "__main__":
    np.random.seed(97)
    charges = []
    for i in range(200):
        charges.append(Charge(np.random.uniform(-5.0,5.0),
        Point(np.random.uniform(-0.5,0.5),np.random.uniform(-0.5,0.5))))
    plot2DChargeAndFields(charges)
    #plot2DChargeAndFields(charges,options={"method":"fmm"})
    print("done")

    
    


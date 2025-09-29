
import numpy as np 
from vector import *
from electrostatics import *
from Utils import *
import math
import plotfield as pfl

class Node:
    def __init__(self,nw,ne,sw,se,wt,cg,l,box):
        self.nw=nw
        self.ne=ne
        self.sw=sw
        self.se=se
        self.l = l 
        self.wt=wt
        self.cg=cg
        self.box=box

    def isleafnode(self):
        return self.ne== None and \
            self.nw == None and \
            self.se == None and \
            self.sw == None
    
    def energy(self):
        u = 0.0
        if(self.ne and not self.ne.isleafnode()):
            u = u+self.ne.energy()
        if(self.nw and not self.nw.isleafnode()):
            u = u+self.nw.energy()
        if(self.se and not self.se.isleafnode()):
            u = u+self.se.energy()
        if(self.sw and not self.sw.isleafnode()):
            u = u+self.sw.energy()

        if self.ne:
            if self.nw: u = u + 0.5*self.ne.wt*self.nw.wt/self.ne.cg.dist(self.nw.cg)
            if self.sw: u = u + 0.5*self.ne.wt*self.sw.wt/self.ne.cg.dist(self.sw.cg)
            if self.se: u = u + 0.5*self.ne.wt*self.se.wt/self.ne.cg.dist(self.se.cg)
        
        if self.nw:
            if self.ne: u = u + 0.5*self.nw.wt*self.ne.wt/self.nw.cg.dist(self.ne.cg)
            if self.sw: u = u + 0.5*self.nw.wt*self.sw.wt/self.nw.cg.dist(self.sw.cg)
            if self.se: u = u + 0.5*self.nw.wt*self.se.wt/self.nw.cg.dist(self.se.cg)

        if self.se:
            if self.nw: u = u + 0.5*self.se.wt*self.nw.wt/self.se.cg.dist(self.nw.cg)
            if self.sw: u = u + 0.5*self.se.wt*self.sw.wt/self.se.cg.dist(self.sw.cg)
            if self.ne: u = u + 0.5*self.se.wt*self.ne.wt/self.se.cg.dist(self.ne.cg)
        
        if self.sw:
            if self.ne: u = u + 0.5*self.sw.wt*self.ne.wt/self.sw.cg.dist(self.ne.cg)
            if self.nw: u = u + 0.5*self.sw.wt*self.nw.wt/self.sw.cg.dist(self.nw.cg)
            if self.se: u = u + 0.5*self.sw.wt*self.se.wt/self.sw.cg.dist(self.se.cg)

        return u

class QuadTree:
    def __init__(self,wts,pts):
        #assumes pts to be in (-1,1)*(1,1)
        self.wts=np.asanyarray(wts)
        self.pts=pts

    @staticmethod
    def getcg(wts,pts):
        wt = np.sum(wts)
        ptx = np.sum(np.asarray([p.x*q for p,q in zip(pts,wts)]))
        pty = np.sum(np.asarray([p.y*q for p,q in zip(pts,wts)]))
        return wt,Point(ptx/wt,pty/wt)
    
    @staticmethod
    def getPotentialEnergy(charges):
        topnode = createTree(np.asarray([c.q for c in charges]),
                              [c.p for c in charges])
        return topnode.energy()
    
def createTree(wts,pts,node=None,box=None, fig=None):
    wt,cg=QuadTree.getcg(wts, pts)
    if node==None:
        level=0
        node = Node(None,None,None,None,wt,cg,level,BoundingBox(-0.5,0.5,-0.5,0.5))
        box = node.box
    else:
        level=node.l
        node = Node(None,None,None,None,wt,cg,level+1,box)
    
    ptsne=[]
    wtsne=[]
    ptsnw=[]
    wtsnw=[]
    ptsse=[]
    wtsse=[]
    ptssw=[]
    wtssw=[]
    side = (box.xmax-box.xmin)/2.0
    for wt,pt in zip(wts,pts):
        if pt.x>box.xmin and pt.x>=box.xmin+side:
            if pt.y>box.ymin and pt.y>=box.ymin+side:
                ptsne.append(pt)
                wtsne.append(wt)
            else:
                ptsse.append(pt)
                wtsse.append(wt)
        else:
            if pt.y>box.ymin and pt.y>=box.ymin+side:
                ptsnw.append(pt)
                wtsnw.append(wt)
            else:
                ptssw.append(pt)
                wtssw.append(wt)
    
    if len(ptsne)>1: 
        node.ne=createTree(wtsne,ptsne,node,BoundingBox(box.xmin+side,box.xmax,box.ymin+side,box.ymax))
    else:
        if len(ptsne)>0:        
            wt,cg=QuadTree.getcg(wtsne, ptsne)
            node.ne=Node(None,None,None,None,wt,cg,level+1,BoundingBox(box.xmin+side,box.xmax,box.ymin+side,box.ymax))

    if len(ptsnw)>1: 
        node.nw=createTree(wtsnw,ptsnw,node,BoundingBox(box.xmin,box.xmin+side,box.ymin+side,box.ymax))
    else:   
        if len(ptsnw)>0:     
            wt,cg=QuadTree.getcg(wtsnw, ptsnw)
            node.nw=Node(None,None,None,None,wt,cg,level+1,BoundingBox(box.xmin,box.xmin+side,box.ymin+side,box.ymax))

    if len(ptsse)>1: 
        node.se=createTree(wtsse,ptsse,node,BoundingBox(box.xmin+side,box.xmax,box.ymin,box.ymin+side))
    else: 
        if len(ptsse)>0:       
            wt,cg=QuadTree.getcg(wtsse, ptsse)
            node.se=Node(None,None,None,None,wt,cg,level+1,BoundingBox(box.xmin+side,box.xmax,box.ymin,box.ymin+side))

    if len(ptssw)>1: 
        node.sw=createTree(wtssw,ptssw,node,BoundingBox(box.xmin,box.xmin+side,box.ymin,box.ymin+side))
    else:   
        if len(ptssw):     
            wt,cg=QuadTree.getcg(wtssw, ptssw)
            node.sw=Node(None,None,None,None,wt,cg,level+1,BoundingBox(box.xmin,box.xmin+side,box.ymin,box.ymin+side))

    return node 


if __name__=="__main__":
    #node = createTree([1.0,1.0,1.0,1.0],[Point(-0.25,0.25),Point(-0.25,-0.25),Point(0.25,-0.25),Point(0.25,0.25)])
    import matplotlib.pyplot as plt
    wts=[1.0,2.0,1.0,2.0,1.0,2.0,1.0,2.0]
    pts=[Point(-0.25,0.25),Point(-0.20,0.20),
        Point(-0.25,-0.25),Point(-0.20,-0.20),
        Point(0.25,-0.25),Point(0.20,-0.20),
        Point(0.25,0.25),Point(0.20,0.20)]
    node = createTree(wts,pts)
    energy = node.energy()
    charges = [Charge(x,p) for x,p in zip(wts,pts)]
    energy2 = Charge.getPotentialEnergy(charges,{})
    print("Energy = {}\tEnergy2 = {}".format(energy,energy2))
    pfl.showTree(node,wts,pts)
    plt.show()
    print("done")
    wts=[]
    pts=[]
    for i in range(200):
        wts.append(np.random.uniform(-5.0,5.0))
        pts.append(Point(np.random.uniform(-0.5,0.5),np.random.uniform(-0.5,0.5)))
    node = createTree(wts,pts)
    energy = node.energy()
    charges = [Charge(x,p) for x,p in zip(wts,pts)]
    energy2 = Charge.getPotentialEnergy(charges,{})
    print("Energy = {}\tEnergy2 = {}".format(energy,energy2))
    pfl.showTree(node,wts,pts)
    plt.show()
    
    
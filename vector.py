# File: vector.py
# Purpose: Support for common vector operations
# Author: Abhikalp Shekhar
import math

def verify(cond,msg):
    if not cond:
        raise Exception(msg)

class BoundingBox:
    def __init__(self,xmin,xmax,ymin,ymax):
        self.xmax=xmax
        self.xmin=xmin
        self.ymax=ymax
        self.ymin=ymin

class Point:
    def __init__(self,x,y,z=0):
        if not isinstance(x,(float,int)):
            raise Exception("Bad x coordinate")
        if not isinstance(y,(float,int)):
            raise Exception("Bad y coordinate")
        self.x=x
        self.y=y
        self.z=z
    def dist(self,point):
        return math.sqrt((self.x - point.x)**2 + (self.y - point.y)**2 + (self.z - point.z)**2)

class vector:
    def __init__(self,xxx,yyy,zzz):
        if not isinstance(xxx,(float,int)):
            raise Exception("Bad x input")
        if not isinstance(yyy,(float,int)):
            raise Exception("Bad y input")
        if not isinstance(zzz,(float,int)):
            raise Exception("Bad z input")
            
        self.x = xxx
        self.y = yyy
        self.z = zzz

    def __str__(self):
        return "{}i + {}j + {}k".format(self.x,self.y,self.z)
        
    def length(self):
        l = math.sqrt(self.x**2 + self.y**2 + self.z**2)
        return l

    def scale(self,d):
        return vector(self.x*d, self.y*d, self.z*d)

    def unitvec(self):
        l = self.length()
        return self.scale(1.0/l)

    def __add__(self,b):
        return vector(self.x+b.x, self.y+b.y, self.z+b.z)

    def dot(self,b):
        return self.x*b.x + self.y*b.y + self.z*b.z

    def cross(self,b):
        return vector(self.y*b.z - self.z*b.y, self.z*b.x - self.x*b.z , self.x*b.y - self.y*b.x)

def polar(r,theta,z=0):
    return vector(r*math.cos(theta),r*math.sin(theta),z)

class ScalarField:
    def __init__(self,fieldFunc,perturb=1e-10):
        self.scalarFunc=fieldFunc
        self.perturb=perturb

    def gradient(self,p):
        pxdash=Point(p.x+self.perturb, p.y,p.z)
        pydash=Point(p.x,p.y+self.perturb,p.z)
        pzdash=Point(p.x,p.y,p.z+self.perturb)
        value=self.scalarFunc(p)
        vxdash=self.scalarFunc(pxdash) 
        vydash=self.scalarFunc(pydash)
        vzdash=self.scalarFunc(pzdash)
        return vector((vxdash-value)/self.perturb,
            (vydash-value)/self.perturb,
            (vzdash-value)/self.perturb)

    def getVal(self,p):
        return self.scalarFunc(p)

class VectorField:
    def __init__(self,fieldFuncs,perturb=1e-10):
        verify(isinstance(fieldFuncs,(list,tuple)) and len(fieldFuncs)==3, "vectorfunc size should be 3")
        self.scalarFuncs=fieldFuncs
        self.perturb=perturb
    
    def getValx(self,p):
        return self.scalarFuncs[0](p)
    
    def getValy(self,p):
        return self.scalarFuncs[1](p)
    
    def getValz(self,p):
        return self.scalarFuncs[2](p)
    
    def getValxDash(self,p,bump=1e-10):
        px=Point(p.x+bump,p.y,p.z)
        return (self.getValx(px)-self.getValx(p))/bump
        
    def getValyDash(self,p,bump=1e-10):
        py=Point(p.x,p.y+bump,p.z)
        return (self.getValy(py)-self.getValy(p))/bump
    
    def getValzDash(self,p,bump=1e-10):
        pz=Point(p.x,p.y,p.z+bump)
        return (self.getValz(pz)-self.getValz(p))/bump
    
    def getVal(self,p):
        return vector(self.getValx(p),self.getValy(p),self.getValz(p))
    
    def divergence(self,p,bump=1e-10):
        return self.getValxDash(p,bump)+self.getValyDash(p,bump)+self.getValzDash(p,bump)

    def curl(self,p,bump=1e-10):
        vx,vy,vz=(self.getValxDash(p,bump),self.getValyDash(p,bump),self.getValzDash(p,bump))
        return vector(vz - vy, vx - vz , vy - vx)        




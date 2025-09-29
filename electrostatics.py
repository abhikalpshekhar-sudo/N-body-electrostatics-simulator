import math
from vector import *
import copy
import numpy as np
import fastmultipole as fmm

class Charge:
    def __init__(self,q,p,units=1.0):
        if not isinstance(q,(float,int)):
            raise Exception("Bad charge")
        if not isinstance(p,Point):
            raise Exception("Bad Point")
        self.q=q
        self.p=p
        self.units=units

    def ElectricField(self,point):
        c = self.units
        d = self.p.dist(point) 
        d2 = d*d
        ef= c*self.q/d2
        unitVect = vector((-self.p.x + point.x)/d, (-self.p.y + point.y)/d, (-self.p.z + point.z)/d)
        return unitVect.scale(ef)

    def ElectricPotential(self,point):
        c = self.units 
        d = self.p.dist(point) 
        return c*self.q/d

    @staticmethod
    def getPotentialEnergy(charges,option):
        if option.get("method","exact")=="exact":
            n=len(charges)
            U = 0.0
            for i in range(n):
                for j in range(i+1,n):
                    dist = charges[i].p.dist(charges[j].p)
                    if abs(dist)>1e-16:
                        U = U + charges[i].q*charges[j].q/dist 
            return U
        else:
            return fmm.QuadTree.getPotentialEnergy(charges)   

    @staticmethod
    def getPotentialAtPoint(charges,p,option=None):
        if option is None:option = {}
        chs=copy.deepcopy(charges)
        chs.append(Charge(1.0,p))
        return Charge.getPotentialEnergy(chs,option) - Charge.getPotentialEnergy(charges,option) 

if __name__ == "__main__":
    
    def singleChargeTest():
        ch = Charge(1.0, Point(0.0,0.0,0.0))
        fld1= ch.ElectricField(Point(1.0,1.0,0.0))
        print(fld1)

    def multipleChargeTest():
        l = [Charge(3.0,Point(4,0)), Charge(3.0,Point(0.0,0.0)), Charge(3.0,Point(0,4))]
        v = vector(0,0,0)
        p = Point(3.0,4.0)
        potential = 0.0
        for c in l:
            fld = c.ElectricField(p)
            v = v + fld
            potential = potential + c.ElectricPotential(p)

        print("Electric Field Vector Sum = {}".format(v))
        print("Electric Potential = {}".format(potential))
        def fldfunc(pnt):
            return Charge.getPotentialAtPoint(l,pnt)
        scalarfld = ScalarField(fldfunc)
        potential2 = scalarfld.getVal(p)
        print("ELectric Potential2 = {}".format(potential2))
        elecfld = scalarfld.gradient(p).scale(-1.0)
        print("Electric Field Gradient of Potential = {}".format(elecfld))

    def dipoleTest():
        l = [Charge(1.0,Point(-0.1,0)),Charge(-1.0,Point(0.1,0)) ]
        v = vector(0,0,0)
        p = Point(0.0,1.0)
        for c in l:
            fld = c.ElectricField(p)
            v = v + fld
        print(v)

    singleChargeTest()
    multipleChargeTest()
    dipoleTest()

import plotfield as plfld
from vector import *
from electrostatics import *

def plotMonoPoleField(q=1.0):
    plfld.plot2DChargeAndFields([Charge(q,Point(0, 0))],options = {'nxpoints':50,'nypoints':50})
def plotDipoleField(q=1.0,d=0.1):
    plfld.plot2DChargeAndFields([Charge(-q,Point(-d/2.0, 0)), Charge(q,Point(d/2.0, 0))],options = {'nxpoints':50,'nypoints':50})
def plotQuadPoleField(q=1.0,d=0.1):
    plfld.plot2DChargeAndFields([Charge(-q,Point(-d/2.0, 0)),Charge(2*q,Point(0, 0)),Charge(-q,Point(d/2.0, 0))],
                                options = {'nxpoints':50,'nypoints':50})
def plotLattice(q=0.1,d=0.1):
    charges = []
    sign = 1
    for i in range(4):
        for j in range(4):
            sign *= -1
            charges.append(Charge(1.0*sign,Point((i+1)*d,(j+1)*d)))
    plfld.plot2DChargeAndFields(charges,options = {'nxpoints':50,'nypoints':50})


if __name__=="__main__":
    plotMonoPoleField()
    plotDipoleField()
    plotQuadPoleField()
    plotLattice(q=0.1,d=0.1)
    print("done")
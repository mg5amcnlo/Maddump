if '__main__' == __name__:
    import sys
    sys.path.append('../../')
    sys.path.append('../../madgraph/various')

import numpy as np
import cmath
import os
import matplotlib.pyplot as plt
from scipy.integrate import quad, dblquad
from meshfitter2D import *
import multiprocessing
import fit2D_card as fit2D
import bremsstrahlung_card as bremss
import models.check_param_card as param_card_mod

class fit2D_z_pt2(CellHistogram):
    
    def __init__(self,pdgcode):

        self.bremss_card = bremss.BremsstrahlungCard(pjoin('../Cards/bremsstrahlung_card.dat'))
        self.fit2D_card = fit2D.Fit2DCard(pjoin('../Cards/fit2D_card.dat'))
        self.param_card = param_card_mod.ParamCard(pjoin('../Cards/param_card.dat'))
        
        self.P = self.bremss_card["pbeam"]
        self.npot = self.bremss_card["npot"]
        self.zmin = self.bremss_card["z_min"]
        self.zmax = self.bremss_card["z_max"]
        self.pt2min = self.bremss_card["pt2_min"] 
        self.pt2max = self.bremss_card["pt2_max"]
        self.mp2 = 0.938**2
        
        self.mV2 = self.param_card['mass'].get(pdgcode).value**2

        self.testplot = True 
        self.data = []
        self.ncores = self.fit2D_card['ncores']
        self.flux_norm = 0.
        
        self.pdgcode = ' '+str(pdgcode)+' ' 
        jobs = []
        k=0
        qlist=[]
        npoints = self.bremss_card["nfit"]
        q = multiprocessing.Queue()
        for i in range(self.ncores):
            #qlist.append(q)
            p = multiprocessing.Process(target= getattr(self, 'get_data'), args=( q,npoints/self.ncores, ) )
            jobs.append(p)
            p.start()
            k +=1

        while k > 0:
            if not (q.empty()):
                self.data.append(q.get())
                k-=1

        data =[]
        for x in self.data:
            for y in x:
                data.append(y)
        
        super(fit2D_z_pt2,self).__init__(Point(self.zmin,self.pt2min), \
                                                self.zmax-self.zmin,self.pt2max-self.pt2min,self.bremss_card["nexit"])
        self.add_pts(data)
        
        
    def get_data(self,q,npoints):
        print(npoints)
        from random import seed,uniform
        seed()
        data = []
#        with open('2dplot.dat','w') as fout:
        for i in range(npoints):
            z = uniform(self.zmin,self.zmax)
            pt2 = uniform(self.pt2min,self.pt2max)
            data.append(WeightedPoint(z,pt2,self.dnV_dzdpt2(z,pt2,self.mV2)))
                
#                fout.write(str(z)+'\t'+str(pt2)+'\t'+str(dnV_dzdpt2(z,pt2,self.mV2))+'\n')
        q.put(data)
        
    def do_fit(self):
        # Generate the 2DMesh, output file: cell_fortran.dat
        self.fit(ncores=self.ncores)
        os.system('mv cell_fortran.dat cell_fortran_z_pt2.dat')
        if self.testplot:
            plot('cell_fortran_z_pt2.dat','mesh2D.png')
        self.flux_norm = self.nV() * self.npot
        self.fit2D_card["flux_norm"] = self.flux_norm
        self.fit2D_card.write('../Cards/fit2D_card.dat',template='../Cards/fit2D_card.dat')

            
    def do_unweight(self, ngen):
        from random import randint,uniform
        z0,pt20,zwidth,pt2width = np.loadtxt('cell_fortran_z_pt2.dat',unpack=True)
        z1=z0+zwidth
        pt21=pt20+pt2width
        
        ncell = len(z0)

        zgen=[]
        pt2gen=[]

        with open('bremsstrahlung.hepmc','w') as fout:
            header = '\nHepMC::Version 2.06.09\nHepMC::IO_GenEvent-START_EVENT_LISTING\n'
            fout.write(header)
            npup=1
            istart=10001                    
            for j in range(ngen):
                i = randint(0,ncell-1)
                zgen = uniform(z0[i],z1[i])
                pt2gen = uniform(pt20[i],pt21[i])
                thetagen = uniform(0,2.*np.pi)
                EV = zgen*self.P + (pt2gen+self.mV2)/(2.*self.P*zgen)
                px = np.sqrt(pt2gen)*np.cos(thetagen)
                py = np.sqrt(pt2gen)*np.sin(thetagen)
                pz = zgen*self.P
                fout.write('E 0 -1 -1.0000000000000000e+00 -1.0000000000000000e+00 -1.0000000000000000e+00 0 0 1 0 0 0 0\nU GEV CM\n')
                fout.write('V -1 0 0 0 0 0 0 '+str(npup)+' 0\n')
                line='P '+str(istart)+self.pdgcode+str(px)+' '+str(py)+' '+\
                    str(pz)+' '+ str(EV)+' '+str(np.sqrt(self.mV2))+' 1 0 0 0 0\n'
                fout.write(line)
                            
            fout.write('HepMC::IO_GenEvent-END_EVENT_LISTING\n')
        
        return 1.

    def wba(self,z,pt2):
        kV = 1.
        H = pt2 + (1.-z)*self.mV2 + z**2*self.mp2
        return kV/(2.*np.pi*H) * ( (1.+(1.-z)**2)/z
                        -2.*z*(1.-z)*( (2.*self.mp2+self.mV2)/H - z**2*2.*self.mp2**2/H**2 )
                        +2.*z*(1.-z)*(z+(1.-z)**2)*self.mp2*self.mV2/H**2
                        +2.*z*(1.-z)**2*self.mV2**2/H**2 )

    def F1p2(self,q2):
        # Returns the F1(Q^2) proton form-factor for time-like Q^2
        # the parametrization and the fit is taken from 

        if(q2<0):
            print('Error: Space-like form factor not implemented!')
            exit(-1)

        mrho = 0.770
        mrho2 = mrho**2
        Garho = 0.150
        grho = 5.03
        f1rhoNN = 3.10117
        f2rhoNN = 20.9418

        mrhop = 1.250
        mrhop2 = mrhop**2
        Garhop = 0.3
        grhop = grho
        f1rhopNN = 1.12272
        f2rhopNN = -22.2148

        mrhopp = 1.450
        mrhopp2 = mrhopp**2
        Garhopp = 0.5
        grhopp = grho
        f1rhoppNN = -1.70888
        f2rhoppNN = 10.6037

        momega = mrho
        momega2 = momega**2
        Gaomega = 0.0085
        gomega = 17.1
        f1omegaNN = 17.301
        f2omegaNN = -2.2966

        momegap = mrhop
        momegap2 = momegap**2
        Gaomegap = Garhop
        gomegap = gomega
        f1omegapNN = -15.0763
        f2omegapNN = 2.64631

        momegapp = mrhopp
        momegapp2 = momegapp**2
        Gaomegapp = Garhopp
        gomegapp = gomega
        f1omegappNN = 6.32533
        f2omegappNN = -1.26315

        F1rho = f1rhoNN/grho * mrho2/(mrho2-1j*mrho*Garho - q2) + f1rhopNN/grhop * mrhop2/(mrhop2-1j*mrhop*Garhop - q2) + f1rhoppNN/grhopp * mrhopp2/(mrhopp2-1j*mrhopp*Garhopp - q2)
        F1omega = f1omegaNN/gomega * momega2/(momega2-1j*momega*Gaomega - q2) + f1omegapNN/gomegap * momegap2/(momegap2-1j*momegap*Gaomegap - q2) + f1omegappNN/gomegapp * momegapp2/(momegapp2-1j*momegapp*Gaomegapp - q2)

        F1p = F1rho+F1omega

        # F2rho = f2rhoNN/grho * mrho2/(mrho2-1j*mrho*Garho - q2) + f2rhopNN/grhop * mrhop2/(mrhop2-1j*mrhop*Garhop - q2) + f2rhoppNN/grhopp * mrhopp2/(mrhopp2-1j*mrhopp*Garhopp - q2)
        # F2omega = f2omegaNN/gomega * momega2/(momega2-1j*momega*Gaomega - q2) + f2omegapNN/gomegap * momegap2/(momegap2-1j*momegap*Gaomegap - q2) + f2omegappNN/gomegapp * momegapp2/(momegapp2-1j*momegapp*Gaomegapp - q2)

        # F2p = F2rho+F2omega
        # G1pM = F1p+F2p

        return (F1p*F1p.conjugate()).real


    def sig_pp(self,s,z,pt2):
        Z = 35.45 
        B = 0.308
        Y1 = 42.53
        Y2 = 33.34
        s0 = 5.38**2
        s1 = 1.
        eta1 = 0.458
        eta2 = 0.545

        return Z + B*(np.log(s/s0))**2 + Y1*np.exp(eta1*np.log(s1/s)) - Y2*np.exp(eta2*np.log(s1/s)) 


    def dnV_dzdpt2(self,z,pt2,mV2):
        Ep = self.P + self.mp2/2./self.P
        return self.sig_pp(2.*np.sqrt(self.mp2)* (Ep-z*self.P-(pt2+mV2)/(2.*z*self.P) ),z,pt2 ) / self.sig_pp(2.*np.sqrt(self.mp2)*Ep,z,pt2) * self.F1p2(mV2) *self.wba(z,pt2)

    def nV(self):
        res = dblquad(lambda pt2,z: self.dnV_dzdpt2(z,pt2,self.mV2), self.zmin, self.zmax, lambda pt2: self.pt2min, lambda pt2: self.pt2max )
        return res[0]


    
# hist2D_z_pt2 = fit2D_z_pt2(400., 0.5, zmin=0.1, zmax=0.9, pt2min=1e-14, pt2max=1., npoints=100000)
# hist2D_z_pt2.do_fit()
# hist2D_z_pt2.do_unweight(100000)

    
#PLOT NV as function of mV 
# t=np.linspace(0.01,3.,100)
# res=np.zeros(len(t))

# for i in range(len(t)):
#     res[i]=nV(t[i]**2)

# plt.xscale("log")

# #plt.plot(t,np.sqrt(F1N2(t)),'-')
# plt.plot(t,1.44e18*1e-12*res/137.,'-')
# plt.show()

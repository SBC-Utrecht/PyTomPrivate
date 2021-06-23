import unittest


class pytom_ScoreTest(unittest.TestCase):

    def setUp(self):
        """set up"""
        from pytom_volume import vol, initSphere
        from pytom.basic.structures import WedgeInfo
        from pytom.simulation.SimpleSubtomogram import simpleSimulation
        from pytom.basic.structures import Wedge, Mask
        from pytom.basic.files import read
        from pytom_numpy import npy2vol

        self.wedge = 30
        self.shift = [-1, 2, 3]
        self.rotation = [0, 0, 0]
        # create sphere
        #self.v = vol(32,32,32)
        #self.mask = vol(32,32,32)
        #initSphere(self.v,10,2,0,15,15,15)
        # there is a slight inconsistency when smoothing > 0 -
        # cleaner implementation would be multipliction with sqrt(mask) in corr function
        #initSphere(self.mask,13,0,0,16,16,16)
        self.wi = WedgeInfo(wedgeAngle=self.wedge, rotation=[10.0,20.0,30.0], 
              cutoffRadius=0.0)
        #self.s = simpleSimulation( volume=self.v, rotation=self.rotation,
              #shiftV=self.shift, wedgeInfo=self.wi, SNR=10.)

        #Ribo template settings
        self.s = read("./testData/ribo.em")
        self.v = vol(100, 100, 100)
        self.mask = Mask("./testData/ribo_mask.em")
    """
    def test_xcfScore(self):
        from pytom.score.score import xcfScore as score
        from pytom_volume import peak

        sc = score()
        # consistency of scoring coefficient and scoring function - difference 
        # due to sub-pixel accuracy for score
        c  = sc.scoringCoefficient( self.s, self.v)
        cf = sc.scoringFunction( self.s, self.v)
        p= peak(cf)
        #print(c)
        #print(cf.getV(16,16,16))
        self.assertAlmostEqual( first=c, second=cf.getV(16,16,16), places=2, 
            msg='Scoring coefficient and scoring function XCF inconsistent')
        self.assertLess( c, cf.getV(p[0],p[1],p[2]), 
            'Scoring coefficient and scoring function XCF inconsistent')
    
    def test_nxcfScore(self):
        
        #test nxcf score
        
        from pytom.score.score import nxcfScore as score
        from pytom_volume import peak

        sc = score()
        # check auto-correlation coefficient
        c  = sc.scoringCoefficient( self.s, self.s)
        self.assertAlmostEqual( first=c,  second=1., places=4, 
             msg='NXCFScore: Auto-correlation not == 1')
        c  = sc.scoringCoefficient( self.s, self.s, self.mask)
        self.assertAlmostEqual( first=c,  second=1., places=4, 
             msg='NXCFScore: Auto-correlation with mask not == 1')
        # check auto-correlation function
        cf  = sc.scoringFunction( self.s, self.s, self.mask)
        self.assertAlmostEqual( first=c, second=cf.getV(16,16,16), places=2, 
            msg='Scoring coefficient and function NXCF inconsistent for auto-corr')
        # consistency of scoring coefficient and scoring function - difference 
        # due to sub-pixel accuracy for score
        c  = sc.scoringCoefficient( self.s, self.v)
        cf = sc.scoringFunction( self.s, self.v)
        p= peak(cf)
        self.assertAlmostEqual( first=c, second=cf.getV(16,16,16), places=2, 
            msg='Scoring coefficient and scoring function NXCF inconsistent')
        self.assertLess( c, cf.getV(p[0],p[1],p[2]), 
            'Scoring coefficient and scoring function NXCF inconsistent')
        # now check mask
        c  = sc.scoringCoefficient( self.s, self.v, self.mask)
        cf = sc.scoringFunction( self.s, self.v, self.mask)
        p= peak(cf)
        pval = cf.getV(p[0],p[1],p[2])
    """
    def test_flcfScore(self):
        """
        test FLCF score
        """
        from pytom.score.score import FLCFScore as score
        from pytom_volume import peak
        
        sc = score()
        # check auto-correlation coefficient
        c  = sc.scoringCoefficient( self.s, self.s)
        self.assertAlmostEqual( first=c,  second=1., places=5, 
             msg='SocScore: Auto-correlation not == 1')
        # consistency of scoring coefficient and scoring function - difference 
        # due to sub-pixel accuracy for score
        c  = sc.scoringCoefficient( self.s, self.v)
        cf = sc.scoringFunction( self.s, self.v)
        p= peak(cf)
        self.assertAlmostEqual( first=c, second=cf.getV(p[0],p[1],p[2]), places=2, 
            msg='Scoring coefficient and scoring function FLCF inconsistent')
    """   
    def test_socScore(self):
        
        #second order correlation score
        
        from pytom.score.score import SOCScore as score
        from pytom_volume import peak

        sc = score()
        # check auto-correlation coefficient
        c  = sc.scoringCoefficient( self.s, self.s)
        self.assertAlmostEqual( first=c,  second=1., places=3, 
             msg='SocScore: Auto-correlation not == 1')
        # consistency of scoring coefficient and scoring function - difference 
        # due to sub-pixel accuracy for score
        c  = sc.scoringCoefficient( self.s, self.v)
        cf = sc.scoringFunction( self.s, self.v)
        p= peak(cf)
        self.assertAlmostEqual( first=c, second=cf.getV(p[0],p[1],p[2]), places=1, 
            msg='Scoring coefficient and scoring function SOC inconsistent')
    """
    def test_pofScore(self):
        """
        Test Phase Only Filter correlation function
        @author: Maria Cristina Trueba
        """
        from pytom.score.score import POFScore as score
        from pytom_volume import peak

        sc = score()
        #Check auto-correlation coefficient
        c = sc.scoringCoefficient(self.s, self.s)
        print("pofScore autocorrelation in cpu is, ", c)
        self.assertAlmostEqual(first=c, second=1, places=5, msg="POFScore: Autocorrelation not == 1")

        # #consistency of scoring coefficient and scoring function - difference due to sub-pixel accuracy for score
        # c = sc.scoringCoefficient(self.v, self.s)
        # cf = sc.scoringFunction(self.v, self.s)
        # p = peak(cf)
        # self.assertAlmostEqual( first = c, second = cf.getV(p[0], p[1], p[2]), places=2, msg = "Scoring coefficient and scoring funtion POF inconsistent")

    def test_mcfScore(self):
        """
        Test Mutual correlation function
        @author: Maria Cristina Trueba
        """
        from pytom.score.score import MCFScore as score
        from pytom_volume import peak

        sc = score()
        #Check auto-correlation coefficient
        c = sc.scoringCoefficient(self.s, self.s)
        print("mcfScore autocorrelation in cpu is, ", c)
        self.assertAlmostEqual(first=c, second=1, places=1, msg="MCFScore: Autocorrelation not == 1")

        #consistency of scoring coefficient and scoring function - difference due to sub-pixel accuracy for score
        c = sc.scoringCoefficient(self.s, self.v)
        cf = sc.scoringFunction(self.s, self.v)
        p = peak(cf)
        self.assertAlmostEqual( first = c, second = cf.getV(p[0], p[1], p[2]), places=1, msg = "Scoring coefficient and scoring funtion MCF inconsistent")
    """
    def test_pofScore_np(self):
        
        #Test Phase only filter function numpy version for gpu
        #@author: Maria Cristina Trueba
        
        import numpy as np
        from pytom.tompy.score import POFScore as score
        from pytom_numpy import vol2npy
        from pytom.tompy.correlation import subPixelMax3D

        s = vol2npy(self.s).copy()
        v = vol2npy(self.v).copy()
        sc = score()
        #Check auto-correlation coefficient
        c = sc.scoringCoefficient(s, s)
        print("pofScore_np autocorrelation is, ", c)
        self.assertAlmostEqual(first=c, second=1, places=5, msg="POFScore: Autocorrelation not == 1")

        #consistency of scoring coefficient and scoring function - difference due to sub-pixel accuracy for score
        c = sc.scoringCoefficient(s, v)
        cf = sc.scoringFunction(s, v)
        p = cf.max()
        self.assertAlmostEqual( first = c, second = p, places=2, msg = "Scoring coefficient and scoring funtion POF inconsistent")

    
    def test_mcfScore_np(self):
        
        #Test Mutual correlation function numpy version for gpu
        #@author: Maria Cristina Trueba
        
        from pytom.tompy.score import MCFScore as score
        from pytom_numpy import vol2npy

        s = vol2npy(self.s).copy()
        v = vol2npy(self.v).copy()
        sc = score()

        #Check auto-correlation coefficient
        c = sc.scoringCoefficient(s, s)
        print("mcfScore_np autocorrelation is, ", c)
        self.assertAlmostEqual(first=c, second=1, places=2, msg= "POFScore: Autocorrelation not == 1")

        #consistency of scoring coefficient and scoring function - difference due to sub-pixel accuracy for score
        c = sc.scoringCoefficient(s, v)
        cf = sc.scoringFunction(s, v)
        p = cf.max()
        self.assertAlmostEqual( first = c, second = p, places=1, msg = "Scoring coefficient and scoring funtion POF inconsistent")
    """
    def test_pof_gpu(self):
        """
        Test POF function in gpu
        @author: Maria Cristina Trueba Sanchez
        """

        from pytom_numpy import vol2npy
        from pytom.localization.extractPeaks import templateMatchingGPU
        from pytom.basic.structures import Wedge, Mask
        from pytom.tompy.correlation import POF as sFunc
        import numpy as np
        from pytom.tompy.io import read

        #Data test
        s = read("./testData/ribo.em")
        v = vol2npy(self.v).copy()
        w = Wedge(self.wedge)
        m = self.mask
        m = Mask("./testData/ribo_mask.em")
        m = m.getVolume()
        m = vol2npy(m).copy()

        #Random object autocorrelation
        # s = vol2npy(self.s)
        # m = vol2npy(self.mask)
        # w = Wedge(self.wedge)

        result = templateMatchingGPU(volume=s, reference=s, rotations=[[0, 0, 0]], scoreFnc=sFunc, mask=m,
                                     maskIsSphere=True, wedgeInfo=w, padding=True, jobid=0, gpuID=1)
        c, a, n, k = result

        #Get the highest peak and compare with == 1
        c=c.max()
        print("pof_gpu autocorrelation is, ", c)
        self.assertAlmostEqual(first=c, second=1, places=4, msg= "POF in GPU Autocorrelation not == 1")

    def test_mcf_gpu(self):
        """
        Test POF function in gpu
        @author: Maria Cristina Trueba Sanchez
        """

        from pytom_numpy import vol2npy
        from pytom.localization.extractPeaks import templateMatchingGPU
        from pytom.basic.structures import Wedge, Mask
        from pytom.tompy.correlation import MCF as sFunc
        from pytom.tompy.io import read

        #Ribosome autocorrelation
        v = vol2npy(self.v).copy()
        s = read("./testData/ribo.em")
        w = Wedge(self.wedge)
        m = Mask("./testData/ribo_mask.em")
        m = m.getVolume()
        m = vol2npy(m).copy()

        #Random object autocorrelation
        # s = vol2npy(self.s)
        # m = vol2npy(self.mask)
        # w = Wedge(self.wedge)

        result = templateMatchingGPU(volume=s, reference=s, rotations=[[0, 0, 0]], scoreFnc=sFunc, mask=m,
                                     maskIsSphere=True, wedgeInfo=w, padding=True, jobid=0, gpuID=1)
        c, a, n, k = result

        #Get the highest peak and compare with == 1
        c=c.max()
        print("mcf_gpu autocorrelation is, ", c)
        self.assertAlmostEqual(first=c, second=1, places=4, msg= "MCF in GPU Autocorrelation not == 1")

    def test_flcf_gpu(self):
        """
        Test POF function in gpu
        @author: Maria Cristina Trueba Sanchez
        """

        from pytom_numpy import vol2npy
        from pytom.localization.extractPeaks import templateMatchingGPU
        from pytom.basic.structures import Wedge, Mask
        from pytom.tompy.correlation import FLCF as sFunc
        from pytom.tompy.io import read

        #Ribo template autocorrelation
        v = vol2npy(self.v).copy()
        s = read("./testData/ribo.em")
        w = Wedge(self.wedge)
        m = Mask("./testData/ribo_mask.em")
        m = m.getVolume()
        m = vol2npy(m).copy()

        #Random Object autocorrelation
        # s = vol2npy(self.s)
        # m = vol2npy(self.mask)
        # w = Wedge(self.wedge)

        result = templateMatchingGPU(volume=s, reference=s, rotations=[[0, 0, 0]], scoreFnc=sFunc, mask=m,
                                     maskIsSphere=True, wedgeInfo=w, padding=True, jobid=0, gpuID=1)
        c, a, n, k = result

        #Get the highest peak and compare with == 1
        c=c.max()
        print("flcf_gpu autocorrelation is, ", c)
        self.assertAlmostEqual(first=c, second=1, places=4, msg= "MCF in GPU Autocorrelation not == 1")

    def RScore_Test(self):
        """
        """
        
    def runTest(self):
        self.test_xcfScore()
        # self.test_nxcfScore()
        self.test_flcfScore()
        # self.test_socScore()
        self.test_pofScore()
        self.test_mcfScore()
        # self.test_pofScore_np()
        # self.test_mcfScore_np()
        self.test_pof_gpu()
        self.test_mcf_gpu()

if __name__ == '__main__':
    unittest.main()

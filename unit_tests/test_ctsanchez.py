import unittest


class pytom_ScoreTest(unittest.TestCase):

    def setUp(self):
        """set up"""
        from pytom_volume import vol, initSphere
        from pytom.basic.structures import WedgeInfo
        from pytom.simulation.SimpleSubtomogram import simpleSimulation
        from pytom.basic.structures import Wedge, Mask
        from pytom.basic.files import read, write_em

        from pytom_numpy import npy2vol

        self.wedge = 30
        self.wedge_zero = 0
        self.shift = [-1, 2, 3]
        self.rotation = [0, 0, 0]

        # create sphere
        self.v = vol(32,32,32)
        self.mask = vol(32,32,32)
        initSphere(self.v,10,2,0,15,15,15)
        # there is a slight inconsistency when smoothing > 0 -
        # cleaner implementation would be multipliction with sqrt(mask) in corr function
        initSphere(self.mask,16,0,0,16,16,16)
        self.wi = WedgeInfo(wedgeAngle=self.wedge,
              cutoffRadius=0.0)
        self.wi_zero = WedgeInfo(wedgeAngle=self.wedge_zero,
              cutoffRadius=0.0)
        self.s_cpu = simpleSimulation( volume=self.v, rotation=self.rotation,
                                       shiftV=self.shift, wedgeInfo=self.wi, SNR=10.)
        self.s_gpu = simpleSimulation( volume=self.v, rotation=self.rotation,
                                       shiftV=self.shift, wedgeInfo=self.wi_zero, SNR=10)

        #write_em("/home/ctsanchez/Desktop/template_random_test.mrc", self.s_gpu)

        #Ribo template settings
        self.s_ribo = read("./testData/ribo.em")
        self.mask_ribo = Mask("./testData/ribo_mask.em")
        self.ribo_cpu = simpleSimulation( volume=self.s_ribo, rotation=self.rotation, shiftV=self.shift, wedgeInfo=self.wi, SNR=10)
        self.ribo_gpu = simpleSimulation(volume=self.s_ribo, rotation=self.rotation, shiftV=self.shift, wedgeInfo=self.wi_zero,
                            SNR=10)

        #write_em("/home/ctsanchez/Desktop/template_ribo_test.mrc", self.ribo_cpu)
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

    def test_flcfScore_random(self):
        """
        test FLCF score with random generated sphere
        """
        from pytom.score.score import FLCFScore as score
        from pytom_volume import peak

        print("Starting FLCF CPU random test...")
        #Test with random sphre
        sc = score()
        # check auto-correlation coefficient
        c  = sc.scoringCoefficient( self.s_cpu, self.s_cpu)
        print("CPU FLCF random sphere score is ", c)
        self.assertAlmostEqual( first=c,  second=1., places=5,
             msg='FLCFScore: Auto-correlation  random sphere not == 1')
        # consistency of scoring coefficient and scoring function - difference
        # due to sub-pixel accuracy for score
        c  = sc.scoringCoefficient( self.v, self.s_cpu)
        cf = sc.scoringFunction( self.v, self.s_cpu)
        p= peak(cf)
        self.assertAlmostEqual( first=c, second=cf.getV(p[0],p[1],p[2]), places=2,
            msg='Scoring coefficient and scoring function FLCF inconsistent')

    def test_flcfScore_ribo(self):
        """
        test FLCF score with ribosome template
        """
        from pytom.score.score import FLCFScore as score
        from pytom_volume import peak

        print("Starting FLCF CPU ribo test...")
        #test with ribo
        sc = score()
        c = sc.scoringCoefficient(self.ribo_cpu, self.ribo_cpu)
        print("CPU FLCF ribo score is ", c)
        self.assertAlmostEqual(first=c, second=1., places=1,
                               msg='FLCFScore: Auto-correlation ribosome not == 1')
        c = sc.scoringCoefficient(self.s_ribo, self.ribo_cpu)
        cf = sc.scoringFunction(self.s_ribo, self.ribo_cpu)
        p = peak(cf)
        self.assertAlmostEqual(first=c, second=cf.getV(p[0], p[1], p[2]), places=2,
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
    def test_pofScore_random(self):
        """
        Test Phase Only Filter correlation function
        @author: Maria Cristina Trueba
        """
        from pytom.score.score import POFScore as score
        from pytom_volume import peak

        print("Starting POF CPU random test...")
        #Random sphere
        sc = score()
        #Check auto-correlation coefficient
        c = sc.scoringCoefficient(self.s_cpu, self.s_cpu)
        print("POF CPU random sphere score is, ", c)
        self.assertAlmostEqual(first=c, second=1, places=1, msg="POFScore: Autocorrelation random sphere not == 1")
        #consistency of scoring coefficient and scoring function - difference due to sub-pixel accuracy for score

        ref = self.wi.apply(self.v)
        c = sc.scoringCoefficient(self.s_cpu, ref)
        cf = sc.scoringFunction(self.s_cpu, ref)
        p = peak(cf)
        print("object and reference POF CPU random sphere", c)
        self.assertAlmostEqual( first = c, second = cf.getV(p[0], p[1], p[2]), places=2, msg = "Scoring coefficient and scoring funtion POF inconsistent")

    def test_pofScore_ribo(self):
        """
        Test Phase Only Filter correlation function for ribosome template
        @author: Maria Cristina Trueba
        """
        from pytom.score.score import POFScore as score
        from pytom_volume import peak

        print("Starting POF CPU ribo test...")
        #Ribo template
        sc = score()
        c = sc.scoringCoefficient(self.ribo_cpu, self.ribo_cpu)
        print("CPU POF ribo score is ", c)
        self.assertAlmostEqual(first=c, second=1., places=1,
                               msg='POFScore: Auto-correlation ribosome not == 1')
        
        # c = sc.scoringCoefficient(self.s_ribo, self.ribo_cpu)
        # cf = sc.scoringFunction(self.s_ribo, self.ribo_cpu)
        # p = peak(cf)
        # self.assertAlmostEqual(first=c, second=cf.getV(p[0], p[1], p[2]), places=2,
        #                        msg='Scoring coefficient and scoring function POF inconsistent')
    def test_mcfScore_random(self):
        """
        Test Mutual correlation function for random generated object
        @author: Maria Cristina Trueba
        """
        from pytom.score.score import MCFScore as score
        from pytom_volume import peak

        print("Starting MCF CPU random test...")
        sc = score()
        #Check auto-correlation coefficient
        c = sc.scoringCoefficient(self.s_cpu, self.s_cpu)
        print("CPU MCF random sphere score is, ", c)
        self.assertAlmostEqual(first=c, second=1, places=1, msg="MCFScore: Autocorrelation random sphere not == 1")

        #consistency of scoring coefficient and scoring function - difference due to sub-pixel accuracy for score
        c = sc.scoringCoefficient(self.s_cpu, self.v)
        cf = sc.scoringFunction(self.s_cpu, self.v)
        p = peak(cf)
        self.assertAlmostEqual( first = c, second = cf.getV(p[0], p[1], p[2]), places=1, msg = "Scoring coefficient and scoring funtion MCF inconsistent")

    def test_mcfScore_ribo(self):
        """
        Test Mutual correlation function for random generated object
        @author: Maria Cristina Trueba
        """
        from pytom.score.score import MCFScore as score
        from pytom_volume import peak

        print("Starting MCF CPU ribo test...")
        # Ribo template
        sc = score()
        c = sc.scoringCoefficient(self.ribo_cpu, self.ribo_cpu)
        print("CPU MCF ribo score is ", c)
        self.assertAlmostEqual(first=c, second=1., places=1,
                               msg='MCF Auto-correlation ribosome not == 1')
        c = sc.scoringCoefficient(self.s_ribo, self.ribo_cpu)
        cf = sc.scoringFunction(self.s_ribo, self.ribo_cpu)
        p = peak(cf)
        self.assertAlmostEqual(first=c, second=cf.getV(p[0], p[1], p[2]), places=2,
                               msg='Scoring coefficient and scoring function FLCF inconsistent')
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
    def test_pof_gpu_random(self):
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

        print("Starting POF GPU random test...")
        # Random object autocorrelation
        s = vol2npy(self.s_cpu)
        m = vol2npy(self.mask)

        w = Wedge(self.wedge)


        result = templateMatchingGPU(volume=s, reference=s.copy(), rotations=[[0, 0, 0]], scoreFnc=sFunc, mask=m,
                                     maskIsSphere=True, wedgeInfo=w, padding=True, jobid=0, gpuID=1)
        c, a, n, k = result
        #Get the highest peak and compare with == 1
        c=c.max()
        print("GPU POF random sphere score is, ", c)
        self.assertAlmostEqual(first=c, second=1, places=1, msg= "POF in GPU Autocorrelation random sphere not == 1")

    def test_pof_gpu_ribo(self):
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

        print("Starting POF GPU ribo test...")
        w = Wedge(self.wedge)

        # Data test
        s_ribo = vol2npy(self.ribo_gpu)
        m_ribo = vol2npy(self.mask_ribo.getVolume()).copy()

        result = templateMatchingGPU(volume=s_ribo, reference=s_ribo, rotations=[[0, 0, 0]], scoreFnc=sFunc, mask=m_ribo,
                                     maskIsSphere=True, wedgeInfo=w, padding=True, jobid=0, gpuID=1)
        c, a, n, k = result
        # Get the highest peak and compare with == 1
        c = c.max()
        print("GPU POF ribo score is, ", c)

        # self.assertAlmostEqual(first=c, second=1, places=1, msg="POF in GPU Autocorrelation ribosome not == 1")
        #
        # result = templateMatchingGPU(volume=s_ribo, reference=vol2npy(self.v).copy(), rotations=[[0, 0, 0]], scoreFnc=sFunc, mask=m_ribo,
        #                              maskIsSphere=True, wedgeInfo=w, padding=True, jobid=0, gpuID=1)[0]
        #
        # print("GPU POF ribo object/volume",result.max())



    def test_mcf_gpu_random(self):
        """
        Test POF function in gpu
        @author: Maria Cristina Trueba Sanchez
        """

        from pytom_numpy import vol2npy
        from pytom.localization.extractPeaks import templateMatchingGPU
        from pytom.basic.structures import Wedge, Mask
        from pytom.tompy.correlation import MCF as sFunc
        from pytom.tompy.io import read

        print("Starting MCF GPU random test...")
        w = Wedge(self.wedge)

        # Random object autocorrelation
        s = vol2npy(self.s_gpu)
        m = vol2npy(self.mask)

        result = templateMatchingGPU(volume=s, reference=s, rotations=[[0, 0, 0]], scoreFnc=sFunc, mask=m,
                                     maskIsSphere=True, wedgeInfo=w, padding=True, jobid=0, gpuID=1)
        c, a, n, k = result
        # Get the highest peak and compare with == 1
        c = c.max()
        print("GPU MCF random sphere score is, ", c)
        self.assertAlmostEqual(first=c, second=1, places=1, msg="MCF in GPU Autocorrelation random sphere not == 1")

    def test_mcf_gpu_ribo(self):
        """
        Test POF function in gpu
        @author: Maria Cristina Trueba Sanchez
        """

        from pytom_numpy import vol2npy
        from pytom.localization.extractPeaks import templateMatchingGPU
        from pytom.basic.structures import Wedge, Mask
        from pytom.tompy.correlation import MCF as sFunc
        from pytom.tompy.io import read

        print("Starting MCF GPU ribo test...")
        w = Wedge(self.wedge)
        # Data test
        s_ribo = vol2npy(self.ribo_gpu)
        m_ribo = vol2npy(self.mask_ribo.getVolume()).copy()

        result = templateMatchingGPU(volume=s_ribo, reference=s_ribo, rotations=[[0, 0, 0]], scoreFnc=sFunc,
                                     mask=m_ribo,
                                     maskIsSphere=True, wedgeInfo=w, padding=True, jobid=0, gpuID=1)
        c, a, n, k = result
        # Get the highest peak and compare with == 1
        c = c.max()
        print("GPU MCF ribo score is, ", c)
        self.assertAlmostEqual(first=c, second=1, places=1, msg="MCF in GPU Autocorrelation ribosome not == 1")

    def test_flcf_gpu_random(self):
        """
        Test POF function in gpu
        @author: Maria Cristina Trueba Sanchez
        """

        from pytom_numpy import vol2npy
        from pytom.localization.extractPeaks import templateMatchingGPU
        from pytom.basic.structures import Wedge, Mask
        from pytom.tompy.correlation import FLCF as sFunc
        from pytom.tompy.io import read

        print("Starting FLCF GPU random test...")
        w = Wedge(self.wedge)

        # Random object autocorrelation
        s = vol2npy(self.s_gpu)
        m = vol2npy(self.mask)

        result = templateMatchingGPU(volume=s, reference=s, rotations=[[0, 0, 0]], scoreFnc=sFunc, mask=m,
                                     maskIsSphere=True, wedgeInfo=w, padding=True, jobid=0, gpuID=1)
        c, a, n, k = result
        # Get the highest peak and compare with == 1
        c = c.max()
        print("GPU FLCF random sphere score is, ", c)
        self.assertAlmostEqual(first=c, second=1, places=1, msg="FLCF in GPU Autocorrelation random sphere not == 1")

    def test_flcf_gpu_ribo(self):
        """
        Test POF function in gpu
        @author: Maria Cristina Trueba Sanchez
        """

        from pytom_numpy import vol2npy
        from pytom.localization.extractPeaks import templateMatchingGPU
        from pytom.basic.structures import Wedge, Mask
        from pytom.tompy.correlation import FLCF as sFunc
        from pytom.tompy.io import read

        print("Starting FLCF GPU ribo test...")
        w = Wedge(self.wedge)
        # Data test
        s_ribo = vol2npy(self.ribo_gpu)
        m_ribo = vol2npy(self.mask_ribo.getVolume()).copy()

        result = templateMatchingGPU(volume=s_ribo, reference=s_ribo, rotations=[[0, 0, 0]], scoreFnc=sFunc,
                                     mask=m_ribo,
                                     maskIsSphere=True, wedgeInfo=w, padding=True, jobid=0, gpuID=1)
        c, a, n, k = result
        # Get the highest peak and compare with == 1
        c = c.max()
        print("GPU FLCF ribo score is, ", c)
        self.assertAlmostEqual(first=c, second=1, places=1, msg="FLCF in GPU Autocorrelation ribosome not == 1")

    def RScore_Test(self):
        """
        """
        
    def runTest(self):
        #self.test_xcfScore()
        # self.test_nxcfScore()
        self.test_flcfScore_random()
        self.test_flcfScore_ribo()
        # self.test_socScore()
        self.test_pofScore_random()
        self.test_pofScore_ribo()
        self.test_mcfScore_random()
        self.test_mcfScore_ribo()
        # self.test_pofScore_np()
        # self.test_mcfScore_np()
        self.test_pof_gpu_random()
        self.test_mcf_gpu_random()
        self.test_flcf_gpu_random()
        self.test_pof_gpu_ribo()
        self.test_mcf_gpu_ribo()
        self.test_flcf_gpu_ribo()

if __name__ == '__main__':
    unittest.main()

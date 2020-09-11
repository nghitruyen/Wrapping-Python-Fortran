from __future__ import print_function
import unittest
import numpy as np
import calcul_volume

class TestExample(unittest.TestCase):

    def setUp(self):
        pass

    def test_auto_diff(self):
        #Setup vallues:
        d1 = calcul_volume.dual_num_auto_diff.DUAL_NUM()
        d1.x_ad_=1 
        d1.xp_ad_=np.array([1,0])
        d2 = calcul_volume.dual_num_auto_diff.DUAL_NUM()
        d2.x_ad_=10
        d2.xp_ad_=np.array([0,1])
        d3 = calcul_volume.mcyldnad.cyldnad(d1, d2)
        # Display
        print("Radius=",d1.x_ad_)
        print("Height=",d2.x_ad_)
        print("--> Volum=",d3.x_ad_)

if __name__ == '__main__':
    unittest.main()
    

# -*- coding: utf-8 -*-
"""
Created on Thu Nov 12 16:30:21 2015

@author: walcob
"""

import unittest
import dagview
import StringIO
import os

class testGFIntermediate(unittest.TestCase):
    def testBlankInit(self):
        self.test_intermediate = dagview.GFIntermediate()
        assert self.test_intermediate.number == 0, 'Number not initialized'
        assert self.test_intermediate.center == (0,0),'Center not initialized'
        assert self.test_intermediate.radius == 0.,'Radius not initialized'
        assert self.test_intermediate.dagfile == '','dagfile not initialized'
        assert self.test_intermediate.iflag == '','iflag not initialized'
        assert self.test_intermediate.state == 0,'State not initialized'
        assert self.test_intermediate.sym == 0.,'sym not initialized'
        assert self.test_intermediate.sas == 0.,'SAS not initialized'
        assert self.test_intermediate.entropy == 0.,'Entropy not initialized'
        assert self.test_intermediate.voids == 0,'Voids not initialized'
        assert self.test_intermediate.hbonds == 0,'Hbonds not initialized'
        assert self.test_intermediate.concentration == 0.,'Concentration not initialized'
        assert self.test_intermediate.barrels == [],'Barrels not initialized'
        
def main():
    GFintermediateSuite = unittest.makeSuite(testGFIntermediate,'test')
    unittest.main()
    
if __name__ == "__main__":
    main()
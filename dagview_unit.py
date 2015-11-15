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
    def setUp(self):
        self.test_intermediate = dagview.GFIntermediate()        
    
    def testBlankInit(self):
        #print '\ntestBlankInit\n'
        assert self.test_intermediate.number == 0, 'Number not initialized'
        assert self.test_intermediate.center == (0,0),'Center not initialized'
        assert self.test_intermediate.radius == 0.,'Radius not initialized'
        assert self.test_intermediate.dagfile == '','dagfile not initialized'
        assert self.test_intermediate.iflag == '','iflag not initialized'
        #assert self.test_intermediate.state == 0,'State not initialized'
        assert self.test_intermediate.sym == 0.,'sym not initialized'
        assert self.test_intermediate.sas == 0.,'SAS not initialized'
        assert self.test_intermediate.entropy == 0.,'Entropy not initialized'
        assert self.test_intermediate.voids == 0,'Voids not initialized'
        assert self.test_intermediate.hbonds == 0,'Hbonds not initialized'
        assert self.test_intermediate.concentration == 0.,'Concentration not initialized'
        assert self.test_intermediate.barrels == [],'Barrels not initialized'
        
    def testInit(self):
        #print '\ntestInit\n'
        try:
            self.test_intermediate = dagview.GFIntermediate(31,(8,90),1.,'emptyfile')
        except AssertionError:            
            pass
        else:
            assert True == False, 'Expected AssertionError'
            
    def testInitBlankDag(self):
        #print '\ntestInitBlankDag\n'
        self.test_intermediate = dagview.GFIntermediate(31,(8,90),1.,'')        
        assert self.test_intermediate.number == 31, 'Number not initialized'
        assert self.test_intermediate.center == (8,90),'Center not initialized'
        assert self.test_intermediate.radius == 1.,'Radius not initialized'
        assert self.test_intermediate.dagfile == '','dagfile not initialized'
        assert self.test_intermediate.iflag == '','iflag not initialized'
        #assert self.test_intermediate.state == 0,'State not initialized'
        assert self.test_intermediate.sym == 0,'sym not initialized'
        assert self.test_intermediate.sas == 0.,'SAS not initialized'
        assert self.test_intermediate.entropy == 0.,'Entropy not initialized'
        assert self.test_intermediate.voids == 0,'Voids not initialized'
        assert self.test_intermediate.hbonds == 0,'Hbonds not initialized'
        assert self.test_intermediate.concentration == 0.,'Concentration not initialized'
        assert self.test_intermediate.barrels == [],'Barrels not initialized'
        
    def testReadDagfileFail(self):
        #print '\ntestReadDagFileFail\n'
        self.success = self.test_intermediate.read_dagfile('')
        assert self.success == False, 'This should fail'
        
    def testFullInitPass(self):
        #print '\ntestFullInitPass\n'
        self.test_intermediate = dagview.GFIntermediate(1,(88,9),8.,'1LMB1_1.dag.out')        
        assert self.test_intermediate.number == 1, 'Number not initialized'
        assert self.test_intermediate.center == (88,9),'Center not initialized'
        assert self.test_intermediate.radius == 8.,'Radius not initialized'
        assert self.test_intermediate.dagfile == '1LMB1_1.dag.out','dagfile not initialized'
        assert self.test_intermediate.iflag == '333333333333333333333333333333333333333333333333333333333333333333333333333333333333333','iflag not initialized'
        #assert self.test_intermediate.state == 0,'State not initialized'
        assert self.test_intermediate.sym == 0,'sym not initialized'
        assert self.test_intermediate.sas == 13734.690,'SAS not initialized'
        assert self.test_intermediate.entropy == 225.270,'Entropy not initialized'
        assert self.test_intermediate.voids == 0,'Voids not initialized'
        assert self.test_intermediate.hbonds == 68,'Hbonds not initialized'
        assert self.test_intermediate.concentration == 0.00000001,'Concentration not initialized'
        assert self.test_intermediate.barrels == [0,0,0,0,0,0,0,0],'Barrels not initialized'
        
        
def main():
    GFintermediateSuite = unittest.makeSuite(testGFIntermediate,'test')
    unittest.main()
    
if __name__ == "__main__":
    main()
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
        
    def testContainsPoint(self):
        self.test_intermediate = dagview.GFIntermediate(1,(88,9),8.,'1LMB1_1.dag.out')
        assert self.test_intermediate.contains_point((89,8)) == True, 'Should be True'
        assert self.test_intermediate.contains_point((0,0)) == False, 'Should be False'
        
class testGFTransition(unittest.TestCase):
    '''Test GFTransition class'''
    def testBlankInit(self):
        self.transition = dagview.GFTransition()
        assert self.transition.number == 0, 'Should be 0'
        assert self.transition.coords == ((0,0),(0,0),(0,0),(0,0)), 'Should initialize to (0,0),(0,0),(0,0),(0,0)'
        assert self.transition.dagfile == '', 'Should be empty string'
        assert self.transition.f == 0, 'Should be 0'
        assert self.transition.u1 == 0, 'Should be 0'
        assert self.transition.u2 == 0, 'Should be 0'
        assert self.transition.entropy == 0., 'Should be 0.'
        assert self.transition.cuttype == '', 'Should be blank'
        assert self.transition.iseam == 0, 'Should be 0'
        assert self.transition.traffic == 0., 'Should be 0.'
        
    def testInitNoFile(self):
        try:
            self.transition = dagview.GFTransition(3,((1,2),(3,4),(5,6),(7,8)),'dagfile')
        except AssertionError:
            pass
        else:
            assert True == False, 'Should have AssertionError from inability to read dagfile'
            
    def testInitPrefill(self):
        self.transition = dagview.GFTransition(3,((1,2),(3,4),(5,6),(7,8)),'')
        assert self.transition.number == 3, 'Should be 3'
        assert self.transition.coords == ((1,2),(3,4),(5,6),(7,8)), 'Should be ((1,2),(3,4),(5,6),(7,8))'
        assert self.transition.dagfile == '', 'Should be empty'
        
    def testFullInitPass(self):
        self.transition = dagview.GFTransition(115,((111,29),(106,32),(111,35),(116,32)),'1LMB1_1.dag.out')
        assert self.transition.number == 115, 'Should be 115'
        assert self.transition.coords == ((111,29),(106,32),(111,35),(116,32)), 'Should be ((111,29),(106,32),(111,35),(116,32))'
        assert self.transition.dagfile == '1LMB1_1.dag.out'
        assert self.transition.f == 1, 'Should be 1'
        assert self.transition.u1 == 97, 'Should be 97'
        assert self.transition.u2 == 9, 'Should be 9'
        assert self.transition.entropy == 0.03, 'Should be 0.03'
        assert self.transition.cuttype == 'p', 'Should be p'
        assert self.transition.traffic == 0.25000323, 'Should be 0.25000323'
        assert self.transition.iseam == 0, 'Should be 0'

    def testContainsPoint(self):
        self.transition = dagview.GFTransition(115,((111,29),(106,32),(111,35),(116,32)),'1LMB1_1.dag.out')
        isFalse = self.transition.contains_point((0,0))
        isTrue = self.transition.contains_point((106,32))
        assert isFalse == False, 'Should be false'
        assert isTrue == True, 'Should be true'
        
class testParse(unittest.TestCase):
    def setUp(self):
        self.goodtxt = '''<html><body><img src="1LMB1_1.dag.png" border=0 usemap="#dagmap">
<map name="dagmap">
<area shape="circle" id="node1" href="../isegment.cgi?iseg=n1&amp;dag=1LMB1_1.dag.out&amp;" title="1" alt="" coords="88,9,8"/>
<area shape="poly" id="node2" href="../isegment.cgi?iseg=tp113&amp;dag=1LMB1_1.dag.out&amp;" alt="" coords="88,24,80,32,88,40,96,32"/>
<area shape="poly" id="node159" href="../isegment.cgi?iseg=tp115&amp;dag=1LMB1_1.dag.out&amp;" alt="" coords="111,29,106,32,111,35,116,32"/>
<area shape="poly" id="node162" href="../isegment.cgi?iseg=tp122&amp;dag=1LMB1_1.dag.out&amp;" alt="" coords="116,119,111,123,116,126,121,123"/>
<area shape="poly" id="node172" href="../isegment.cgi?iseg=tp127&amp;dag=1LMB1_1.dag.out&amp;" alt="" coords="181,74,176,77,181,80,186,77"/>
<area shape="circle" id="node3" href="../isegment.cgi?iseg=n2&amp;dag=1LMB1_1.dag.out&amp;" alt="" coords="9,281,9"/>'''
        self.goodOut = open('testgood.map','w+')
        self.goodOut.write(self.goodtxt)
        self.goodOut.close()
    
    def tearDown(self):
        os.remove('testgood.map')
        
    def testParseFiles(self):
        intermediates,transitions = dagview.parseImgMap('testgood.map')
        assert len(intermediates) == 2, 'Should be 2 -->%s'%(len(intermediates))
        assert len(transitions) == 4, 'Should be 4'
        intermediate = intermediates[0]
        assert intermediate.number == 1, 'Should be 1'
        assert intermediate.center == (88,9), 'Should be (88,9)'
        assert intermediate.radius == 8, 'Should be 8'
        assert intermediate.dagfile == '1LMB1_1.dag.out', 'Should be 1LMB1_1.dag.out'
        transition = transitions[0]
        assert transition.number == 113, 'Should be 113, but is %d'%(transition.number)
        assert transition.coords == ((88,24),(80,32),(88,40),(96,32)), 'Should be ((88,24),(80,32),(88,40),(96,32))'
        assert transition.dagfile == '1LMB1_1.dag.out', 'Should be 1LMB1_1.dag.out'
    
def main():
    GFIntermediateSuite = unittest.makeSuite(testGFIntermediate,'test')
    GFTransitionSuite = unittest.makeSuite(testGFTransition,'test')
    testParseSuite = unittest.makeSuite(testParse,'test')
    unittest.main()
    
if __name__ == "__main__":
    main()
import os
import sys
sys.path.insert(1,"../")
import unittest

import numpy as np
import pandas as pd
from pandas.util.testing import assert_frame_equal


from sgrsea import *


class Test_around(unittest.TestCase):
  def setUp(self):
    self.n = 0
    self.data = pd.Series([0.0,1.2,1.8,"4.5"])
	
  def test_around(self):
    expect = pd.Series([0,1,2,5])
    result = self.data.apply(lambda x: Normalization.around(x))
    self.assertTrue(expect.equals(result))

class Test_norm(unittest.TestCase):
	def setUp(self):
#		self.testfile = pd.read_table("testNormalization.txt",sep="\t")
#		self.method = ["upperquartile"]#,"total"]
#	
#	def test_norm_big(self):
#		for method in self.method:
#			output = Normalization.norm(self.testfile, method, "single")
#
#			if method == "upperquartile":
#				resultfile = pd.read_table("testNormalization_output.txt",sep="\t")
#				#self.assertEqual
#			#try:
#			assert_frame_equal(output.sort(axis=1),resultfile.sort(axis=1))
#			#	return True
#			#except:
#			#	return False

if __name__ == '__main__':
	unittest.main()

import sys
import os

sys.path.append(f'{os.path.dirname(os.path.realpath(__file__))}../../dependencies/OpenCOPTER')
sys.path.append(f'{os.path.dirname(os.path.realpath(__file__))}../../dependencies/wopwopd')

from libopencopter import *
from libwopwopd import *
import numpy as np


def compute_acoustics(case_list, rank):
	exec_wopwop3(case_list, './', 1, f'cases_{rank}')

def l2_norm(_x1, _x2):

	x1 = np.asarray(_x1)
	x2 = np.asarray(_x2)

	return np.sqrt(np.sum((x1 - x2)**2.0)/float(len(_x1)))

def rms(_x):
	x = np.asarray(_x)
	return np.sqrt(np.sum(x**2.0)/float(len(_x)))

from rpy2.rinterface import RRuntimeError
from rpy2.robjects.packages import importr
utils = importr('utils')

def importr_tryhard(packname, contriburl):
    try:
        rpack = importr(packname)
    except RRuntimeError:
		try:
			utils.install_packages(packname, contriburl = contriburl)
			rpack = importr(packname)
		except RRuntimeError:
			print 'no pack'
			rpack = 'none'
    return rpack

packname = 'rgl'
contriburl = 'http://cran.stat.ucla.edu/'
importr_tryhard(packname, contriburl)
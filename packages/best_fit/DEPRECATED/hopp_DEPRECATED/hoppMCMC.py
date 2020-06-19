"""
Source: https://github.com/kerguler/hoppMCMC (01/12/18)

adaptive basin-hopping Markov-chain Monte Carlo for Bayesian optimisation

This is the python (v2.7) implementation of the hoppMCMC algorithm aiming to
identify and sample from the high-probability regions of a posterior
distribution. The algorithm combines three strategies: (i) parallel MCMC,
(ii) adaptive Gibbs sampling and (iii) simulated annealing.

Overall, hoppMCMC resembles the basin-hopping algorithm implemented in the
optimize module of scipy, but it is developed for a wide range of modelling
approaches including stochastic models with or without time-delay.

"""

import os
# import sys
import numpy as np
from struct import pack, unpack
from scipy.stats import ttest_1samp as ttest

MPI_MASTER = 0
try:
    from mpi4py import MPI
    MPI_SIZE = MPI.COMM_WORLD.Get_size()
    MPI_RANK = MPI.COMM_WORLD.Get_rank()

    def Abort(str):
        print("ERROR: " + str)
        MPI.COMM_WORLD.Abort(1)
except:
    MPI_SIZE = 1
    MPI_RANK = 0

    def Abort(str):
        raise errorMCMC(str)

EPS_PULSE_VAR_MIN = 1e-12
EPS_VARMAT_MIN = 1e-7
EPS_VARMAT_MAX = 1e7


class errorMCMC(Exception):
    def __init__(self, value):
        self.value = value

    def __str__(self):
        return repr(self.value)


class binfile():

    def __init__(self, fname, mode, rowsize=1):
        self.headsize = 4
        self.bitsize = 8
        self.fname = fname
        self.mode = mode
        self.rowsize = rowsize
        if self.mode == 'r':
            try:
                self.f = open(self.fname, "rb")
            except IOError:
                Abort("File not found: " + self.fname)
            self.rowsize = unpack('<i', self.f.read(self.headsize))[0]
        elif self.mode == 'w':
            try:
                self.f = open(self.fname, "w+b")
            except IOError:
                Abort("File not found: " + self.fname)
            tmp = self.f.read(self.headsize)
            if not tmp:
                self.f.write(pack('<i', self.rowsize))
            else:
                self.rowsize = unpack('<i', tmp)[0]
        else:
            Abort("Wrong i/o mode: " + self.mode)
        self.fmt = '<' + 'd' * self.rowsize
        self.size = self.bitsize * self.rowsize
        print("Success: %s opened for %s size %d double rows" % (
            self.fname, "reading" if self.mode == 'r' else
            "writing", self.rowsize))

    def writeRow(self, row):
        self.f.write(pack(self.fmt, *row))
        self.f.flush()

    def readRows(self):
        self.f.seek(0, os.SEEK_END)
        filesize = self.f.tell()
        self.f.seek(self.headsize)
        tmp = self.f.read(filesize - self.headsize)
        ret = np.array(
            unpack('<' + 'd' * int(np.floor((
                filesize - self.headsize) / self.bitsize)), tmp),
            dtype=np.float64).reshape((int(np.floor(((
                filesize - self.headsize) / self.bitsize) / self.rowsize)),
                self.rowsize))
        return ret

    def close(self):
        self.f.close()
        print("Success: %s closed" % (self.fname if self.fname else "file"))


def readFile(filename):
    """
    Reads a binary output file and returns all rows/columns

    Parameters
    ----------

    filename:
         name of the output file

    Returns
    -------

    a np array with all the rows/columns

    """
    a = binfile(filename, "r")
    b = a.readRows()
    a.close()
    return b


def parMin(parmat):
    prm = min(parmat[:, 0])
    prmi = np.where(parmat[:, 0] == prm)[0][0]
    return {'i': prmi, 'f': prm}


def diagdot(mat, vec):
    for n in np.arange(mat.shape[0]):
        mat[n, n] *= vec[n]


def covariance(mat):
    return np.cov(mat, rowvar=False)


def coefVar(mat):
    mean0 = 1.0 / np.mean(mat, 0)
    return mean0 * (np.cov(mat, rowvar=False).T * mean0).T


def sensVar(mat):
    mean0 = np.mean(mat, 0)
    return mean0 * (np.linalg.inv(np.cov(mat, rowvar=False)).T * mean0).T


cov = covariance
rnorm = np.random.multivariate_normal
determinant = np.linalg.det


def join(a, s):
    return s.join(["%.16g" % (x) for x in a])


def logsumexp(x):
    a = np.max(x)
    return a + np.log(np.sum(np.exp(x - a)))


def compareAUCs(parmats, groups, tol=1.0):
    from scipy.stats import gaussian_kde
    ids = np.unique(groups)
    parmats = parmats[:, np.var(parmats, axis=0) != 0]
    parmat_list = [parmats[groups == n, :] for n in ids]
    counts = [np.sum(n == groups) for n in ids]
    try:
        kde = gaussian_kde(parmats[:, 1:].T)
    except:
        print("Warning: Problem encountered in compareAUC!")
        return {ids[g]: 0 for g in range(len(ids))}

    wt_list = np.array([
        [kde.evaluate(pr) for pr in parmat[:, 1:]] for parmat in parmat_list])
    mn = np.min(wt_list)
    wt_list = np.log(wt_list / mn)
    # Importance sampling for Monte Carlo integration:
    favg_exp = np.array([
        [-parmat_list[n][m, 0] / tol - wt_list[n][m] for m in
         range(parmat_list[n].shape[0])] for n in range(len(parmat_list))])
    favg_logsums = np.array([logsumexp(x) for x in favg_exp]) - np.log(counts)
    favg_list = np.exp(favg_logsums - logsumexp(favg_logsums))
    return {ids[g]: favg_list[g] for g in range(len(ids))}


def compareAUC(parmat0, parmat1, T):
    from scipy.stats import gaussian_kde
    parmats = np.vstack((parmat0, parmat1))
    parmat0cp = parmat0[:, np.var(parmats, axis=0) != 0].copy()
    parmat1cp = parmat1[:, np.var(parmats, axis=0) != 0].copy()
    parmats = parmats[:, np.var(parmats, axis=0) != 0]
    try:
        kde = gaussian_kde(parmats[:, 1:].T)
    except:
        print("Warning: Problem encountered in compareAUC!")
        return {'acc': 0, 'favg0': 0, 'favg1': 0}
    wt0 = np.array([kde.evaluate(pr) for pr in parmat0cp[:, 1:]])
    wt1 = np.array([kde.evaluate(pr) for pr in parmat1cp[:, 1:]])
    mn = np.min([wt0, wt1])
    wt0 /= mn
    wt1 /= mn
    # Importance sampling for Monte Carlo integration:
    favg0 = np.mean([
        np.exp(-parmat0cp[m, 0] / T) / wt0[m] for m in
        range(parmat0cp.shape[0])])
    favg1 = np.mean([
        np.exp(-parmat1cp[m, 0] / T) / wt1[m] for m in
        range(parmat1cp.shape[0])])
    acc = (not np.isnan(favg0) and not np.isnan(favg1) and
           favg0 > 0 and favg1 >= 0 and
           (favg1 >= favg0 or
           (np.random.uniform() < (favg1 / favg0))))
    return {'acc': acc, 'favg0': favg0, 'favg1': favg1}


def anneal_exp(y0, y1, steps):
    return y0 * np.exp(-(np.arange(
        steps, dtype=np.float64) / (steps - 1.)) * np.log(np.float64(y0) / y1))


def anneal_linear(y0, y1, steps):
    return np.append(np.arange(
        y0, y1, (y1 - y0) / (steps - 1), dtype=np.float64), y1)


def anneal_sigma(y0, y1, steps):
    return y0 + (y1 - y0) * (1. - 1. / (1. + np.exp(-(
        np.arange(steps) - (0.5 * steps)))))


def anneal_sigmasoft(y0, y1, steps):
    return y0 + (y1 - y0) * (1. - 1. / (1. + np.exp(-12.5 * (
        np.arange(steps) - (0.5 * steps)) / steps)))


def ldet(mat):
    try:
        det = determinant(mat)
    except:
        return -np.Inf
    if det == 0:
        return -np.Inf
    else:
        return np.log(det)


def finalTest(fitFun, param, testnum=10):
    for n in range(testnum):
        f = fitFun(param)
        if not np.isnan(f) and not np.isinf(f):
            return f
    Abort("Incompatible parameter set (%g): %s" % (f, join(param, ",")))


class hoppMCMC:
    def __init__(
        self, fitFun, param, varmat, inferpar=None, gibbs=True, num_hopp=3,
        num_adapt=25, num_chain=12, chain_length=50, rangeT=None,
            model_comp=1000.0, outfilename=''):
        """
        Adaptive Basin-Hopping MCMC Algorithm

        Parameters
        ----------

        fitFun:
              fitFun(x) - objective function which takes a np array as the
              only argument

        param:
              initial parameter vector

        varmat:
              2-dimensional array of initial covariance matrix

        inferpar:
              an array of indexes of parameter dimensions to be inferred
              (all parameters are inferred by default)

        gibbs:
              indicates the type of chain iteration
              True - Gibbs iteration where each parameter dimension has its
              own univariate Gaussian proposal distribution (default)
              False - Metropolis-Hastings iteration where there is a single
              multivariate Gaussian proposal distribution

        num_hopp:
              number of hopp-steps (default=3)

        num_adapt:
              number of adaptation steps (default=25)

        num_chain:
              number of MCMC chains (default=12)

        chain_length:
              size of each chain (default=50)

        rangeT:
              [min,max] - range of annealing temperatures for each hopp-step
              (default=[1,1000])
              min should be as low as possible but not lower
              max should be sufficiently permissive to be able to jump between
              posterior modes

        model_comp:
              tolerance for accepting subsequent hopp-steps (default=1000)
              this should ideally be equal to or higher than max(rangeT)

        outfilename:
              name of the output file (default='')
              use this option for a detailed account of the results

              outfilename.final lists information on hopp-steps:
                   hopp-step
                   acceptance (0/1)
                   weighted average of exp{-f0/model_comp}
                   weighted average of exp{-f1/model_comp}
              outfilename.parmat lists chain status at the end of each
              adaptation step:
                   adaptation step
                   chain id
                   annealing temperature
                   score (f)
                   parameter values
              both files can be read using the readFile function

        Returns
        -------

        A hoppMCMC object with

              parmat: an array of num_chain x (1+len(param))
                      current score (f) and parameter values for each chain
              varmat: an array of len(inferpar) x len(inferpar)
                      latest proposal distribution
              parmats: a list of parameter values (i.e. parmat) accepted at the
              end of hopp-steps

        See also
        --------

        the documentation (doc/hoppMCMC_manual.pdf) for more information and
        examples
        """
        self.multi = {'cov': covariance,
                      'rnorm': np.random.multivariate_normal,
                      'det': np.linalg.det}
        self.single = {'cov': covariance,
                       'rnorm': np.random.normal,
                       'det': abs}
        self.stat = self.multi
        self.gibbs = gibbs
        self.num_hopp = num_hopp
        self.num_adapt = num_adapt
        self.num_chain = num_chain
        self.chain_length = chain_length
        self.rangeT = np.sort([1.0, 1000.0] if rangeT is None else rangeT)
        self.model_comp = model_comp
        # ---
        self.fitFun = fitFun
        self.param = np.array(param, dtype=np.float64, ndmin=1)
        f0 = finalTest(self.fitFun, self.param)
        self.parmat = np.array([
            [f0] + self.param.tolist() for n in range(self.num_chain)],
            dtype=np.float64)
        self.varmat = np.array(varmat, dtype=np.float64, ndmin=2)
        self.parmats = []
        # ---
        if inferpar is None:
            self.inferpar = np.arange(len(self.param), dtype=np.int32)
        else:
            self.inferpar = np.array(inferpar, dtype=np.int32)
        print("Parameters to infer: %s" % (join(self.inferpar, ",")))
        # ---
        self.rank_indices = [
            np.arange(i, self.num_chain, MPI_SIZE) for i in range(MPI_SIZE)]
        self.worker_indices = np.delete(range(MPI_SIZE), MPI_MASTER)
        # ---
        self.outfilename = outfilename
        self.outparmat = None
        self.outfinal = None
        if MPI_RANK == MPI_MASTER:
            if self.outfilename:
                self.outparmat = binfile(
                    self.outfilename + '.parmat', 'w',
                    self.parmat.shape[1] + 3)
                self.outfinal = binfile(self.outfilename + '.final', 'w', 4)
        # ---
        for hopp_step in range(self.num_hopp):
            self.anneal = anneal_sigmasoft(
                self.rangeT[0], self.rangeT[1], self.num_adapt)
            for adapt_step in range(self.num_adapt):
                self.runAdaptStep(hopp_step * self.num_adapt + adapt_step)
            if MPI_RANK == MPI_MASTER:
                test = {'acc': True, 'favg0': np.nan, 'favg1': np.nan} if\
                    len(self.parmats) == 0 else compareAUC(
                        self.parmats[-1][
                            :, [0] + (1 + self.inferpar).tolist()],
                        self.parmat[
                            :, [0] + (1 + self.inferpar).tolist()],
                        self.model_comp)
                if test['acc']:
                    self.parmats.append(self.parmat)
                else:
                    self.parmat = self.parmats[-1].copy()
                    self.param = self.parmat[
                        parMin(self.parmat)['i'], 1:].copy()
                # ---
                if self.outfinal:
                    self.outfinal.writeRow([
                        hopp_step, test['acc'], test['favg0'], test['favg1']])
                else:
                    print("parMatAcc.final: %d,%s" % (
                        hopp_step,
                        join([test['acc'], test['favg0'], test['favg1']],
                             ",")))
                # ---
            if MPI_SIZE > 1:
                self.parmat = MPI.COMM_WORLD.bcast(
                    self.parmat, root=MPI_MASTER)
                self.param = MPI.COMM_WORLD.bcast(self.param, root=MPI_MASTER)
        # ---
        if MPI_SIZE > 1:
            self.parmats = MPI.COMM_WORLD.bcast(self.parmats, root=MPI_MASTER)
            if self.outparmat:
                self.outparmat.close()
            if self.outfinal:
                self.outfinal.close()

    def runAdaptStep(self, adapt_step):
        if MPI_RANK == MPI_MASTER:
            pm = parMin(self.parmat)
            self.param = self.parmat[pm['i'], 1:].copy()
            self.parmat = np.array([
                self.parmat[pm['i'], :].tolist() for n in
                range(self.num_chain)], dtype=np.float64)
        if MPI_SIZE > 1:
            self.param = MPI.COMM_WORLD.bcast(self.param, root=MPI_MASTER)
            self.parmat = MPI.COMM_WORLD.bcast(self.parmat, root=MPI_MASTER)
        # ---
        for chain_id in self.rank_indices[MPI_RANK]:
            # ---
            mcmc = chainMCMC(self.fitFun,
                             self.param,
                             self.varmat,
                             gibbs=self.gibbs,
                             chain_id=chain_id,
                             pulsevar=1.0,
                             anneal=self.anneal[0],
                             accthr=0.5,
                             inferpar=self.inferpar,
                             varmat_change=0,
                             pulse_change=10,
                             pulse_change_ratio=2,
                             print_iter=0)
            for m in range(self.chain_length):
                mcmc.iterate()
            self.parmat[chain_id, :] = mcmc.getParam()
            # ---
        if MPI_RANK == MPI_MASTER:
            for worker in self.worker_indices:
                parmat = MPI.COMM_WORLD.recv(source=worker, tag=1)
                for chain_id in self.rank_indices[worker]:
                    self.parmat[chain_id, :] = parmat[chain_id, :]
            # ---
            self.varmat = np.array(
                self.stat['cov'](self.parmat[:, 1 + self.inferpar]), ndmin=2)
            self.varmat[np.abs(self.varmat) < EPS_VARMAT_MIN] = 1.0
            # ---
            for chain_id in range(self.num_chain):
                if self.outparmat:
                    tmp = [adapt_step, chain_id, self.anneal[0]] +\
                        self.parmat[chain_id, :].tolist()
                    self.outparmat.writeRow(tmp)
                else:
                    print("param.mat.step: %d,%d,%g,%s" % (
                        adapt_step, chain_id, self.anneal[0],
                        join(self.parmat[chain_id, :], ",")))
            # ---
            if len(self.anneal) > 1:
                self.anneal = self.anneal[1:]
        else:
            MPI.COMM_WORLD.send(self.parmat, dest=MPI_MASTER, tag=1)

        if MPI_SIZE > 1:
            self.parmat = MPI.COMM_WORLD.bcast(self.parmat, root=MPI_MASTER)
            self.varmat = MPI.COMM_WORLD.bcast(self.varmat, root=MPI_MASTER)
            self.anneal = MPI.COMM_WORLD.bcast(self.anneal, root=MPI_MASTER)


class chainMCMC:
    def __init__(self,
                 fitFun,
                 param,
                 varmat,
                 inferpar=None,
                 gibbs=True,
                 chain_id=0,
                 pulsevar=1.0,
                 anneal=1,
                 accthr=0.5,
                 varmat_change=0,
                 pulse_change=10,
                 pulse_change_ratio=2,
                 pulse_allow_decrease=True,
                 pulse_allow_increase=True,
                 pulse_min=1e-7,
                 pulse_max=1e7,
                 print_iter=0):
        """
                 
        MCMC Chain with Adaptive Proposal Distribution

        Usage
        -----

        Once created, a chainMCMC is iterated using the iterate method.
        Depending on the value of gibbs, this method calls either iterateMulti or iterateSingle.
        
        Parameters
        ----------
        
        fitFun:
              fitFun(x) - objective function which takes a np array as the only argument

        param:
              initial parameter vector

        varmat:
              2-dimensional array of initial covariance matrix

        inferpar:
              an array of indexes of parameter dimensions to be inferred
              (all parameters are inferred by default)
              
        gibbs:
              indicates the type of chain iteration
              True - Gibbs iteration where each parameter dimension has its own univariate
                     Gaussian proposal distribution (default)
              False - Metropolis-Hastings iteration where there is a single multivariate
                     Gaussian proposal distribution

        chain_id:
              a chain identifier (default=0)

        pulsevar:
              scaling factor for the variance of the proposal distribution (default=1)

        anneal:
              annealing temperature (default=1)

        accthr:
              desired acceptance rate (default=0.5)

        varmat_change:
              how often variance should be updated? (default=0)
              varmat_change=0 - fixed variance
              varmat_change=n - variance is updated at each nth step

        pulse_change:
              how often pulsevar should be updated? (default=10)
              pulse_change=0 - fixed pulsevar
              pulse_change=n - pulsevar is updated at each nth step

        pulse_change_ratio:
              how should pulsevar be updated? (default=2)
              (pulsevar *= pulse_change_ratio)

        pulse_allow_increase:
             allow pulse to increase (default=True)

        pulse_allow_decrease:
             allow pulse to decrease (default=True)

        pulse_min:
             minimum value of pulse (default=1e-7)
             
        pulse_max:
             maximum value of pulse (default=1e7)
             
        print_iter:
              how often chain status should be printed? (default=0)
              print_iter=0 - do not print status
              print_iter=n - print status at each nth step
              (default=0)
              
        Returns
        -------

        A chainMCMC object with 

              getParam: a method for obtaining the latest iteration (f + parameter values)

              getVarmat: a method for obtaining the latest proposal distribution (varmat * pulsevar)

        See also
        --------

        the documentation (doc/hoppMCMC_manual.pdf) for more information and examples
                            
        """
        self.multi = {'cov': covariance,
                      'rnorm': np.random.multivariate_normal,
                      'det': np.linalg.det}
        self.single = {'cov': covariance,
                       'rnorm': np.random.normal,
                       'det': abs}
        # ---
        self.chain_id = chain_id;
        self.fitFun = fitFun
        self.parmat = np.array(param,dtype=np.float64,ndmin=1)
        self.varmat = np.array(varmat,dtype=np.float64,ndmin=2)
        if inferpar is None:
            self.inferpar = np.arange(len(param),dtype=np.int32)
        else:
            self.inferpar = np.array(inferpar,dtype=np.int32)
        self.anneal = np.array(anneal,dtype=np.float64,ndmin=1)
        self.accthr = np.array(accthr,dtype=np.float64,ndmin=1)
        # ---
        self.pulse_change_ratio = pulse_change_ratio
        self.pulse_nochange = np.float64(1)
        self.pulse_increase = np.float64(self.pulse_change_ratio)
        self.pulse_decrease = np.float64(1.0/self.pulse_change_ratio)
        self.pulse_change = pulse_change
        self.pulse_collect = max(1,self.pulse_change)
        self.allow_pincr = pulse_allow_increase
        self.allow_pdecr = pulse_allow_decrease
        self.pulse_min = pulse_min
        self.pulse_max = pulse_max
        self.varmat_change = varmat_change
        self.varmat_collect = max(1,self.varmat_change)
        # ---
        if self.parmat.ndim==1:
            f0 = finalTest(self.fitFun,self.parmat)
            self.parmat = np.array([[f0]+self.parmat.tolist() for i in range(self.varmat_collect)])
        elif self.parmat.shape[0]!=self.varmat_collect:
            Abort("Dimension mismatch in chainMCMC! parmat.shape[0]=%d collect=%d" %(self.parmat.shape[0],self.varmat_collect))
        if self.varmat.shape and self.inferpar.shape[0] != self.varmat.shape[0]:
            Abort("Dimension mismatch in chainMCMC! inferpar.shape[0]=%d varmat.shape[0]=%d" %(self.inferpar.shape[0],self.varmat.shape[0]))
        # ---
        self.gibbs = gibbs
        if self.gibbs:
            self.pulsevar = np.array(np.repeat(pulsevar,len(self.inferpar)),dtype=np.float64)
            self.acc_vecs = [np.repeat(False,self.pulse_collect) for n in range(len(self.inferpar))]
            self.iterate = self.iterateSingle
        else:
            self.pulsevar = pulsevar
            self.acc_vec = np.repeat(False,self.pulse_collect)
            self.iterate = self.iterateMulti
        self.pulsevar0 = self.pulsevar
        if not self.gibbs and (self.parmat.shape[1]==2 or self.inferpar.shape[0]==1):
            Abort("Please set gibbs=True!")
        self.varmat = self.varmat*self.pulsevar
        # ---
        self.halfa = 0.025
        if self.pulse_change<25:
            self.halfa = 0.05
        self.print_iter = print_iter
        self.step = 0
        self.index = 0
        self.index_acc = 0

    def getParam(self):
        return self.parmat[self.index,:].copy()

    def getVarmat(self):
        return (self.varmat*self.pulsevar).copy()

    def getVarPar(self):
        return ldet(self.multi['cov'](self.parmat[:,1+self.inferpar]))

    def getVarVar(self):
        return ldet(self.varmat)

    def getAcc(self):
        if self.gibbs:
            return np.array([np.mean(acc_vec) for acc_vec in self.acc_vecs])
        else:
            return np.mean(self.acc_vec)

    def setParam(self,parmat):
        self.parmat[self.index,:] = np.array(parmat,dtype=np.float64,ndmin=1).copy()

    def newParamSingle(self,param,param_id):
        try:
            param1 = self.single['rnorm'](param,
                                          self.varmat[param_id,param_id]*self.pulsevar[param_id])
        except:
            print("Warning: Failed to generate a new parameter set")
            param1 = np.copy(param)
        return param1

    def newParamMulti(self):
        try:
            param1 = self.multi['rnorm'](self.parmat[self.index,1:][self.inferpar],self.varmat*self.pulsevar)
        except np.linalg.linalg.LinAlgError:
            print("Warning: Failed to generate a new parameter set")
            param1 = np.copy(self.parmat[self.index,1:][self.inferpar])
        return param1

    def checkMove(self,f0,f1):
        acc = (not np.isnan(f1) and not np.isinf(f1) and
            f1 >= 0 and
            (f1 <= f0 or
            (np.log(np.random.uniform()) < (f0-f1)/self.anneal[0])))
        # --- f0 = 0.5*SS_0
        # --- f = 0.5*SS
        # --- sqrt(anneal) == st.dev.
        # --- 0.5*x^2/(T*s^2)
        # --- return exp(-0.5*SS/anneal)/exp(-0.5*SS_0/anneal)
        return(acc)

    def pulsevarUpdate(self,acc_vec):
        # --- Test if mean(acc_vec) is equal to accthr
        try:
            r = ttest(acc_vec,self.accthr)
        except ZeroDivisionError:
            if all(acc_vec)<=0 and self.allow_pdecr:
                return self.pulse_decrease
            elif all(acc_vec)>=0 and self.allow_pincr:
                return self.pulse_increase
            else:
                return self.pulse_nochange
        # ---
        if r[1]>=self.halfa: return self.pulse_nochange
        if r[0]>0 and self.allow_pincr: return self.pulse_increase
        if r[0]<0 and self.allow_pdecr: return self.pulse_decrease
        # --- Return default
        return self.pulse_nochange

    def iterateMulti(self):
        self.step += 1
        # ---
        acc = False
        f0 = self.parmat[self.index,0]
        param1 = np.copy(self.parmat[self.index,1:])
        param1[self.inferpar] = self.newParamMulti()
        f1 = self.fitFun(param1)
        acc = self.checkMove(f0,f1)
        # ---
        self.index_acc = (self.index_acc+1)%self.pulse_collect
        if acc:
            self.index = (self.index+1)%self.varmat_collect
        # ---
        self.acc_vec[self.index_acc] = acc
        if acc:
            self.parmat[self.index,0] = f1
            self.parmat[self.index,1:] = param1
        # ---
        if self.print_iter and (self.step%self.print_iter)==0:
            print("param.mat.chain: %d,%d,%s" %(self.step,self.chain_id,join(self.parmat[self.index,:],",")))
        # ---
        if self.step>1:
            # --- 
            if self.pulse_change and (self.step%self.pulse_change)==0:
                self.pulsevar = max(1e-7,self.pulsevar*self.pulsevarUpdate(self.acc_vec))
            # --- 
            if self.varmat_change and (self.step%self.varmat_change)==0:
                self.varmat = np.array(self.multi['cov'](self.parmat[:,1+self.inferpar]),ndmin=2)
                a = np.diag(self.varmat)<EPS_VARMAT_MIN
                self.varmat[a,a] = EPS_VARMAT_MIN
        # --- 
        if self.print_iter and (self.step%self.print_iter)==0:
            print("parMatAcc.chain: %s" %(join([self.step,self.chain_id,ldet(self.multi['cov'](self.parmat[:,1+self.inferpar])),ldet(self.varmat),np.mean(self.acc_vec),self.pulsevar],",")))

    def iterateSingle(self):
        self.step += 1
        self.index_acc = (self.index_acc+1)%self.pulse_collect
        # ---
        acc_steps = False
        f0 = self.parmat[self.index,0]
        param0 = self.parmat[self.index,1:].copy()
        for param_id in np.arange(len(self.inferpar)):
            param1 = param0.copy()
            param1[self.inferpar[param_id]] = self.newParamSingle(param1[self.inferpar[param_id]],param_id)
            f1 = self.fitFun(param1)
            acc = self.checkMove(f0,f1)
            if acc:
                acc_steps = True
                f0 = f1
                param0[self.inferpar[param_id]] = np.copy(param1[self.inferpar[param_id]])
            self.acc_vecs[param_id][self.index_acc] = acc
            # ---
        if acc_steps:
            self.index = (self.index+1)%self.varmat_collect
            self.parmat[self.index,0] = f0
            self.parmat[self.index,1:] = param0
            if np.isnan(f0) or np.isinf(f0):
                Abort("Iterate single failed with %g: %s" %(f0,join(param0,",")))
        # ---
        if self.print_iter and (self.step%self.print_iter)==0:
            print("param.mat.chain: %d,%d,%s" %(self.step,self.chain_id,join(self.parmat[self.index,:],",")))
        # ---
        if self.step>1:
            # --- 
            if self.pulse_change and (self.step%self.pulse_change)==0:
                for param_id in np.arange(len(self.inferpar)):
                    tmp = min(self.pulse_max,max(self.pulse_min,self.pulsevar[param_id]*self.pulsevarUpdate(self.acc_vecs[param_id])))
                    if np.abs(self.varmat[param_id,param_id]*tmp) >= EPS_PULSE_VAR_MIN:
                        self.pulsevar[param_id] = tmp
            # --- 
            if self.varmat_change and (self.step%self.varmat_change)==0:
                for param_id in np.arange(len(self.inferpar)):
                    tmp = max(EPS_VARMAT_MIN,self.single['cov'](self.parmat[:,1+self.inferpar[param_id]]))
                    if np.abs(tmp*self.pulsevar[param_id]) >= EPS_PULSE_VAR_MIN:
                        self.varmat[param_id,param_id] = tmp
        # --- 
        if self.print_iter and (self.step%self.print_iter)==0:
            print("parMatAcc.chain: %s" %(join([self.step,self.chain_id,ldet(self.multi['cov'](self.parmat[:,1+self.inferpar])),ldet(self.varmat)],",")))
            print("parMatAcc.chain.accs: %d,%d,%s" %(self.step,self.chain_id,join([np.mean(acc_vec) for acc_vec in self.acc_vecs],",")))
            print("parMatAcc.chain.pulses: %d,%d,%s" %(self.step,self.chain_id,join(self.pulsevar,",")))

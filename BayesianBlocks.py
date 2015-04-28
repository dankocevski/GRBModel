"""
@brief Pure python implementation of the Bayesian Blocks algorithm
described by Jackson, Scargle et al. 2005, IEEE Signal Processing
Letters, 12, 105. (http://arxiv.org/abs/math/0309285)

@author J. Chiang <jchiang@slac.stanford.edu>
"""
import sys
import copy
import numpy as num

def gammln(xx):
    cof = [76.18009172947146, -86.50532032941677,
           24.01409824083091, -1.231739572450155,
           0.1208650973866179e-2, -0.5395239384953e-5]
    y = xx
    x = xx
    tmp = x + 5.5
    tmp -= (x + 0.5)*num.log(tmp)
    ser = 1.000000000190015
    for j in range(6):
        y += 1
        ser += cof[j]/y
    return -tmp + num.log(2.5066282746310005*ser/x)

class BayesianBlocks(object):
    """ 
    Unbinned mode:
    >>> bb = BayesianBlocks(arrival_times)

    Binned:
    >>> bb = BayesianBlocks(bin_content, bin_sizes, start_time)

    Obtaining the piecewise constant light curve:
    >>> time, rate = bb.globalOpt(ncp_prior=1)
    """
    def __init__(self, *argv):
        if len(argv) == 1:
            events = list(argv[0])
            events.sort()
            events = num.array(events)
            self.cellContent = num.ones(len(argv[0]))
            self.cellSizes = self._generateCells(events)
            self.binned = False
        else:
            self.cellContent = copy.deepcopy(argv[0])
            self.cellSizes = copy.deepcopy(argv[1])
            self.tstart = argv[2]
            self.binned = True
        self._rescaleCells()
        
    def globalOpt(self, ncp_prior=1):
        opt, last = [], []
        opt.append(self.blockCost(0, 0) - ncp_prior)
        last.append(0)
        npts = len(self.cellContent)
        for nn in range(1, npts):
            max_opt = self.blockCost(0, nn) - ncp_prior
            jmax = 0
            for j in range(1, nn+1):
                my_opt = opt[j-1] + self.blockCost(j, nn) - ncp_prior
                if my_opt > max_opt:
                    max_opt = my_opt
                    jmax = j
            opt.append(max_opt)
            last.append(jmax)
        changePoints = []
        indx = last[-1]
        while indx > 0:
            changePoints.insert(0, indx)
            indx = last[indx-1]
        changePoints.insert(0, 0)
        changePoints.append(npts)
        return self._lightCurve(changePoints)
        
    def _lightCurve(self, changePoints):
        xx = []
        yy = []
        cell_sizes = self.cellSizes/self.cellScale
        for imin, imax in zip(changePoints[:-1], changePoints[1:]):
            xx.extend([self.tstart + sum(cell_sizes[:imin]),
                       self.tstart + sum(cell_sizes[:imax])])
            yval = sum(self.cellContent[imin:imax])/sum(cell_sizes[imin:imax])
            yy.extend([yval, yval])
        return xx, yy
        
    def blockCost(self, imin, imax):
        size = self.blockSize(imin, imax)
        content = self.blockContent(imin, imax)
        arg = size - content
        if arg > 0:
            my_cost = gammln(content + 1) + gammln(arg + 1) - gammln(size + 2)
            return my_cost
        return -num.log(size)
        
    def blockSize(self, imin, imax):
        return sum(self.cellSizes[imin:imax+1])
        
    def blockContent(self, imin, imax):
        return sum(self.cellContent[imin:imax+1])
        
    def _generateCells(self, events):
        self.tstart = (3*events[0] - events[1])/2.
        bounds = ((events[1:] + events[:-1])/2.).tolist()
        bounds.insert(0, self.tstart)
        bounds.append((3*events[-1] - events[-2])/2.)
        bounds = num.array(bounds)
        return bounds[1:] - bounds[:-1]
        
    def _rescaleCells(self):
        self.tbounds = [self.tstart]
        for item in self.cellSizes:
            self.tbounds.append(self.tbounds[-1] + item)
        smallest_cell = min(self.cellSizes)
        self.cellScale = 2./smallest_cell
        if self.binned:
            self.cellScale *= sum(self.cellContent)
        self.cellSizes *= self.cellScale

if __name__ == '__main__':
#    import hippoplotter as plot
#    import distributions as dist
#    nsamp = 200
#    events = dist.sample(dist.stepFunction(0.5, 0.7, amp=0.7), nsamp)
#
#    output = open('events_bb.dat', 'w')
#    for event in events:
#        output.write("%12.4e\n" % event)
#    output.close()

    class Histogram(object):
        def __init__(self, xmin, xmax, nx):
            self.xmin = xmin
            self.dx = (xmax - xmin)/float(nx)
            self.binContent = num.zeros(nx)
            self.binSizes = self.dx*num.ones(nx)
        def add(self, xx, wt=1):
            indx = int((xx - self.xmin)/self.dx)
            self.binContent[indx] += wt

    events = [float(x.strip()) for x in open(sys.argv[1], 'r')]

#    hist = Histogram(0, 1001, 10)
#    for event in events:
#        hist.add(event)

#    bb = BayesianBlocks(events)
#    xx, yy = bb.globalOpt(ncp_prior=1)
#    print xx
#    print yy
    
    
    binsize = num.ones(len(events))
    	
    bb2 = BayesianBlocks(events, binsize, 0)
    time, rate = bb2.globalOpt(ncp_prior=0.1)
    
#    print events
#    print time
#    print rate

#    output = open('/Users/Kocevski/Temp/idl_tmp_bb.txt', 'w')
    output = open(sys.argv[1], 'w')

    for timei, ratei in zip(time, rate):
    	output.write("%12.2f %12.4f\n" % (timei, ratei) )
    output.close()




#    plot.histogram(events)
#    plot.scatter(xx, yy, oplot=1, pointRep='Line', color='red', autoscale=1)
#    plot.scatter(xx2, yy2, oplot=1, pointRep='Line', color='blue')

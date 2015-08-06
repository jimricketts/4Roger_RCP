#support for regression of a single gcm, with optional regridding 
import numpy as np
import scipy
import os
import sys
from scipy import stats
from scipy.stats import t
import math
import string
import glob
import subprocess as sp
#import Nio
#import Ngl
#import gcm
import globalw
import numpy.ma as ma
import copy
SVNRevision="$Revision: 325 $"

#import pdb 

class regressException(Exception):
  pass
#define a varable to control print diagnostics

printRegressDiags=False

statkeylist=["sse","ssb","ssa","sx" ,"sxx","sxy","sy" ,"syy" ,"ssx","ssy"]

def regress(data, xs, stats=[], asDict=False):
    '''
    This routine computes linear regressions on a gridded basis
    It returns slopes and then offsets
    The optional third parameter specifies some descriptive stats to return also returned in order specified.
    see http://en.wikipedia.org/wiki/Simple_linear_regression
    Values are any of 
        "sse" : sum squared error
        "ssb" : sum square on beta
        "ssa" : sum squared on alpha
        "sx"  : sum of x
        "sxx" : sum of squares of x
        "sxy" : sum of x * y
        "sy"  : sum of y
        "syy  : sum of squares of y
        "ssx" : sum squared of x
        "ssy" : sum squared of y
    '''
    
    if len(xs) != len(data):
        raise regressException("error: regress.regress. rank of xs and ys needs to agree")
    rank=len(np.shape(data))
    #print np.shape(data), len(xs)
    if rank == 1:
        n = np.shape(data)[0]
        sy = np.sum(data)
        syy = np.sum(data * data)
        sxy = np.sum(data * xs)
        sx = np.sum(xs)
        sxx = np.sum(xs * xs)
        ssx = sx * sx
    else:
        
        (n,x,y)=np.shape(data)
        xa=np.array(xs)
        sx = np.sum(xa)
        sxx= np.sum(xa * xa)
        #sxx=np.zeros((x,y))
        sxy=np.zeros((x,y))
        sy=np.zeros((x,y))
        syy=np.zeros((x,y))
        
        for i in range(n):
            dat=data[i][:]
            sy += dat
            syy += dat * dat
            sxy += dat * xa[i]
            #sx += xs[i]
            #sxx += xs[i] * xs[i]
        ssx = sx * sx
    
    beta = (n * sxy - sx * sy)/(n * sxx - ssx)
    alpha = (sy-beta*sx)/n
    
    if stats != []:
        if asDict:
            returns=dict()
        else:
            returns = []
        #then we set up to return extra statistics
        statset = set(stats)
        nullset = set([])
        if statset & set(["sse", "ssy", "ssb", "ssa"]) != nullset:
            ssy=sy*sy
        if statset & set(["sse", "ssb", "ssa"]) != nullset:
            if n > 2:
                sse= (1.0/(n * (n -2)))*(n * syy - ssy - beta * beta *(n * sxx - ssx))
            else:
                sse = np.NaN
        if statset & set(["ssb", "ssa"]) != nullset:
            ssb= (n * sse)/(n*sxx-ssx)
        if "ssa" in statset :
            ssa= ssb * sxx/n
        #print stats
        if not asDict:
            for k in stats:
                if k == "sse": returns.append((k,sse))
                if k == "ssb": returns.append((k,ssb))
                if k == "ssa": returns.append((k,ssa))
                if k == "sx":  returns.append((k,sx))
                if k == "sxx": returns.append((k,sxx))
                if k == "sxy": returns.append((k,sxy))
                if k == "sy" : returns.append((k,sy))
                if k == "syy": returns.append((k,syy))
                if k == "ssx": returns.append((k,ssx))
                if k == "ssy": returns.append((k,ssy))
                #print k, len(returns)
        else:
            for k in stats:
                if k == "sse": returns[k] =sse
                if k == "ssb": returns[k] =ssb
                if k == "ssa": returns[k] =ssa
                if k == "sx":  returns[k] =sx
                if k == "sxx": returns[k] =sxx
                if k == "sxy": returns[k] =sxy
                if k == "sy" : returns[k] =sy
                if k == "syy": returns[k] =syy
                if k == "ssx": returns[k] =ssx
                if k == "ssy": returns[k] =ssy
        
        return beta, alpha, returns
    else:
        return beta, alpha

def xs_matched_to_data( data, xs):
    xa = ma.array(np.zeros(np.shape(data)), mask=data.mask, keep_mask=True)
    for i in range(len(xs)):
        xa[i]=xs[i]
    xa = ma.array(xa, mask=data.mask, keep_mask=True)
    #print "SHAPE", np.shape(xa) , xa   
#    xa = ma.array(xs, copy=True, keep_mask=True)
    return xa
    
def masked_regress(data, xs, stats=[], asDict=False):
    '''
    This routine computes linear regressions on a gridded basis of Masked Arrays!@!
    It returns slopes and then offsets
    The optional third parameter specifies some descriptive stats to return also returned in order specified.
    see http://en.wikipedia.org/wiki/Simple_linear_regression
    Values are any of 
        "sse" : sum squared error
        "ssb" : sum square on beta
        "ssa" : sum squared on alpha
        "sx"  : sum of x
        "sxx" : sum of squares of x
        "sxy" : sum of x * y
        "sy"  : sum of y
        "syy  : sum of squares of y
        "ssx" : sum squared of x
        "ssy" : sum squared of y
    '''
    #data is expected to be a masked array enforce this 
    n=data.count(0) #count the non masked 
    
    xa=xs_matched_to_data( data, xs)
    
    rank=len(np.shape(data))
    #print np.shape(data), len(xs)
    if rank == 1:
        #n = ma.shape(data)[0]
        sy = ma.sum(data)
        syy = ma.sum(data * data)
        sxy = ma.sum(data * xa)
        sx = ma.sum(xa)
        sxx = ma.sum(xs * xa)
        ssx = sx * sx
    else:
        
        #(n,x,y)=ma.shape(data)  ###OOOPS needless and incorrect over write on n
        #print "masked_regress.shape(data)", ma.shape(data)
        #print "XA", xa.count(0)
        #xa=ma.array(xa)
        sx = ma.sum(xa, axis=0)
        sxx= ma.sum(xa * xa, axis=0)
        ssx = sx * sx
        sy = ma.sum(data, axis=0)
        syy = ma.sum(data * data, axis=0)
        sxy= ma.sum(data * xa, axis = 0)
        
#    print "shapes",ma.shape(sx),ma.shape(xa), ma.shape(sy), ma.shape(syy), ma.shape(ssx), ma.shape(sxy)
    beta = (n * sxy - sx * sy)/(n * sxx - ssx)
    alpha = (sy-beta*sx)/n
    if stats != []:
        if asDict:
            returns=dict()
        else:
            returns = []
        #then we set up to return extra statistics
        statset = set(stats)
        nullset = set([])
        if statset & set(["sse", "ssy", "ssb", "ssa"]) != nullset:
            ssy=sy*sy
        if statset & set(["sse", "ssb", "ssa"]) != nullset:
            sse= (1.0/(n * (n -2)))*(n * syy - ssy - beta * beta *(n * sxx - ssx))
        if statset & set(["ssb", "ssa"]) != nullset:
            ssb= (n * sse)/(n*sxx-ssx)
        if "ssa" in statset :
            ssa= ssb * sxx/n
        #print stats
        if not asDict:
            for k in stats:
                if k == "sse": returns.append((k,sse))
                if k == "ssb": returns.append((k,ssb))
                if k == "ssa": returns.append((k,ssa))
                if k == "sx":  returns.append((k,sx))
                if k == "sxx": returns.append((k,sxx))
                if k == "sxy": returns.append((k,sxy))
                if k == "sy" : returns.append((k,sy))
                if k == "syy": returns.append((k,syy))
                if k == "ssx": returns.append((k,ssx))
                if k == "ssy": returns.append((k,ssy))
                #print k, len(returns)
        else:
            for k in stats:
                if k == "sse": returns[k] =sse
                if k == "ssb": returns[k] =ssb
                if k == "ssa": returns[k] =ssa
                if k == "sx":  returns[k] =sx
                if k == "sxx": returns[k] =sxx
                if k == "sxy": returns[k] =sxy
                if k == "sy" : returns[k] =sy
                if k == "syy": returns[k] =syy
                if k == "ssx": returns[k] =ssx
                if k == "ssy": returns[k] =ssy
        
        return beta, alpha, returns
    else:
        return beta, alpha
 
def analysed_regress(data, xs):
    #print type(data), type(ma.masked_array([1]))
    if isinstance(data, ma.masked_array):
        beta, alpha, stats = masked_regress(data, xs, statkeylist, True)
        n = data.count(0)
        #print "CALLED ",beta, alpha
#        print "masked"
    else:
        beta, alpha, stats = regress(data, xs, statkeylist, True)
        n = len(xs)
    #print beta, alpha        
    stats["n"]=n
    #print 'stats["n"]', stats["n"]
    Sxx = (n*stats["sxx"] - stats["sx"]*stats["sx"])/n
    #print "Sxx",Sxx
    Syy = (n*stats["syy"] - stats["sy"]*stats["sy"])/n
    #print 'stats["syy"]', np.shape(stats["syy"]), stats["syy"], 'stats["sy"]', stats["sy"]
    #print "Syy", Syy
    Sxy=  (n*stats["sxy"] - stats["sx"] * stats["sy"])/n
    SSR = beta*Sxy
    SSE= Syy - SSR
    stats["SSR"]=SSR
    stats["SSE"]=SSE
    sigma=np.sqrt(SSE/(n-2.))
    stats["sigma"]=sigma
    t95=t.ppf(0.975, n-2)
    stats["t95"]=t95
    sqrtsxx=np.sqrt(Sxx)
    beta_conf=t95 * sigma /sqrtsxx
    stats["beta_conf"]=beta_conf
    stats["beta"]=beta
    stats["alpha"]=alpha
    stats["t"]=sqrtsxx*(beta)/sigma
    stats["prob"]=t.pdf(stats["t"], n-2)
    y_conf=np.array([t95*sigma*np.sqrt(1./n + ((x-stats["sx"]/n) * (x-stats["sx"]/n))/Sxx) for x in xs])
    stats["y_conf"]=y_conf
    stats["mse"]= SSE/n
    stats["stderr"]=stats["mse"]/sqrtsxx
    stats["rsq"]=(Sxy * Sxy)/(Sxx * Syy)
    return stats
    
def residuals(data, xs, stats,dump=False):
    '''
    given an analyis compute residuals, does not work for masked arrays.
    '''
   #if data is masked then ensure that xs are masked
    if isinstance(data, np.ma.masked_array) and not isinstance(xs, np.ma.masked_array):
        yhat=xs_matched_to_data( data, xs)
        yhat=stats["beta"] * yhat  + stats["alpha"] 
        #print "HERE"
    else:
        #print "THRERE"
        yhat=np.array([stats["beta"] * x  + stats["alpha"] for x in xs])
    resid = data-yhat
    if dump:
        print "data"
        print data
        print "yhat"
        print yhat
        print "xs"
        print xs
        print "stats"
        print stats
    #    sys.exit()
    return yhat, resid
    
def multi_process_regress(data, xs):
    '''
    an implementation of the ricketts multi-process regression technique
    Given the data, determine the set of regression equations that explain the data
    as being explained by independent models 
    '''
    stats = analysed_regress(data, xs)
    Yhat, resid = residuals(data, xs, stats)
    #ACCUMULATE the statistics of the data points
    #df=[[t.fit(resid[:, y1, x1]) for x1 in range(np.shape(data)[2])] for y1 in range(np.shape(data)[1])]  
    
    
#if __name__ == "__main__":
#    import sys
##AN INTERPRETER   
#
#    xs=globalw.QCCCEgw('csiro_mk3_5')
#    x = list()
#    for k in sorted(xs.keys())[:100] :
#        x.append(xs[k])
#        
#    print x
#    nc = Nio.open_file('~/ozclim/pr_regressions/pr/csiro_mk3_5/run1/pr_A1_2001_2100.nc','r')
#    data= nc.variables['pr'][100:]
#    for i in [11]:#range(12):
#        data=nc.variables['pr'][i:1200+i:12]
#        data =ma.masked_array(data)
#        stats=analysed_regress(data, x)
#        yhat, r=residuals(data, x, stats)
#        print r
#        #print stats
#        print np.shape(data), np.shape(np.array(x)), np.shape(stats["alpha"]), np.shape(r)
#        #sys.exit()
#        df=[[t.fit(r.data[:, y1, x1]) for x1 in range(np.shape(data)[2])] for y1 in range(np.shape(data)[1])] 
#        print np.shape(df)
#        #sys.exit()
#        
#        for m in range(len(x)):
#            print x[m], data[m][50][50], yhat[m][50][50], r[m][50][50], yhat[m][50][50] - stats["y_conf"][m][50][50], 
#            print yhat[m][50][50] + stats["y_conf"][m][50][50]  , 
#            print t.cdf(r[m][50][50], df[50][50][0], df[50][50][1], df[50][50][2]), 
#            print t.pdf(r[m][50][50], df[50][50][0], df[50][50][1], df[50][50][2]), 
#            print t.isf(0.05, df[50][50][0], df[50][50][1], df[50][50][2]), stats["t"][50][50], stats["mse"][50][50]
#        #print t.fit(r[50][50])#, len(x)-2)    
#        
#    sys.exit()
#    
#    print np.shape(nc.variables['pr']), np.shape(data)
#    print type(data)
#    xi=np.array([1.47 ,1.50,1.52,1.55,1.57 ,1.60,1.63,1.65,1.68,1.70,1.73,1.75,1.78,1.80,1.83])
#    yis=np.array([[[ 52.21,53.12,54.48,55.84,57.20,58.57,59.93,61.29,63.11,64.47,66.28,68.10,69.92,72.19,74.46], [ 52.21,53.12,54.48,55.84,57.20,58.57,59.93,61.29,63.11,64.47,66.28,68.10,69.92,72.19,74.46]],[[ 52.21,53.12,54.48,55.84,57.20,58.57,59.93,61.29,63.11,64.47,66.28,68.10,69.92,72.19,74.46], [ 52.21,53.12,54.48,55.84,57.20,58.57,59.93,61.29,63.11,64.47,66.28,68.10,69.92,72.19,74.46]]])
#    yiz=[ 52.21,53.12,54.48,55.84,57.20,58.57,59.93,61.29,63.11,64.47,66.28,68.10,69.92,72.19,74.46]
#    yy=[]
#    for r in yiz:
#        yy.append([[r,r],[r,r*2]])
#    yi = np.array(yy)
#    #print regress(yi, xi,statkeylist), np.shape(yi), np.shape(xi)
#    print multi_process_regress(yi, xi)
#
#    xi=np.array([1.5, 1.8, 2.4, 3.0, 3.5, 3.9, 4.4, 4.8, 5.0])
#    yi=np.array([4.8, 5.7, 7.0, 8.3, 10.9, 12.4, 13.1, 13.6, 15.3])
#    stats= multi_process_regress(yi, xi)
#    yhat, r=residuals(yi, xi, stats)
#    
#    print "=================="
#    print stats
#    for m in range(len(yi)):
#        print xi[m], yi[m],yhat[m], r[m]
#    sys.exit()
#    
#    d=[
#    [16.10435577	,9.37e-05],
#    [16.19011795	,6.82e-05],
#    [16.1262555	,1.32e-04],
#    [15.91938073	,6.94e-05],
#    [15.92240431	,1.91e-04],
#    [16.05840871	,1.49e-05],
#    [16.2221622	,2.20e-05],
#    [16.28947209	,3.56e-05],
#    [16.21970378	,1.36e-05],
#    [16.28653329	,7.81e-06],
#    [16.42200079	,6.60e-05],
#    [16.63373589	,1.77e-05],
#    [16.77101189	,1.04e-04],
#    [16.37523427	,3.41e-05],
#    [16.41202582	,2.10e-05],
#    [16.58295677	,6.35e-05],
#    [16.80195409	,4.07e-05],
#    [16.45737946	,1.21e-04],
#    [16.58236335	,2.17e-05],
#    [16.75806986	,2.47e-05],
#    [17.01801262	,2.42e-05],
#    [16.74461918	,4.38e-05],
#    [16.75900236	,2.54e-05],
#    [16.98998096	,3.10e-05],
#    [16.96983321	,1.35e-04],
#    [16.5922253	,3.91e-05],
#    [16.7225499	,3.82e-05],
#    [17.0836553	,1.46e-04],
#    [17.01467821	,2.07e-05],
#    [16.88791408	,7.71e-05],
#    [16.97681286	,7.26e-05],
#    [17.03205671	,4.95e-06],
#    [16.93160051	,7.88e-05],
#    [17.03977106	,1.17e-04],
#    [17.30143755	,2.09e-05],
#    [17.16667648	,1.21e-04],
#    [17.23455153	,3.73e-05],
#    [17.3806722	,8.56e-05],
#    [17.40590634	,4.78e-05],
#    [17.31322102	,9.16e-05],
#    [17.39465977	,4.65e-05],
#    [17.54682759	,5.81e-05],
#    [17.19651664	,1.98e-05],
#    [17.38310236	,4.45e-05],
#    [17.50596693	,4.56e-05],
#    [17.40062215	,2.91e-05],
#    [17.50949914	,1.13e-04],
#    [17.77277632	,3.94e-05],
#    [17.69354167	,1.65e-04],
#    [17.52077397	,5.24e-06],
#    [17.75237424	,1.72e-05],
#    [17.79891471	,3.35e-05],
#    [17.78159272	,4.03e-05],
#    [17.92714824	,8.81e-05],
#    [17.94413113	,2.90e-05],
#    [17.73632386	,1.11e-05],
#    [17.87458888	,6.01e-05],
#    [18.05357328	,3.28e-05],
#    [17.8662246	,1.77e-04],
#    [18.01731862	,8.89e-06],
#    [18.23837876	,6.09e-05],
#    [17.94500712	,1.58e-04],
#    [18.0575011	,1.49e-05],
#    [18.1814112	,8.58e-05],
#    [18.19455104	,6.05e-05],
#    [18.1923752	,6.92e-05],
#    [18.21404887	,8.35e-05],
#    [18.27019696	,2.76e-05],
#    [18.4290336	,1.33e-04],
#    [18.287745	,4.54e-05],
#    [18.33671563	,4.92e-05],
#    [18.53726067	,5.02e-05],
#    [18.3968198	,4.24e-05],
#    [18.60434449	,2.05e-06],
#    [18.74063148	,3.68e-05],
#    [18.62994599	,2.71e-05],
#    [18.50103427	,2.93e-05],
#    [18.50312534	,2.60e-05],
#    [18.58671168	,5.94e-05],
#    [18.71565165	,1.99e-04],
#    [18.67501705	,3.73e-05],
#    [18.70895457	,9.15e-05],
#    [18.56760946	,1.46e-05],
#    [18.62296633	,3.90e-05],
#    [18.74780894	,5.06e-06],
#    [18.61268052	,2.78e-05],
#    [18.72805679	,9.72e-05]
#    ]
#
#    yi = d[1]
#    xi = d[0:1]
#    stats= multi_process_regress(yi, xi)
#    yhat, r=residuals(yi, xi, stats)
#    
#    print "=================="
#    print stats
#    for m in range(len(yi)):
#        print xi[m], yi[m],yhat[m], r[m]
#

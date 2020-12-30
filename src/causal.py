from sklearn.linear_model import LinearRegression
from scipy import stats
from statsmodels.tsa.stattools import grangercausalitytests
import os
os.environ['R_HOME'] = '/global/homes/z/zhangtao/.conda/envs/luffy/lib/R'
from rpy2.robjects.packages import importr
import rpy2.robjects as ro
import rpy2.robjects.numpy2ri
rpy2.robjects.numpy2ri.activate()
pcalg = importr('pcalg')

sys.path.append("/global/homes/z/zhangtao/climate_causal/notears/notears/")
import linear
import utils
import nonlinear
from sklearn import preprocessing

def calc_lag_corr(x,y,max_lag=1):
    if np.ma.is_masked(y) or np.ma.is_masked(x):
        corr = y[0] + x[0]
        return corr, corr
    else:
        max_corr = 0
        max_corr_p = 0
        ntime = len(x)
        for i in range(max_lag):
            corr, corr_p = stats.pearsonr(x[:ntime-i], y[i:])
            if abs(corr) > abs(max_corr):
                max_corr = corr
                max_corr_p = corr_p
        
        return max_corr, max_corr_p

def calc_corr(x,y):
    if np.ma.is_masked(y) or np.ma.is_masked(x):
        corr = y[0] + x[0]
        return corr, corr
    else:
        corr, corr_p = stats.pearsonr(x, y)
        
        return corr, corr_p

def calc_reg(x,y):
    if np.ma.is_masked(y) or np.ma.is_masked(x):
        reg = y[0] + x[0]
        return reg, reg
    else:
        X = sm.add_constant(x)
        #X = np.column_stack((x, x))
        reg_mod = sm.OLS(y,X).fit()
        reg = reg_mod.params[1]
        #reg = reg_mod.params[1]
        #reg_p = reg_mod.pvalues[0] + reg_mod.pvalues[0]
        reg_p = reg_mod.pvalues[1]
        return reg, reg_p
            
def calc_lingam(x,y):
    if np.ma.is_masked(y) or np.ma.is_masked(x):
        lingam = y[0] + x[0]
        return lingam,lingam
    else:
        data = np.transpose(np.array([x, y]))
        ro.globalenv['aa'] = data
        ro.r('res <- lingam(aa)')
        ro.r('resmat <- as(res, "amat")')
        resmat = np.array(ro.r['resmat'])
        lingam = resmat[0,1]
        return resmat[0,1], resmat[1,0]
    
def calc_pc(x,y):
    if np.ma.is_masked(y) or np.ma.is_masked(x):
        lingam = y[0] + x[0]
        return lingam,lingam
    else:
        data = np.transpose(np.array([x, y]))
        ro.globalenv['aa'] = data
        ro.r('suffStat <- list(C = cor(aa), n = nrow(aa))')
        ro.r('varname <- c("a","b")')
        ro.r('pc.fit = pc(suffStat, indepTest=gaussCItest, labels=varname, alpha=0.05)')
        ro.r('pcmat <- as(pc.fit@graph, "matrix")')
        pcmat = np.array(ro.r['pcmat'])
        #print(pcmat)
        return pcmat[0,1], pcmat[1,0]
    
def calc_granger(x,y):
    if np.ma.is_masked(y) or np.ma.is_masked(x):
        granger = y[0] + x[0]
        return granger, granger
    else:
        nlag = 7
        data = np.transpose(np.array([y, x]))
        granger_model = grangercausalitytests(data, nlag, verbose=False)
        granger = granger_model[1][0]['ssr_chi2test'][1]
        
        data = np.transpose(np.array([x, y]))
        granger_model = grangercausalitytests(data, nlag, verbose=False)
        granger_r = granger_model[1][0]['ssr_chi2test'][1]
                
        return granger, granger_r 

def causality_test(var1,var2):
    data1= pd.Series(var1, name='Var1')
    data2= pd.Series(var2, name='Var2')
    mdata=pd.concat([data1, data2], axis=1)
    #data = np.log(mdata).diff().dropna()
    data=mdata

    model = VAR(data)
    results=model.fit(7)
    #print(results.summary())

    foo=crit=results.test_causality('Var2', ['Var1'], kind='f')
    crit=foo.crit_value
    stat=foo.test_statistic
    if(stat > crit):
        cause=1
    else:
        cause=0

    return cause, cause
    
def calc_notears_block(x,y):
    block_x = y.shape[1]
    blocl_y = y.shape[2]
    
    notears = np.zeros_like(y[0,:,:])
    
    for i in range(block_x):
        for j in range(block_y):    
            if np.ma.is_masked(y) or np.ma.is_masked(x):
                notears = y[0,0,0] + x[0]
            else:
                data = np.transpose(np.array([x[:], y[:,i,j]]))         
                aa = linear.notears_linear(data, lambda1=0.1, loss_type='l2', max_iter=200)
                notears[i,j] = aa[0,1]

    return notears

def calc_notears(x,y):
    if np.ma.is_masked(y) or np.ma.is_masked(x):
        notears = y[0] + x[0]
        return notears, notears
    else:
        data = np.transpose(np.array([x, y]))         
        aa = linear.notears_linear(data, lambda1=0.1, loss_type='l2', max_iter=200, w_threshold=0.05)
        notears = aa[0,1]
        return aa[0,1], aa[1,0]

def calc_corr_wrapper(args):
    #return calc_corr(*args)
    return calc_lag_corr(*args)

def calc_reg_wrapper(args):
    return calc_reg(*args)

def calc_notears_wrapper(args):
    return calc_notears(*args)

def calc_lingam_wrapper(args):
    return calc_lingam(*args)

def calc_pc_wrapper(args):
    return calc_pc(*args)

def calc_granger_wrapper(args):
    return calc_granger(*args)
    #return causality_test(*args)

def causal_algs_2d_to_2d(data,nlat,nlon):
    results={}

    pool = mp.Pool(10)
    logger.info("Calc correlation")
    r = pool.map(calc_corr_wrapper,data)
    a = np.ma.masked_array([c[0] for c in r]).reshape(nlat,nlon)
    p = np.ma.masked_array([c[1] for c in r]).reshape(nlat,nlon)
    results['correlation'] = a
    results['correlation_p'] = p
    
#     logger.info("Calc regression")
#     r = pool.map(calc_reg_wrapper,data)
#     a = np.ma.masked_array([c[0] for c in r]).reshape(nlat,nlon)
#     p = np.ma.masked_array([c[1] for c in r]).reshape(nlat,nlon)
#     results['reg'] = a
#     results['reg_p'] = p

    logger.info("Calc notears linear")
    r = pool.map(calc_notears_wrapper,data)
    a = np.ma.masked_array([c[0] for c in r]).reshape(nlat,nlon)
    b = np.ma.masked_array([c[1] for c in r]).reshape(nlat,nlon)
    results['notears_linear'] = a
    results['notears_linear_r'] = b 

    logger.info("Calc lingam")
    r = pool.map(calc_lingam_wrapper,data)
    a = np.ma.masked_array([c[0] for c in r]).reshape(nlat,nlon)
    b = np.ma.masked_array([c[1] for c in r]).reshape(nlat,nlon)
    results['lingam'] = a
    results['lingam_r'] = b
    
    logger.info("Calc Granger")
    r = pool.map(calc_granger_wrapper,data)
    a = np.ma.masked_array([c[0] for c in r]).reshape(nlat,nlon)
    b = np.ma.masked_array([c[1] for c in r]).reshape(nlat,nlon)
    results['granger'] = a
    results['granger_r'] = b
    
    logger.info("Calc pc")
    r = pool.map(calc_pc_wrapper,data)
    a = np.ma.masked_array([c[0] for c in r]).reshape(nlat,nlon)
    b = np.ma.masked_array([c[1] for c in r]).reshape(nlat,nlon)
    results['pc'] = a
    results['pc_r'] = b
        
    return results

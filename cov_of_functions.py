import numpy as np




def jacfun(fun,xd,params):
    npar = len(params)
    valf = fun(xd,*params)
    dpar = 1e-5*params
    wzero = np.where(dpar == 0)
    dpar[wzero]=1e-5
    dgi = []
    for k in range(npar):
        pars = params.copy()
        pars[k]+=dpar[k]
        dgi.append((fun(xd,*pars)-valf)/dpar[k])
    return np.array(dgi).transpose()



def hesfun(fun,xd,params):
    jacv = jacfun(fun,xd,params)
    npar = len(params)
    dpar = 1e-5*params
    wzero = np.where(dpar == 0)
    dpar[wzero]=1e-5
    ddgi = []
    for j in range(npar):
        pars = params.copy()
        pars[j]+=dpar[j]
        ddgi.append(((jacfun(fun,xd,pars)-jacv)/dpar[j]).transpose())
    return np.array(ddgi)



def covma(fun,xd,yval,stds,params):
    j = jacfun(fun,xd,params)
    h = hesfun(fun,xd,params)
    jj = np.dot(j.transpose(),j)
    f = fun(xd,*params)
    r = f - yval
    dchi = jj + (r*h).sum(axis=2)
    cov = np.linalg.inv(dchi)*stds**2
    return cov
    
    

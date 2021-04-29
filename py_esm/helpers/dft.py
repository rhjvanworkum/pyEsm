import numpy as np

def LSDA(rho):
    clda = (3 / 4) * (3.0 / np.pi) ** (1 / 3)
    return -clda * rho ** (4 / 3)

def xs(rho, alpha=2/3):
    fac = -2.25 * alpha * np.power(0.75/np.pi, 1/3)
    rho3 = np.power(rho, 1/3)
    fx = fac * rho * rho3
    dfxdna = (4/3) * fac * rho3
    return fx, dfxdna

def get_xc(grid, P):
    rho = grid.get_rho(P)

    # print(rho)

    fx, dfxa = xs(rho)
    fc, dfca, dfcb = cvwn5(rho, rho)
    # fx = dfxa = LSDA(rho)
    # fc, dfca, dfcb = 0, 0, 0

    w = grid.points[:, 3]

    vxc = np.einsum('g,g,gI,gJ->IJ', w, dfxa+dfca, grid.funcs, grid.funcs)
    # print('VXCs: ')
    # print(vxc)
    # print(grid.integrate(dfxa, grid.funcs, grid.funcs))

    exc = np.dot(w, 2 * fx + fc)
    # exc = excfunction(rho, np)

    return exc, vxc

def zero_low_density(rho,cut=1e-10):
    rho[rho<cut]=0
    return rho

def cvwn5(rhoa,rhob,tol=1e-10):
    rhoa = zero_low_density(rhoa)
    rhob = zero_low_density(rhob)

    ecs = []
    vcrhoas = []
    vcrhobs = []
    for na,nb in zip(rhoa,rhob):
        rho = na+nb
        ec = vcrhoa = vcrhob = 0
        if rho>tol:
            zeta=(na-nb)/rho
            x = pow(3./4./np.pi/rho,1/6.)
            epsp = vwn_epsp(x)
            epsf = vwn_epsf(x)
            g = vwn_g(zeta)
            eps = epsp + g*(epsf-epsp)
            ec = eps*rho

            depsp = vwn_depsp(x)
            depsf = vwn_depsf(x)
            dg = vwn_dg(zeta)
            deps_dx = depsp + g*(depsf-depsp)
            deps_dg = (epsf-epsp)*dg
            vcrhoa = eps - (x/6.)*deps_dx + deps_dg*(1-zeta)
            vcrhob = eps - (x/6.)*deps_dx - deps_dg*(1+zeta)
        ecs.append(ec)
        vcrhoas.append(vcrhoa)
        vcrhobs.append(vcrhob)
    return np.array(ecs),np.array(vcrhoas),np.array(vcrhobs)


def vwn_depsp(x): return vwn_deps(x,0.0310907,-0.10498,3.72744,12.9352)
#def vwn_depsf(x): return vwn_deps(x,0.01554535,-0.32500,7.06042,13.0045)
def vwn_depsf(x): return vwn_deps(x,0.01554535,-0.32500,7.06042,18.0578)
def vwn_deps(x,a,x0,b,c):
    q = np.sqrt(4*c-b*b)
    deps = a*(2/x - (2*x+b)/vwn_xx(x,b,c)
              - 4*b/(np.power(2*x+b,2)+q*q) - (b*x0/vwn_xx(x0,b,c))
              * (2/(x-x0)-(2*x+b)/vwn_xx(x,b,c)-4*(2*x0+b)/(np.power(2*x+b,2)+q*q)))
    return deps

def vwn_g(z): return 1.125*(np.power(1+z,4./3.)+np.power(1-z,4./3.)-2)
def vwn_dg(z): return 1.5*(np.power(1+z,1./3.)-np.power(1-z,1./3.))

def vwn_xx(x,b,c): return x*x+b*x+c
def vwn_epsp(x): return vwn_eps(x,0.0310907,-0.10498,3.72744,12.9352)
#def vwn_epsf(x): return vwn_eps(x,0.01554535,-0.32500,7.06042,13.0045)
def vwn_epsf(x): return vwn_eps(x,0.01554535,-0.32500,7.06042,18.0578)

def vwn_eps(x,a,x0,b,c):
    Q = np.sqrt(4*c-b*b)
    eps = a*(np.log(x*x/vwn_xx(x,b,c))
             - b*(x0/vwn_xx(x0,b,c))*np.log(np.power(x-x0,2)/vwn_xx(x,b,c))
             + (2*b/Q)*(1-(x0*(2*x0+b)/vwn_xx(x0,b,c))) * np.arctan(Q/(2*x+b)))
    #eps = a*(np.log(x*x/vwn_xx(x,b,c)) + (2*b/Q)*np.arctan(Q/(2*x+b))
    #         - (b*x0/vwn_xx(x0,b,c))*np.log(np.power(x-x0,2)/vwn_xx(x,b,c))
    #         + (2*(b+2*x0)/Q)*np.arctan(Q/(2*x+b)))
    return eps

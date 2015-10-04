

def f(y,dy,mu):
    return dy, -y -mu*(y**2 - 1)*dy

mu = 1.350

def get_K1(y_n,dy_n,mu,h,f):
    return h*f(y_n, dy_n, mu)

def get_K2(y_n,dy_n,mu,h,f):
    K1=get_K1(y_n, dy_n, mu, h, f)
    return h*f(y_n + K1[0]/2., dy_n + K1[1]/2., mu)

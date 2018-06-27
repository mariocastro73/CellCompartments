# 2018. Compartment with division, death and differentiation
# Study fate of a cell: how many cells exit? How long does it take?
# Generate realisations and plot histograms

# print statement syntax is python2

import numpy as np
from math import factorial as fa
from pylab import figure,hist,plot,show,xlim,xlabel,ylabel,xticks,yticks,ylim,legend,title,rc,text

mu = 0.05  # death rate
lamb = 0.20 # division rate
nu = 0.2 # exit rate

mut = mu/(mu+lamb+nu)
nut = nu/(mu+lamb+nu)
lat = lamb/(mu+lamb+nu)
ett = nu/(mu+lamb+nu)

def myprint(mylist):
    try:
#        [print("{:.3f}".format(x),end=' ') for x in mylist]
        for p in mylist: print p
#        [print(x) for x in mylist]
    except:
        for p in mylist: print p

def f(x):
    return (1-np.sqrt(1-4*x))

def ft(t):
    delta = mu+nu-lamb
    return delta**2*(mu+nu)*np.exp(delta*t)/(lamb-(mu+nu)*np.exp(delta*t))**2

beta = lat/np.sqrt(1-4*mut*lat)
delta = np.sqrt((mu+lamb+nu)**2-4*mu*lamb)

# predicted 
q0 = 0.5/lat*(1-np.sqrt(1-4*lat*mut))
q1 = nut/np.sqrt(1-4*mut*lat)
q2 = lat*q1*q1/(1-f(lat*mut))
q3 = lat*2*q2*q1/(1-f(lat*mut))
qpred = [q0,q1]
n=2
while n < 40:
    thisq = 0
    for i in range(1,n):
        thisq += qpred[i]*qpred[n-i]
    qpred.append(thisq*lat/(1-f(lat*mut)))
    n += 1
print('prediction')
myprint(qpred[:10])

def onecellsfate(mu,nu,lamb):
    ''' Gillespie realisation, starting with one cell.
        returns number of cells exiting, dying, dividing,
        time of extinction and (conditional) time of first passage'''
    n = 1
    x,y,z = 0,0,0
    t = 0
    s = tmax*0.9
    while n > 0:
        sumrates = mu+nu+lamb
        t += -np.log(np.random.random())/sumrates
        urv = np.random.uniform()
        if urv < mu/sumrates:  # death
            n -= 1
            y += 1
        elif urv < (mu+nu)/sumrates:   #exit
            n -= 1
            x += 1
            if x == 1:
                s = t
        else:   # division
            n += 1
            z += 1
    return x,y,z,t,s

tmax = 50.0
fates=[]
for i in range(100000):
    fates.append(onecellsfate(mu,nu,lamb))
print'last realisation:',fates[-1]
print'mean number of cells exiting =',np.mean([f[0] for f in fates])
print'mean-square number of cells exiting =',np.mean([f[0]**2 for f in fates])

f = figure()
rc('font', family='serif', serif='cm10')
rc('text', usetex=True)

mybins = np.arange(100)-0.5

ax = f.add_axes([0.1,0.55,0.375,0.40])
qlist = hist([thisfate[0] for thisfate in fates],bins=mybins,normed=True)
plot(qpred,'ro',label='exact')
xlim([-0.5,19.5])
ylim([0,0.6])
yticks([0.1,0.3,0.5])
title('number of cells exiting')
ylabel('$P$')
text(5,0.5,'$\lambda = $'+str(lamb)+'$\quad\mu = $'+str(mu)+' $\quad\\nu = $'+str(nu))

ax = f.add_axes([0.55,0.55,0.375,0.40],yscale='log')
qlist = hist([thisfate[0] for thisfate in fates],bins=mybins,normed=True,label='numerical')
plot(qpred,'ro',label='exact')
xlim([-0.5,29.5])
ylim([0.0003,0.6])
yticks([0.001,0.01,0.1])
title('number of cells exiting')
#ylabel('$P$')
legend()

ax = f.add_axes([0.55,0.30,0.375,0.16])
numans=hist([thisfate[1] for thisfate in fates],bins=mybins,normed=True)
xlim([-0.5,10.5])
ylim([0,0.7])
yticks([0.1,0.3,0.5])
xticks(range(10))
title('total number of cells dying')

ax = f.add_axes([0.1,0.30,0.375,0.16])
ax.hist([thisfate[2] for thisfate in fates],bins=mybins,normed=True)
xlim([-0.5,10.5])
ylim([0,0.7])
yticks([0.1,0.3,0.5])
xticks(range(10))
title('number of rounds of division')
ylabel('$P$')

print [thisfate[1] for thisfate in fates[-10:]]
print [thisfate[2] for thisfate in fates[-10:]]

ax = f.add_axes([0.1,0.04,0.375,0.16])
tbins = np.arange(100)*0.01*tmax
ax.hist([thisfate[3] for thisfate in fates],bins=tbins,normed=True)
ax.plot(tbins,ft(tbins),'r')
xlim([0,tmax*0.6])
ylabel('$P$')
title('extinction time')

ax = f.add_axes([0.55,0.04,0.375,0.16])
ax.hist([thisfate[4] for thisfate in fates],bins=tbins,normed=True)
xlim([0,tmax*0.6])
title('first passage time')

#myn,myp,myc = [],[],[]

qlist = qlist[0]
print('numerical:')
myprint(qlist[:10])
predm = nu/(mu+nu-lamb)
print'predicted mean, mean square =',predm,2*predm + 2*lamb*predm*predm/(mu+nu-lamb)

show()


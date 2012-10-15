Z
=

Z-analysis data handling

from pylab import*
import numpy as np
params={'legend.fontsize':8}#,'text.usetex': True}
rcParams.update(params)
def calc_coeff(num_points, pol_degree, diff_order=0):

    """ calculates filter coefficients for symmetric savitzky-golay filter.
        see: http://www.nrbook.com/a/bookcpdf/c14-8.pdf

        num_points   means that 2*num_points+1 values contribute to the
                     smoother.

        pol_degree   is degree of fitting polynomial

        diff_order   is degree of implicit differentiation.
                     0 means that filter results in smoothing of function
                     1 means that filter results in smoothing the first 
                                                 derivative of function.
                     and so on ...

    """

    # setup interpolation matrix
    # ... you might use other interpolation points
    # and maybe other functions than monomials ....

    x = arange(-num_points, num_points+1, dtype=int)
    monom = lambda x, deg : pow(x, deg)

    A = zeros((2*num_points+1, pol_degree+1), float)
    for i in range(2*num_points+1):
        for j in range(pol_degree+1):
            A[i,j] = monom(x[i], j)
        
    # calculate diff_order-th row of inv(A^T A)
    ATA = dot(A.transpose(), A)
    rhs = zeros((pol_degree+1,), float)
    rhs[diff_order] = (-1)**diff_order
    wvec = linalg.solve(ATA, rhs)

    # calculate filter-coefficients
    coeff = dot(A, wvec)

    return coeff
def smooth(signal, coeff):
    
    """ applies coefficients calculated by calc_coeff()
        to signal """
    
    N = (size(coeff)-1)/2
    res = convolve(signal, coeff)
    return res[N:-N]
def peakdet(v, delta, x = None):
   
    maxtab = []
    mintab = []
       
    if x is None:
        x = arange(len(v))
    
    v = asarray(v)
    
    if len(v) != len(x):
        sys.exit('Input vectors v and x must have same length')
    
    if not isscalar(delta):
        sys.exit('Input argument delta must be a scalar')
    
    if delta <= 0:
        sys.exit('Input argument delta must be positive')
    
    mn, mx = Inf, -Inf
    mnpos, mxpos = NaN, NaN
    
    lookformax = True
    
    for i in arange(len(v)):
        this = v[i]
        if this > mx:
            mx = this
            mxpos = x[i]
        if this < mn:
            mn = this
            mnpos = x[i]
        
        if lookformax:
            if this < mx-delta:
                maxtab.append((mxpos))#, mx))
                mn = this
                mnpos = x[i]
                lookformax = False
        else:
            if this > mn+delta:
                mintab.append((mnpos, mn))
                mx = this
                mxpos = x[i]
                lookformax = True

    return maxtab#, mintab
T60=loadtxt('/home/pharnedy/Data_/Z-analysis/DashSuzieComparison_61668/60.000000.dat')
T80=loadtxt('/home/pharnedy/Data_/Z-analysis/DashSuzieComparison_61668/80.000000.dat')  
T100=loadtxt('/home/pharnedy/Data_/Z-analysis/DashSuzieComparison_61668/100.000000.dat')	
T120=loadtxt('/home/pharnedy/Data_/Z-analysis/DashSuzieComparison_61668/120.000000.dat')	
T140=loadtxt('/home/pharnedy/Data_/Z-analysis/DashSuzieComparison_61668/140.000000.dat')	
T160=loadtxt('/home/pharnedy/Data_/Z-analysis/DashSuzieComparison_61668/160.000000.dat')	
T180=loadtxt('/home/pharnedy/Data_/Z-analysis/DashSuzieComparison_61668/180.000000.dat')	
T200=loadtxt('/home/pharnedy/Data_/Z-analysis/DashSuzieComparison_61668/200.000000.dat')	
T220=loadtxt('/home/pharnedy/Data_/Z-analysis/DashSuzieComparison_61668/220.000000.dat')
T240=loadtxt('/home/pharnedy/Data_/Z-analysis/DashSuzieComparison_61668/240.000000.dat')
T260=loadtxt('/home/pharnedy/Data_/Z-analysis/DashSuzieComparison_61668/260.000000.dat')
T280=loadtxt('/home/pharnedy/Data_/Z-analysis/DashSuzieComparison_61668/280.000000.dat')
T300=loadtxt('/home/pharnedy/Data_/Z-analysis/DashSuzieComparison_61668/300.000000.dat')
T320=loadtxt('/home/pharnedy/Data_/Z-analysis/DashSuzieComparison_61668/320.000000.dat')
T340=loadtxt('/home/pharnedy/Data_/Z-analysis/DashSuzieComparison_61668/340.000000.dat')
T350=loadtxt('/home/pharnedy/Data_/Z-analysis/DashSuzieComparison_61668/350.000000.dat')
T330=loadtxt('/home/pharnedy/Data_/Z-analysis/DashSuzieComparison_61668/330.000000.dat')
T310=loadtxt('/home/pharnedy/Data_/Z-analysis/DashSuzieComparison_61668/310.000000.dat')
T290=loadtxt('/home/pharnedy/Data_/Z-analysis/DashSuzieComparison_61668/290.000000.dat')
T270=loadtxt('/home/pharnedy/Data_/Z-analysis/DashSuzieComparison_61668/270.000000.dat')
T250=loadtxt('/home/pharnedy/Data_/Z-analysis/DashSuzieComparison_61668/250.000000.dat')
T230=loadtxt('/home/pharnedy/Data_/Z-analysis/DashSuzieComparison_61668/230.000000.dat')
T210=loadtxt('/home/pharnedy/Data_/Z-analysis/DashSuzieComparison_61668/210.000000.dat')
T190=loadtxt('/home/pharnedy/Data_/Z-analysis/DashSuzieComparison_61668/190.000000.dat')
T170=loadtxt('/home/pharnedy/Data_/Z-analysis/DashSuzieComparison_61668/170.000000.dat')
T150=loadtxt('/home/pharnedy/Data_/Z-analysis/DashSuzieComparison_61668/150.000000.dat')
T130=loadtxt('/home/pharnedy/Data_/Z-analysis/DashSuzieComparison_61668/130.000000.dat')
T110=loadtxt('/home/pharnedy/Data_/Z-analysis/DashSuzieComparison_61668/110.000000.dat')
T90=loadtxt('/home/pharnedy/Data_/Z-analysis/DashSuzieComparison_61668/90.000000.dat')
T70=loadtxt('/home/pharnedy/Data_/Z-analysis/DashSuzieComparison_61668/70.000000.dat')
T50=loadtxt('/home/pharnedy/Data_/Z-analysis/DashSuzieComparison_61668/50.000000.dat')	


std50=std(T50[:,4])
std60=std(T60[:,4])
std70=std(T70[:,4])
std80=std(T80[:,4])
std90=std(T90[:,4])
std100=std(T100[:,4])
std110=std(T110[:,4])
std120=std(T120[:,4])
std130=std(T130[:,4])
std140=std(T140[:,4])
std150=std(T150[:,4])
std160=std(T160[:,4])
std170=std(T170[:,4])
std180=std(T180[:,4])
#std190=std(T190[:,4])
std200=std(T200[:,4])
std210=std(T210[:,4])
std220=std(T220[:,4])
std230=std(T230[:,4])
std240=std(T240[:,4])
std250=std(T250[:,4])
std260=std(T260[:,4])
std270=std(T270[:,4])
std280=std(T280[:,4])
std290=std(T290[:,4])
std300=std(T300[:,4])
std310=std(T310[:,4])
std320=std(T320[:,4])
std330=std(T330[:,4])
std340=std(T340[:,4])
std350=std(T350[:,4])

Tu=(50,70,90,110,130,150,170,210,230,250,270,290,310,330,350)
Td=(60,80,100,120,140,160,180,200,220,240,260,280,300,320,340)
stdu=(std50,std70,std90,std110,std130,std150,std170,std210,std230,std250,std270,std290,std310,std330,std350)
stdd=(std60,std80,std100,std120,std140,std160,std180,std200,std220,std240,std260,std280,std300,std320,std340)
 
plot(Tu,stdu,'o''b',label='up')
plot(Td,stdd,'o''r',label='down') 
legend(loc='lower right', ncol=1, shadow=True, numpoints=1)
show() 
plot(T60[:,4])
plot(T80[:,4])
plot(T100[:,4])
plot(T120[:,4])
plot(T140[:,4])
plot(T160[:,4])
plot(T180[:,4])
plot(T200[:,4])
plot(T220[:,4])
plot(T240[:,4])
plot(T260[:,4])
plot(T280[:,4])
plot(T300[:,4])
plot(T320[:,4])
plot(T340[:,4])
plot(T350[:,4])
plot(T330[:,4])
plot(T310[:,4])
plot(T290[:,4])
plot(T270[:,4])
plot(T250[:,4])
plot(T230[:,4])
plot(T210[:,4])
#plot(T190[:,4])
plot(T170[:,4])
plot(T150[:,4])
plot(T130[:,4])
plot(T110[:,4])
plot(T90[:,4])
plot(T70[:,4])
plot(T50[:,4])
show() 
 
 
                      ###### Data arrays #########
#0=V, 1=I, 2=F, 3=SE, 4=T, 5=St.Dev.(F), 6=St.Dev.(SE), 7=St.Dev.(I)# 
#ax1.plot(T60[:,1],T60[:,2],label=('T=60K'))
plot(T60[:,1],T60[:,2],'--',label=('T=60K'))
plot(T80[:,1],T80[:,2],'--',label=('T=80K'))
plot(T100[:,1],T100[:,2],'--',label=('T=100K'))
plot(T120[:,1],T120[:,2],'--',label=('T=120K'))
plot(T140[:,1],T140[:,2],'--',label=('T=140K'))
plot(T160[:,1],T160[:,2],'--',label=('T=160K'))
plot(T180[:,1],T180[:,2],'--',label=('T=180K'))
plot(T200[:,1],T200[:,2],'--',label=('T=200K'))
plot(T220[:,1],T220[:,2],'--',label=('T=220K'))
plot(T240[:,1],T240[:,2],'--',label=('T=240K'))
plot(T260[:,1],T260[:,2],'--',label=('T=260K'))
plot(T280[:,1],T280[:,2],'--',label=('T=280K'))
plot(T300[:,1],T300[:,2],'--',label=('T=300K'))
plot(T320[:,1],T320[:,2],'--',label=('T=320K'))
plot(T340[:,1],T340[:,2],'--',label=('T=340K'))
plot(T350[:,1],T350[:,2],label=('T=350K'))
plot(T330[:,1],T330[:,2],label=('T=330K'))
plot(T310[:,1],T310[:,2],label=('T=310K'))
plot(T290[:,1],T290[:,2],label=('T=290K'))
plot(T270[:,1],T270[:,2],label=('T=270K'))
plot(T250[:,1],T250[:,2],label=('T=250K'))
plot(T230[:,1],T230[:,2],label=('T=230K'))
plot(T210[:,1],T210[:,2],label=('T=210K'))
plot(T190[:,1],T190[:,2],label=('T=190K'))
plot(T170[:,1],T170[:,2],label=('T=170K'))
plot(T150[:,1],T150[:,2],label=('T=150K'))
plot(T130[:,1],T130[:,2],label=('T=130K'))
plot(T110[:,1],T110[:,2],label=('T=110K'))
plot(T90[:,1],T90[:,2],label=('T=90K'))
plot(T70[:,1],T70[:,2],label=('T=70K'))
plot(T50[:,1],T50[:,2],label=('T=50K'))
xlabel('Current (mA)')
ylabel('Light (a.u.)')
title('L-I Characteristic (Facet)')
legend(loc='lower right', ncol=1, shadow=True, numpoints=1)
show()
###################################################################
plot(T60[:,1],T60[:,3],'--',label=('T=60K'))
plot(T80[:,1],T80[:,3],'--',label=('T=80K'))
plot(T100[:,1],T100[:,3],'--',label=('T=100K'))
plot(T120[:,1],T120[:,3],'--',label=('T=120K'))
plot(T140[:,1],T140[:,3],'--',label=('T=140K'))
plot(T160[:,1],T160[:,3],'--',label=('T=160K'))
plot(T180[:,1],T180[:,3],'--',label=('T=180K'))
plot(T200[:,1],T200[:,3],'--',label=('T=200K'))
plot(T220[:,1],T220[:,3],'--',label=('T=220K'))
plot(T240[:,1],T240[:,3],'--',label=('T=240K'))
plot(T260[:,1],T260[:,3],'--',label=('T=260K'))
plot(T280[:,1],T280[:,3],'--',label=('T=280K'))
plot(T300[:,1],T300[:,3],'--',label=('T=300K'))
plot(T320[:,1],T320[:,3],'--',label=('T=320K'))
plot(T340[:,1],T340[:,3],'--',label=('T=340K'))
plot(T350[:,1],T350[:,3],label=('T=350K'))
plot(T330[:,1],T330[:,3],label=('T=330K'))
plot(T310[:,1],T310[:,3],label=('T=310K'))
plot(T290[:,1],T290[:,3],label=('T=290K'))
plot(T270[:,1],T270[:,3],label=('T=270K'))
plot(T250[:,1],T250[:,3],label=('T=250K'))
plot(T230[:,1],T230[:,3],label=('T=230K'))
plot(T210[:,1],T210[:,3],label=('T=210K'))
plot(T190[:,1],T190[:,3],label=('T=190K'))
plot(T170[:,1],T170[:,3],label=('T=170K'))
plot(T150[:,1],T150[:,3],label=('T=150K'))
plot(T130[:,1],T130[:,3],label=('T=130K'))
plot(T110[:,1],T110[:,3],label=('T=110K'))
plot(T90[:,1],T90[:,3],label=('T=90K'))
plot(T70[:,1],T70[:,3],label=('T=70K'))
plot(T50[:,1],T50[:,3],label=('T=50K'))
xlabel('Current (mA)')
ylabel('SE (a.u.)')
title('L-I Characteristic (window)')
legend(loc='lower right', ncol=1, shadow=True, numpoints=1)
show()
#######################################################################
plot(T60[:,0],T60[:,1],label=('T=60K'))
plot(T80[:,0],T80[:,1],label=('T=80K'))
plot(T100[:,0],T100[:,1],label=('T=100K'))
plot(T120[:,0],T120[:,1],label=('T=120K'))
plot(T140[:,0],T140[:,1],label=('T=140K'))
plot(T160[:,0],T160[:,1],label=('T=160K'))
plot(T180[:,0],T180[:,1],label=('T=180K'))
plot(T200[:,0],T200[:,1],label=('T=200K'))
plot(T220[:,0],T220[:,1],label=('T=230K'))
plot(T240[:,0],T240[:,1],label=('T=240K'))
plot(T260[:,0],T260[:,1],label=('T=260K'))
plot(T280[:,0],T280[:,1],label=('T=280K'))
plot(T300[:,0],T300[:,1],label=('T=300K'))
plot(T320[:,0],T320[:,1],label=('T=330K'))
plot(T340[:,0],T340[:,1],label=('T=340K'))
plot(T350[:,0],T350[:,1],label=('T=350K'))
plot(T330[:,0],T330[:,1],label=('T=330K'))
plot(T310[:,0],T310[:,1],label=('T=310K'))
plot(T290[:,0],T290[:,1],label=('T=290K'))
plot(T270[:,0],T270[:,1],label=('T=270K'))
plot(T250[:,0],T250[:,1],label=('T=250K'))
plot(T230[:,0],T230[:,1],label=('T=230K'))
plot(T210[:,0],T210[:,1],label=('T=210K'))
plot(T190[:,0],T190[:,1],label=('T=190K'))
plot(T170[:,0],T170[:,1],label=('T=170K'))
plot(T150[:,0],T150[:,1],label=('T=150K'))
plot(T130[:,0],T130[:,1],label=('T=130K'))
plot(T110[:,0],T110[:,1],label=('T=110K'))
plot(T90[:,0],T90[:,1],label=('T=90K'))
plot(T70[:,0],T70[:,1],label=('T=70K'))
plot(T50[:,0],T50[:,1],label=('T=50K'))
xlabel('Voltage (V)')
ylabel('Current (mA)')
title('I-V Charactarstic')
legend(loc='lower right', ncol=1, shadow=True, numpoints=1)
xlim([0,2.5])
show()


#				SWITCH ON VOLTAGE
############################################################
############################################################
############################################################
############################################################
x=T50[:,0]
y=T50[:,1]
C1=[]
for n in range(2,len(x)-4):
	c1=((y[n+1]-2*y[n]+y[n-1])/(x[n]-x[n-1])*(x[n]-x[n-1]))
	C1.append(c1)
C1=array(C1)
x1=[]
for n in range(len(C1)):
	a=T50[n,0]
	x1.append(a)
x1=array(x1)
SO50=peakdet(C1,.0000000000000002,x1)
#print 'T=50 swith on voltage',(SO50)

############################################################
x=T60[:,0]
y=T60[:,1]
C1=[]
for n in range(2,len(x)-4):
	c1=((y[n+1]-2*y[n]+y[n-1])/(x[n]-x[n-1])*(x[n]-x[n-1]))
	C1.append(c1)
C1=array(C1)
x1=[]
for n in range(len(C1)):
	a=T60[n,0]
	x1.append(a)
x1=array(x1)
SO60=peakdet(C1,.0000000000000002,x1)
#print 'T=60 swith on voltage',(SO60)

############################################################
x=T70[:,0]
y=T70[:,1]
C1=[]
for n in range(2,len(x)-4):
	c1=((y[n+1]-2*y[n]+y[n-1])/(x[n]-x[n-1])*(x[n]-x[n-1]))
	C1.append(c1)
C1=array(C1)
x1=[]
for n in range(len(C1)):
	a=T70[n,0]
	x1.append(a)
x1=array(x1)
SO70=peakdet(C1,.0000000000000002,x1)
#print 'T=70 swith on voltage',(SO70)

############################################################
x=T80[:,0]
y=T80[:,1]
C1=[]
for n in range(2,len(x)-4):
	c1=((y[n+1]-2*y[n]+y[n-1])/(x[n]-x[n-1])*(x[n]-x[n-1]))
	C1.append(c1)
C1=array(C1)
x1=[]
for n in range(len(C1)):
	a=T80[n,0]
	x1.append(a)
x1=array(x1)
SO80=peakdet(C1,.0000000000000002,x1)
#print 'T=80 swith on voltage',(SO80)

############################################################
x=T90[:,0]
y=T90[:,1]
C1=[]
for n in range(2,len(x)-4):
	c1=((y[n+1]-2*y[n]+y[n-1])/(x[n]-x[n-1])*(x[n]-x[n-1]))
	C1.append(c1)
C1=array(C1)
x1=[]
for n in range(len(C1)):
	a=T90[n,0]
	x1.append(a)
x1=array(x1)
SO90=peakdet(C1,.0000000000000002,x1)
#print 'T=90 swith on voltage',(SO90)

############################################################
x=T100[:,0]
y=T100[:,1]
C1=[]
for n in range(2,len(x)-4):
	c1=((y[n+1]-2*y[n]+y[n-1])/(x[n]-x[n-1])*(x[n]-x[n-1]))
	C1.append(c1)
C1=array(C1)
x1=[]
for n in range(len(C1)):
	a=T100[n,0]
	x1.append(a)
x1=array(x1)
SO100=peakdet(C1,.0000000000000002,x1)
#print 'T=100 swith on voltage',(SO100)

############################################################
x=T110[:,0]
y=T110[:,1]
C1=[]
for n in range(2,len(x)-4):
	c1=((y[n+1]-2*y[n]+y[n-1])/(x[n]-x[n-1])*(x[n]-x[n-1]))
	C1.append(c1)
C1=array(C1)
x1=[]
for n in range(len(C1)):
	a=T110[n,0]
	x1.append(a)
x1=array(x1)
SO110=peakdet(C1,.0000000000000002,x1)
#print 'T=110 swith on voltage',(SO110)

############################################################
x=T120[:,0]
y=T120[:,1]
C1=[]
for n in range(2,len(x)-4):
	c1=((y[n+1]-2*y[n]+y[n-1])/(x[n]-x[n-1])*(x[n]-x[n-1]))
	C1.append(c1)
C1=array(C1)
x1=[]
for n in range(len(C1)):
	a=T120[n,0]
	x1.append(a)
x1=array(x1)
SO120=peakdet(C1,.0000000000000002,x1)
#print 'T=120 swith on voltage',(SO120)

############################################################
x=T130[:,0]
y=T130[:,1]
C1=[]
for n in range(2,len(x)-4):
	c1=((y[n+1]-2*y[n]+y[n-1])/(x[n]-x[n-1])*(x[n]-x[n-1]))
	C1.append(c1)
C1=array(C1)
x1=[]
for n in range(len(C1)):
	a=T130[n,0]
	x1.append(a)
x1=array(x1)
SO130=peakdet(C1,.0000000000000002,x1)
#print 'T=130 swith on voltage',(SO130)

############################################################
x=T140[:,0]
y=T140[:,1]
C1=[]
for n in range(2,len(x)-4):
	c1=((y[n+1]-2*y[n]+y[n-1])/(x[n]-x[n-1])*(x[n]-x[n-1]))
	C1.append(c1)
C1=array(C1)
x1=[]
for n in range(len(C1)):
	a=T140[n,0]
	x1.append(a)
x1=array(x1)
SO140=peakdet(C1,.0000000000000002,x1)
#print 'T=140 swith on voltage',(SO140)

############################################################
x=T150[:,0]
y=T150[:,1]
C1=[]
for n in range(2,len(x)-4):
	c1=((y[n+1]-2*y[n]+y[n-1])/(x[n]-x[n-1])*(x[n]-x[n-1]))
	C1.append(c1)
C1=array(C1)
x1=[]
for n in range(len(C1)):
	a=T150[n,0]
	x1.append(a)
x1=array(x1)
SO150=peakdet(C1,.0000000000000002,x1)
#print 'T=150 swith on voltage',(SO150)

############################################################
x=T160[:,0]
y=T160[:,1]
C1=[]
for n in range(2,len(x)-4):
	c1=((y[n+1]-2*y[n]+y[n-1])/(x[n]-x[n-1])*(x[n]-x[n-1]))
	C1.append(c1)
C1=array(C1)
x1=[]
for n in range(len(C1)):
	a=T160[n,0]
	x1.append(a)
x1=array(x1)
SO160=peakdet(C1,.0000000000000002,x1)
#print 'T=160 swith on voltage',(SO160)

############################################################
x=T170[:,0]
y=T170[:,1]
C1=[]
for n in range(2,len(x)-4):
	c1=((y[n+1]-2*y[n]+y[n-1])/(x[n]-x[n-1])*(x[n]-x[n-1]))
	C1.append(c1)
C1=array(C1)
x1=[]
for n in range(len(C1)):
	a=T170[n,0]
	x1.append(a)
x1=array(x1)
SO170=peakdet(C1,.0000000000000002,x1)
#print 'T=170 swith on voltage',(SO170)

############################################################
x=T180[:,0]
y=T180[:,1]
C1=[]
for n in range(2,len(x)-4):
	c1=((y[n+1]-2*y[n]+y[n-1])/(x[n]-x[n-1])*(x[n]-x[n-1]))
	C1.append(c1)
C1=array(C1)
x1=[]
for n in range(len(C1)):
	a=T180[n,0]
	x1.append(a)
x1=array(x1)
SO180=peakdet(C1,.0000000000000002,x1)
#print 'T=180 swith on voltage',(SO180)

############################################################
x=T190[:,0]
y=T190[:,1]
C1=[]
for n in range(2,len(x)-4):
	c1=((y[n+1]-2*y[n]+y[n-1])/(x[n]-x[n-1])*(x[n]-x[n-1]))
	C1.append(c1)
C1=array(C1)
x1=[]
for n in range(len(C1)):
	a=T190[n,0]
	x1.append(a)
x1=array(x1)
SO190=peakdet(C1,.0000000000000002,x1)
#print 'T=190 swith on voltage',(SO190)

############################################################
x=T200[:,0]
y=T200[:,1]
C1=[]
for n in range(2,len(x)-4):
	c1=((y[n+1]-2*y[n]+y[n-1])/(x[n]-x[n-1])*(x[n]-x[n-1]))
	C1.append(c1)
C1=array(C1)
x1=[]
for n in range(len(C1)):
	a=T200[n,0]
	x1.append(a)
x1=array(x1)
SO200=peakdet(C1,.0000000000000002,x1)
#print 'T=200 swith on voltage',(SO200)

############################################################
x=T210[:,0]
y=T210[:,1]
C1=[]
for n in range(2,len(x)-4):
	c1=((y[n+1]-2*y[n]+y[n-1])/(x[n]-x[n-1])*(x[n]-x[n-1]))
	C1.append(c1)
C1=array(C1)
x1=[]
for n in range(len(C1)):
	a=T210[n,0]
	x1.append(a)
x1=array(x1)
SO210=peakdet(C1,.0000000000000002,x1)
#print 'T=210 swith on voltage',(SO210)

############################################################
x=T220[:,0]
y=T220[:,1]
C1=[]
for n in range(2,len(x)-4):
	c1=((y[n+1]-2*y[n]+y[n-1])/(x[n]-x[n-1])*(x[n]-x[n-1]))
	C1.append(c1)
C1=array(C1)
x1=[]
for n in range(len(C1)):
	a=T220[n,0]
	x1.append(a)
x1=array(x1)
SO220=peakdet(C1,.0000000000000002,x1)
#print 'T=220 swith on voltage',(SO220)

############################################################
x=T230[:,0]
y=T230[:,1]
C1=[]
for n in range(2,len(x)-4):
	c1=((y[n+1]-2*y[n]+y[n-1])/(x[n]-x[n-1])*(x[n]-x[n-1]))
	C1.append(c1)
C1=array(C1)
x1=[]
for n in range(len(C1)):
	a=T230[n,0]
	x1.append(a)
x1=array(x1)
SO230=peakdet(C1,.0000000000000002,x1)
#print 'T=230 swith on voltage',(SO230)

############################################################
x=T240[:,0]
y=T240[:,1]
C1=[]
for n in range(2,len(x)-4):
	c1=((y[n+1]-2*y[n]+y[n-1])/(x[n]-x[n-1])*(x[n]-x[n-1]))
	C1.append(c1)
C1=array(C1)
x1=[]
for n in range(len(C1)):
	a=T240[n,0]
	x1.append(a)
x1=array(x1)
SO240=peakdet(C1,.0000000000000002,x1)
#print 'T=240 swith on voltage',(SO240)

############################################################
x=T250[:,0]
y=T250[:,1]
C1=[]
for n in range(2,len(x)-4):
	c1=((y[n+1]-2*y[n]+y[n-1])/(x[n]-x[n-1])*(x[n]-x[n-1]))
	C1.append(c1)
C1=array(C1)
x1=[]
for n in range(len(C1)):
	a=T250[n,0]
	x1.append(a)
x1=array(x1)
SO250=peakdet(C1,.0000000000000002,x1)
#print 'T=250 swith on voltage',(SO250)

############################################################
x=T260[:,0]
y=T260[:,1]
C1=[]
for n in range(2,len(x)-4):
	c1=((y[n+1]-2*y[n]+y[n-1])/(x[n]-x[n-1])*(x[n]-x[n-1]))
	C1.append(c1)
C1=array(C1)
x1=[]
for n in range(len(C1)):
	a=T260[n,0]
	x1.append(a)
x1=array(x1)
SO260=peakdet(C1,.0000000000000002,x1)
#print 'T=260 swith on voltage',(SO260)

############################################################
x=T270[:,0]
y=T270[:,1]
C1=[]
for n in range(2,len(x)-4):
	c1=((y[n+1]-2*y[n]+y[n-1])/(x[n]-x[n-1])*(x[n]-x[n-1]))
	C1.append(c1)
C1=array(C1)
x1=[]
for n in range(len(C1)):
	a=T270[n,0]
	x1.append(a)
x1=array(x1)
SO270=peakdet(C1,.0000000000000002,x1)
#print 'T=270 swith on voltage',(SO270)

############################################################
x=T280[:,0]
y=T280[:,1]
C1=[]
for n in range(2,len(x)-4):
	c1=((y[n+1]-2*y[n]+y[n-1])/(x[n]-x[n-1])*(x[n]-x[n-1]))
	C1.append(c1)
C1=array(C1)
x1=[]
for n in range(len(C1)):
	a=T280[n,0]
	x1.append(a)
x1=array(x1)
SO280=peakdet(C1,.0000000000000002,x1)
#print 'T=280 swith on voltage',(SO280)

############################################################
x=T290[:,0]
y=T290[:,1]
C1=[]
for n in range(2,len(x)-4):
	c1=((y[n+1]-2*y[n]+y[n-1])/(x[n]-x[n-1])*(x[n]-x[n-1]))
	C1.append(c1)
C1=array(C1)
x1=[]
for n in range(len(C1)):
	a=T290[n,0]
	x1.append(a)
x1=array(x1)
SO290=peakdet(C1,.0000000000000002,x1)
#print 'T=290 swith on voltage',(SO290)

############################################################
x=T300[:,0]
y=T300[:,1]
C1=[]
for n in range(2,len(x)-4):
	c1=((y[n+1]-2*y[n]+y[n-1])/(x[n]-x[n-1])*(x[n]-x[n-1]))
	C1.append(c1)
C1=array(C1)
x1=[]
for n in range(len(C1)):
	a=T300[n,0]
	x1.append(a)
x1=array(x1)
SO300=peakdet(C1,.0000000000000002,x1)
#print 'T=300 swith on voltage',(SO300)

############################################################
x=T310[:,0]
y=T310[:,1]
C1=[]
for n in range(2,len(x)-4):
	c1=((y[n+1]-2*y[n]+y[n-1])/(x[n]-x[n-1])*(x[n]-x[n-1]))
	C1.append(c1)
C1=array(C1)
x1=[]
for n in range(len(C1)):
	a=T310[n,0]
	x1.append(a)
x1=array(x1)
SO310=peakdet(C1,.0000000000000002,x1)
#print 'T=310 swith on voltage',(SO310)

############################################################
x=T320[:,0]
y=T320[:,1]
C1=[]
for n in range(2,len(x)-4):
	c1=((y[n+1]-2*y[n]+y[n-1])/(x[n]-x[n-1])*(x[n]-x[n-1]))
	C1.append(c1)
C1=array(C1)
x1=[]
for n in range(len(C1)):
	a=T320[n,0]
	x1.append(a)
x1=array(x1)
SO320=peakdet(C1,.0000000000000002,x1)
#print 'T=320 swith on voltage',(SO320)

############################################################
x=T330[:,0]
y=T330[:,1]
C1=[]
for n in range(2,len(x)-4):
	c1=((y[n+1]-2*y[n]+y[n-1])/(x[n]-x[n-1])*(x[n]-x[n-1]))
	C1.append(c1)
C1=array(C1)
x1=[]
for n in range(len(C1)):
	a=T330[n,0]
	x1.append(a)
x1=array(x1)
SO330=peakdet(C1,.0000000000000002,x1)
#print 'T=330 swith on voltage',(SO330)

############################################################
x=T340[:,0]
y=T340[:,1]
C1=[]
for n in range(2,len(x)-4):
	c1=((y[n+1]-2*y[n]+y[n-1])/(x[n]-x[n-1])*(x[n]-x[n-1]))
	C1.append(c1)
C1=array(C1)
x1=[]
for n in range(len(C1)):
	a=T340[n,0]
	x1.append(a)
x1=array(x1)
SO340=peakdet(C1,.0000000000000002,x1)
#print 'T=340 swith on voltage',(SO340)

############################################################
x=T350[:,0]
y=T350[:,1]
C1=[]
for n in range(2,len(x)-4):
	c1=((y[n+1]-2*y[n]+y[n-1])/(x[n]-x[n-1])*(x[n]-x[n-1]))
	C1.append(c1)
C1=array(C1)
x1=[]
for n in range(len(C1)):
	a=T350[n,0]
	x1.append(a)
x1=array(x1)
SO350=peakdet(C1,.0000000000000002,x1)
#print 'T=350 swith on voltage',(SO350)



#             THRESHOLD
########################################################################
########################################################################
########################################################################
coeff=calc_coeff(3,2,0)
x=T50[:,1]
y=T50[:,2]
x2=smooth(x,coeff)
y2=smooth(y,coeff)
#plot(x,y)
C=[]
for n in range(2,len(y)-4):
	c=((y2[n+1]-y2[n-1])/(x2[n+1]-x2[n-1]))
	C.append(c)
C=array(C)
C1=[]
for n in range(2,len(x)-4):
	c1=((y2[n+1]-2*y2[n]+y2[n-1])/(x2[n]-x2[n-1])*(x2[n]-x2[n-1]))
	C1.append(c1)
C1=array(C1)
x1=[]
for n in range(len(C)):
	a=T50[n,1]
	x1.append(a)
x1=array(x1)
x2=[]
for n in range(len(C1)):
	b=T50[n,1]
	x2.append(b)
x2=array(x2)
th50=peakdet(C1,50,x2)

d=x[2]-x[1]
thPoint=int(((th50[0])/d)+1)
########## Nd(ext)-fitting ##############
x_1=[]
y_1=[]
for n in range(thPoint+5,len(x)):
		b=x[n]
		c=y[n]
		x_1.append(b)
		y_1.append(c)
x_1=array(x_1)
y_1=array(y_1)
(ar,br)=polyfit(x_1,y_1,1)
xr=polyval([ar,br], x_1)
Z=diff(xr)
X=diff(x_1)
Ne50=(Z[0]/X[0])

mx=sqrt(max((y_1-xr)**2))
err=mx
Z1=diff(((xr[len(xr)-2]-err),(xr[1]+err)))
X1=diff((x_1[len(x_1)-2],x_1[1]))
err1=(Z1[0]/X1[0])
Nerr50=Ne50-err1
	
Irad50=T50[thPoint,3]

#show()


########################################################################
########################################################################


########################################################################
########################################################################
########################################################################
coeff=calc_coeff(4,2,0)
x=T60[:,1]
y=T60[:,2]
x2=smooth(x,coeff)
y2=smooth(y,coeff)
#plot(x,y)
C=[]
for n in range(2,len(y)-4):
	c=((y2[n+1]-y2[n-1])/(x2[n+1]-x2[n-1]))
	C.append(c)
C=array(C)
C1=[]
for n in range(2,len(x)-4):
	c1=((y2[n+1]-2*y2[n]+y2[n-1])/(x2[n]-x2[n-1])*(x2[n]-x2[n-1]))
	C1.append(c1)
C1=array(C1)
x1=[]
for n in range(len(C)):
	a=T60[n,1]
	x1.append(a)
x1=array(x1)
x2=[]
for n in range(len(C1)):
	b=T60[n,1]
	x2.append(b)
x2=array(x2)
th60=peakdet(C1,50,x2)

d=x[2]-x[1]
thPoint=int((th60[0]/d)+1)
########## Nd(ext)-fitting ##############
x_1=[]
y_1=[]
for n in range(thPoint+5,len(x)):
		b=x[n]
		c=y[n]
		x_1.append(b)
		y_1.append(c)
x_1=array(x_1)
y_1=array(y_1)
(ar,br)=polyfit(x_1,y_1,1)
xr=polyval([ar,br], x_1)
Z=diff(xr)
X=diff(x_1)
Ne60=(Z[0]/X[0])


mx=sqrt(max((y_1-xr)**2))
err=mx
Z1=diff(((xr[len(xr)-2]-err),(xr[1]+err)))
X1=diff((x_1[len(x_1)-2],x_1[1]))
err1=(Z1[0]/X1[0])
Nerr60=Ne60-err1

Irad60=T60[thPoint,3]

#show()
########################################################################
########################################################################

########################################################################
########################################################################
########################################################################
x=T70[:,1]
y=T70[:,2]
x2=smooth(x,coeff)
y2=smooth(y,coeff)
#plot(x,y)
C=[]
for n in range(2,len(y)-4):
	c=((y2[n+1]-y2[n-1])/(x2[n+1]-x2[n-1]))
	C.append(c)
C=array(C)
C1=[]
for n in range(2,len(x)-4):
	c1=((y2[n+1]-2*y2[n]+y2[n-1])/(x2[n]-x2[n-1])*(x2[n]-x2[n-1]))
	C1.append(c1)
C1=array(C1)

x1=[]
for n in range(len(C)):
	a=T70[n,1]
	x1.append(a)
x1=array(x1)
x2=[]
for n in range(len(C1)):
	b=T70[n,1]
	x2.append(b)
x2=array(x2)
th70=peakdet(C1,50,x2)

d=x[2]-x[1]
thPoint=int((th70[0]/d)+1)
########## Nd(ext)-fitting ##############
x_1=[]
y_1=[]
for n in range(thPoint+5,len(x)):
		b=x[n]
		c=y[n]
		x_1.append(b)
		y_1.append(c)
x_1=array(x_1)
y_1=array(y_1)
(ar,br)=polyfit(x_1,y_1,1)
xr=polyval([ar,br], x_1)
Z=diff(xr)
X=diff(x_1)
Ne70=(Z[0]/X[0])

mx=sqrt(max((y_1-xr)**2))
err=mx
Z1=diff(((xr[len(xr)-2]-err),(xr[1]+err)))
X1=diff((x_1[len(x_1)-2],x_1[1]))
err1=(Z1[0]/X1[0])
Nerr70=Ne70-err1

Irad70=T70[thPoint,3]
	
#show()
########################################################################
########################################################################

########################################################################
########################################################################
########################################################################
x=T80[:,1]
y=T80[:,2]
x2=smooth(x,coeff)
y2=smooth(y,coeff)
#plot(x,y)
C=[]
for n in range(2,len(y)-4):
	c=((y2[n+1]-y2[n-1])/(x2[n+1]-x2[n-1]))
	C.append(c)
C=array(C)
C1=[]
for n in range(2,len(x)-4):
	c1=((y2[n+1]-2*y2[n]+y2[n-1])/(x2[n]-x2[n-1])*(x2[n]-x2[n-1]))
	C1.append(c1)
C1=array(C1)

x1=[]
for n in range(len(C)):
	a=T80[n,1]
	x1.append(a)
x1=array(x1)
x2=[]
for n in range(len(C1)):
	b=T80[n,1]
	x2.append(b)
x2=array(x2)
th80=peakdet(C1,50,x2)

d=x[2]-x[1]
thPoint=int((th80[0]/d)+1)
########## Nd(ext)-fitting ##############
x_1=[]
y_1=[]
for n in range(thPoint+5,len(x)):
		b=x[n]
		c=y[n]
		x_1.append(b)
		y_1.append(c)
x_1=array(x_1)
y_1=array(y_1)
(ar,br)=polyfit(x_1,y_1,1)
xr=polyval([ar,br], x_1)
Z=diff(xr)
X=diff(x_1)
Ne80=(Z[0]/X[0])

mx=sqrt(max((y_1-xr)**2))
err=mx
Z1=diff(((xr[len(xr)-2]-err),(xr[1]+err)))
X1=diff((x_1[len(x_1)-2],x_1[1]))
err1=(Z1[0]/X1[0])
Nerr80=Ne80-err1

Irad80=T80[thPoint,3]
	
#show()
########################################################################
########################################################################

########################################################################
########################################################################
########################################################################
x=T90[:,1]
y=T90[:,2]
x2=smooth(x,coeff)
y2=smooth(y,coeff)
#plot(x,y)
C=[]
for n in range(2,len(y)-4):
	c=((y2[n+1]-y2[n-1])/(x2[n+1]-x2[n-1]))
	C.append(c)
C=array(C)
C1=[]
for n in range(2,len(x)-4):
	c1=((y2[n+1]-2*y2[n]+y2[n-1])/(x2[n]-x2[n-1])*(x2[n]-x2[n-1]))
	C1.append(c1)
C1=array(C1)

x1=[]
for n in range(len(C)):
	a=T90[n,1]
	x1.append(a)
x1=array(x1)
x2=[]
for n in range(len(C1)):
	b=T90[n,1]
	x2.append(b)
x2=array(x2)
th90=peakdet(C1,50,x2)

d=x[2]-x[1]
thPoint=int((th90[0]/d)+1)
########## Nd(ext)-fitting ##############
x_1=[]
y_1=[]
for n in range(thPoint+5,len(x)):
		b=x[n]
		c=y[n]
		x_1.append(b)
		y_1.append(c)
x_1=array(x_1)
y_1=array(y_1)
(ar,br)=polyfit(x_1,y_1,1)
xr=polyval([ar,br], x_1)
Z=diff(xr)
X=diff(x_1)
Ne90=(Z[0]/X[0])

mx=sqrt(max((y_1-xr)**2))
err=mx
Z1=diff(((xr[len(xr)-2]-err),(xr[1]+err)))
X1=diff((x_1[len(x_1)-2],x_1[1]))
err1=(Z1[0]/X1[0])
Nerr90=Ne90-err1

Irad90=T90[thPoint,3]
	
#show()
########################################################################
########################################################################

########################################################################
########################################################################
########################################################################
x=T100[:,1]
y=T100[:,2]
x2=smooth(x,coeff)
y2=smooth(y,coeff)
#plot(x,y)
C=[]
for n in range(2,len(y)-4):
	c=((y2[n+1]-y2[n-1])/(x2[n+1]-x2[n-1]))
	C.append(c)
C=array(C)
C1=[]
for n in range(2,len(x)-4):
	c1=((y2[n+1]-2*y2[n]+y2[n-1])/(x2[n]-x2[n-1])*(x2[n]-x2[n-1]))
	C1.append(c1)
C1=array(C1)

x1=[]
for n in range(len(C)):
	a=T100[n,1]
	x1.append(a)
x1=array(x1)
x2=[]
for n in range(len(C1)):
	b=T100[n,1]
	x2.append(b)
x2=array(x2)
th100=peakdet(C1,40,x2)

d=x[2]-x[1]
thPoint=int((th100[0]/d)+1)
########## Nd(ext)-fitting ##############
x_1=[]
y_1=[]
for n in range(thPoint+5,len(x)):
		b=x[n]
		c=y[n]
		x_1.append(b)
		y_1.append(c)
x_1=array(x_1)
y_1=array(y_1)
(ar,br)=polyfit(x_1,y_1,1)
xr=polyval([ar,br], x_1)
Z=diff(xr)
X=diff(x_1)
Ne100=(Z[0]/X[0])

mx=sqrt(max((y_1-xr)**2))
err=mx
Z1=diff(((xr[len(xr)-2]-err),(xr[1]+err)))
X1=diff((x_1[len(x_1)-2],x_1[1]))
err1=(Z1[0]/X1[0])
Nerr100=Ne100-err1

Irad100=T100[thPoint,3]

#show()
########################################################################
########################################################################

########################################################################
########################################################################
########################################################################
x=T110[:,1]
y=T110[:,2]
x2=smooth(x,coeff)
y2=smooth(y,coeff)
#plot(x,y)
C=[]
for n in range(2,len(y)-4):
	c=((y2[n+1]-y2[n-1])/(x2[n+1]-x2[n-1]))
	C.append(c)
C=array(C)
C1=[]
for n in range(2,len(x)-4):
	c1=((y2[n+1]-2*y2[n]+y2[n-1])/(x2[n]-x2[n-1])*(x2[n]-x2[n-1]))
	C1.append(c1)
C1=array(C1)

x1=[]
for n in range(len(C)):
	a=T110[n,1]
	x1.append(a)
x1=array(x1)
x2=[]
for n in range(len(C1)):
	b=T110[n,1]
	x2.append(b)
x2=array(x2)
th110=peakdet(C1,50,x2)

d=x[2]-x[1]
thPoint=int((th110[0]/d)+1)
########## Nd(ext)-fitting ##############
x_1=[]
y_1=[]
for n in range(thPoint+5,len(x)):
		b=x[n]
		c=y[n]
		x_1.append(b)
		y_1.append(c)
x_1=array(x_1)
y_1=array(y_1)
(ar,br)=polyfit(x_1,y_1,1)
xr=polyval([ar,br], x_1)
Z=diff(xr)
X=diff(x_1)
Ne110=(Z[0]/X[0])

mx=sqrt(max((y_1-xr)**2))
err=mx
Z1=diff(((xr[len(xr)-2]-err),(xr[1]+err)))
X1=diff((x_1[len(x_1)-2],x_1[1]))
err1=(Z1[0]/X1[0])
Nerr110=Ne110-err1

Irad110=T110[thPoint,3]
	
#show()
########################################################################
########################################################################

########################################################################
########################################################################
########################################################################
x=T120[:,1]
y=T120[:,2]
x2=smooth(x,coeff)
y2=smooth(y,coeff)
#plot(x,y)
C=[]
for n in range(2,len(y)-4):
	c=((y2[n+1]-y2[n-1])/(x2[n+1]-x2[n-1]))
	C.append(c)
C=array(C)
C1=[]
for n in range(2,len(x)-4):
	c1=((y2[n+1]-2*y2[n]+y2[n-1])/(x2[n]-x2[n-1])*(x2[n]-x2[n-1]))
	C1.append(c1)
C1=array(C1)

x1=[]
for n in range(len(C)):
	a=T120[n,1]
	x1.append(a)
x1=array(x1)
x2=[]
for n in range(len(C1)):
	b=T120[n,1]
	x2.append(b)
x2=array(x2)
th120=peakdet(C1,50,x2)

d=x[2]-x[1]
thPoint=int((th120[0]/d)+1)
########## Nd(ext)-fitting ##############
x_1=[]
y_1=[]
for n in range(thPoint+5,len(x)):
		b=x[n]
		c=y[n]
		x_1.append(b)
		y_1.append(c)
x_1=array(x_1)
y_1=array(y_1)
(ar,br)=polyfit(x_1,y_1,1)
xr=polyval([ar,br], x_1)
Z=diff(xr)
X=diff(x_1)
Ne120=(Z[0]/X[0])

mx=sqrt(max((y_1-xr)**2))
err=mx
Z1=diff(((xr[len(xr)-2]-err),(xr[1]+err)))
X1=diff((x_1[len(x_1)-2],x_1[1]))
err1=(Z1[0]/X1[0])
Nerr120=Ne120-err1

Irad120=T120[thPoint,3]
	
#show()
########################################################################
########################################################################

########################################################################
########################################################################
########################################################################
x=T130[:,1]
y=T130[:,2]
x2=smooth(x,coeff)
y2=smooth(y,coeff)
#plot(x,y)
C=[]
for n in range(2,len(y)-4):
	c=((y2[n+1]-y2[n-1])/(x2[n+1]-x2[n-1]))
	C.append(c)
C=array(C)
C1=[]
for n in range(2,len(x)-4):
	c1=((y2[n+1]-2*y2[n]+y2[n-1])/(x2[n]-x2[n-1])*(x2[n]-x2[n-1]))
	C1.append(c1)
C1=array(C1)

x1=[]
for n in range(len(C)):
	a=T130[n,1]
	x1.append(a)
x1=array(x1)
x2=[]
for n in range(len(C1)):
	b=T130[n,1]
	x2.append(b)
x2=array(x2)
th130=peakdet(C1,50,x2)

d=x[2]-x[1]
thPoint=int((th130[0]/d)+1)
########## Nd(ext)-fitting ##############
x_1=[]
y_1=[]
for n in range(thPoint+5,len(x)):
		b=x[n]
		c=y[n]
		x_1.append(b)
		y_1.append(c)
x_1=array(x_1)
y_1=array(y_1)
(ar,br)=polyfit(x_1,y_1,1)
xr=polyval([ar,br], x_1)
Z=diff(xr)
X=diff(x_1)
Ne130=(Z[0]/X[0])

mx=sqrt(max((y_1-xr)**2))
err=mx
Z1=diff(((xr[len(xr)-2]-err),(xr[1]+err)))
X1=diff((x_1[len(x_1)-2],x_1[1]))
err1=(Z1[0]/X1[0])
Nerr130=Ne130-err1

Irad130=T130[thPoint,3]

#show()
########################################################################
########################################################################

########################################################################
########################################################################
########################################################################
x=T140[:,1]
y=T140[:,2]
x2=smooth(x,coeff)
y2=smooth(y,coeff)
#plot(x,y)
C=[]
for n in range(2,len(y)-4):
	c=((y2[n+1]-y2[n-1])/(x2[n+1]-x2[n-1]))
	C.append(c)
C=array(C)
C1=[]
for n in range(2,len(x)-4):
	c1=((y2[n+1]-2*y2[n]+y2[n-1])/(x2[n]-x2[n-1])*(x2[n]-x2[n-1]))
	C1.append(c1)
C1=array(C1)

x1=[]
for n in range(len(C)):
	a=T140[n,1]
	x1.append(a)
x1=array(x1)
x2=[]
for n in range(len(C1)):
	b=T140[n,1]
	x2.append(b)
x2=array(x2)
th140=peakdet(C1,50,x2)

d=x[2]-x[1]
thPoint=int((th140[0]/d)+1)
########## Nd(ext)-fitting ##############
x_1=[]
y_1=[]
for n in range(thPoint+5,len(x)):
		b=x[n]
		c=y[n]
		x_1.append(b)
		y_1.append(c)
x_1=array(x_1)
y_1=array(y_1)
(ar,br)=polyfit(x_1,y_1,1)
xr=polyval([ar,br], x_1)
Z=diff(xr)
X=diff(x_1)
Ne140=(Z[0]/X[0])

mx=sqrt(max((y_1-xr)**2))
err=mx
Z1=diff(((xr[len(xr)-2]-err),(xr[1]+err)))
X1=diff((x_1[len(x_1)-2],x_1[1]))
err1=(Z1[0]/X1[0])
Nerr140=Ne140-err1

Irad140=T140[thPoint,3]
	
#show()
########################################################################
########################################################################

########################################################################
########################################################################
########################################################################
x=T150[:,1]
y=T150[:,2]
x2=smooth(x,coeff)
y2=smooth(y,coeff)
#plot(x,y)
C=[]
for n in range(2,len(y)-4):
	c=((y2[n+1]-y2[n-1])/(x2[n+1]-x2[n-1]))
	C.append(c)
C=array(C)
C1=[]
for n in range(2,len(x)-4):
	c1=((y2[n+1]-2*y2[n]+y2[n-1])/(x2[n]-x2[n-1])*(x2[n]-x2[n-1]))
	C1.append(c1)
C1=array(C1)

x1=[]
for n in range(len(C)):
	a=T150[n,1]
	x1.append(a)
x1=array(x1)
x2=[]
for n in range(len(C1)):
	b=T150[n,1]
	x2.append(b)
x2=array(x2)
th150=peakdet(C1,50,x2)

d=x[2]-x[1]
thPoint=int((th150[0]/d)+1)
########## Nd(ext)-fitting ##############
x_1=[]
y_1=[]
for n in range(thPoint+5,len(x)):
		b=x[n]
		c=y[n]
		x_1.append(b)
		y_1.append(c)
x_1=array(x_1)
y_1=array(y_1)
(ar,br)=polyfit(x_1,y_1,1)
xr=polyval([ar,br], x_1)
Z=diff(xr)
X=diff(x_1)
Ne150=(Z[0]/X[0])

mx=sqrt(max((y_1-xr)**2))
err=mx
Z1=diff(((xr[len(xr)-2]-err),(xr[1]+err)))
X1=diff((x_1[len(x_1)-2],x_1[1]))
err1=(Z1[0]/X1[0])
Nerr150=Ne150-err1

Irad150=T150[thPoint,3]
	
#show()
########################################################################
########################################################################

########################################################################
########################################################################
########################################################################
x=T160[:,1]
y=T160[:,2]
x2=smooth(x,coeff)
y2=smooth(y,coeff)
#plot(x,y)
C=[]
for n in range(2,len(y)-4):
	c=((y2[n+1]-y2[n-1])/(x2[n+1]-x2[n-1]))
	C.append(c)
C=array(C)
C1=[]
for n in range(2,len(x)-4):
	c1=((y2[n+1]-2*y2[n]+y2[n-1])/(x2[n]-x2[n-1])*(x2[n]-x2[n-1]))
	C1.append(c1)
C1=array(C1)

x1=[]
for n in range(len(C)):
	a=T160[n,1]
	x1.append(a)
x1=array(x1)
x2=[]
for n in range(len(C1)):
	b=T160[n,1]
	x2.append(b)
x2=array(x2)
th160=peakdet(C1,50,x2)

d=x[2]-x[1]
thPoint=int((th160[0]/d)+1)
########## Nd(ext)-fitting ##############
x_1=[]
y_1=[]
for n in range(thPoint+5,len(x)):
		b=x[n]
		c=y[n]
		x_1.append(b)
		y_1.append(c)
x_1=array(x_1)
y_1=array(y_1)
(ar,br)=polyfit(x_1,y_1,1)
xr=polyval([ar,br], x_1)
Z=diff(xr)
X=diff(x_1)
Ne160=(Z[0]/X[0])

mx=sqrt(max((y_1-xr)**2))
err=mx
Z1=diff(((xr[len(xr)-2]-err),(xr[1]+err)))
X1=diff((x_1[len(x_1)-2],x_1[1]))
err1=(Z1[0]/X1[0])
Nerr160=Ne160-err1

Irad160=T160[thPoint,3]
	
#show()
########################################################################
########################################################################

########################################################################
########################################################################
########################################################################
x=T170[:,1]
y=T170[:,2]
x2=smooth(x,coeff)
y2=smooth(y,coeff)
#plot(x,y)
C=[]
for n in range(2,len(y)-4):
	c=((y2[n+1]-y2[n-1])/(x2[n+1]-x2[n-1]))
	C.append(c)
C=array(C)
C1=[]
for n in range(2,len(x)-4):
	c1=((y2[n+1]-2*y2[n]+y2[n-1])/(x2[n]-x2[n-1])*(x2[n]-x2[n-1]))
	C1.append(c1)
C1=array(C1)

x1=[]
for n in range(len(C)):
	a=T170[n,1]
	x1.append(a)
x1=array(x1)
x2=[]
for n in range(len(C1)):
	b=T170[n,1]
	x2.append(b)
x2=array(x2)
th170=peakdet(C1,50,x2)

d=x[2]-x[1]
thPoint=int((th170[0]/d)+1)
########## Nd(ext)-fitting ##############
x_1=[]
y_1=[]
for n in range(thPoint+5,len(x)):
		b=x[n]
		c=y[n]
		x_1.append(b)
		y_1.append(c)
x_1=array(x_1)
y_1=array(y_1)
(ar,br)=polyfit(x_1,y_1,1)
xr=polyval([ar,br], x_1)
Z=diff(xr)
X=diff(x_1)
Ne170=(Z[0]/X[0])

mx=sqrt(max((y_1-xr)**2))
err=mx
Z1=diff(((xr[len(xr)-2]-err),(xr[1]+err)))
X1=diff((x_1[len(x_1)-2],x_1[1]))
err1=(Z1[0]/X1[0])
Nerr170=Ne170-err1
	
Irad170=T170[thPoint,3]

#show()
########################################################################
########################################################################

########################################################################
########################################################################
########################################################################
x=T180[:,1]
y=T180[:,2]
x2=smooth(x,coeff)
y2=smooth(y,coeff)
#plot(x,y)
C=[]
for n in range(2,len(y)-4):
	c=((y2[n+1]-y2[n-1])/(x2[n+1]-x2[n-1]))
	C.append(c)
C=array(C)
C1=[]
for n in range(2,len(x)-4):
	c1=((y2[n+1]-2*y2[n]+y2[n-1])/(x2[n]-x2[n-1])*(x2[n]-x2[n-1]))
	C1.append(c1)
C1=array(C1)

x1=[]
for n in range(len(C)):
	a=T180[n,1]
	x1.append(a)
x1=array(x1)
x2=[]
for n in range(len(C1)):
	b=T180[n,1]
	x2.append(b)
x2=array(x2)
th180=peakdet(C1,50,x2)

d=x[2]-x[1]
thPoint=int((th180[0]/d)+1)
########## Nd(ext)-fitting ##############
x_1=[]
y_1=[]
for n in range(thPoint+5,len(x)):
		b=x[n]
		c=y[n]
		x_1.append(b)
		y_1.append(c)
x_1=array(x_1)
y_1=array(y_1)
(ar,br)=polyfit(x_1,y_1,1)
xr=polyval([ar,br], x_1)
Z=diff(xr)
X=diff(x_1)
Ne180=(Z[0]/X[0])

mx=sqrt(max((y_1-xr)**2))
err=mx
Z1=diff(((xr[len(xr)-2]-err),(xr[1]+err)))
X1=diff((x_1[len(x_1)-2],x_1[1]))
err1=(Z1[0]/X1[0])
Nerr180=Ne180-err1

Irad180=T180[thPoint,3]

#show()
########################################################################
########################################################################

########################################################################
########################################################################
########################################################################
x=T190[:,1]
y=T190[:,2]
x2=smooth(x,coeff)
y2=smooth(y,coeff)
#plot(x,y)
C=[]
for n in range(2,len(y)-4):
	c=((y2[n+1]-y2[n-1])/(x2[n+1]-x2[n-1]))
	C.append(c)
C=array(C)
C1=[]
for n in range(2,len(x)-4):
	c1=((y2[n+1]-2*y2[n]+y2[n-1])/(x2[n]-x2[n-1])*(x2[n]-x2[n-1]))
	C1.append(c1)
C1=array(C1)

x1=[]
for n in range(len(C)):
	a=T190[n,1]
	x1.append(a)
x1=array(x1)
x2=[]
for n in range(len(C1)):
	b=T190[n,1]
	x2.append(b)
x2=array(x2)
th190=peakdet(C1,50,x2)

d=x[2]-x[1]
thPoint=int((th190[0]/d)+1)
########## Nd(ext)-fitting ##############
x_1=[]
y_1=[]
for n in range(thPoint+5,len(x)):
		b=x[n]
		c=y[n]
		x_1.append(b)
		y_1.append(c)
x_1=array(x_1)
y_1=array(y_1)
(ar,br)=polyfit(x_1,y_1,1)
xr=polyval([ar,br], x_1)
Z=diff(xr)
X=diff(x_1)
Ne190=(Z[0]/X[0])

mx=sqrt(max((y_1-xr)**2))
err=mx
Z1=diff(((xr[len(xr)-2]-err),(xr[1]+err)))
X1=diff((x_1[len(x_1)-2],x_1[1]))
err1=(Z1[0]/X1[0])
Nerr190=Ne190-err1

Irad190=T190[thPoint,3]

#show()
########################################################################
########################################################################

########################################################################
########################################################################
########################################################################
x=T200[:,1]
y=T200[:,2]
x2=smooth(x,coeff)
y2=smooth(y,coeff)
#plot(x,y)
C=[]
for n in range(2,len(y)-4):
	c=((y2[n+1]-y2[n-1])/(x2[n+1]-x2[n-1]))
	C.append(c)
C=array(C)
C1=[]
for n in range(2,len(x)-4):
	c1=((y2[n+1]-2*y2[n]+y2[n-1])/(x2[n]-x2[n-1])*(x2[n]-x2[n-1]))
	C1.append(c1)
C1=array(C1)

x1=[]
for n in range(len(C)):
	a=T200[n,1]
	x1.append(a)
x1=array(x1)
x2=[]
for n in range(len(C1)):
	b=T200[n,1]
	x2.append(b)
x2=array(x2)
th200=peakdet(C1,50,x2)

d=x[2]-x[1]
thPoint=int((th200[0]/d)+1)
########## Nd(ext)-fitting ##############
x_1=[]
y_1=[]
for n in range(thPoint+5,len(x)):
		b=x[n]
		c=y[n]
		x_1.append(b)
		y_1.append(c)
x_1=array(x_1)
y_1=array(y_1)
(ar,br)=polyfit(x_1,y_1,1)
xr=polyval([ar,br], x_1)
Z=diff(xr)
X=diff(x_1)
Ne200=(Z[0]/X[0])

mx=sqrt(max((y_1-xr)**2))
err=mx
Z1=diff(((xr[len(xr)-2]-err),(xr[1]+err)))
X1=diff((x_1[len(x_1)-2],x_1[1]))
err1=(Z1[0]/X1[0])
Nerr200=Ne200-err1

Irad200=T200[thPoint,3]
	
#show()
########################################################################
########################################################################

########################################################################
########################################################################
########################################################################
x=T210[:,1]
y=T210[:,2]
x2=smooth(x,coeff)
y2=smooth(y,coeff)
#plot(x,y)
C=[]
for n in range(2,len(y)-4):
	c=((y2[n+1]-y2[n-1])/(x2[n+1]-x2[n-1]))
	C.append(c)
C=array(C)
C1=[]
for n in range(2,len(x)-4):
	c1=((y2[n+1]-2*y2[n]+y2[n-1])/(x2[n]-x2[n-1])*(x2[n]-x2[n-1]))
	C1.append(c1)
C1=array(C1)

x1=[]
for n in range(len(C)):
	a=T210[n,1]
	x1.append(a)
x1=array(x1)
x2=[]
for n in range(len(C1)):
	b=T210[n,1]
	x2.append(b)
x2=array(x2)
th210=peakdet(C1,50,x2)

d=x[2]-x[1]
thPoint=int((th210[0]/d)+1)
########## Nd(ext)-fitting ##############
x_1=[]
y_1=[]
for n in range(thPoint+5,len(x)):
		b=x[n]
		c=y[n]
		x_1.append(b)
		y_1.append(c)
x_1=array(x_1)
y_1=array(y_1)
(ar,br)=polyfit(x_1,y_1,1)
xr=polyval([ar,br], x_1)
Z=diff(xr)
X=diff(x_1)
Ne210=(Z[0]/X[0])

mx=sqrt(max((y_1-xr)**2))
err=mx
Z1=diff(((xr[len(xr)-2]-err),(xr[1]+err)))
X1=diff((x_1[len(x_1)-2],x_1[1]))
err1=(Z1[0]/X1[0])
Nerr210=Ne210-err1

Irad210=T210[thPoint,3]

#show()
########################################################################
########################################################################

########################################################################
########################################################################
########################################################################
x=T220[:,1]
y=T220[:,2]
x2=smooth(x,coeff)
y2=smooth(y,coeff)
#plot(x,y)
C=[]
for n in range(2,len(y)-4):
	c=((y2[n+1]-y2[n-1])/(x2[n+1]-x2[n-1]))
	C.append(c)
C=array(C)
C1=[]
for n in range(2,len(x)-4):
	c1=((y2[n+1]-2*y2[n]+y2[n-1])/(x2[n]-x2[n-1])*(x2[n]-x2[n-1]))
	C1.append(c1)
C1=array(C1)

x1=[]
for n in range(len(C)):
	a=T220[n,1]
	x1.append(a)
x1=array(x1)
x2=[]
for n in range(len(C1)):
	b=T220[n,1]
	x2.append(b)
x2=array(x2)
th220=peakdet(C1,50,x2)

d=x[2]-x[1]
thPoint=int((th220[0]/d)+1)
########## Nd(ext)-fitting ##############
x_1=[]
y_1=[]
for n in range(thPoint+5,len(x)):
		b=x[n]
		c=y[n]
		x_1.append(b)
		y_1.append(c)
x_1=array(x_1)
y_1=array(y_1)
(ar,br)=polyfit(x_1,y_1,1)
xr=polyval([ar,br], x_1)
Z=diff(xr)
X=diff(x_1)
Ne220=(Z[0]/X[0])

mx=sqrt(max((y_1-xr)**2))
err=mx
Z1=diff(((xr[len(xr)-2]-err),(xr[1]+err)))
X1=diff((x_1[len(x_1)-2],x_1[1]))
err1=(Z1[0]/X1[0])
Nerr220=Ne220-err1

Irad220=T220[thPoint,3]

#show()
########################################################################
########################################################################

########################################################################
########################################################################
########################################################################
x=T230[:,1]
y=T230[:,2]
x2=smooth(x,coeff)
y2=smooth(y,coeff)
#plot(x,y)
C=[]
for n in range(2,len(y)-4):
	c=((y2[n+1]-y2[n-1])/(x2[n+1]-x2[n-1]))
	C.append(c)
C=array(C)
C1=[]
for n in range(2,len(x)-4):
	c1=((y2[n+1]-2*y2[n]+y2[n-1])/(x2[n]-x2[n-1])*(x2[n]-x2[n-1]))
	C1.append(c1)
C1=array(C1)

x1=[]
for n in range(len(C)):
	a=T230[n,1]
	x1.append(a)
x1=array(x1)
x2=[]
for n in range(len(C1)):
	b=T230[n,1]
	x2.append(b)
x2=array(x2)
th230=peakdet(C1,50,x2)

d=x[2]-x[1]
thPoint=int((th230[0]/d)+1)
########## Nd(ext)-fitting ##############
x_1=[]
y_1=[]
for n in range(thPoint+5,len(x)):
		b=x[n]
		c=y[n]
		x_1.append(b)
		y_1.append(c)
x_1=array(x_1)
y_1=array(y_1)
(ar,br)=polyfit(x_1,y_1,1)
xr=polyval([ar,br], x_1)
Z=diff(xr)
X=diff(x_1)
Ne230=(Z[0]/X[0])

mx=sqrt(max((y_1-xr)**2))
err=mx
Z1=diff(((xr[len(xr)-2]-err),(xr[1]+err)))
X1=diff((x_1[len(x_1)-2],x_1[1]))
err1=(Z1[0]/X1[0])
Nerr230=Ne230-err1

Irad230=T230[thPoint,3]

#show()
########################################################################
########################################################################

########################################################################
########################################################################
########################################################################
x=T240[:,1]
y=T240[:,2]
x2=smooth(x,coeff)
y2=smooth(y,coeff)
#plot(x,y)
C=[]
for n in range(2,len(y)-4):
	c=((y2[n+1]-y2[n-1])/(x2[n+1]-x2[n-1]))
	C.append(c)
C=array(C)
C1=[]
for n in range(2,len(x)-4):
	c1=((y2[n+1]-2*y2[n]+y2[n-1])/(x2[n]-x2[n-1])*(x2[n]-x2[n-1]))
	C1.append(c1)
C1=array(C1)

x1=[]
for n in range(len(C)):
	a=T240[n,1]
	x1.append(a)
x1=array(x1)
x2=[]
for n in range(len(C1)):
	b=T240[n,1]
	x2.append(b)
x2=array(x2)
th240=peakdet(C1,50,x2)

d=x[2]-x[1]
thPoint=int((th240[0]/d)+1)
########## Nd(ext)-fitting ##############
x_1=[]
y_1=[]
for n in range(thPoint+5,len(x)):
		b=x[n]
		c=y[n]
		x_1.append(b)
		y_1.append(c)
x_1=array(x_1)
y_1=array(y_1)
(ar,br)=polyfit(x_1,y_1,1)
xr=polyval([ar,br], x_1)
Z=diff(xr)
X=diff(x_1)
Ne240=(Z[0]/X[0])

mx=sqrt(max((y_1-xr)**2))
err=mx
Z1=diff(((xr[len(xr)-2]-err),(xr[1]+err)))
X1=diff((x_1[len(x_1)-2],x_1[1]))
err1=(Z1[0]/X1[0])
Nerr240=Ne240-err1

Irad240=T240[thPoint,3]

#show()
########################################################################
########################################################################

########################################################################
########################################################################
########################################################################
x=T250[:,1]
y=T250[:,2]
x2=smooth(x,coeff)
y2=smooth(y,coeff)
#plot(x,y)
C=[]
for n in range(2,len(y)-4):
	c=((y2[n+1]-y2[n-1])/(x2[n+1]-x2[n-1]))
	C.append(c)
C=array(C)
C1=[]
for n in range(2,len(x)-4):
	c1=((y2[n+1]-2*y2[n]+y2[n-1])/(x2[n]-x2[n-1])*(x2[n]-x2[n-1]))
	C1.append(c1)
C1=array(C1)

x1=[]
for n in range(len(C)):
	a=T250[n,1]
	x1.append(a)
x1=array(x1)
x2=[]
for n in range(len(C1)):
	b=T250[n,1]
	x2.append(b)
x2=array(x2)
th250=peakdet(C1,50,x2)

d=x[2]-x[1]
thPoint=int((th250[0]/d)+1)
########## Nd(ext)-fitting ##############
x_1=[]
y_1=[]
for n in range(thPoint+5,len(x)):
		b=x[n]
		c=y[n]
		x_1.append(b)
		y_1.append(c)
x_1=array(x_1)
y_1=array(y_1)
(ar,br)=polyfit(x_1,y_1,1)
xr=polyval([ar,br], x_1)
Z=diff(xr)
X=diff(x_1)
Ne250=(Z[0]/X[0])

mx=sqrt(max((y_1-xr)**2))
err=mx
Z1=diff(((xr[len(xr)-2]-err),(xr[1]+err)))
X1=diff((x_1[len(x_1)-2],x_1[1]))
err1=(Z1[0]/X1[0])
Nerr250=Ne250-err1

Irad250=T250[thPoint,3]

#show()
########################################################################
########################################################################

########################################################################
########################################################################
########################################################################
x=T260[:,1]
y=T260[:,2]
x2=smooth(x,coeff)
y2=smooth(y,coeff)
#plot(x,y)
C=[]
for n in range(2,len(y)-4):
	c=((y2[n+1]-y2[n-1])/(x2[n+1]-x2[n-1]))
	C.append(c)
C=array(C)
C1=[]
for n in range(2,len(x)-4):
	c1=((y2[n+1]-2*y2[n]+y2[n-1])/(x2[n]-x2[n-1])*(x2[n]-x2[n-1]))
	C1.append(c1)
C1=array(C1)

x1=[]
for n in range(len(C)):
	a=T260[n,1]
	x1.append(a)
x1=array(x1)
x2=[]
for n in range(len(C1)):
	b=T260[n,1]
	x2.append(b)
x2=array(x2)
th260=peakdet(C1,50,x2)

d=x[2]-x[1]
thPoint=int((th260[0]/d)+1)
########## Nd(ext)-fitting ##############
x_1=[]
y_1=[]
for n in range(thPoint+5,len(x)):
		b=x[n]
		c=y[n]
		x_1.append(b)
		y_1.append(c)
x_1=array(x_1)
y_1=array(y_1)
(ar,br)=polyfit(x_1,y_1,1)
xr=polyval([ar,br], x_1)
Z=diff(xr)
X=diff(x_1)
Ne260=(Z[0]/X[0])

mx=sqrt(max((y_1-xr)**2))
err=mx
Z1=diff(((xr[len(xr)-2]-err),(xr[1]+err)))
X1=diff((x_1[len(x_1)-2],x_1[1]))
err1=(Z1[0]/X1[0])
Nerr260=Ne260-err1

Irad260=T260[thPoint,3]
	
#show()
########################################################################
########################################################################

########################################################################
########################################################################
########################################################################
x=T270[:,1]
y=T270[:,2]
x2=smooth(x,coeff)
y2=smooth(y,coeff)
#plot(x,y)
C=[]
for n in range(2,len(y)-4):
	c=((y2[n+1]-y2[n-1])/(x2[n+1]-x2[n-1]))
	C.append(c)
C=array(C)
C1=[]
for n in range(2,len(x)-4):
	c1=((y2[n+1]-2*y2[n]+y2[n-1])/(x2[n]-x2[n-1])*(x2[n]-x2[n-1]))
	C1.append(c1)
C1=array(C1)

x1=[]
for n in range(len(C)):
	a=T270[n,1]
	x1.append(a)
x1=array(x1)
x2=[]
for n in range(len(C1)):
	b=T270[n,1]
	x2.append(b)
x2=array(x2)
th270=peakdet(C1,50,x2)

d=x[2]-x[1]
thPoint=int((th270[0]/d)+1)
########## Nd(ext)-fitting ##############
x_1=[]
y_1=[]
for n in range(thPoint+5,len(x)):
		b=x[n]
		c=y[n]
		x_1.append(b)
		y_1.append(c)
x_1=array(x_1)
y_1=array(y_1)
(ar,br)=polyfit(x_1,y_1,1)
xr=polyval([ar,br], x_1)
Z=diff(xr)
X=diff(x_1)
Ne270=(Z[0]/X[0])

mx=sqrt(max((y_1-xr)**2))
err=mx
Z1=diff(((xr[len(xr)-2]-err),(xr[1]+err)))
X1=diff((x_1[len(x_1)-2],x_1[1]))
err1=(Z1[0]/X1[0])
Nerr270=Ne270-err1

Irad270=T270[thPoint,3]

#show()
########################################################################
########################################################################

########################################################################
########################################################################
########################################################################
x=T280[:,1]
y=T280[:,2]
x2=smooth(x,coeff)
y2=smooth(y,coeff)
#plot(x,y)
C=[]
for n in range(2,len(y)-4):
	c=((y2[n+1]-y2[n-1])/(x2[n+1]-x2[n-1]))
	C.append(c)
C=array(C)
C1=[]
for n in range(2,len(x)-4):
	c1=((y2[n+1]-2*y2[n]+y2[n-1])/(x2[n]-x2[n-1])*(x2[n]-x2[n-1]))
	C1.append(c1)
C1=array(C1)

x1=[]
for n in range(len(C)):
	a=T280[n,1]
	x1.append(a)
x1=array(x1)
x2=[]
for n in range(len(C1)):
	b=T280[n,1]
	x2.append(b)
x2=array(x2)
th280=peakdet(C1,50,x2)

d=x[2]-x[1]
thPoint=int((th280[0]/d)+1)
########## Nd(ext)-fitting ##############
x_1=[]
y_1=[]
for n in range(thPoint+5,len(x)):
		b=x[n]
		c=y[n]
		x_1.append(b)
		y_1.append(c)
x_1=array(x_1)
y_1=array(y_1)
(ar,br)=polyfit(x_1,y_1,1)
xr=polyval([ar,br], x_1)
Z=diff(xr)
X=diff(x_1)
Ne280=(Z[0]/X[0])	

mx=sqrt(max((y_1-xr)**2))
err=mx
Z1=diff(((xr[len(xr)-2]-err),(xr[1]+err)))
X1=diff((x_1[len(x_1)-2],x_1[1]))
err1=(Z1[0]/X1[0])
Nerr280=Ne280-err1

Irad280=T280[thPoint,3]

#show()
########################################################################
########################################################################

########################################################################
########################################################################
########################################################################
x=T290[:,1]
y=T290[:,2]
x2=smooth(x,coeff)
y2=smooth(y,coeff)
#plot(x,y)
C=[]
for n in range(2,len(y)-4):
	c=((y2[n+1]-y2[n-1])/(x2[n+1]-x2[n-1]))
	C.append(c)
C=array(C)
C1=[]
for n in range(2,len(x)-4):
	c1=((y2[n+1]-2*y2[n]+y2[n-1])/(x2[n]-x2[n-1])*(x2[n]-x2[n-1]))
	C1.append(c1)
C1=array(C1)

x1=[]
for n in range(len(C)):
	a=T290[n,1]
	x1.append(a)
x1=array(x1)
x2=[]
for n in range(len(C1)):
	b=T290[n,1]
	x2.append(b)
x2=array(x2)
th290=peakdet(C1,50,x2)

d=x[2]-x[1]
thPoint=int((th290[0]/d)+1)
########## Nd(ext)-fitting ##############
x_1=[]
y_1=[]
for n in range(thPoint+5,len(x)):
		b=x[n]
		c=y[n]
		x_1.append(b)
		y_1.append(c)
x_1=array(x_1)
y_1=array(y_1)
(ar,br)=polyfit(x_1,y_1,1)
xr=polyval([ar,br], x_1)
Z=diff(xr)
X=diff(x_1)
Ne290=(Z[0]/X[0])	

mx=sqrt(max((y_1-xr)**2))
err=mx
Z1=diff(((xr[len(xr)-2]-err),(xr[1]+err)))
X1=diff((x_1[len(x_1)-2],x_1[1]))
err1=(Z1[0]/X1[0])
Nerr290=Ne290-err1

Irad290=T290[thPoint,3]

#show()
########################################################################
########################################################################

########################################################################
########################################################################
########################################################################
x=T300[:,1]
y=T300[:,2]
x2=smooth(x,coeff)
y2=smooth(y,coeff)
#plot(x,y)
C=[]
for n in range(2,len(y)-4):
	c=((y2[n+1]-y2[n-1])/(x2[n+1]-x2[n-1]))
	C.append(c)
C=array(C)
C1=[]
for n in range(2,len(x)-4):
	c1=((y2[n+1]-2*y2[n]+y2[n-1])/(x2[n]-x2[n-1])*(x2[n]-x2[n-1]))
	C1.append(c1)
C1=array(C1)

x1=[]
for n in range(len(C)):
	a=T300[n,1]
	x1.append(a)
x1=array(x1)
x2=[]
for n in range(len(C1)):
	b=T300[n,1]
	x2.append(b)
x2=array(x2)
th300=peakdet(C1,50,x2)

d=x[2]-x[1]
thPoint=int((th300[0]/d)+1)
########## Nd(ext)-fitting ##############
x_1=[]
y_1=[]
for n in range(thPoint+5,len(x)):
		b=x[n]
		c=y[n]
		x_1.append(b)
		y_1.append(c)
x_1=array(x_1)
y_1=array(y_1)
(ar,br)=polyfit(x_1,y_1,1)
xr=polyval([ar,br], x_1)
Z=diff(xr)
X=diff(x_1)
Ne300=(Z[0]/X[0])

mx=sqrt(max((y_1-xr)**2))
err=mx
Z1=diff(((xr[len(xr)-2]-err),(xr[1]+err)))
X1=diff((x_1[len(x_1)-2],x_1[1]))
err1=(Z1[0]/X1[0])
Nerr300=Ne300-err1

Irad300=T300[thPoint,3]
	
#show()
########################################################################
########################################################################

########################################################################
########################################################################
########################################################################
x=T310[:,1]
y=T310[:,2]
x2=smooth(x,coeff)
y2=smooth(y,coeff)
#plot(x,y)
C=[]
for n in range(2,len(y)-4):
	c=((y2[n+1]-y2[n-1])/(x2[n+1]-x2[n-1]))
	C.append(c)
C=array(C)
C1=[]
for n in range(2,len(x)-4):
	c1=((y2[n+1]-2*y2[n]+y2[n-1])/(x2[n]-x2[n-1])*(x2[n]-x2[n-1]))
	C1.append(c1)
C1=array(C1)

x1=[]
for n in range(len(C)):
	a=T310[n,1]
	x1.append(a)
x1=array(x1)
x2=[]
for n in range(len(C1)):
	b=T310[n,1]
	x2.append(b)
x2=array(x2)
th310=peakdet(C1,50,x2)

d=x[2]-x[1]
thPoint=int((th310[0]/d)+1)
########## Nd(ext)-fitting ##############
x_1=[]
y_1=[]
for n in range(thPoint+5,len(x)):
		b=x[n]
		c=y[n]
		x_1.append(b)
		y_1.append(c)
x_1=array(x_1)
y_1=array(y_1)
(ar,br)=polyfit(x_1,y_1,1)
xr=polyval([ar,br], x_1)
Z=diff(xr)
X=diff(x_1)
Ne310=(Z[0]/X[0])

mx=sqrt(max((y_1-xr)**2))
err=mx
Z1=diff(((xr[len(xr)-2]-err),(xr[1]+err)))
X1=diff((x_1[len(x_1)-2],x_1[1]))
err1=(Z1[0]/X1[0])
Nerr310=Ne310-err1

Irad310=T310[thPoint,3]
	
#show()
########################################################################
########################################################################

########################################################################
########################################################################
########################################################################
x=T320[:,1]
y=T320[:,2]
x2=smooth(x,coeff)
y2=smooth(y,coeff)
#plot(x,y)
C=[]
for n in range(2,len(y)-4):
	c=((y2[n+1]-y2[n-1])/(x2[n+1]-x2[n-1]))
	C.append(c)
C=array(C)
C1=[]
for n in range(2,len(x)-4):
	c1=((y2[n+1]-2*y2[n]+y2[n-1])/(x2[n]-x2[n-1])*(x2[n]-x2[n-1]))
	C1.append(c1)
C1=array(C1)

x1=[]
for n in range(len(C)):
	a=T320[n,1]
	x1.append(a)
x1=array(x1)
x2=[]
for n in range(len(C1)):
	b=T320[n,1]
	x2.append(b)
x2=array(x2)
th320=peakdet(C1,50,x2)

d=x[2]-x[1]
thPoint=int((th320[0]/d)+1)
########## Nd(ext)-fitting ##############
x_1=[]
y_1=[]
for n in range(thPoint+5,len(x)):
		b=x[n]
		c=y[n]
		x_1.append(b)
		y_1.append(c)
x_1=array(x_1)
y_1=array(y_1)
(ar,br)=polyfit(x_1,y_1,1)
xr=polyval([ar,br], x_1)
Z=diff(xr)
X=diff(x_1)
Ne320=(Z[0]/X[0])

mx=sqrt(max((y_1-xr)**2))
err=mx
Z1=diff(((xr[len(xr)-2]-err),(xr[1]+err)))
X1=diff((x_1[len(x_1)-2],x_1[1]))
err1=(Z1[0]/X1[0])
Nerr320=Ne320-err1

Irad320=T320[thPoint,3]

#show()
########################################################################
########################################################################

########################################################################
########################################################################
########################################################################
x=T330[:,1]
y=T330[:,2]
x2=smooth(x,coeff)
y2=smooth(y,coeff)
#plot(x,y)
C=[]
for n in range(2,len(y)-4):
	c=((y2[n+1]-y2[n-1])/(x2[n+1]-x2[n-1]))
	C.append(c)
C=array(C)
C1=[]
for n in range(2,len(x)-4):
	c1=((y2[n+1]-2*y2[n]+y2[n-1])/(x2[n]-x2[n-1])*(x2[n]-x2[n-1]))
	C1.append(c1)
C1=array(C1)

x1=[]
for n in range(len(C)):
	a=T330[n,1]
	x1.append(a)
x1=array(x1)
x2=[]
for n in range(len(C1)):
	b=T330[n,1]
	x2.append(b)
x2=array(x2)
th330=peakdet(C1,50,x2)

d=x[2]-x[1]
thPoint=int((th330[0]/d)+1)
########## Nd(ext)-fitting ##############
x_1=[]
y_1=[]
for n in range(thPoint+5,len(x)):
		b=x[n]
		c=y[n]
		x_1.append(b)
		y_1.append(c)
x_1=array(x_1)
y_1=array(y_1)
(ar,br)=polyfit(x_1,y_1,1)
xr=polyval([ar,br], x_1)
Z=diff(xr)
X=diff(x_1)
Ne330=(Z[0]/X[0])

mx=sqrt(max((y_1-xr)**2))
err=mx
Z1=diff(((xr[len(xr)-2]-err),(xr[1]+err)))
X1=diff((x_1[len(x_1)-2],x_1[1]))
err1=(Z1[0]/X1[0])
Nerr330=Ne330-err1

Irad330=T330[thPoint,3]
	
#show()
########################################################################
########################################################################

########################################################################
########################################################################
########################################################################
x=T340[:,1]
y=T340[:,2]
x2=smooth(x,coeff)
y2=smooth(y,coeff)
#plot(x,y)
C=[]
for n in range(2,len(y)-4):
	c=((y2[n+1]-y2[n-1])/(x2[n+1]-x2[n-1]))
	C.append(c)
C=array(C)
C1=[]
for n in range(2,len(x)-4):
	c1=((y2[n+1]-2*y2[n]+y2[n-1])/(x2[n]-x2[n-1])*(x2[n]-x2[n-1]))
	C1.append(c1)
C1=array(C1)

x1=[]
for n in range(len(C)):
	a=T340[n,1]
	x1.append(a)
x1=array(x1)
x2=[]
for n in range(len(C1)):
	b=T340[n,1]
	x2.append(b)
x2=array(x2)


th340=peakdet(C1,50,x2)
d=x[2]-x[1]
thPoint=int((th340[0]/d)+1)

########## Nd(ext)-fitting ##############
x_1=[]
y_1=[]
for n in range(thPoint+5,len(x)):
		b=x[n]
		c=y[n]
		x_1.append(b)
		y_1.append(c)
x_1=array(x_1)
y_1=array(y_1)
(ar,br)=polyfit(x_1,y_1,1)
xr=polyval([ar,br], x_1)
Z=diff(xr)
X=diff(x_1)
Ne340=(Z[0]/X[0])

mx=sqrt(max((y_1-xr)**2))
err=mx
Z1=diff(((xr[len(xr)-2]-err),(xr[1]+err)))
X1=diff((x_1[len(x_1)-2],x_1[1]))
err1=(Z1[0]/X1[0])
Nerr340=Ne340-err1

Irad340=T340[thPoint,3]	
show()
########################################################################
########################################################################

########################################################################
########################################################################
########################################################################











#				Z PARAMETER
#################################
coeff1=calc_coeff(1,1,0)
########################################################################
########################################################################
########################################################################
x=(log(sqrt(T50[:,3])))
y=(log(T50[:,1]))
x2=T50[:,3]
y2=T50[:,1]
x3=x#smooth(x,coeff1)
y3=y#smooth(x,coeff1)
d=y2[2]-y2[1]
thPoint=int((th50[0]*.8/d)-1)
thPoint1=int(((th50[0]*0.3)/d)+1)
x50=diff(x3)
y50=diff(y3)

########## Z-fitting ##############
x_150=[]
y_1=[]
for n in range(thPoint1,thPoint):
		b=x3[n]
		c=y3[n]
		x_150.append(b)
		y_1.append(c)
x_150=array(x_150)
y_1=array(y_1)
(ar,br)=polyfit(x_150,y_1,1)
xr50=polyval([ar,br], x_150)

Z=diff(xr50)
X=diff(x_150)
Z50=(Z[0]/X[0])


mx=sqrt(max((y_1-xr50)**2))
err=mx
Z1=diff(((xr50[len(xr50)-2]-err),(xr50[1]+err)))
X1=diff((x_150[len(x_150)-2]*-1,x_150[1]*-1))
err1=(Z1[0]/-X1[0])
err50=Z50-err1
###################################	


########################################################################
########################################################################
########################################################################
x=(log(sqrt(T60[:,3])))
y=(log(T60[:,1]))
x2=T60[:,3]
y2=T60[:,1]
x3=x#smooth(x,coeff1)
y3=y#smooth(x,coeff1)
d=y2[2]-y2[1]
thPoint=int((th60[0]*.8/d)-1)
thPoint1=int(((th60[0]*0.3)/d)+1)
########## Z-fitting ##############
x_160=[]
y_1=[]
for n in range(thPoint1,thPoint):
		b=x3[n]
		c=y3[n]
		x_160.append(b)
		y_1.append(c)
x_160=array(x_160)
y_1=array(y_1)
(ar,br)=polyfit(x_160,y_1,1)
xr60=polyval([ar,br], x_160)

plot(x,y,label='SE @ 60K')
plot(x_160,xr60,'o',label='Linear fit between Ith/3 and Ith')
Z=diff(xr60)
X=diff(x_160)
Z60=(Z[0]/X[0])

mx=sqrt(max((y_1-xr60)**2))
err=mx
Z1=diff(((xr60[len(xr60)-2]-err),(xr60[1]+err)))
X1=diff((x_160[len(x_160)-2]*-1,x_160[1]*-1))
err1=(Z1[0]/-X1[0])
err60=Z60-err1

xlabel('ln sqrt SE')
ylabel('ln I')
title('Linear fitting')
legend(loc='lower right', ncol=1, shadow=True, numpoints=1)
show()
plot(x_160,xr60,label='60')


### Z Vs I ####
x60=diff(x_160)
y60=diff(y_1)
x1_60=[]
for n in range(thPoint1,thPoint-1):
	x1=y2[n]
	x1_60.append(x1)
x1_60=array(x1_60)
###################################	


########################################################################
########################################################################
########################################################################
x=(log(sqrt(T70[:,3])))
y=(log(T70[:,1]))
x2=T70[:,3]
y2=T70[:,1]
x3=x#smooth(x,coeff1)
y3=y#smooth(x,coeff1)
C=[]
for n in range(2,len(y)-1):
	c=((y3[n+1]-y3[n-1])/(x3[n+1]-x3[n-1]))
	C.append(c)
C=array(C)
x1=[]
for n in range(len(C)):
	a=x[n]
	x1.append(a)
x1=array(x1)
d=y2[2]-y2[1]
thPoint=int((th70[0]*.8/d)-1)
thPoint1=int(((th70[0]*0.3)/d)+1)

########## Z-fitting ##############
x_170=[]
y_1=[]
for n in range(thPoint1,thPoint):
		b=x3[n]
		c=y3[n]
		x_170.append(b)
		y_1.append(c)
x_170=array(x_170)
y_1=array(y_1)
(ar,br)=polyfit(x_170,y_1,1)
xr70=polyval([ar,br], x_170)


plot(x_170,xr70,label='70')
Z=diff(xr70)
X=diff(x_170)
Z70=(Z[0]/X[0])

mx=sqrt(max((y_1-xr70)**2))
err=mx
Z1=diff(((xr70[len(xr70)-2]-err),(xr70[1]+err)))
X1=diff((x_170[len(x_170)-2]*-1,x_170[1]*-1))
err1=(Z1[0]/-X1[0])
err70=Z70-err1
###################################	

########################################################################
########################################################################
########################################################################
x=(log(sqrt(T80[:,3])))
y=(log(T80[:,1]))
x2=T80[:,3]
y2=T80[:,1]
x3=x#smooth(x,coeff1)
y3=y#smooth(x,coeff1)
C=[]
for n in range(2,len(y)-1):
	c=((y3[n+1]-y3[n-1])/(x3[n+1]-x3[n-1]))
	C.append(c)
C=array(C)
x1=[]
for n in range(len(C)):
	a=x[n]
	x1.append(a)
x1=array(x1)
d=y2[2]-y2[1]
thPoint=int((th80[0]*.8/d)-1)
thPoint1=int(((th80[0]*0.3)/d)+1)

########## Z-fitting ##############
x_180=[]
y_1=[]
for n in range(thPoint1,thPoint):
		b=x3[n]
		c=y3[n]
		x_180.append(b)
		y_1.append(c)
x_180=array(x_180)
y_1=array(y_1)
(ar,br)=polyfit(x_180,y_1,1)
xr80=polyval([ar,br], x_180)

plot(x_180,xr80,label='80')
Z=diff(xr80)
X=diff(x_180)
Z80=(Z[0]/X[0])

mx=sqrt(max((y_1-xr80)**2))
err=mx
Z1=diff(((xr80[len(xr80)-2]-err),(xr80[1]+err)))
X1=diff((x_180[len(x_180)-2]*-1,x_180[1]*-1))
err1=(Z1[0]/-X1[0])
err80=Z80-err1
###################################	


########################################################################
########################################################################
########################################################################
x=(log(sqrt(T90[:,3])))
y=(log(T90[:,1]))
x2=T90[:,3]
y2=T90[:,1]
x3=x#smooth(x,coeff1)
y3=y#smooth(x,coeff1)
C=[]
for n in range(2,len(y)-1):
	c=((y3[n+1]-y3[n-1])/(x3[n+1]-x3[n-1]))
	C.append(c)
C=array(C)
x1=[]
for n in range(len(C)):
	a=x[n]
	x1.append(a)
x1=array(x1)
d=y2[2]-y2[1]
thPoint=int((th90[0]*.8/d)-1)
thPoint1=int(((th90[0]*0.3)/d)+1)

########## Z-fitting ##############
x_190=[]
y_1=[]
for n in range(thPoint1,thPoint):
		b=x3[n]
		c=y3[n]
		x_190.append(b)
		y_1.append(c)
x_190=array(x_190)
y_1=array(y_1)
(ar,br)=polyfit(x_190,y_1,1)
xr90=polyval([ar,br], x_190)

plot(x_190,xr90,label='90')
Z=diff(xr90)
X=diff(x_190)
Z90=(Z[0]/X[0])

mx=sqrt(max((y_1-xr90)**2))
err=mx
Z1=diff(((xr90[len(xr90)-2]-err),(xr90[1]+err)))
X1=diff((x_190[len(x_190)-2]*-1,x_190[1]*-1))
err1=(Z1[0]/-X1[0])
err90=Z90-err1
###################################	


########################################################################
########################################################################
########################################################################
x=(log(sqrt(T100[:,3])))
y=(log(T100[:,1]))
x2=T100[:,3]
y2=T100[:,1]
x3=x#smooth(x,coeff1)
y3=y#smooth(x,coeff1)
d=y2[2]-y2[1]
thPoint=int((th100[0]*.8/d)-1)
thPoint1=int(((th100[0]*0.3)/d)+1)


########## Z-fitting ##############
x_1100=[]
y_1=[]
for n in range(thPoint1,thPoint):
		b=x3[n]
		c=y3[n]
		x_1100.append(b)
		y_1.append(c)
x_1100=array(x_1100)
y_1=array(y_1)
(ar,br)=polyfit(x_1100,y_1,1)
xr100=polyval([ar,br], x_1100)


plot(x_1100,xr100,label='100')
Z=diff(xr100)
X=diff(x_1100)
Z100=(Z[0]/X[0])

mx=sqrt(max((y_1-xr100)**2))
err=mx
Z1=diff(((xr100[len(xr100)-2]-err),(xr100[1]+err)))
X1=diff((x_1100[len(x_1100)-2]*-1,x_1100[1]*-1))
err1=(Z1[0]/-X1[0])
err100=Z100-err1


### Z Vs I ####
x100=diff(x_1100)
y100=diff(y_1)
x1_100=[]
for n in range(thPoint1,thPoint-1):
	x1=y2[n]
	x1_100.append(x1)
x1_100=array(x1_100)
###################################	


########################################################################
########################################################################
########################################################################
x=(log(sqrt(T110[:,3])))
y=(log(T110[:,1]))
x2=T110[:,3]
y2=T110[:,1]
x3=x#smooth(x,coeff1)
y3=y#smooth(x,coeff1)
C=[]
for n in range(2,len(y)-1):
	c=((y3[n+1]-y3[n-1])/(x3[n+1]-x3[n-1]))
	C.append(c)
C=array(C)
x1=[]
for n in range(len(C)):
	a=x[n]
	x1.append(a)
x1=array(x1)
d=y2[2]-y2[1]
thPoint=int((th110[0]*.8/d)-1)
thPoint1=int(((th110[0]*0.3)/d)+1)

########## Z-fitting ##############
x_1110=[]
y_1=[]
for n in range(thPoint1,thPoint):
		b=x3[n]
		c=y3[n]
		x_1110.append(b)
		y_1.append(c)
x_1110=array(x_1110)
y_1=array(y_1)
(ar,br)=polyfit(x_1110,y_1,1)
xr110=polyval([ar,br], x_1110)

plot(x_1110,xr110,label='110')
Z=diff(xr110)
X=diff(x_1110)
Z110=(Z[0]/X[0])

mx=sqrt(max((y_1-xr110)**2))
err=mx
Z1=diff(((xr110[len(xr110)-2]-err),(xr110[1]+err)))
X1=diff((x_1110[len(x_1110)-2]*-1,x_1110[1]*-1))
err1=(Z1[0]/-X1[0])
err110=Z110-err1
###################################	


########################################################################
########################################################################
########################################################################
x=(log(sqrt(T120[:,3])))
y=(log(T120[:,1]))
x2=T120[:,3]
y2=T120[:,1]
x3=x#smooth(x,coeff1)
y3=y#smooth(x,coeff1)
C=[]
for n in range(2,len(y)-1):
	c=((y3[n+1]-y3[n-1])/(x3[n+1]-x3[n-1]))
	C.append(c)
C=array(C)
x1=[]
for n in range(len(C)):
	a=x[n]
	x1.append(a)
x1=array(x1)
d=y2[2]-y2[1]
thPoint=int((th120[0]*.8/d)-1)
thPoint1=int(((th120[0]*0.3)/d)+1)

########## Z-fitting ##############
x_1120=[]
y_1=[]
for n in range(thPoint1,thPoint):
		b=x3[n]
		c=y3[n]
		x_1120.append(b)
		y_1.append(c)
x_1120=array(x_1120)
y_1=array(y_1)
(ar,br)=polyfit(x_1120,y_1,1)
xr120=polyval([ar,br], x_1120)

plot(x_1120,xr120,label='120')
Z=diff(xr120)
X=diff(x_1120)
Z120=(Z[0]/X[0])

mx=sqrt(max((y_1-xr120)**2))
err=mx
Z1=diff(((xr120[len(xr120)-2]-err),(xr120[1]+err)))
X1=diff((x_1120[len(x_1120)-2]*-1,x_1120[1]*-1))
err1=(Z1[0]/-X1[0])
err120=Z120-err1
###################################	


########################################################################
########################################################################
########################################################################
x=(log(sqrt(T130[:,3])))
y=(log(T130[:,1]))
x2=T130[:,3]
y2=T130[:,1]
x3=x#smooth(x,coeff1)
y3=y#smooth(x,coeff1)
C=[]
for n in range(2,len(y)-1):
	c=((y3[n+1]-y3[n-1])/(x3[n+1]-x3[n-1]))
	C.append(c)
C=array(C)
x1=[]
for n in range(len(C)):
	a=x[n]
	x1.append(a)
x1=array(x1)
d=y2[2]-y2[1]
thPoint=int((th130[0]*.8/d)-1)
thPoint1=int(((th130[0]*0.3)/d)+1)

########## Z-fitting ##############
x_1130=[]
y_1=[]
for n in range(thPoint1,thPoint):
		b=x3[n]
		c=y3[n]
		x_1130.append(b)
		y_1.append(c)
x_1130=array(x_1130)
y_1=array(y_1)
(ar,br)=polyfit(x_1130,y_1,1)
xr130=polyval([ar,br], x_1130)

plot(x_1130,xr130,label='130')
Z=diff(xr130)
X=diff(x_1130)
Z130=(Z[0]/X[0])

mx=sqrt(max((y_1-xr130)**2))
err=mx
Z1=diff(((xr130[len(xr130)-2]-err),(xr130[1]+err)))
X1=diff((x_1130[len(x_1130)-2]*-1,x_1130[1]*-1))
err1=(Z1[0]/-X1[0])
err130=Z130-err1
###################################	


########################################################################
########################################################################
########################################################################
x=(log(sqrt(T140[:,3])))
y=(log(T140[:,1]))
x2=T140[:,3]
y2=T140[:,1]
x3=x#smooth(x,coeff1)
y3=y#smooth(x,coeff1)
C=[]
for n in range(2,len(y)-1):
	c=((y3[n+1]-y3[n-1])/(x3[n+1]-x3[n-1]))
	C.append(c)
C=array(C)
x1=[]
for n in range(len(C)):
	a=x[n]
	x1.append(a)
x1=array(x1)
d=y2[2]-y2[1]
thPoint=int((th140[0]*.8/d)-1)
thPoint1=int(((th140[0]*0.3)/d)+1)

########## Z-fitting ##############
x_1140=[]
y_1=[]
for n in range(thPoint1,thPoint):
		b=x3[n]
		c=y3[n]
		x_1140.append(b)
		y_1.append(c)
x_1140=array(x_1140)
y_1=array(y_1)
(ar,br)=polyfit(x_1140,y_1,1)
xr140=polyval([ar,br], x_1140)

plot(x_1140,xr140,label='140')
Z=diff(xr140)
X=diff(x_1140)
Z140=(Z[0]/X[0])

mx=sqrt(max((y_1-xr140)**2))
err=mx
Z1=diff(((xr140[len(xr140)-2]-err),(xr140[1]+err)))
X1=diff((x_1140[len(x_1140)-2]*-1,x_1140[1]*-1))
err1=(Z1[0]/-X1[0])
err140=Z140-err1
###################################	


########################################################################
########################################################################
########################################################################
x=(log(sqrt(T150[:,3])))
y=(log(T150[:,1]))
x2=T150[:,3]
y2=T150[:,1]
x3=x#smooth(x,coeff1)
y3=y#smooth(x,coeff1)
C=[]
for n in range(2,len(y)-1):
	c=((y3[n+1]-y3[n-1])/(x3[n+1]-x3[n-1]))
	C.append(c)
C=array(C)
x1=[]
for n in range(len(C)):
	a=x[n]
	x1.append(a)
x1=array(x1)
d=y2[2]-y2[1]
thPoint=int((th150[0]*.8/d)-1)
thPoint1=int(((th150[0]*0.3)/d)+1)

########## Z-fitting ##############
x_1150=[]
y_1=[]
for n in range(thPoint1,thPoint):
		b=x3[n]
		c=y3[n]
		x_1150.append(b)
		y_1.append(c)
x_1150=array(x_1150)
y_1=array(y_1)
(ar,br)=polyfit(x_1150,y_1,1)
xr150=polyval([ar,br], x_1150)

plot(x_1150,xr150,label='150')
Z=diff(xr150)
X=diff(x_1150)
Z150=(Z[0]/X[0])

mx=sqrt(max((y_1-xr150)**2))
err=mx
Z1=diff(((xr150[len(xr150)-2]-err),(xr150[1]+err)))
X1=diff((x_1150[len(x_1150)-2]*-1,x_1150[1]*-1))
err1=(Z1[0]/-X1[0])
err150=Z150-err1
###################################	


########################################################################
########################################################################
########################################################################
x=(log(sqrt(T160[:,3])))
y=(log(T160[:,1]))
x2=T160[:,3]
y2=T160[:,1]
x3=x#smooth(x,coeff1)
y3=y#smooth(x,coeff1)
C=[]
for n in range(2,len(y)-1):
	c=((y3[n+1]-y3[n-1])/(x3[n+1]-x3[n-1]))
	C.append(c)
C=array(C)
x1=[]
for n in range(len(C)):
	a=x[n]
	x1.append(a)
x1=array(x1)
d=y2[2]-y2[1]
thPoint=int((th160[0]*.8/d)-1)
thPoint1=int(((th160[0]*0.3)/d)+1)


########## Z-fitting ##############
x_1160=[]
y_1=[]
for n in range(thPoint1,thPoint):
		b=x3[n]
		c=y3[n]
		x_1160.append(b)
		y_1.append(c)
x_1160=array(x_1160)
y_1=array(y_1)
(ar,br)=polyfit(x_1160,y_1,1)
xr160=polyval([ar,br], x_1160)

plot(x_1160,xr160,label='160')
Z=diff(xr160)
X=diff(x_1160)
Z160=(Z[0]/X[0])

mx=sqrt(max((y_1-xr160)**2))
err=mx
Z1=diff(((xr160[len(xr160)-2]-err),(xr160[1]+err)))
X1=diff((x_1160[len(x_1160)-2]*-1,x_1160[1]*-1))
err1=(Z1[0]/-X1[0])
err160=Z160-err1
###################################	


########################################################################
########################################################################
########################################################################
x=(log(sqrt(T170[:,3])))
y=(log(T170[:,1]))
x2=T170[:,3]
y2=T170[:,1]
x3=x#smooth(x,coeff1)
y3=y#smooth(x,coeff1)
C=[]
for n in range(2,len(y)-1):
	c=((y3[n+1]-y3[n-1])/(x3[n+1]-x3[n-1]))
	C.append(c)
C=array(C)
x1=[]
for n in range(len(C)):
	a=x[n]
	x1.append(a)
x1=array(x1)
d=y2[2]-y2[1]
thPoint=int((th170[0]*.8/d)-1)
thPoint1=int(((th170[0]*0.3)/d)+1)

########## Z-fitting ##############
x_1170=[]
y_1=[]
for n in range(thPoint1,thPoint):
		b=x3[n]
		c=y3[n]
		x_1170.append(b)
		y_1.append(c)
x_1170=array(x_1170)
y_1=array(y_1)
(ar,br)=polyfit(x_1170,y_1,1)
xr170=polyval([ar,br], x_1170)

plot(x_1170,xr170,label='170')
Z=diff(xr170)
X=diff(x_1170)
Z170=(Z[0]/X[0])

mx=sqrt(max((y_1-xr170)**2))
err=mx
Z1=diff(((xr170[len(xr170)-2]-err),(xr170[1]+err)))
X1=diff((x_1170[len(x_1170)-2]*-1,x_1170[1]*-1))
err1=(Z1[0]/-X1[0])
err170=Z170-err1
###################################	


########################################################################
########################################################################
########################################################################
x=(log(sqrt(T180[:,3])))
y=(log(T180[:,1]))
x2=T180[:,3]
y2=T180[:,1]
x3=x#smooth(x,coeff1)
y3=y#smooth(x,coeff1)
C=[]
for n in range(2,len(y)-1):
	c=((y3[n+1]-y3[n-1])/(x3[n+1]-x3[n-1]))
	C.append(c)
C=array(C)
x1=[]
for n in range(len(C)):
	a=x[n]
	x1.append(a)
x1=array(x1)
d=y2[2]-y2[1]
thPoint=int((th180[0]*.8/d)-1)
thPoint1=int(((th180[0]*0.3)/d)+1)

########## Z-fitting ##############
x_1180=[]
y_1=[]
for n in range(thPoint1,thPoint):
		b=x3[n]
		c=y3[n]
		x_1180.append(b)
		y_1.append(c)
x_1180=array(x_1180)
y_1=array(y_1)
(ar,br)=polyfit(x_1180,y_1,1)
xr180=polyval([ar,br], x_1180)

plot(x_1180,xr180,label='180')
Z=diff(xr180)
X=diff(x_1180)
Z180=(Z[0]/X[0])

mx=sqrt(max((y_1-xr180)**2))
err=mx
Z1=diff(((xr180[len(xr180)-2]-err),(xr180[1]+err)))
X1=diff((x_1180[len(x_1180)-2]*-1,x_1180[1]*-1))
err1=(Z1[0]/-X1[0])
err180=Z180-err1
###################################	


########################################################################
########################################################################
########################################################################
x=(log(sqrt(T190[:,3])))
y=(log(T190[:,1]))
x2=T190[:,3]
y2=T190[:,1]
x3=x#smooth(x,coeff1)
y3=y#smooth(x,coeff1)
C=[]
for n in range(2,len(y)-1):
	c=((y3[n+1]-y3[n-1])/(x3[n+1]-x3[n-1]))
	C.append(c)
C=array(C)
x1=[]
for n in range(len(C)):
	a=x[n]
	x1.append(a)
x1=array(x1)
d=y2[2]-y2[1]
thPoint=int((th190[0]*.8/d)-1)
thPoint1=int(((th190[0]*0.3)/d)+1)

########## Z-fitting ##############
x_1190=[]
y_1=[]
for n in range(thPoint1,thPoint):
		b=x3[n]
		c=y3[n]
		x_1190.append(b)
		y_1.append(c)
x_1190=array(x_1190)
y_1=array(y_1)
(ar,br)=polyfit(x_1190,y_1,1)
xr190=polyval([ar,br], x_1190)

plot(x_1190,xr190,label='190')
Z=diff(xr190)
X=diff(x_1190)
Z190=(Z[0]/X[0])

mx=sqrt(max((y_1-xr190)**2))
err=mx
Z1=diff(((xr190[len(xr190)-2]-err),(xr190[1]+err)))
X1=diff((x_1190[len(x_1190)-2]*-1,x_1190[1]*-1))
err1=(Z1[0]/-X1[0])
err190=Z190-err1
###################################
	

########################################################################
########################################################################
########################################################################
x=(log(sqrt(T200[:,3])))
y=(log(T200[:,1]))
x2=T200[:,3]
y2=T200[:,1]
x3=x#smooth(x,coeff1)
y3=y#smooth(x,coeff1)
d=y2[2]-y2[1]
thPoint=int((th200[0]*.8/d)-1)
thPoint1=int(((th200[0]*0.3)/d)+1)

########## Z-fitting ##############
x_1200=[]
y_1=[]
for n in range(thPoint1,thPoint):
		b=x3[n]
		c=y3[n]
		x_1200.append(b)
		y_1.append(c)
x_1200=array(x_1200)
y_1=array(y_1)
(ar,br)=polyfit(x_1200,y_1,1)
xr200=polyval([ar,br], x_1200)

plot(x_1200,xr200,label='200')
Z=diff(xr200)
X=diff(x_1200)
Z200=(Z[0]/X[0])

mx=sqrt(max((y_1-xr200)**2))
err=mx
Z1=diff(((xr200[len(xr200)-2]-err),(xr200[1]+err)))
X1=diff((x_1200[len(x_1200)-2]*-1,x_1200[1]*-1))
err1=(Z1[0]/-X1[0])
err200=Z200-err1

### Z Vs I ####
x200=diff(x_1200)
y200=diff(y_1)
x1_200=[]
for n in range(thPoint1,thPoint-1):
	x1=y2[n]
	x1_200.append(x1)
x1_200=array(x1_200)
###################################	


########################################################################
########################################################################
########################################################################
x=(log(sqrt(T210[:,3])))
y=(log(T210[:,1]))
x2=T210[:,3]
y2=T210[:,1]
x3=x#smooth(x,coeff1)
y3=y#smooth(x,coeff1)
C=[]
for n in range(2,len(y)-1):
	c=((y3[n+1]-y3[n-1])/(x3[n+1]-x3[n-1]))
	C.append(c)
C=array(C)
x1=[]
for n in range(len(C)):
	a=x[n]
	x1.append(a)
x1=array(x1)
d=y2[2]-y2[1]
thPoint=int((th210[0]*.8/d)-1)
thPoint1=int(((th210[0]*0.3)/d)+1)

########## Z-fitting ##############
x_1210=[]
y_1=[]
for n in range(thPoint1,thPoint):
		b=x3[n]
		c=y3[n]
		x_1210.append(b)
		y_1.append(c)
x_1210=array(x_1210)
y_1=array(y_1)
(ar,br)=polyfit(x_1210,y_1,1)
xr210=polyval([ar,br], x_1210)

plot(x_1210,xr210,label='210')
Z=diff(xr210)
X=diff(x_1210)
Z210=(Z[0]/X[0])

mx=sqrt(max((y_1-xr210)**2))
err=mx
Z1=diff(((xr210[len(xr210)-2]-err),(xr210[1]+err)))
X1=diff((x_1210[len(x_1210)-2]*-1,x_1210[1]*-1))
err1=(Z1[0]/-X1[0])
err210=Z210-err1
###################################	

########################################################################
########################################################################
########################################################################
x=(log(sqrt(T220[:,3])))
y=(log(T220[:,1]))
x2=T220[:,3]
y2=T220[:,1]
x3=x#smooth(x,coeff1)
y3=y#smooth(x,coeff1)
C=[]
for n in range(2,len(y)-1):
	c=((y3[n+1]-y3[n-1])/(x3[n+1]-x3[n-1]))
	C.append(c)
C=array(C)
x1=[]
for n in range(len(C)):
	a=x[n]
	x1.append(a)
x1=array(x1)
d=y2[2]-y2[1]
thPoint=int((th220[0]*.8/d)-1)
thPoint1=int(((th220[0]*0.3)/d)+1)

########## Z-fitting ##############
x_1220=[]
y_1=[]
for n in range(thPoint1,thPoint):
		b=x3[n]
		c=y3[n]
		x_1220.append(b)
		y_1.append(c)
x_1220=array(x_1220)
y_1=array(y_1)
(ar,br)=polyfit(x_1220,y_1,1)
xr220=polyval([ar,br], x_1220)

plot(x_1220,xr220,label='220')
Z=diff(xr220)
X=diff(x_1220)
Z220=(Z[0]/X[0])

mx=sqrt(max((y_1-xr220)**2))
err=mx
Z1=diff(((xr220[len(xr220)-2]-err),(xr220[1]+err)))
X1=diff((x_1220[len(x_1220)-2]*-1,x_1220[1]*-1))
err1=(Z1[0]/-X1[0])
err220=Z220-err1
###################################	


########################################################################
########################################################################
########################################################################
x=(log(sqrt(T230[:,3])))
y=(log(T230[:,1]))
x2=T230[:,3]
y2=T230[:,1]
x3=x#smooth(x,coeff1)
y3=y#smooth(x,coeff1)
C=[]
for n in range(2,len(y)-1):
	c=((y3[n+1]-y3[n-1])/(x3[n+1]-x3[n-1]))
	C.append(c)
C=array(C)
x1=[]
for n in range(len(C)):
	a=x[n]
	x1.append(a)
x1=array(x1)
d=y2[2]-y2[1]
thPoint=int((th230[0]*.8/d)-1)
thPoint1=int(((th230[0]*0.3)/d)+1)

########## Z-fitting ##############
x_1230=[]
y_1=[]
for n in range(thPoint1,thPoint):
		b=x3[n]
		c=y3[n]
		x_1230.append(b)
		y_1.append(c)
x_1230=array(x_1230)
y_1=array(y_1)
(ar,br)=polyfit(x_1230,y_1,1)
xr230=polyval([ar,br], x_1230)

plot(x_1230,xr230,label='230')
Z=diff(xr230)
X=diff(x_1230)
Z230=(Z[0]/X[0])

mx=sqrt(max((y_1-xr230)**2))
err=mx
Z1=diff(((xr230[len(xr230)-2]-err),(xr230[1]+err)))
X1=diff((x_1230[len(x_1230)-2]*-1,x_1230[1]*-1))
err1=(Z1[0]/-X1[0])
err230=Z230-err1
###################################	

########################################################################
########################################################################
########################################################################
x=(log(sqrt(T240[:,3])))
y=(log(T240[:,1]))
x2=T240[:,3]
y2=T240[:,1]
x3=x#smooth(x,coeff1)
y3=y#smooth(x,coeff1)
C=[]
for n in range(2,len(y)-1):
	c=((y3[n+1]-y3[n-1])/(x3[n+1]-x3[n-1]))
	C.append(c)
C=array(C)
x1=[]
for n in range(len(C)):
	a=x[n]
	x1.append(a)
x1=array(x1)
d=y2[2]-y2[1]
thPoint=int((th240[0]*.8/d)-1)
thPoint1=int(((th240[0]*0.3)/d)+1)

########## Z-fitting ##############
x_1240=[]
y_1=[]
for n in range(thPoint1,thPoint):
		b=x3[n]
		c=y3[n]
		x_1240.append(b)
		y_1.append(c)
x_1240=array(x_1240)
y_1=array(y_1)
(ar,br)=polyfit(x_1240,y_1,1)
xr240=polyval([ar,br], x_1240)

plot(x_1240,xr240,label='240')
Z=diff(xr240)
X=diff(x_1240)
Z240=(Z[0]/X[0])

mx=sqrt(max((y_1-xr240)**2))
err=mx
Z1=diff(((xr240[len(xr240)-2]-err),(xr240[1]+err)))
X1=diff((x_1240[len(x_1240)-2]*-1,x_1240[1]*-1))
err1=(Z1[0]/-X1[0])
err240=Z240-err1
###################################	


########################################################################
########################################################################
########################################################################
x=(log(sqrt(T250[:,3])))
y=(log(T250[:,1]))
x2=T250[:,3]
y2=T250[:,1]
x3=x#smooth(x,coeff1)
y3=y#smooth(x,coeff1)
C=[]
for n in range(2,len(y)-1):
	c=((y3[n+1]-y3[n-1])/(x3[n+1]-x3[n-1]))
	C.append(c)
C=array(C)
x1=[]
for n in range(len(C)):
	a=x[n]
	x1.append(a)
x1=array(x1)
d=y2[2]-y2[1]
thPoint=int((th250[0]*.8/d)-1)
thPoint1=int(((th250[0]*0.3)/d)+1)

########## Z-fitting ##############
x_1250=[]
y_1=[]
for n in range(thPoint1,thPoint):
		b=x3[n]
		c=y3[n]
		x_1250.append(b)
		y_1.append(c)
x_1250=array(x_1250)
y_1=array(y_1)
(ar,br)=polyfit(x_1250,y_1,1)
xr250=polyval([ar,br], x_1250)

plot(x_1250,xr250,label='250')
Z=diff(xr250)
X=diff(x_1250)
Z250=(Z[0]/X[0])

mx=sqrt(max((y_1-xr250)**2))
err=mx
Z1=diff(((xr250[len(xr250)-2]-err),(xr250[1]+err)))
X1=diff((x_1250[len(x_1250)-2]*-1,x_1250[1]*-1))
err1=(Z1[0]/-X1[0])
err250=Z250-err1
###################################	


########################################################################
########################################################################
########################################################################
x=(log(sqrt(T260[:,3])))
y=(log(T260[:,1]))
x2=T260[:,3]
y2=T260[:,1]
x3=x#smooth(x,coeff1)
y3=y#smooth(x,coeff1)
C=[]
for n in range(2,len(y)-1):
	c=((y3[n+1]-y3[n-1])/(x3[n+1]-x3[n-1]))
	C.append(c)
C=array(C)
x1=[]
for n in range(len(C)):
	a=x[n]
	x1.append(a)
x1=array(x1)
d=y2[2]-y2[1]
thPoint=int((th260[0]*.8/d)-1)
thPoint1=int(((th260[0]*0.3)/d)+1)

########## Z-fitting ##############
x_1260=[]
y_1=[]
for n in range(thPoint1,thPoint):
		b=x3[n]
		c=y3[n]
		x_1260.append(b)
		y_1.append(c)
x_1260=array(x_1260)
y_1=array(y_1)
(ar,br)=polyfit(x_1260,y_1,1)
xr260=polyval([ar,br], x_1260)

plot(x_1260,xr260,label='260')
Z=diff(xr260)
X=diff(x_1260)
Z260=(Z[0]/X[0])

mx=sqrt(max((y_1-xr260)**2))
err=mx
Z1=diff(((xr260[len(xr260)-2]-err),(xr260[1]+err)))
X1=diff((x_1260[len(x_1260)-2]*-1,x_1260[1]*-1))
err1=(Z1[0]/-X1[0])
err260=Z260-err1
###################################	


########################################################################
########################################################################
########################################################################
x=(log(sqrt(T270[:,3])))
y=(log(T270[:,1]))
x2=T270[:,3]
y2=T270[:,1]
x3=x#smooth(x,coeff1)
y3=y#smooth(x,coeff1)
C=[]
for n in range(2,len(y)-1):
	c=((y3[n+1]-y3[n-1])/(x3[n+1]-x3[n-1]))
	C.append(c)
C=array(C)
x1=[]
for n in range(len(C)):
	a=x[n]
	x1.append(a)
x1=array(x1)
d=y2[2]-y2[1]
thPoint=int((th270[0]*.8/d)-1)
thPoint1=int(((th270[0]*0.3)/d)+1)

########## Z-fitting ##############
x_1270=[]
y_1=[]
for n in range(thPoint1,thPoint):
		b=x3[n]
		c=y3[n]
		x_1270.append(b)
		y_1.append(c)
x_1270=array(x_1270)
y_1=array(y_1)
(ar,br)=polyfit(x_1270,y_1,1)
xr270=polyval([ar,br], x_1270)

plot(x_1270,xr270,label='270')
Z=diff(xr270)
X=diff(x_1270)
Z270=(Z[0]/X[0])

mx=sqrt(max((y_1-xr270)**2))
err=mx
Z1=diff(((xr270[len(xr270)-2]-err),(xr270[1]+err)))
X1=diff((x_1270[len(x_1270)-2]*-1,x_1270[1]*-1))
err1=(Z1[0]/-X1[0])
err270=Z270-err1
###################################	


########################################################################
########################################################################
########################################################################
x=(log(sqrt(T280[:,3])))
y=(log(T280[:,1]))
x2=T280[:,3]
y2=T280[:,1]
x3=x#smooth(x,coeff1)
y3=y#smooth(x,coeff1)
C=[]
for n in range(2,len(y)-1):
	c=((y3[n+1]-y3[n-1])/(x3[n+1]-x3[n-1]))
	C.append(c)
C=array(C)
x1=[]
for n in range(len(C)):
	a=x[n]
	x1.append(a)
x1=array(x1)
d=y2[2]-y2[1]
thPoint=int((th280[0]*.8/d)-1)
thPoint1=int(((th280[0]*0.3)/d)+1)

########## Z-fitting ##############
x_1280=[]
y_1=[]
for n in range(thPoint1,thPoint):
		b=x3[n]
		c=y3[n]
		x_1280.append(b)
		y_1.append(c)
x_1280=array(x_1280)
y_1=array(y_1)
(ar,br)=polyfit(x_1280,y_1,1)
xr280=polyval([ar,br], x_1280)

plot(x_1280,xr280,label='280')
Z=diff(xr280)
X=diff(x_1280)
Z280=(Z[0]/X[0])

mx=sqrt(max((y_1-xr280)**2))
err=mx
Z1=diff(((xr280[len(xr280)-2]-err),(xr280[1]+err)))
X1=diff((x_1280[len(x_1280)-2]*-1,x_1280[1]*-1))
err1=(Z1[0]/-X1[0])
err280=Z280-err1
###################################	


########################################################################
########################################################################
########################################################################
x=(log(sqrt(T290[:,3])))
y=(log(T290[:,1]))
x2=T290[:,3]
y2=T290[:,1]
x3=x#smooth(x,coeff1)
y3=y#smooth(x,coeff1)
C=[]
for n in range(2,len(y)-1):
	c=((y3[n+1]-y3[n-1])/(x3[n+1]-x3[n-1]))
	C.append(c)
C=array(C)
x1=[]
for n in range(len(C)):
	a=x[n]
	x1.append(a)
x1=array(x1)
d=y2[2]-y2[1]
thPoint=int((th290[0]*.8/d)-1)
thPoint1=int(((th290[0]*0.3)/d)+1)

########## Z-fitting ##############
x_1290=[]
y_1=[]
for n in range(thPoint1,thPoint):
		b=x3[n]
		c=y3[n]
		x_1290.append(b)
		y_1.append(c)
x_1290=array(x_1290)
y_1=array(y_1)
(ar,br)=polyfit(x_1290,y_1,1)
xr290=polyval([ar,br], x_1290)

plot(x_1290,xr290,label='290')
Z=diff(xr290)
X=diff(x_1290)
Z290=(Z[0]/X[0])

mx=sqrt(max((y_1-xr290)**2))
err=mx
Z1=diff(((xr290[len(xr290)-2]-err),(xr290[1]+err)))
X1=diff((x_1290[len(x_1290)-2]*-1,x_1290[1]*-1))
err1=(Z1[0]/-X1[0])
err290=Z290-err1
###################################	


########################################################################
########################################################################
########################################################################
x=(log(sqrt(T300[:,3])))
y=(log(T300[:,1]))
x2=T300[:,3]
y2=T300[:,1]
x3=x#smooth(x,coeff1)
y3=y#smooth(x,coeff1)
d=y2[2]-y2[1]
thPoint=int((th300[0]*.8/d)-1)
thPoint1=int(((th300[0]*0.3)/d)+1)

########## Z-fitting ##############
x_1300=[]
y_1=[]
for n in range(thPoint1,thPoint):
		b=x3[n]
		c=y3[n]
		x_1300.append(b)
		y_1.append(c)
x_1300=array(x_1300)
y_1=array(y_1)
(ar,br)=polyfit(x_1300,y_1,1)
xr300=polyval([ar,br], x_1300)

plot(x_1300,xr300,label='300')
Z=diff(xr300)
X=diff(x_1300)
Z300=(Z[0]/X[0])

mx=sqrt(max((y_1-xr300)**2))
err=mx
Z1=diff(((xr300[len(xr300)-2]-err),(xr300[1]+err)))
X1=diff((x_1300[len(x_1300)-2]*-1,x_1300[1]*-1))
err1=(Z1[0]/-X1[0])
err300=Z300-err1
###################################	


########################################################################
########################################################################
########################################################################
x=(log(sqrt(T310[:,3])))
y=(log(T310[:,1]))
x2=T310[:,3]
y2=T310[:,1]
x3=x#smooth(x,coeff1)
y3=y#smooth(x,coeff1)
C=[]
for n in range(2,len(y)-1):
	c=((y3[n+1]-y3[n-1])/(x3[n+1]-x3[n-1]))
	C.append(c)
C=array(C)
x1=[]
for n in range(len(C)):
	a=x[n]
	x1.append(a)
x1=array(x1)
d=y2[2]-y2[1]
thPoint=int((th310[0]*.8/d)-1)
thPoint1=int(((th310[0]*0.3)/d)+1)

########## Z-fitting ##############
x_1310=[]
y_1=[]
for n in range(thPoint1,thPoint):
		b=x3[n]
		c=y3[n]
		x_1310.append(b)
		y_1.append(c)
x_1310=array(x_1310)
y_1=array(y_1)
(ar,br)=polyfit(x_1310,y_1,1)
xr310=polyval([ar,br], x_1310)


plot(x_1310,xr310,label='310')
Z=diff(xr310)
X=diff(x_1310)
Z310=(Z[0]/X[0])

mx=sqrt(max((y_1-xr310)**2))
err=mx
Z1=diff(((xr310[len(xr310)-2]-err),(xr310[1]+err)))
X1=diff((x_1310[len(x_1310)-2]*-1,x_1310[1]*-1))
err1=(Z1[0]/-X1[0])
err310=Z310-err1
###################################	


########################################################################
########################################################################
########################################################################
x=(log(sqrt(T320[:,3])))
y=(log(T320[:,1]))
x2=T320[:,3]
y2=T320[:,1]
x3=x#smooth(x,coeff1)
y3=y#smooth(x,coeff1)
C=[]
for n in range(2,len(y)-1):
	c=((y3[n+1]-y3[n-1])/(x3[n+1]-x3[n-1]))
	C.append(c)
C=array(C)
x1=[]
for n in range(len(C)):
	a=x[n]
	x1.append(a)
x1=array(x1)
d=y2[2]-y2[1]
thPoint=int((th320[0]*.8/d)-1)
thPoint1=int(((th320[0]*0.3)/d)+1)

########## Z-fitting ##############
x_1320=[]
y_1=[]
for n in range(thPoint1,thPoint):
		b=x3[n]
		c=y3[n]
		x_1320.append(b)
		y_1.append(c)
x_1320=array(x_1320)
y_1=array(y_1)
(ar,br)=polyfit(x_1320,y_1,1)
xr320=polyval([ar,br], x_1320)

plot(x_1320,xr320,label='320')
Z=diff(xr320)
X=diff(x_1320)
Z320=(Z[0]/X[0])
mx=sqrt(max((y_1-xr320)**2))
err=mx
Z1=diff(((xr320[len(xr320)-2]-err),(xr320[1]+err)))
X1=diff((x_1320[len(x_1320)-2]*-1,x_1320[1]*-1))
err1=(Z1[0]/-X1[0])
err320=Z320-err1
###################################	


########################################################################
########################################################################
########################################################################

x=(log(sqrt(T330[:,3])))
y=(log(T330[:,1]))
x2=T330[:,3]
y2=T330[:,1]
x3=x#smooth(x,coeff1)
y3=y#smooth(x,coeff1)
C=[]
for n in range(2,len(y)-1):
	c=((y3[n+1]-y3[n-1])/(x3[n+1]-x3[n-1]))
	C.append(c)
C=array(C)
x1=[]
for n in range(len(C)):
	a=x[n]
	x1.append(a)
x1=array(x1)
d=y2[2]-y2[1]
thPoint=int((th330[0]*.8/d)-1)
thPoint1=int(((th330[0]*0.3)/d)+1)

########## Z-fitting ##############
x_1330=[]
y_1=[]
for n in range(thPoint1,thPoint):
		b=x3[n]
		c=y3[n]
		x_1330.append(b)
		y_1.append(c)
x_1330=array(x_1330)
y_1=array(y_1)
(ar,br)=polyfit(x_1330,y_1,1)
xr330=polyval([ar,br], x_1330)



plot(x_1330,xr330,label='330')
Z=diff(xr330)
X=diff(x_1330)
Z330=(Z[0]/X[0])

mx=sqrt(max((y_1-xr330)**2))
err=mx
Z1=diff(((xr330[len(xr330)-2]-err),(xr330[1]+err)))
X1=diff((x_1330[len(x_1330)-2]*-1,x_1330[1]*-1))
err1=(Z1[0]/-X1[0])
err330=Z330-err1
###################################	


########################################################################
########################################################################
########################################################################
x=(log(sqrt(T340[:,3])))
y=(log(T340[:,1]))
x2=T340[:,3]
y2=T340[:,1]
x3=x#smooth(x,coeff1)
y3=y#smooth(x,coeff1)
d=y2[2]-y2[1]
thPoint=int((th340[0]*.8/d)-1)
thPoint1=int(((th340[0]*0.3)/d)+1)
				

########## Z-fitting ##############
x_1340=[]
y_1=[]
for n in range(thPoint1,thPoint):
		b=x3[n]
		c=y3[n]
		x_1340.append(b)
		y_1.append(c)
x_1340=array(x_1340)
y_1=array(y_1)
(ar,br)=polyfit(x_1340,y_1,1)
xr340=polyval([ar,br], x_1340)

plot(x_1340,xr340,label='340')

Z=diff(xr340)
X=diff(x_1340)
Z340=(Z[0]/X[0])

mx=sqrt(max((y_1-xr340)**2))
err=mx
Z1=diff(((xr340[len(xr340)-2]-err),(xr340[1]+err)))
X1=diff((x_1340[len(x_1340)-2]*-1,x_1340[1]*-1))
err1=(Z1[0]/-X1[0])
err340=Z340-err1

xlabel('ln sqrt SE')
ylabel('ln I')
title('Z = gradient')
legend(loc='upper right', ncol=1, shadow=True, numpoints=1)
show()
###################################		

########################################################################
########################################################################
########################################################################

plot(x1_60, (y60/x60),'o',label='Z Vs I @ 60K')
xlabel('Current I')
ylabel('Z')
title('Z Vs I')
ylim([.5,4])
legend(loc='lower right', ncol=1, shadow=True, numpoints=1)
show()
plot(x1_100, (y100/x100),'o',label='Z Vs I @ 100K')
xlabel('Current I')
ylabel('Z')
title('Z Vs I')
ylim([.5,4])
legend(loc='lower right', ncol=1, shadow=True, numpoints=1)
show()
plot(x1_200, (y200/x200),'o',label='Z Vs I @ 200K')
xlabel('Current I')
ylabel('Z')
title('Z Vs I')
ylim([.5,4])
legend(loc='lower right', ncol=1, shadow=True, numpoints=1)
show()

subplot(1,3,1)
ithu=[th60[0],th80[0],th100[0],th120[0],th140[0],th160[0],th180[0],th200[0],th220[0],th240[0],th260[0],th280[0],th300[0],th320[0],th340[0]]
Tu=[60,80,100,120,140,160,180,200,220,240,260,280,300,320, 340]
ithd=[th50[0],th70[0],th90[0],th110[0],th130[0],th150[0],th170[0],th190[0],th210[0],th230[0],th250[0],th270[0],th290[0],th310[0],th330[0]]
Td=[50,70,90,110,130,150,170,190,210,230,250,270,290,310,330]
norm=(th90[0]/Irad90)
Iradu=[Irad60*norm,Irad80*norm,Irad100*norm,Irad120*norm,Irad140*norm,Irad160*norm,Irad180*norm,Irad200*norm,Irad220*norm,Irad240*norm,Irad260*norm,Irad280*norm,Irad300*norm,Irad320*norm,Irad340*norm]
Iradd=[Irad50*norm,Irad70*norm,Irad90*norm,Irad110*norm,Irad130*norm,Irad150*norm,Irad170*norm,Irad190*norm,Irad210*norm,Irad230*norm,Irad250*norm,Irad270*norm,Irad290*norm,Irad310*norm,Irad330*norm]
semilogy(Tu,ithu,'o''r',label='I(th) up')
semilogy(Td,ithd,'o''b',label='I(th) down')
semilogy(Tu,Iradu,'^''y',label='I(rad) up')
semilogy(Td,Iradd,'^''k',label='I(rad) down')
xlabel('$Temperature$ $(K)$')
ylabel('$Threshold$ $current$ $(mA)$')
title('$Ith$ $Vs$ $T$')
legend(loc='lower right', ncol=1, shadow=True, numpoints=1)
T0u=1/((log(ithu[14])-log(ithu[0]))/(Tu[14]-Tu[0]))
T0d=1/((log(ithd[14])-log(ithd[0]))/(Td[14]-Td[0]))
print 'T0u =',T0u,'K',        'T0d =', T0d,'K'
xlim([0,380])
xticks([100,200,300])
#show()
subplot(1,3,2)
Tu=[60,80,100,120,140,160,180,200,220,240,260,280,300,320, 340]
Td=[50,70,90,110,130,150,170,190,210,230,250,270,290,310,330]
Neu=[Ne60,Ne80,Ne100,Ne120,Ne140,Ne160,Ne180,Ne200,Ne220,Ne240,Ne260,Ne280,Ne300,Ne320,Ne340]
Ned=[Ne50,Ne70,Ne90,Ne110,Ne130,Ne150,Ne170,Ne190,Ne210,Ne230,Ne250,Ne270,Ne290,Ne310,Ne330]
mx=max(Ned)
ErrNu=[Nerr60,Nerr80,Nerr100,Nerr120,Nerr140,Nerr160,Nerr180,Nerr200,Nerr220,Nerr240,Nerr260,Nerr280,Nerr300,Nerr320,Nerr340]
ErrNd=[Nerr50,Nerr70,Nerr90,Nerr110,Nerr130,Nerr150,Nerr170,Nerr190,Nerr210,Nerr230,Nerr250,Nerr270,Nerr290,Nerr310,Nerr330]
plot(Tu,Neu*1/mx,'o''r', label='up')
errorbar(Tu,Neu*1/mx,yerr=ErrNu*1/mx,xerr=None,fmt='.''r')
plot(Td,Ned*1/mx,'o''b', label='down')
errorbar(Td,Ned*1/mx,yerr=ErrNd*1/mx,xerr=None,fmt='.''b')
xlabel('$Temperature$ $(K)$')
ylabel('$\eta^{d}_{ext}$')
#yticks([])
title('$\eta^{d}_{ext}$ $Vs$ $T$')
legend(loc='upper right', ncol=1, shadow=True, numpoints=1)
xlim([0,380])
xticks([100,200,300])
#show()

subplot(1,3,3)
Tu=[60,80,100,120,140,160,180,200,220,240,260,280,300,320, 340]
Td=[50,70,90,110,130,150,170,190,210,230,250,270,290,310,330]
Zu=[Z60,Z80,Z100,Z120,Z140,Z160,Z180,Z200,Z220,Z240,Z260,Z280,Z300,Z320, Z340]
Zd=[Z50,Z70,Z90,Z110,Z130,Z150,Z170,Z190,Z210,Z230,Z250,Z270,Z290,Z310,Z330]
ErrZu=[err60,err80,err100,err120,err140,err160,err180,err200,err220,err240,err260,err280,err300,err320,err340]
ErrZd=[err50,err70,err90,err110,err130,err150,err170,err190,err210,err230,err250,err270,err290,err310,err330]
plot(Tu,Zu,'o''r', label='up')
errorbar(Tu,Zu,yerr=ErrZu,xerr=None,fmt='.''r')
plot(Td,Zd,'o''b', label='down')
errorbar(Td,Zd,yerr=ErrZd,xerr=None,fmt='.''b')
y=(2,2)
y1=(3,3)
x=(0,380)
plot(x,y,'--''k')
plot(x,y1,'--''k')
xlabel('$Temperature$ $(K)$')
ylabel('$Z$')
title('$Z$ $Vs$ $T$')
legend(loc='lower right', ncol=1, shadow=True, numpoints=1)
xlim([0,380])
ylim([1,4])
xticks([100,200,300])
show()

T01=1/((log(ithd[0])-log(ithu[0]))/(Td[0]-Tu[0]))
T02=1/((log(ithd[1])-log(ithu[1]))/(Td[1]-Tu[1]))
T03=1/((log(ithd[2])-log(ithu[2]))/(Td[2]-Tu[2]))
T04=1/((log(ithd[3])-log(ithu[3]))/(Td[3]-Tu[3]))
T05=1/((log(ithd[4])-log(ithu[4]))/(Td[4]-Tu[4]))
T06=1/((log(ithd[5])-log(ithu[5]))/(Td[5]-Tu[5]))
T07=1/((log(ithd[6])-log(ithu[6]))/(Td[6]-Tu[6]))
T08=1/((log(ithd[7])-log(ithu[7]))/(Td[7]-Tu[7]))
T09=1/((log(ithd[8])-log(ithu[8]))/(Td[8]-Tu[8]))
T010=1/((log(ithd[9])-log(ithu[9]))/(Td[9]-Tu[9]))
T011=1/((log(ithd[10])-log(ithu[10]))/(Td[10]-Tu[10]))
T012=1/((log(ithd[11])-log(ithu[11]))/(Td[11]-Tu[11]))
T013=1/((log(ithd[12])-log(ithu[12]))/(Td[12]-Tu[12]))
T014=1/((log(ithd[13])-log(ithu[13]))/(Td[13]-Tu[13]))

T=[55,75,95,115,135,155,175,195,215,235,255,275,295,315]
T0=[T01,T02,T03,T04,T05,T06,T07,T08,T09,T010,T011,T012,T013,T014]
plot(T,T0,'o',label='T0 Vs T')
xlabel('Temperature (K)')
ylabel('T0')
title('T0 Vs T')
legend(loc='lower right', ncol=1, shadow=True, numpoints=1)
show()





T=[50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200,210,220,230,240,250,260,270,280,290,300,310,320,330,340,350]
SOV=[SO50[0],SO60[0],SO70[0],SO80[0],SO90[0],SO100[0],SO110[0],SO120[0],SO130[0],SO140[0],SO150[0],SO160[0],SO170[0],SO180[0],SO190[0],SO200[0],SO210[0],SO220[0],SO230[0],SO240[0],SO250[0],SO260[0],SO270[0],SO280[0],SO290[0],SO300[0],SO310[0],SO320[0],SO330[0],SO340[0],SO350[0]]
plot(T,SOV,'o',label='V Vs T')
xlabel('Temperature (K)')
ylabel('V')
title('switch on voltage Vs Temperature')
legend(loc='lower right', ncol=1, shadow=True, numpoints=1)
show()
thp50=(th50[0]/(T50[1,1]-T50[0,1]))+1
thp60=(th60[0]/(T60[1,1]-T60[0,1]))+1
thp70=(th70[0]/(T70[1,1]-T70[0,1]))+1
thp80=(th80[0]/(T80[1,1]-T80[0,1]))+1
thp90=(th90[0]/(T90[1,1]-T90[0,1]))+1
thp100=(th100[0]/(T100[1,1]-T100[0,1]))+1
thp110=(th110[0]/(T110[1,1]-T110[0,1]))+1
thp120=(th120[0]/(T120[1,1]-T120[0,1]))+1
thp130=(th130[0]/(T130[1,1]-T130[0,1]))+1
thp140=(th140[0]/(T140[1,1]-T140[0,1]))+1
thp150=(th150[0]/(T150[1,1]-T150[0,1]))+1
thp160=(th160[0]/(T160[1,1]-T160[0,1]))+1
thp170=(th170[0]/(T170[1,1]-T170[0,1]))+1
thp180=(th180[0]/(T180[1,1]-T180[0,1]))+1
thp190=(th190[0]/(T190[1,1]-T190[0,1]))+1
thp200=(th200[0]/(T200[1,1]-T200[0,1]))+1
thp210=(th210[0]/(T210[1,1]-T210[0,1]))+1
thp220=(th220[0]/(T220[1,1]-T220[0,1]))+1
thp230=(th230[0]/(T230[1,1]-T230[0,1]))+1
thp240=(th240[0]/(T240[1,1]-T240[0,1]))+1
thp250=(th250[0]/(T250[1,1]-T250[0,1]))+1
thp260=(th260[0]/(T260[1,1]-T260[0,1]))+1
thp270=(th270[0]/(T270[1,1]-T270[0,1]))+1
thp280=(th280[0]/(T280[1,1]-T280[0,1]))+1
thp290=(th290[0]/(T290[1,1]-T290[0,1]))+1
thp300=(th300[0]/(T300[1,1]-T300[0,1]))+1
thp310=(th310[0]/(T310[1,1]-T310[0,1]))+1
thp320=(th320[0]/(T320[1,1]-T320[0,1]))+1
thp330=(th330[0]/(T330[1,1]-T330[0,1]))+1
thp340=(th340[0]/(T340[1,1]-T340[0,1]))+1
subplot(2,1,1)
plot(T60[:,1],T60[:,2],label=('T=60K'))
plot(T80[:,1],T80[:,2],label=('T=80K'))
plot(T100[:,1],T100[:,2],label=('T=100K'))
plot(T120[:,1],T120[:,2],label=('T=120K'))
plot(T140[:,1],T140[:,2],label=('T=140K'))
plot(T160[:,1],T160[:,2],label=('T=160K'))
plot(T180[:,1],T180[:,2],label=('T=180K'))
plot(T200[:,1],T200[:,2],label=('T=200K'))
plot(T220[:,1],T220[:,2],label=('T=220K'))
plot(T240[:,1],T240[:,2],label=('T=240K'))
plot(T260[:,1],T260[:,2],label=('T=260K'))
plot(T280[:,1],T280[:,2],label=('T=280K'))
plot(T300[:,1],T300[:,2],label=('T=300K'))
plot(T320[:,1],T320[:,2],label=('T=320K'))
plot(T340[:,1],T340[:,2],label=('T=340K'))
plot(T350[:,1],T350[:,2],label=('T=350K'))
plot(T330[:,1],T330[:,2],label=('T=330K'))
plot(T310[:,1],T310[:,2],label=('T=310K'))
plot(T290[:,1],T290[:,2],label=('T=290K'))
plot(T270[:,1],T270[:,2],label=('T=270K'))
plot(T250[:,1],T250[:,2],label=('T=250K'))
plot(T230[:,1],T230[:,2],label=('T=230K'))
plot(T210[:,1],T210[:,2],label=('T=210K'))
plot(T190[:,1],T190[:,2],label=('T=190K'))
plot(T170[:,1],T170[:,2],label=('T=170K'))
plot(T150[:,1],T150[:,2],label=('T=150K'))
plot(T130[:,1],T130[:,2],label=('T=130K'))
plot(T110[:,1],T110[:,2],label=('T=110K'))
plot(T90[:,1],T90[:,2],label=('T=90K'))
plot(T70[:,1],T70[:,2],label=('T=70K'))
plot(T50[:,1],T50[:,2],label=('T=50K'))
plot(th50[0],T50[thp50,2],'o')
plot(th60[0],T60[thp60,2],'o')
plot(th70[0],T70[thp70,2],'o')
plot(th80[0],T80[thp80,2],'o')
plot(th90[0],T90[thp90,2],'o')
plot(th100[0],T100[thp100,2],'o')
plot(th110[0],T110[thp110,2],'o')
plot(th120[0],T120[thp120,2],'o')
plot(th130[0],T130[thp130,2],'o')
plot(th140[0],T140[thp140,2],'o')
plot(th150[0],T150[thp150,2],'o')
plot(th160[0],T160[thp160,2],'o')
plot(th170[0],T170[thp170,2],'o')
plot(th180[0],T180[thp180,2],'o')
plot(th190[0],T190[thp190,2],'o')
plot(th200[0],T200[thp200,2],'o')
plot(th210[0],T210[thp210,2],'o')
plot(th220[0],T220[thp220,2],'o')
plot(th230[0],T230[thp230,2],'o')
plot(th240[0],T240[thp240,2],'o')
plot(th250[0],T250[thp250,2],'o')
plot(th260[0],T260[thp260,2],'o')
plot(th270[0],T270[thp270,2],'o')
plot(th280[0],T280[thp280,2],'o')
plot(th290[0],T290[thp290,2],'o')
plot(th300[0],T300[thp300,2],'o')
plot(th310[0],T310[thp310,2],'o')
plot(th320[0],T320[thp320,2],'o')
plot(th330[0],T330[thp330,2],'o')
plot(th340[0],T340[thp340,2],'o')
xlabel('Current (mA)')
ylabel('Light (a.u.)')
title('L-I Characteristic (Facet)')
legend(loc='lower right', ncol=1, shadow=True, numpoints=1)
#show()
###################################################################
subplot(2,1,2)
plot(T60[:,1],T60[:,3],'--',label=('T=60K'))
plot(T80[:,1],T80[:,3],'--',label=('T=80K'))
plot(T100[:,1],T100[:,3],'--',label=('T=100K'))
plot(T120[:,1],T120[:,3],'--',label=('T=120K'))
plot(T140[:,1],T140[:,3],'--',label=('T=140K'))
plot(T160[:,1],T160[:,3],'--',label=('T=160K'))
plot(T180[:,1],T180[:,3],'--',label=('T=180K'))
plot(T200[:,1],T200[:,3],'--',label=('T=200K'))
plot(T220[:,1],T220[:,3],'--',label=('T=220K'))
plot(T240[:,1],T240[:,3],'--',label=('T=240K'))
plot(T260[:,1],T260[:,3],'--',label=('T=260K'))
plot(T280[:,1],T280[:,3],'--',label=('T=280K'))
plot(T300[:,1],T300[:,3],'--',label=('T=300K'))
plot(T320[:,1],T320[:,3],'--',label=('T=320K'))
plot(T340[:,1],T340[:,3],'--',label=('T=340K'))
plot(T350[:,1],T350[:,3],label=('T=350K'))
plot(T330[:,1],T330[:,3],label=('T=330K'))
plot(T310[:,1],T310[:,3],label=('T=310K'))
plot(T290[:,1],T290[:,3],label=('T=290K'))
plot(T270[:,1],T270[:,3],label=('T=270K'))
plot(T250[:,1],T250[:,3],label=('T=250K'))
plot(T230[:,1],T230[:,3],label=('T=230K'))
plot(T210[:,1],T210[:,3],label=('T=210K'))
plot(T190[:,1],T190[:,3],label=('T=190K'))
plot(T170[:,1],T170[:,3],label=('T=170K'))
plot(T150[:,1],T150[:,3],label=('T=150K'))
plot(T130[:,1],T130[:,3],label=('T=130K'))
plot(T110[:,1],T110[:,3],label=('T=110K'))
plot(T90[:,1],T90[:,3],label=('T=90K'))
plot(T70[:,1],T70[:,3],label=('T=70K'))
plot(T50[:,1],T50[:,3],label=('T=50K'))
plot(th50[0],T50[thp50,3],'o')
plot(th60[0],T60[thp60,3],'o')
plot(th70[0],T70[thp70,3],'o')
plot(th80[0],T80[thp80,3],'o')
plot(th90[0],T90[thp90,3],'o')
plot(th100[0],T100[thp100,3],'o')
plot(th110[0],T110[thp110,3],'o')
plot(th120[0],T120[thp120,3],'o')
plot(th130[0],T130[thp130,3],'o')
plot(th140[0],T140[thp140,3],'o')
plot(th150[0],T150[thp150,3],'o')
plot(th160[0],T160[thp160,3],'o')
plot(th170[0],T170[thp170,3],'o')
plot(th180[0],T180[thp180,3],'o')
plot(th190[0],T190[thp190,3],'o')
plot(th200[0],T200[thp200,3],'o')
plot(th210[0],T210[thp210,3],'o')
plot(th220[0],T220[thp220,3],'o')
plot(th230[0],T230[thp230,3],'o')
plot(th240[0],T240[thp240,3],'o')
plot(th250[0],T250[thp250,3],'o')
plot(th260[0],T260[thp260,3],'o')
plot(th270[0],T270[thp270,3],'o')
plot(th280[0],T280[thp280,3],'o')
plot(th290[0],T290[thp290,3],'o')
plot(th300[0],T300[thp300,3],'o')
plot(th310[0],T310[thp310,3],'o')
plot(th320[0],T320[thp320,3],'o')
plot(th330[0],T330[thp330,3],'o')
plot(th340[0],T340[thp340,3],'o')
xlabel('Current (mA)')
ylabel('SE (a.u.)')
title('L-I Characteristic (window)')
legend(loc='lower right', ncol=1, shadow=True, numpoints=1)
show()






















#### 3D plot ######
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt
mpl.rcParams['legend.fontsize'] = 10
z50=[]
for n in range(len(x_150)):
	z=50
	z50.append(z)
z50=array(z50)

z60=[]
for n in range(len(x_160)):
	z=60
	z60.append(z)
z60=array(z60)

z70=[]
for n in range(len(x_170)):
	z=70
	z70.append(z)
z70=array(z70)

z80=[]
for n in range(len(x_180)):
	z=80
	z80.append(z)
z80=array(z80)

z90=[]
for n in range(len(x_190)):
	z=90
	z90.append(z)
z90=array(z90)

z100=[]
for n in range(len(x_1100)):
	z=100
	z100.append(z)
z100=array(z100)

z110=[]
for n in range(len(x_1110)):
	z=110
	z110.append(z)
z110=array(z110)

z120=[]
for n in range(len(x_1120)):
	z=120
	z120.append(z)
z120=array(z120)

z130=[]
for n in range(len(x_1130)):
	z=130
	z130.append(z)
z130=array(z130)

z140=[]
for n in range(len(x_1140)):
	z=140
	z140.append(z)
z140=array(z140)

z150=[]
for n in range(len(x_1150)):
	z=150
	z150.append(z)
z150=array(z150)

z160=[]
for n in range(len(x_1160)):
	z=160
	z160.append(z)
z160=array(z160)

z170=[]
for n in range(len(x_1170)):
	z=170
	z170.append(z)
z170=array(z170)

z180=[]
for n in range(len(x_1180)):
	z=180
	z180.append(z)
z180=array(z180)

z190=[]
for n in range(len(x_1190)):
	z=190
	z190.append(z)
z190=array(z190)

z200=[]
for n in range(len(x_1200)):
	z=200
	z200.append(z)
z200=array(z200)

z210=[]
for n in range(len(x_1210)):
	z=210
	z210.append(z)
z210=array(z210)

z220=[]
for n in range(len(x_1220)):
	z=220
	z220.append(z)
z220=array(z220)

z230=[]
for n in range(len(x_1230)):
	z=230
	z230.append(z)
z230=array(z230)

z240=[]
for n in range(len(x_1240)):
	z=240
	z240.append(z)
z240=array(z240)

z250=[]
for n in range(len(x_1250)):
	z=250
	z250.append(z)
z250=array(z250)

z260=[]
for n in range(len(x_1260)):
	z=260
	z260.append(z)
z260=array(z260)

z270=[]
for n in range(len(x_1270)):
	z=270
	z270.append(z)
z270=array(z270)

z280=[]
for n in range(len(x_1280)):
	z=280
	z280.append(z)
z280=array(z280)

z290=[]
for n in range(len(x_1290)):
	z=290
	z290.append(z)
z290=array(z290)

z300=[]
for n in range(len(x_1300)):
	z=300
	z300.append(z)
z300=array(z300)

z310=[]
for n in range(len(x_1310)):
	z=310
	z310.append(z)
z310=array(z310)

z320=[]
for n in range(len(x_1320)):
	z=320
	z320.append(z)
z320=array(z320)

z330=[]
for n in range(len(x_1330)):
	z=330
	z330.append(z)
z330=array(z330)

z340=[]
for n in range(len(x_1340)):
	z=340
	z340.append(z)
z340=array(z340)

fig = plt.figure()
ax = fig.gca(projection='3d')
x=np.concatenate((x_150,x_160,x_170,x_180,x_190,x_1100,x_1110,x_1120,x_1130,x_1140,x_1150,x_1160,x_1170,x_1180,x_1190,x_1200,x_1210,x_1220,x_1230,x_1240,x_1250,x_1260,x_1270,x_1280,x_1290,x_1300,x_1310,x_1320,x_1330,x_1340), axis=0)
z=np.concatenate((xr50,xr60,xr70,xr80,xr90,xr100,xr110,xr120,xr130,xr140,xr150,xr160,xr170,xr180,xr190,xr200,xr210,xr220,xr230,xr240,xr250,xr260,xr270,xr280,xr290,xr300,xr310,xr320,xr330,xr340),axis=0)
y=np.concatenate((z50,z60,z70,z80,z90,z100,z110,z120,z130,z140,z150,z160,z170,z180,z190,z200,z210,z220,z230,z240,z250,z260,z270,z280,z290,z300,z310,z320,z330,z340),axis=0)


x1=[x_150[len(x_150)-1],x_160[len(x_160)-1],x_170[len(x_170)-1],x_180[len(x_180)-1],x_190[len(x_190)-1],x_1100[len(x_1100)-1],x_1110[len(x_1110)-1],x_1120[len(x_1120)-1],x_1130[len(x_1130)-1],x_1140[len(x_1140)-1],x_1150[len(x_1150)-1],x_1160[len(x_1160)-1],x_1170[len(x_1170)-1],x_1180[len(x_1180)-1],x_1190[len(x_1190)-1],x_1200[len(x_1200)-1],x_1210[len(x_1210)-1],x_1220[len(x_1220)-1],x_1230[len(x_1230)-1],x_1240[len(x_1240)-1],x_1250[len(x_1250)-1],x_1260[len(x_1260)-1],x_1270[len(x_1270)-1],x_1280[len(x_1280)-1],x_1290[len(x_1290)-1],x_1300[len(x_1300)-1],x_1310[len(x_1310)-1],x_1320[len(x_1320)-1],x_1330[len(x_1330)-1],x_1340[len(x_1340)-1]]
z1=[xr50[len(xr50)-1],xr60[len(xr60)-1],xr70[len(xr70)-1],xr80[len(xr80)-1],xr90[len(xr90)-1],xr100[len(xr100)-1],xr110[len(xr110)-1],xr120[len(xr120)-1],xr130[len(xr130)-1],xr140[len(xr140)-1],xr150[len(xr150)-1],xr160[len(xr160)-1],xr170[len(xr170)-1],xr180[len(xr180)-1],xr190[len(xr190)-1],xr200[len(xr200)-1],xr210[len(xr210)-1],xr220[len(xr220)-1],xr230[len(xr230)-1],xr240[len(xr240)-1],xr250[len(xr250)-1],xr260[len(xr260)-1],xr270[len(xr270)-1],xr280[len(xr280)-1],xr290[len(xr290)-1],xr300[len(xr300)-1],xr310[len(xr310)-1],xr320[len(xr320)-1],xr330[len(xr330)-1],xr340[len(xr340)-1]]
y1=[z50[len(z50)-1],z60[len(z60)-1],z70[len(z70)-1],z80[len(z80)-1],z90[len(z90)-1],z100[len(z100)-1],z110[len(z110)-1],z120[len(z120)-1],z130[len(z130)-1],z140[len(z140)-1],z150[len(z150)-1],z160[len(z160)-1],z170[len(z170)-1],z180[len(z180)-1],z190[len(z190)-1],z200[len(z200)-1],z210[len(z210)-1],z220[len(z220)-1],z230[len(z230)-1],z240[len(z240)-1],z250[len(z250)-1],z260[len(z260)-1],z270[len(z270)-1],z280[len(z280)-1],z290[len(z290)-1],z300[len(z300)-1],z310[len(z310)-1],z320[len(z320)-1],z330[len(z330)-1],z340[len(z340)-1]]

ci = (np.linspace(min(z), max(z),4))/340

 
ax.plot(x, y, z, '.')
ax.plot(x1, y1, z1,'o''r',label='I(th)')
ax.plot(x1, y1, z1,'b')
ax.set_xlabel('ln sqrt(SE)')
ax.set_zlabel('ln (I)')
ax.set_ylabel('Temperature (K)')
ax.legend()

plt.show()


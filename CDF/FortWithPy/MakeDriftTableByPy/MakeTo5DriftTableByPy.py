def opt3Cfunc(r,a,b,er):
    from scipy.optimize import minimize

    def func(x):
        return x[0]**2 \
            + (x[1]-1.0)**2*(1.0-r) + x[1]**2*r \
           + (x[2]-1.0)**2*r + x[2]**2*r

    cons = ({'type': 'eq',
         'fun' : lambda x: -x[0]+x[2]-r},
        {'type': 'eq',
         'fun' : lambda x: x[0]+x[1]+x[2]-1.0},
         {'type': 'eq',
         'fun' : lambda x: a*x[0]+x[1]+b*x[2]-er})

    res = minimize(func, [-0.0001, 1.000, 0.0001],
                  constraints=cons)
    dmy=0.0
    print(res.x[0],res.x[1],res.x[2],dmy,dmy)
    fout.write(str(res.x[0])+" ")
    fout.write(str(res.x[1])+" ")
    fout.write(str(res.x[2])+" ")
    fout.write(str(dmy)+" ")
    fout.write(str(dmy)+"\n")
    
def opt5Cfunc(r,a,b,c,d,er):
    from scipy.optimize import minimize

    def func(x):
        return x[0]**2 + x[1]**2 \
            + (x[2]-1.0)**2*(1.0-r) + x[2]**2*r \
           + (x[3]-1.0)**2*r + x[3]**2*r + x[4]**2

    cons = ({'type': 'eq',
         'fun' : lambda x: -2.0*x[0]-x[1]+x[3]+2.0*x[4]-r},
        {'type': 'eq',
         'fun' : lambda x: x[0]+x[1]+x[2]+x[3]+x[4]-1.0},
         {'type': 'eq',
         'fun' : lambda x: a*x[0]+b*x[1]+x[2]+c*x[3]+d*x[4]-er})

    res = minimize(func, [-0.0001, -0.0001, 1.000, 0.0001, 0.0001],
                  constraints=cons)

    print(res.x[0],res.x[1],res.x[2],res.x[3],res.x[4])
    fout.write(str(res.x[0])+" ")
    fout.write(str(res.x[1])+" ")
    fout.write(str(res.x[2])+" ")
    fout.write(str(res.x[3])+" ")
    fout.write(str(res.x[4])+"\n")    

##### Main #############################################################
import math
fout = open("PreDriftTable.txt", "w")

f = open("CoefForOpt.txt", "r")

Ef = float( input('Ef [kV/cm] >> '))
Ef = Ef*1000.0     #   1kV/cm = 1000V/cm
line = f.readline()
data = line.split()
dkxUnit = float(data[0])
rUnit = float(data[1])
hbar = float(data[2])
mdos = float(data[3])
alpha = float(data[4])
nkx = int(data[5])
i_cell_max = int(data[6])

dkx = dkxUnit*Ef
r = rUnit*Ef
print (dkxUnit,rUnit,dkx,r,nkx,i_cell_max)
fout.write(str(Ef/1000.0)+"\n")
fout.write(str(nkx)+"  ")
fout.write(str(i_cell_max)+"\n")
     
for i in range(1,i_cell_max+1):   
    line = f.readline()
    data = line.split()
    i_cell=int(data[0])
    ikx=int(data[1])
    kx=float(data[2])
    ky=float(data[3])
    kz=float(data[4])
    eee=float(data[5])
    a=float(data[6])/eee
    b=float(data[7])/eee
    c=float(data[8])/eee
    d=float(data[9])/eee
    print (i_cell,ikx,r,a,b,c,d)
    
    if ikx == 1 :
        p0=1.0-r
        p1=r
        p2=0.0
        p3=0.0
        p4=0.0
        print (i_cell)
        print (p0,p1,p2,p3,p4)
        fout.write(str(i_cell)+" ")
        fout.write(str(ikx)+" ")
        fout.write(str(p0)+" ")
        fout.write(str(p1)+" ")
        fout.write(str(p2)+" ")
        fout.write(str(p3)+" ")
        fout.write(str(p4)+"\n")
    elif ikx == 2 :
        print (i_cell)
        fout.write(str(i_cell)+" ")
        fout.write(str(ikx)+" ")
        gk = (hbar*hbar/(2.0*mdos))*( (kx+dkx)**2 + ky**2 + kz**2 )
        er = ( 2.0*gk/(math.sqrt(1.0+4.0*alpha*gk)+1.0) )/eee
#        print (er)
        opt3Cfunc(r,a,b,er)
    elif 2 < ikx < nkx-1 :
        print (i_cell)
        fout.write(str(i_cell)+" ")
        fout.write(str(ikx)+" ")
        gk = (hbar*hbar/(2.0*mdos))*( (kx+dkx)**2 + ky**2 + kz**2 )
        er = ( 2.0*gk/(math.sqrt(1.0+4.0*alpha*gk)+1.0) )/eee
#        print (er)
        opt5Cfunc(r,a,b,c,d,er)
    elif ikx == nkx-1 :
        print (i_cell)
        fout.write(str(i_cell)+" ")
        fout.write(str(ikx)+" ")
        gk = (hbar*hbar/(2.0*mdos))*( (kx+dkx)**2 + ky**2 + kz**2 )
        er = ( 2.0*gk/(math.sqrt(1.0+4.0*alpha*gk)+1.0) )/eee
#        print (er)
        opt3Cfunc(r,a,b,er)        
    elif ikx == nkx :
        p0=1.0
        p1=0.0
        p2=0.0
        p3=0.0
        p4=0.0
        print (i_cell)
        print (p0,p1,p2,p3,p4)
        fout.write(str(i_cell)+" ")
        fout.write(str(ikx)+" ")
        fout.write(str(p0)+" ")
        fout.write(str(p1)+" ")
        fout.write(str(p2)+" ")
        fout.write(str(p3)+" ")
        fout.write(str(p4)+"\n")
        
f.close()  

fout.close() # result.txtファイルを閉じる
#!/usr/bin/python
def fatphonon():
    f=open("band.yaml")
    nq=[]
    na=[]
    fff=[]
    qpoint=[]
    distance=[]
    frequency=[]
    eigen=[]
    eigenvx=[]
    eigenvy=[]
    eigenvz=[]
    ieigenvx=[]
    ieigenvy=[]
    ieigenvz=[]
    i=0
    for line in f:
        fff.append(line.split())
        if 'nqpoint' in line:
            nq.append(line.split())
        if 'natom' in line:
            na.append(line.split())
        if 'q-position' in line:
            qpoint.append(line.split())
        if 'distance' in line:
            distance.append(line.split())
        if 'frequency' in line:
            frequency.append(line.split())
        if 'eigenvector' in line:
            eigen.append(i)
        i=i+1
    f.close()
    #print nq[0][1]
    #print na[0][1]
    #print i
    for ii in range(int(nq[0][1])):
        qpoint[ii][3]=qpoint[ii][3].replace(",","")
        qpoint[ii][4]=qpoint[ii][4].replace(",","")
        qpoint[ii][5]=qpoint[ii][5].replace(",","")
    #    print qpoint[ii][3],qpoint[ii][4],qpoint[ii][5]
    #    print distance[ii][1]
    #for iii in range(int(na[0][1])*3*int(nq[0][1])):
    #    print frequency[iii][1]
    for aa in range(len(eigen)):
        for aaa in range(int(na[0][1])):
            eigenvx.append(fff[eigen[aa]+aaa*4+2][2])
            eigenvy.append(fff[eigen[aa]+aaa*4+3][2])
            eigenvz.append(fff[eigen[aa]+aaa*4+4][2])
            ieigenvx.append(fff[eigen[aa]+aaa*4+2][3])
            ieigenvy.append(fff[eigen[aa]+aaa*4+3][3])
            ieigenvz.append(fff[eigen[aa]+aaa*4+4][3])
    for bb in range(len(eigenvx)):
        eigenvx[bb]=eigenvx[bb].replace(",","")
        eigenvy[bb]=eigenvy[bb].replace(",","")
        eigenvz[bb]=eigenvz[bb].replace(",","")
        ieigenvx[bb]=ieigenvx[bb].replace(",","")
        ieigenvy[bb]=ieigenvy[bb].replace(",","")
        ieigenvz[bb]=ieigenvz[bb].replace(",","")
    #    print eigenvx[bb],eigenvy[bb],eigenvz[bb]
    #    print ieigenvx[bb],ieigenvy[bb],ieigenvz[bb]
    for iii in range(int(na[0][1])):
        filename = 'fat-band%s' %iii
        f=open(filename,'w')
        for i in range(int(na[0][1])*3):
            for ii in range(int(nq[0][1])):
                j=ii*(int(na[0][1])*3*int(na[0][1]))+i*int(na[0][1])+iii
    #            print iii,i,ii,j
                ampx=(float(eigenvx[j])**2+float(ieigenvx[j])**2)**0.5
                ampy=(float(eigenvy[j])**2+float(ieigenvy[j])**2)**0.5
                ampz=(float(eigenvz[j])**2+float(ieigenvz[j])**2)**0.5
    #            print eigenvx[j],ieigenvx[j],eigenvy[j],ieigenvy[j],eigenvz[j],ieigenvz[j]
                ampl=(ampx**2+ampy**2+ampz**2)**0.5
                print >> f, distance[ii][1],frequency[ii*int(na[0][1])*3+i][1],ampl
            print >> f, ''
        f.close()

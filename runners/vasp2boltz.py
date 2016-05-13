from numpy  import *
from numpy.linalg import *
from ase.calculators.vasp import *
from ase import units
import os
import sys
try:
 from pyspglib import spglib
 has_spglib=True
except:
 has_spglib=False

######################################################################          
def get_vasp_bandstructure(pathname='./', filename_eigenval='EIGENVAL', pathname_sc=None, filename_outcar_sc='OUTCAR'):
    # return a dictionary with band structure information, keywords are:
    # 'unit_E': energy unit ('eV', 'Hartree', 'Rydberg')
    # 'E_Fermi': Fermi energy in above unit
    # 'type': 'lines' or 'BZ'
    # 'nkp': number of k-points
    # 'nband': number of bands
    # 'nspin': number of spins
    # 'e_n_k': list with actual band structure;
    # e_n_k[ikp] is a dictionary with 'kpoint' (k-vector) and 'energies' (e_n(k)) for k-point #ikp

    bandstructure={}
    parentdir=os.getcwd()
    # get the Fermi energy from a selfconsistent OUTCAR file (not from bs calculation!)
    if pathname_sc is None:
        pathname_sc=pathname
    if filename_outcar_sc is not None:
        outcarfile=os.path.join(pathname_sc,filename_outcar_sc)
        if os.path.isfile(outcarfile):
            os.chdir(pathname_sc)
            calc=Vasp()
            ef=calc.read_fermi()
            nelect=calc.read_number_of_electrons()
            os.chdir(parentdir)
        else:
            print 'get_vasp_bandstructure(): could not find OUTCAR file <%s>' % outcarfile
            return None
    else:
        ef=None
    bandstructure['unit_E']='eV'
    bandstructure['E_Fermi']=ef    
    bandstructure['n_electrons']=nelect    
    # evaluate eigenval file
    fin=None
    if (pathname is not None) and (filename_eigenval is not None):
        eigenvalfile=os.path.join(pathname,filename_eigenval)
        if os.path.isfile(eigenvalfile):
            fin=open(eigenvalfile,"r")
    else:
        eigenvalfile=None
    if fin is None:
        print 'get_vasp_bandstructure(): could not find EIGENVAL file <%s>' % eigenvalfile
        return None
    #read header
    line=fin.readline() #??/??/??/nspin
    sline=line.split()
    nspin=int(sline[3])
    bandstructure['nspin']=nspin
    line=fin.readline() #unknown data
    line=fin.readline() #unknown data
    line=fin.readline() #unknown data
    line=fin.readline() #unknown data
    line=fin.readline() #??/nkp/nband
    sline=line.split()
    nband=int(sline[2])
    nkp=int(sline[1])
    bandstructure['nkp']=nkp
    bandstructure['nband']=nband
    bandstructure['e_n_k']=[]
    for ikp in range(nkp):
        ek={}
        line=fin.readline() #empty line
        line=fin.readline() #kx ky kz wk
        sline=line.split()
        ek['kpoint']=[float(sline[0]),float(sline[1]),float(sline[2])]
        energies_up=[]
        energies_dn=[]
        for iband in range(nband):
            line=fin.readline() #no. energy
            sline=line.split()
            energies_up.append(float(sline[1]))
            if nspin==2:
                energies_dn.append(float(sline[2]))
        ek['energies_up']=energies_up
        if nspin==2:
            ek['energies_dn']=energies_dn
        bandstructure['e_n_k'].append(ek)
    fin.close()
    return bandstructure



def write_bandstructure_boltztrap(bandstructure, unit='Ry', yscale=None, E_Fermi_zero=True, pathname='./', filename="energies.boltztrap", ktrafo=None, runBoltz=False):
    # output of bandstructure in a format which can be used by boltztrap code
    # input: bandstrcuture as dictionary like in get_vasp_bandstructure
    # default is unit conversion to Ry; if yscale is given this is taken as conversion factor
    if bandstructure['nspin']!=1:
        print 'write_bandstructure_boltztrap: No idea what to do for nspin=%d' % bandstructure['nspin']
        return False
    if  (E_Fermi_zero==True) and (bandstructure['E_Fermi'] is None):
        print 'write_bandstructure_boltztrap: Fermi energy not known, nothing done'
        return False
    #    
    eb=bandstructure['e_n_k']
    if yscale is None:
        if (unit is not None) and (bandstructure['unit_E'] is not None):
            yscale=get_conversion(bandstructure['unit_E'],unit)
        if yscale is None:
            print 'write_bandstructure_boltztrap: unit conversion failed, nothing done!'
            return False
    outfile=os.path.join(pathname,filename)
    fout=open(outfile,"w")
    fout.write('HTE output'+'\n') # title
    fout.write(str(len(eb))+'\n') # no. of k-points
    print '***************'
    print ktrafo
    for ikp in range(len(eb)):
        enk=eb[ikp]
        #print ikp,enk['kpoint']
        if ktrafo is None:
            kp=enk['kpoint']
        else:
            kp=[0.0,0.0,0.0]
            for i in range(3):
                for j in range(3):
                    kp[i]=kp[i]+ktrafo[j][i]*enk['kpoint'][j]
                    #kp[i]=kp[i]+ktrafo[i][j]*enk['kpoint'][j]
        #print ikp,enk['kpoint'],' -> ',kp
        #line=''
        #for i in range(3):
        #    line=line+' '+str(kp[i])
        #line=line+' '+str(len(enk['energies_up']))+'\n'
        #fout.write(line)
        fout.write("%12.8f %12.8f %12.8f %d\n" %(kp[0],kp[1],kp[2],len(enk['energies_up'])))
        line=''
        for i in range(len(enk['energies_up'])):
            e_up=enk['energies_up'][i]
            if (E_Fermi_zero==True):
                e_up=e_up-bandstructure['E_Fermi']
            e_up=e_up*yscale
            fout.write("%18.8f\n" %e_up)
    fout.close()
    if runBoltz==True:
        parentdir=os.getcwd()
        os.chdir(pathname)
        exitcode, out=commands.getstatusoutput('/home/users/opahlivs/bin/BoltzTraP BoltzTraP.def')
        os.chdir(parentdir)
        print exitcode, out
    return 0
    
def write_structure_boltztrap(ao, pathname='./', filename="hte.struct"):
    if not os.path.isdir(pathname):
        os.path.mkdir(pathname)
    fname=os.path.join(pathname,filename)
    fout=open(fname,"w")
    fout.write('HTE output'+'\n') # title
    latt=ao.get_cell()/0.5291772083
    for i in range(3):
        line=''
        for j in range(3):
            line=line+"%12.5f"%latt[i][j]
        fout.write(line+'\n')
    krot=get_kspace_operations(ao)
    fout.write(str(len(krot))+'\n')
    for iop in range(len(krot)):
        for i in range(3):
            for j in range(3):
                fout.write(str(krot[iop][i][j])+' ')
            fout.write('\n')
    fout.close()
    return 0




def get_transport_boltztrap(pathname='./', filename='*.trace', Tmin=299.0, Tmax=301.0):
    tracefile=os.path.join(pathname,filename)
    if os.path.isfile(tracefile):
        fin=open(tracefile,"r")
        #read header line
        line=fin.readline()
        line=fin.readline()
        transport={}
        transport['E_Fermi']=[]
        transport['T']=[]
        transport['N']=[]
        transport['DOS']=[]
        transport['Seebeck']=[]
        transport['sigotau']=[]
        transport['specific_heat']=[]
        transport['R_H']=[]
        while line:
            property={}
            sline=line.split()
            Ef=float(sline[0]) #*27.2113834/2.0 #use eV
            T=float(sline[1])
            N=float(sline[2])
            DOS=float(sline[3])
            S=float(sline[4])
            sigotau=float(sline[5]) 
            R_H=float(sline[6])
            kap_0=float(sline[7])
            c=float(sline[8])
            chi=float(sline[9])
            if (T>Tmin) and (T<Tmax):
                transport['E_Fermi'].append(Ef)
                transport['T'].append(T)
                transport['N'].append(N)
                transport['DOS'].append(DOS)
                transport['Seebeck'].append(S)
                transport['sigotau'].append(sigotau)
                transport['specific_heat'].append(c)
                transport['R_H'].append(R_H)
            line=fin.readline()
        fin.close()
        return transport
    else:
        print 'Could not open ',tracefile
    return None




def read_genstruct_boltztrap(filename):
    latt=None
    kops=None
    if os.path.isfile(filename):
        fin=open(filename,"r")
        line=fin.readline()
        #lattice vectors
        latt=[]
        for i in range(3):
            sline=fin.readline().split()
            vec=[]
            for j in range(3):
                vec.append(float(sline[j]))
            latt.append(vec)
        print latt
        #kops
        nkops=int(fin.readline())
        kops=[]
        for iops in range(nkops):
            mat=[]
            for i in range(3):
                sline=fin.readline().split()
                vec=[]
                for j in range(3):
                    vec.append(float(sline[j]))
                mat.append(vec)
            kops.append(mat)
        print kops
    return latt,kops



def write_intrans_boltztrap(boltzdir='./',filename='hte.intrans', boltz_generic=True, E_Fermi=0.0, n_electrons=1.0):
    parentdir=os.getcwd()
    if not os.path.isdir(boltzdir):
        os.mkdir(boltzdir)
    os.chdir(boltzdir)
    fout=open(filename,'w')
    if boltz_generic:
        fout.write("GENE          # use generic interface\n")
    else:
        fout.write("WIEN          # use wien interface\n")
    fout.write("0 0 0 0.0         # iskip (not presently used) idebug setgap shiftgap \n")
    fout.write("%7.5f 0.0005 0.4 %6.1f     # Fermilevel (Ry), energygrid, energy span around Fermilevel, number of electrons\n"%(E_Fermi,n_electrons))
    fout.write("CALC                    # CALC (calculate expansion coeff), NOCALC read from file\n")
    fout.write("5                         # lpfac, number of latt-points per k-point\n")
    fout.write("BOLTZ                     # run mode (only BOLTZ is supported)\n")
    fout.write(".15                       # (efcut) energy range of chemical potential\n")
    fout.write("300. 10.                  # Tmax, temperature grid\n")
    fout.write("-1.                       # energyrange of bands given individual DOS output sig_xxx and dos_xxx (xxx is band number)\n")
    fout.write("HISTO\n")
    fout.close()
    os.chdir(parentdir)
    return True


def get_kspace_operations(ao, methods=['atoms_info','spglib'], symprec_spglib=1e-5):
    # returns k-space operations for the atoms object ao
    kops=None
    for method in methods:
        if method=='atoms_info':
            # get operations from ao.info
            if ('spacegroup' in ao.info) and (ao.info['spacegroup'] is not None):
                if ('unit_cell' in ao.info):
                    if (ao.info['unit_cell']=='conventional'):
                        primitive_cell=False
                    else:
                        primitive_cell=True
                else:
                    primitive_cell=True
                    print 'get_kspace_operations(): Warning, assuming primitive cell'
                rot=ao.info['spacegroup'].rotations
                p2c_dir=ao.info['spacegroup'].scaled_primitive_cell
                p2c_rec=ao.info['spacegroup'].reciprocal_cell 
                kops=[]
                for iop in range(len(rot)):
                    if primitive_cell:
                        mat=dot(p2c_rec.T, dot(rot[iop],p2c_dir))
                    else:
                        mat=rot[iop]
                    newop=True
                    for i in range(len(kops)):
                        if cmp_mat(kops[i],mat):
                            newop=False
                            break
                    if newop:
                        kops.append(mat)
                return kops
            else:
                print 'get_kspace_operations(): atoms object has no space group information'
        if method=='spglib':
            if has_spglib:
                rot=spglib.get_symmetry_dataset(ao, symprec=symprec_spglib)['rotations']
                kops=[]
                for iop in range(len(rot)):
                    mat=rot[iop]
                    newop=True
                    for i in range(len(kops)):
                        if cmp_mat(kops[i],mat):
                            newop=False
                            break
                    if newop:
                        kops.append(mat)
                return kops
            else:
                print 'get_kspace_operations(): spglib not found'
        

def get_conversion(unitin,unitout):
    conv = {'Ry': units.Ry, 'eV': 1.0, 'J': units._e}
    try:
        return conv[unitin] / conv[unitout]
    except KeyError:
        print 'Error: Unknown unit.'
        sys.exit()

def cmp_mat(mat1, mat2):
    absdiff = abs(mat1 - mat2)
    value = sum(sum(absdiff))
    if value > 1.0e-10:
        return False
    else:
        return True

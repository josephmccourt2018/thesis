import numpy as np
import MDAnalysis as mda
import matplotlib.pyplot as plt
plt.style.use('presentation')
from tqdm import tqdm
import pandas as pd


def dot_prod(a,b):
    return a[0]*b[0]+a[1]*b[1]+a[2]*b[2]

def mag(a):
    return np.sqrt(a[0]**2+a[1]**2+a[2]**2)

#finding normal vector (n) and d for plane defined as ax + by + cz = d with
#crete_A matrix containing list of coordinates of three atoms defining plane
#see https://www.had2know.com/academics/equation-plane-through-3-points.html
def nd(A):
    t1=A[1]-A[0]
    t2=A[2]-A[0]
    n=np.cross(t1,t2)
    d=n[0]*A[0][0]+n[1]*A[0][1]+n[2]*A[0][2]
    return n,d
    
def dihedral(t1,t2,t3,t4):
    #input vectors defined by the 3 points for the plane
    n1=np.cross(t1,t2)
    n2=np.cross(t3,t4)
    #calculate angle between two normal vectors
    numer = dot_prod(n1,n2)
    denom1=mag(n1)
    denom2=mag(n2)
    theta = np.arccos(numer/(denom1*denom2))
    theta_deg=theta*180/np.pi
    return theta_deg

# create list of positions for 3 atoms defining a plane from AtomGroups
def create_A(atom1, atom2, atom3):
    return np.array([[atom1.position[0], atom1.position[1], atom1.position[2]], [atom2.position[0], atom2.position[1], atom2.position[2]], [atom3.position[0], atom3.position[1], atom3.position[2]]])
    

#first need to make trajectory corrected for periodic boundaries and molecules
#jumping across box
#
# gmx trjconv -f combined.xtc -s npt.tpr -o combined_whole.xtc -pbc whole
#

charge=[0,50,100]

for c in charge:
    u=mda.Universe(str(c)+'/npt.tpr',str(c)+'/combined_whole.xtc')

    #dihedral 1: CG CB CA N
    #dihedral 2: NT  C CA N
    #(long) dihedral 3: NZ CD CA N
    #dihedral 4: NT C CA CB
    #dihedral 5: NZ CE CD CG
    #dihedral 6: CE CD CG CB
    #dihedral 7: CD CG CB CA
    z_box = u.dimensions[2]

    CAs_t=u.select_atoms('resname LYS or resname LSN and name CA and prop z >='+str(z_box/2))
    CGs_t=u.select_atoms('resname LYS or resname LSN and name CG and prop z >='+str(z_box/2))
    CBs_t=u.select_atoms('resname LYS or resname LSN and name CB and prop z >='+str(z_box/2))
    CDs_t=u.select_atoms('resname LYS or resname LSN and name CD and prop z >='+str(z_box/2))
    CEs_t=u.select_atoms('resname LYS or resname LSN and name CE and prop z >='+str(z_box/2))
    Cs_t=u.select_atoms('resname LYS or resname LSN and name C and prop z >='+str(z_box/2))
    Ns_t=u.select_atoms('resname LYS or resname LSN and name N and prop z >='+str(z_box/2))
    NTs_t=u.select_atoms('resname LYS or resname LSN and name NT and prop z >='+str(z_box/2))
    NZs_t=u.select_atoms('resname LYS or resname LSN and name NZ and prop z >='+str(z_box/2))
    
    # dihedral calculations
    dih_1_all_frames=[]
    dih_2_all_frames=[]
    dih_3_all_frames=[]
    dih_4_all_frames=[]
    dih_5_all_frames=[]
    dih_6_all_frames=[]
    dih_7_all_frames=[]
    for ts in tqdm(u.trajectory):
        dih_1=[]
        dih_2=[]
        dih_3=[]
        dih_4=[]
        dih_5=[]
        dih_6=[]
        dih_7=[]
        for i in range(0,len(CAs_t)):
            #######################################################################
            # dihedral 1
            # solving for normal vector of plane defined by CG CB CA
            # ax+by+cz=d
            #A1 = np.array([[CGs_t[i].position[0], CGs_t[i].position[1], CGs_t[i].position[2]], [CBs_t[i].position[0], CBs_t[i].position[1], CBs_t[i].position[2]], [CAs_t[i].position[0], CAs_t[i].position[1], CAs_t[i].position[2]]])
            A1=create_A(CGs_t[i],CBs_t[i],CAs_t[i])
            #b = np.array([1, 1, 1])
            #n = np.linalg.solve(A, b)
            temp1=A1[1]-A1[0]
            temp2=A1[2]-A1[0]
            n1 = np.cross(temp1,temp2)
            d1 = n1[0]*A1[0][0]+n1[1]*A1[0][1]+n1[2]*A1[0][2]
            # solving for normal vector of plane defined by CB CA N
            # ax+by+cz=d
            A2=create_A(CBs_t[i],CAs_t[i],Ns_t[i])
            temp3=A2[1]-A2[0]
            temp4=A2[2]-A2[0]
            n2 = np.cross(temp3,temp4)
            d2 = n2[0]*A2[0][0]+n2[1]*A2[0][1]+n2[2]*A2[0][2]
            
            dih_1.append(dihedral(temp1,temp2,temp3,temp4))
            
            #######################################################################
            # dihedral 2
            # solving for normal vector of plane defined by NT C CA
            # ax+by+cz=d
            A1=create_A(NTs_t[i],Cs_t[i],CAs_t[i])
            temp1=A1[1]-A1[0]
            temp2=A1[2]-A1[0]
            
            # solving for normal vector of plane defined by C CA N
            # ax+by+cz=d
            A2 = create_A(Cs_t[i],CAs_t[i],Ns_t[i])
            temp3=A2[1]-A2[0]
            temp4=A2[2]-A2[0]
            
            dih_2.append(dihedral(temp1,temp2,temp3,temp4))
            
            #######################################################################
            # dihedral 3
            # solving for normal vector of plane defined by NZ CD CA
            # ax+by+cz=d
            A1 = create_A(NZs_t[i],CDs_t[i],CAs_t[i])
            temp1=A1[1]-A1[0]
            temp2=A1[2]-A1[0]
            
            # solving for normal vector of plane defined by CD CA N
            # ax+by+cz=d
            A2 = create_A(CDs_t[i],CAs_t[i],Ns_t[i])
            temp3=A2[1]-A2[0]
            temp4=A2[2]-A2[0]
            
            dih_3.append(dihedral(temp1,temp2,temp3,temp4))
            
            
            #######################################################################
            # dihedral 4
            # solving for normal vector of plane defined by NT C CA
            # ax+by+cz=d
            A1 = create_A(NTs_t[i],Cs_t[i],CAs_t[i])
            temp1=A1[1]-A1[0]
            temp2=A1[2]-A1[0]
            
            # solving for normal vector of plane defined by C CA CB
            # ax+by+cz=d
            A2 = create_A(Cs_t[i],CAs_t[i],CBs_t[i])
            temp3=A2[1]-A2[0]
            temp4=A2[2]-A2[0]
            
            dih_4.append(dihedral(temp1,temp2,temp3,temp4))
            
            
            #######################################################################
            # dihedral 5
            # solving for normal vector of plane defined by NZ CE CD
            # ax+by+cz=d
            A1 = create_A(NZs_t[i],CEs_t[i],CDs_t[i])
            temp1=A1[1]-A1[0]
            temp2=A1[2]-A1[0]
            
            # solving for normal vector of plane defined by CE CD CG
            # ax+by+cz=d
            A2 = create_A(CEs_t[i],CDs_t[i],CGs_t[i])
            temp3=A2[1]-A2[0]
            temp4=A2[2]-A2[0]
            
            dih_5.append(dihedral(temp1,temp2,temp3,temp4))

            #######################################################################
            # dihedral 6
            # solving for normal vector of plane defined by CE CD CG
            # ax+by+cz=d
            A1 = create_A(CEs_t[i],CDs_t[i],CGs_t[i])
            temp1=A1[1]-A1[0]
            temp2=A1[2]-A1[0]
            
            # solving for normal vector of plane defined by CD CG CB
            # ax+by+cz=d
            A2 = create_A(CDs_t[i],CGs_t[i],CBs_t[i])
            temp3=A2[1]-A2[0]
            temp4=A2[2]-A2[0]
            
            dih_6.append(dihedral(temp1,temp2,temp3,temp4))

            #######################################################################
            # dihedral 7
            # solving for normal vector of plane defined by CD CG CB
            # ax+by+cz=d
            A1 = create_A(CDs_t[i],CGs_t[i],CBs_t[i])
            temp1=A1[1]-A1[0]
            temp2=A1[2]-A1[0]
            
            # solving for normal vector of plane defined by CG CB CA
            # ax+by+cz=d
            A2 = create_A(CGs_t[i],CBs_t[i],CAs_t[i])
            temp3=A2[1]-A2[0]
            temp4=A2[2]-A2[0]
            
            dih_7.append(dihedral(temp1,temp2,temp3,temp4))
            
        dih_1_all_frames.append(dih_1)
        dih_2_all_frames.append(dih_2)
        dih_3_all_frames.append(dih_3)
        dih_4_all_frames.append(dih_4)
        dih_5_all_frames.append(dih_5)
        dih_6_all_frames.append(dih_6)
        dih_7_all_frames.append(dih_7)
        
    #dihedral 1: CG CB CA N
    #dihedral 2: NT  C CA N
    #(long) dihedral 3: NZ CD CA N
    #dihedral 4: NT C CA CB
    #dihedral 5: NZ CE CD CG
    #dihedral 6: CE CD CG CB
    #dihedral 7: CD CG CB CA
    flat_dih_1  = [val for sublist in dih_1_all_frames for val in sublist]
    plt.figure()
    plt.hist(flat_dih_1,bins=100)
    plt.title("CG CB CA N Dihedral 1, Charge "+str(c)+"%")
    plt.xlabel("Theta (deg)")
    plt.ylabel("Counts")
    plt.savefig(str(c)+"/dih1_c"+str(c)+".png")
    #plt.show()
    
    flat_dih_2  = [val for sublist in dih_2_all_frames for val in sublist]
    plt.figure()
    plt.hist(flat_dih_2,bins=100)
    plt.title("NT C CA N Dihedral 2, Charge "+str(c)+"%")
    plt.xlabel("Theta (deg)")
    plt.ylabel("Counts")
    plt.savefig(str(c)+"/dih2_c"+str(c)+".png")
    #plt.show()
    
    flat_dih_3  = [val for sublist in dih_3_all_frames for val in sublist]
    plt.figure()
    plt.hist(flat_dih_3,bins=100)
    plt.title("NZ CD CA N Dihedral 3, Charge "+str(c)+"%")
    plt.xlabel("Theta (deg)")
    plt.ylabel("Counts")
    plt.savefig(str(c)+"/dih3_c"+str(c)+".png")
    #plt.show()
    
    flat_dih_4  = [val for sublist in dih_4_all_frames for val in sublist]
    plt.figure()
    plt.hist(flat_dih_4,bins=100)
    plt.title("NT C CA CB Dihedral 4, Charge "+str(c)+"%")
    plt.xlabel("Theta (deg)")
    plt.ylabel("Counts")
    plt.savefig(str(c)+"/dih4_c"+str(c)+".png")
    #plt.show()
    
    flat_dih_5  = [val for sublist in dih_5_all_frames for val in sublist]
    plt.figure()
    plt.hist(flat_dih_5,bins=100)
    plt.title("NZ CE CD CG Dihedral 5, Charge "+str(c)+"%")
    plt.xlabel("Theta (deg)")
    plt.ylabel("Counts")
    plt.savefig(str(c)+"/dih5_c"+str(c)+".png")
    #plt.show()

    flat_dih_6  = [val for sublist in dih_6_all_frames for val in sublist]
    plt.figure()
    plt.hist(flat_dih_6,bins=100)
    plt.title("CE CD CG CB Dihedral 6, Charge "+str(c)+"%")
    plt.xlabel("Theta (deg)")
    plt.ylabel("Counts")
    plt.savefig(str(c)+"/dih6_c"+str(c)+".png")
    #plt.show()

    flat_dih_7  = [val for sublist in dih_7_all_frames for val in sublist]
    plt.figure()
    plt.hist(flat_dih_7,bins=100)
    plt.title("CD CG CB CA Dihedral 7, Charge "+str(c)+"%")
    plt.xlabel("Theta (deg)")
    plt.ylabel("Counts")
    plt.savefig(str(c)+"/dih7_c"+str(c)+".png")
    #plt.show()
    
    dictionary={'dih1':flat_dih_1,'dih2':flat_dih_2,'dih3':flat_dih_3,'dih4':flat_dih_4,'dih5':flat_dih_5, 'dih6':flat_dih_6, 'dih7':flat_dih_7}
    df=pd.DataFrame(dictionary)
    df.to_csv(str(c)+"/dihedrals_c"+str(c)+".csv")

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
plt.style.use('presentation')
from tqdm import tqdm
import MDAnalysis as mda
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from scipy.fftpack import fft, fft2, ifft2
from scipy.interpolate import UnivariateSpline as us
from scipy.optimize import curve_fit
from PIL import Image
import re
import os

#origin o, voxel size
def voxel(o,size=(1,1,1)):
    #defining a voxel by size planes of a cube (4 points each)
    X = [[[0, 1, 0], [0, 0, 0], [1, 0, 0], [1, 1, 0]],
         [[0, 0, 0], [0, 0, 1], [1, 0, 1], [1, 0, 0]],
         [[1, 0, 1], [1, 0, 0], [1, 1, 0], [1, 1, 1]],
         [[0, 0, 1], [0, 0, 0], [0, 1, 0], [0, 1, 1]],
         [[0, 1, 0], [0, 1, 1], [1, 1, 1], [1, 1, 0]],
         [[0, 1, 1], [0, 0, 1], [1, 0, 1], [1, 1, 1]]]
    X = np.array(X).astype(float)
    for i in range(3):
        X[:,:,i] *= size[i]
    X += np.array(o)
    return X

def in_voxel(v,atom):
    #max and min boundaries are set by voxel planes
    bmax=np.max(v,0)[0]
    bmin=np.min(v,0)[0]
    count=0
    #loop over all coordinates for each atom position (x,y,z)
    for i in range(0,3):
        if bmin[i]<atom.position[i]<bmax[i]:
            count+=1
        else:
            continue;
    if count==3:
        return True
    else:
        return False
    
#voxelize a box of size Lx,Ly,Lz into n voxels
def voxelize(Lx,Ly,Lz,n):
    v_list=[voxel((i*Lx/n,j*Ly/n,k*Lz/n),size=(Lx/n,Ly/n,Lz/n)) for k in range(0,n) for j in range(0,n)for i in range(0,n)]
    return v_list

def voxelize_prism(Lx,Ly,Lz,n):
    v_list=[voxel((i*Lx/n,j*Ly/n,0),size=(Lx/n,Ly/n,Lz)) for j in range(0,n) for i in range(0,n)]
    return v_list

#https://stackoverflow.com/questions/42611342/representing-voxels-with-matplotlib

# plotting voxels

def plotCubeAt(positions,sizes=(1,1,1),colors=None, **kwargs):
    if not isinstance(colors,(list,np.ndarray)): colors=["C0"]*len(positions)
    if not isinstance(sizes,(list,np.ndarray)): sizes=[sizes]*len(positions)
    g = []
    for p,s,c in zip(positions,sizes,colors):
        g.append( voxel(p, size=s) )
    return Poly3DCollection(np.concatenate(g),  
                            facecolors=np.repeat(colors,6, axis=0),alpha=0.2, **kwargs)
def plotPrisms(Lx,Ly,Lz,n,v_list):
    #Lx=69
    #Ly=60
    #Lz=60
    #n=8
    N1 = int(Lx/n)
    N2 = int(Ly/n)
    N3 = int(Lz)
    #v_list=voxelize(Lx,Ly,Lz,n)
    #v_list=voxelize_prism(Lx,Ly,Lz,n)
    pos=[]
    counter=0
    for v in v_list:
        pos.append(np.min(v,0)[0])
        counter+=1
    colors= np.random.rand(len(pos),3)
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    pc = plotCubeAt(pos,sizes=(N1,N2,N3),color='r',edgecolor="k")#,colors=colors,edgecolor="k")
    ax.add_collection3d(pc)
    ax.set_xlim([0,Lx])
    ax.set_ylim([0,Ly])
    ax.set_zlim([0,Lz])
    plt.show()


# Main script
    
#charge=["0","10","20","30","40","50","60","70","80","90","100"]
#basedir_list=["/home/joey/Documents/quest_b1021/CK1/C12/charge_variation/10nm_box/","/home/joey/Documents/quest_b1021/CK1/C14/charge_variation/10nm_box/","/home/joey/Documents/quest_b1021/CK1/L_C16K_charge_variation/10nm_box/"]
#basedir_list=["/home/joey/Documents/quest_b1021/CK1/L_C16K_charge_variation/10nm_box/"]
charge=["40","60","70"]#,"100"]#"10","20","30","40","50","60","70","80","90"]
basedir_list=["/home/joey/Documents/quest_b1021/CK1/MARTINI_2/bilayers/charge_variation/"]
#tail_length=["3"]




for basedir in basedir_list:
    q_avg_list=[]
    hq2_avg_list=[]
    print(basedir)
    for c in charge:
        print(c)
        d_list=[]
        file1=basedir+c+"/5/bilayer/4x_bilayer/md_2.tpr"
        file2=basedir+c+"/5/bilayer/4x_bilayer/md_2_skip10.xtc"

        u = mda.Universe(file1,file2)
        #bilayer = u.select_atoms("resname CK8 or resname CK8+")#, 2
        #bilayer = u.select_atoms("resname CK12 or resname C12+")#, 3
        #bilayer = u.select_atoms("resname CK1 or resname CK1+")#, 4
        bilayer = u.select_atoms("resname C20K or resname C20+")#, 5
        
        # voxelizing the simulation box with rectangular prisms
        count=0
        heights=[]
        frames=[]
        n=16
        Lx=u.dimensions[0]
        Ly=u.dimensions[1]
        Lz=u.dimensions[2]
        v_list=voxelize_prism(Lx,Ly,Lz,n)
        os.mkdir(basedir +c+ "/5/bilayer/4x_bilayer/images"+str(n)+"/")

        #frames
        f1=1
        f2=500



        
        for ts in tqdm(u.trajectory):
            #if count%8==0:
            if f1<=count<f2:
            #if count==0:
                #print(count//5)
                #Lx=u.dimensions[0]
                #Ly=u.dimensions[1]
                #Lz=u.dimensions[2]
                # create voxels of box for every frame in trajectory (since box size changes a little bit)
                #v_list=voxelize_prism(Lx,Ly,Lz,n)
                j=0
                voxel_heights=[]
                #for v in tqdm(v_list):
                for v in v_list:
                    atom_list=[]
                    for atom in bilayer.atoms:
                        if in_voxel(v,atom):
                            atom_list.append(atom)
                        else:
                            continue

                    v_xyz=((np.min(v,0)[0]+np.max(v,0)[0])/2)

                    total=0
                    for atom in atom_list:
                        total+=atom.position[2] #z coordinate sum for all atoms

                    if len(atom_list)!=0:
                        voxel_heights.append({"Voxel {} Location,Height".format(j):[v_xyz,total/len(atom_list)]})
                    else:
                        voxel_heights.append({"Voxel {} Location,Height".format(j):[v_xyz,None]})

                    j+=1
                heights.append({"Frame {}".format(count):voxel_heights})
                frames.append(count)
                count+=1
            else:
                count+=1

        #for prisms
        #fig = plt.figure()
        #ax = fig.add_subplot(111, projection='3d')
        h0=[]
        #colors= np.random.rand(n**3,3) #random color for every voxel (don't really need this)
        v_layer_x=[]
        v_layer_y=[]
        v_layer_h=[]
        for i in range(0,len(heights)):
            #color=colors[i]
            for j in range(0,n**2): #n**2 for prisms (plane defines space)
                #heights[i]['Frame {}'.format(frame)][voxel_num]['Voxel {} Location,Height'.format(voxel_num)]
                if heights[i]['Frame {}'.format(frames[i])][j]['Voxel {} Location,Height'.format(j)][1]!=None:
                    #ax.scatter(heights[i]['Frame {}'.format(frames[i])][j]['Voxel {} Location,Height'.format(j)][0][0],heights[i]['Frame {}'.format(frames[i])][j]['Voxel {} Location,Height'.format(j)][0][1],heights[i]['Frame {}'.format(frames[i])][j]['Voxel {} Location,Height'.format(j)][1],color=color) 
                    v_layer_h.append(heights[i]['Frame {}'.format(frames[i])][j]['Voxel {} Location,Height'.format(j)][1])
                    v_layer_x.append(heights[i]['Frame {}'.format(frames[i])][j]['Voxel {} Location,Height'.format(j)][0][0])
                    v_layer_y.append(heights[i]['Frame {}'.format(frames[i])][j]['Voxel {} Location,Height'.format(j)][0][1])
        #ax.set_xlim(0,60)
        #ax.set_ylim(0,60)
        #ax.set_zlim(30,90)
        #plt.show()
        d_list.append({c:[v_layer_h,v_layer_x,v_layer_y]})

        df=pd.DataFrame({"x":d_list[0][c][1],"y":d_list[0][c][2],"h(x,y)":d_list[0][c][0]})

        h_xy_snapshots=[]
        for k in range(0,len(df)//n**2):
            h_xy_grid_2 = [[0 for i in range(n)] for j in range(n)]
            m=0
            for i in range(0,n):
                for j in range(m*n,(m+1)*n):
                    h_xy_grid_2[j-m*n][i]=df[df.columns[2]][j+k*n**2]
                m+=1
            #h_xy_grid_2=h_xy_grid_2-np.mean(h_xy_grid_2)
            h_xy_snapshots.append(h_xy_grid_2)

        hq2_snapshots=[]
        q_snapshots=[]
        for i in range(0,len(h_xy_snapshots)):
            """CREATING REAL AND MOMENTUM SPACES GRIDS"""
            N_x, N_y = n,n
            range_x, range_y = np.arange(N_x), np.arange(N_y)
            #dx, dy = 0.005,0.005
            #x0=v_layer_x[1]-v_layer_x[0]
            #y0=v_layer_y[n]-v_layer_y[0]
            #dx, dy = x0,y0
            # real space grid vectors
            #xv, yv = dx * (range_x - 0.5 * N_x), dy * (range_y - 0.5 * N_y)
            xv = df["x"][i*n:n*(i+1)].values
            yv=list(set(df["y"][i*n**2:n**2*(i+1)])) #remove duplicates from v_layer_y
            yv.sort() #hopefully this is ok

            #xv,yv=xv-0.5 * N_x*dx,yv-0.5 * N_y*dy

            dk_x, dk_y = np.pi / np.max(xv), np.pi / np.max(yv)
            # momentum space grid vectors, shifted to center for zero frequency
            k_xv, k_yv = dk_x * np.append(range_x[:N_x//2], -range_x[N_x//2:0:-1]), \
                 dk_y * np.append(range_y[:N_y//2], -range_y[N_y//2:0:-1])

            # create real and momentum spaces grids
            x, y = np.meshgrid(xv, yv, sparse=False, indexing='ij')
            kx, ky = np.meshgrid(k_xv, k_yv, sparse=False, indexing='ij')

            f = h_xy_snapshots[i]
            F = fft2(f)
            
            # for T=300
            kT=4.1*10**(-12)*10**(-9) # 4.1 pN * nM
            hq = []
            q=[]
            m=0
            for i in range(0,n):
                for j in range(0,n):
                    hq.append(np.abs(F[i][j])**2)
                    q.append(np.sqrt(kx[i][j]**2+ky[i][j]**2))
                m+=1
            hq2_snapshots.append(hq)
            q_snapshots.append(q)
        hq2_avg=np.mean(hq2_snapshots,axis=0)
        q_avg=np.mean(q_snapshots,axis=0)
        hq2_avg_list.append(hq2_avg)
        q_avg_list.append(q_avg)
        print(c)

        ###########################################################################################################
        # Visualization
        ###########################################################################################################

        #
#         """CREATING REAL AND MOMENTUM SPACES GRIDS"""
        for i in range(0,len(h_xy_snapshots)):
             N_x, N_y = n,n
             range_x, range_y = np.arange(N_x), np.arange(N_y)
             #dx, dy = 0.005,0.005
             #dx, dy = x0,y0
             # real space grid vectors
             #xv, yv = dx * (range_x - 0.5 * N_x), dy * (range_y - 0.5 * N_y)
             xv = v_layer_x[:n]
             yv=list(set(v_layer_y[:n**2])) #remove duplicates from v_layer_y
             yv.sort() #hopefully this is ok

             #xv,yv=xv-0.5 * N_x*dx,yv-0.5 * N_y*dy

             dk_x, dk_y = np.pi / np.max(xv), np.pi / np.max(yv)
             # momentum space grid vectors, shifted to center for zero frequency
             k_xv, k_yv = dk_x * np.append(range_x[:N_x//2], -range_x[N_x//2:0:-1]), \
                      dk_y * np.append(range_y[:N_y//2], -range_y[N_y//2:0:-1])

             # create real and momentum spaces grids
             x, y = np.meshgrid(xv, yv, sparse=False, indexing='ij')
             kx, ky = np.meshgrid(k_xv, k_yv, sparse=False, indexing='ij')

#             """HEIGHT MAP"""
             #f = h_xy_grid
             #f = h_xy_grid_2
             f = h_xy_snapshots[i]
             F = fft2(f)
#             """PLOTTING"""
             fig = plt.figure()
             ax = Axes3D(fig)
             surf = ax.plot_surface(x, y, np.abs(f),cmap='BuGn')

             # for other plots I changed to
             #fig2 = plt.figure()
             #ax2 =Axes3D(fig2)
             #surf = ax2.plot_surface(kx, ky, np.abs(F))#*dx*dy,cmap='viridis')
             #ax.set_xlim(0,20)
             #ax.set_ylim(0,20)
             ax.set_zlim(20,160)
             #plt.show()
             filename=basedir +c+ "/5/bilayer/4x_bilayer/images"+str(n)+"/100_height_map{}.png".format(i*n)
             plt.savefig(filename,dpi=75)
             plt.close(fig)

        png_count = len(h_xy_snapshots)
        files = []
        for i in range(png_count):
             seq = str(i*n)
             file_names = basedir + c+"/5/bilayer/4x_bilayer/images"+str(n)+"/100_height_map"+ seq +".png"
             files.append(file_names)

        # Create the frames
        frames = []
        for i in files:
             new_frame = Image.open(i)
             frames.append(new_frame)

        # Save into a GIF file that loops forever   
        frames[0].save(basedir + c+"/5/bilayer/4x_bilayer/images"+str(n)+"/3d_vis.gif", format='GIF',
                        append_images=frames[1:],
                        save_all=True,
                        duration=40, loop=0)

        #NEED COMMAND THAT WORKS TO DELETE FILES
        os.chdir(basedir+c+"/5/bilayer/4x_bilayer/images"+str(n)+"/")
        for file in os.listdir("."):
            if os.path.isfile(file) and file.endswith("png"):
                try:
                    os.remove(file)
                except Exception:
                    print(Exception)

#         ###########################################################################################################
#         # End Visualization
#         ###########################################################################################################

        for i in range(0,len(q_avg_list)):
            plt.scatter(q_avg_list[i][1:],hq2_avg_list[i][1:])

        for j in range(0,len(q_avg_list)):
            x=[]
            for i in range(0,len(q_avg_list[j])):
                x.append([q_avg_list[j][i],hq2_avg_list[j][i]])
            x.sort()
            temp1=[]
            temp2=[]
            for i in range(0,len(x)):
                temp1.append(x[i][0])
                temp2.append(x[i][1])

            #guess=[0.02,10]
            #popt,pcov = curve_fit(fit,temp1[1:],temp2[1:],bounds=[0.000001,20],p0=guess)#,sigma=temp2[1:])
            #plt.scatter(temp1[1:],temp2[1:])
            #plt.plot(temp1[1:],[fit(t,*popt) for t in temp1[1:]],label=charge[j])
            #plt.xlabel('q')
            #plt.ylabel('$\langle |h_{q}|^{2} \\rangle $')
            #plt.ylim(0,125)
            #plt.legend()
            df1=pd.DataFrame({"q_avg":temp1,"<hq2>":temp2})
            df1.to_csv(basedir+c+"/5/bilayer/4x_bilayer/bending_fluctuations_{}_n{}_f{}_to_f{}.csv".format(c,n,f1,f2))
            df2=pd.DataFrame({"h_xy":h_xy_snapshots})
            df2.to_csv(basedir+c+"/5/bilayer/4x_bilayer/height_map_{}_n{}_f{}_to_f{}.csv".format(c,n,f1,f2))
        #print("kappa{}={}".format(charge[j],popt[0]))
    #print("n[n*n grid] = {}".format(n))

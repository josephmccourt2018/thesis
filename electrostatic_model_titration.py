from scipy.optimize import fsolve
from scipy.special import kv,kn
from scipy import optimize
from scipy import integrate
from latex2sympy2 import latex2sympy, latex2latex
from titration_fit_funcs import *


####################
# take TeXForm in Mathematica of function and convert into python:
#
# latex2sympy(input) gives psiS_prime (then have to manually remove some pieces like "\text{}" or some variable names
#
#input=r"\frac{2 \left(1-e^{-0.3 a}\right) \left(1-\tanh \left(\frac{P}{4}\right)\right)\left(-\frac{\left(\tanh \left(\frac{P}{4}\right)+1\right) \tanh\left(\frac{P}{4}\right)}{\left(1-\tanh\left(\frac{P}{4}\right)\right)^2}-\frac{\tanh \left(\frac{P}{4}\right)}{1-\tanh\left(\frac{P}{4}\right)}\right)}{\tanh \left(\frac{P}{4}\right)+1}-\frac{2\left(1-e^{-0.3 a}\right) \log \left(\frac{\tanh\left(\frac{P}{4}\right)+1}{1-\tanh \left(\frac{P}{4}\right)}\right)}{a}+\frac{2e^{-0.3 a} \left(1-\tanh \left(\frac{P}{4}\right)\right) \left(\frac{-\frac{P\text{sech}^2\left(\frac{P}{4}\right)}{4 a}-\tanh \left(\frac{P}{4}\right)}{1-\tanh\left(\frac{P}{4}\right)}-\frac{\left(\tanh \left(\frac{P}{4}\right)+1\right)\left(\frac{P \text{sech}^2\left(\frac{P}{4}\right)}{4 a}+\tanh\left(\frac{P}{4}\right)\right)}{\left(1-\tanh\left(\frac{P}{4}\right)\right)^2}\right)}{\tanh \left(\frac{P}{4}\right)+1}"
####################

#BILAYER
thickness=3.45
height_bil=1.35
#CYLINDER
#5mM
radius_cyl=1.54
height_cyl=0.90
#250mM
radius_cyl=1.64
height_cyl=1.0
#MICELLE BEFORE TRANSITION
#pH 6.3 sample142
radius_mic=1.88
height_mic=0.78
#pH 8.0 sample146
radius_mic_2=2.19
height_mic_2=0.71


def main():
    counter_ions=0.008 #conc of PA * 2 [M]
    nacls=[0,0.005,0.05,0.25,0.5] #nacl concentrations [M]
    #nacls=np.linspace(0,0.025,25)
    #nacls=np.append(nacls,[0.05,0.25,0.5])
    nacls_exp=nacls #experimental nacls
    nacls=[nacl+counter_ions for nacl in nacls] # to account for counter ions in some way
    print(nacls)
    titration_files=["c16k2/Chengrui.dat","c16k2/L_C16KK_3.dat",\
                     "c16k2/L_C16KK_2.dat","c16k2/L_C16KK_4_edit.dat","c16k2/L_C16KK.dat"]#titration data corresponding to "nacls"
    #titration_files=["c16k2/L_C16KK_3.dat","c16k2/L_C16KK_3.dat",\
    #                 "c16k2/L_C16KK_2.dat","c16k2/L_C16KK_4.dat","c16k2/L_C16KK.dat"]#titration data corresponding to "nacls"
    vols=[6400,5000,5800,6556,3500] #volumes for titration
    #vols=[5000,5000,5800,6556,3500] #volumes for titration
    pKa_HH=[10.5 for i in range(0,len(nacls))]

    #with rho_tails=270, rho_h=500
    #Rs=[0.11*(c*1000)**(0.22)+1.85 for c in nacls]#radii of micelle [nm] for 0, 5, 50, 250 and 500 + counter_ions mM salt using power law
    #hs=[-0.45*(c*1000)**0.08+1.13 for c in nacls]

    #with rho_tails=205, rho_h=500
    #Rs=[0.15*(c*1000)**(0.21)+1.58 for c in nacls]
    #hs=[-0.12*(c*1000)**(0.22)+1.15 for c in nacls]

    #with rho_tails=220, rho_h=500
    #Rs=[(0.15*(c*1000)**(0.21)+1.62) for c in nacls]
    #hs=[-0.16*(c*1000)**(0.19)+1.13 for c in nacls]

    #with rho_tail=270, rho_h=450
    #c=[NACL]
    Rs=[0.14*(c*1000)**(0.21)+1.71 for c in nacls_exp]
    hs=[-0.24*(c*1000)**(0.14)+1.1 for c in nacls_exp]
    #c=IONIC STRENGTH ([NACL]+2*[PA])
    Rs=[0.52*(c*1000)**(0.1)+1.24 for c in nacls]
    hs=[-0.06*(c*1000)**(0.28)+0.86 for c in nacls]

    
    Rs_2=[radius_mic_2 for R in Rs]
    hs_2=[height_mic_2 for h in hs]
    

    plot_count=0
    #agg_nums=[20,22,25,28,30] # from MD

    
    rho_tails=270#electron density used in fits (around 80% rho_water = 334 e /nm^3)
    n_electrons=16*6+32*1+1*7+1*8 #number of electrons in tail: 16 carbons + 32 hydrogens + 1 nitrogen + 1 oxygen
    agg_nums=[rho_tails*(4/3)*np.pi*Rs[i]**3/n_electrons for i in range(0,len(Rs))]
    n_electrons_heads=12*6+27*1+4*7+2*8 #number of electrons in head groups: 12 carbons + 27 hydrogens + 4 nitrogen + 2 oxygen


    
    S0s=[4*np.pi*Rs[i]**2/(2*agg_nums[i]) for i in range(0,len(Rs))] #apl [nm^(2)](per 1 charge)
    
    print(agg_nums)


    #MICELLE
    S0s=[3*n_electrons/Rs[i]/rho_tails/2 for i in range(0,len(Rs))]
    print("area per charge (mic)={}".format(S0s))
    ##CYLINDRICAL MICELLE
    #S0s_cyl=[2*n_electrons/radius_cyl/rho_tails/2 for i in range(0,len(Rs))]
    #print("area per charge (cyl)={}".format(S0s_cyl))
    ##BILAYER
    #S0s_planar=[4*n_electrons/thickness/rho_tails/2 for i in range(0,len(Rs))]
    #print("area per charge (bil)={}".format(S0s_planar))

    S0s_2=[3*n_electrons/Rs_2[i]/rho_tails/2 for i in range(0,len(Rs))]


    for i in range(0,len(nacls)):

        v1=False
        v2=False
        v3=False
        v4=False
        pHrange=np.arange(5.,11.5,0.05)
       
        pKa_iso=10.5
        plt.figure()
        r1=run(pHrange,pKa_iso,nacls[i],Rs[i],S0s[i],hs[i],shape='sphere',verbose=v1)
        #r1_2=run(pHrange,pKa_iso,nacls[i],Rs_2[i],S0s_2[i],hs_2[i],shape='sphere',verbose=v1)
        r2=run(pHrange,pKa_iso,nacls[i],Rs[i],S0s[i],hs[i],shape='cyl',verbose=v2)
        r3=run(pHrange,pKa_iso,nacls[i],Rs[i],S0s[i],hs[i],shape='planar',verbose=v3)
        r4=run(pHrange,pKa_iso,nacls[i],Rs[i],S0s[i],hs[i],shape='isolated',verbose=v4)
        plt.figure()
        plt.scatter(-np.log10(r1[0]),r1[1:],label='sphere')
        #plt.scatter(-np.log10(r1_2[0]),r1_2[1:],label='sphere')#,c='#E41A1C')
        plt.scatter(-np.log10(r2[0]),r2[1:],label='cyl')
        plt.scatter(-np.log10(r3[0]),r3[1:],label='planar')
        plt.scatter(-np.log10(r4[0]),r4[1:],label='isolated')

        ################################################################
        #to plot shading for alpha plot
        #xs=-np.log10(r1[0])
        #ys=r4[1:]
        #ys2=r1[1:]
        #plt.fill_between(xs[1:110], ys[1:110], ys2[1:110], where=xs[1:110]<=8.1,  color='k',alpha=0.2)

        #xs=-np.log10(r2[0])
        #ys=r1[1:]
        #ys2=r2[1:]
        #plt.fill_between(xs[49:110], ys[49:110], ys2[49:110], where=xs[49:110]<=8.4,  color='k',alpha=0.2)

        #xs=-np.log10(r2[0])
        #ys=r2[1:]
        #ys2=r3[1:]
        #plt.fill_between(xs[55:116], ys[55:116], ys2[55:116], where=xs[55:116]<=12,  color='k',alpha=0.2)
        
        #plt.xlabel("pH",fontsize=28)
        #plt.ylabel("$\\alpha$, degree of ionization",fontsize=28)
        #plt.plot(pHrange,[hh(pH,pKa_HH[i]) for pH in pHrange],color='black',linestyle='dotted', label='HH, pKa={}'.format(pKa_iso))
        #plt.legend(fontsize=20)
        #plt.show()

        ################################################################

        
        ################################################################
        
        plt.figure()
        df_r1=titration2(r1[0],r1[1:], titration_files[i],vols[i],label="sphere",color='#E41A1c')
        #df_r1_2=titration2(r1_2[0],r1_2[1:], titration_files[i],vols[i],label="sphere",color='#E41A1C',alpha=0.5)
        df_r2=titration2(r2[0],r2[1:], titration_files[i],vols[i],label="cyl",color='#377EB8')
        df_r3=titration2(r3[0],r3[1:], titration_files[i],vols[i],label="planar",color='#4DAF4A')
        df_r4=titration2(r4[0],r4[1:], titration_files[i],vols[i],label="isolated",color='#F884BC')

        #df_r1.to_csv('sphere1.csv')
        #df_r1_2.to_csv('sphere2.csv')
        #df_r2.to_csv('cyl.csv')
        #df_r3.to_csv('bil.csv')
        #df_r4.to_csv('iso.csv')

        ###################################################

        #for i in range(0,len(df_r1[1])-85):
        #    def func(variables):
        #        (a,b)=variables
        #        Vb=df_r1[1]['Vb_red'][i]
        #        pH=df_r1[1]['pH'][i]
        #        eqn1=assembly_vol_frac_fit(Vb,a,b,df_r1,df_r2)-pH
        #        eqn2=a+b-1
        #        return [eqn1,eqn2]
        #    result = fsolve(func,(0.8,0.2))
        #    print(result)

            
            
        #plt.title("Exp={}mM NaCl, Model={}mM NaCl".format(nacls_exp[i]*1000,np.round(nacls[i]*1000,1)))
        plt.title("Exp={}mM NaCl".format(nacls_exp[i]*1000))
        plt.xlabel('V, 0.1 M NaOH / V, Equivalence Point',fontsize=20)
        plt.ylabel('pH',fontsize=20)
        plt.legend()
        plt.show()
        plot_count+=1

        def assembly_vol_frac_fit(Vb,a,b,dfa,dfb):
            return a*dfa[0](Vb)+b*dfb[0](Vb)

        #plt.figure()
        #for j in range(0,len(df_r1[1])-85):
        #    Vb=df_r1[1]['Vb_red'][j]
        #    pH=df_r1[1]['pH'][j]
        #    plt.scatter(Vb,pH,color='#E41A1C')
        #    Vb=df_r4[1]['Vb_red'][j]
        #    pH=df_r4[1]['pH'][j]
        #    plt.scatter(Vb,pH,color='#984EA3')
        #    model=assembly_vol_frac_fit(Vb,0.2,0.8,df_r1,df_r4)
        #    plt.scatter(Vb,model,color='#A65628')
        #plt.show()


#DLVO stability of colloidal systems, attractive and repulsive
#scaled by kT/e
#h is distance between particles
#R is radius
#a is kappa*R
#P is surface potential Psi0
##http://lib.tkk.fi/Diss/2007/isbn9789512288861/isbn9789512288861.pdf
#H is Haymaker constant
def V_att_sphere(h,R,H):
    return (-H/6)*((2*R**2/(h**2+4*R*h))+(2*R**2/(h+2*R)**2)+np.log((h**2+4*R*h)/(h+2*R)**2))
# t is thickness of charged layers
def V_att_planar(h,t,H):
    return (-H/(12*np.pi))*((1/h**2)-2/(h+t)**2+1/(h+2*t)**2)
#need cylinder, this just does evenly weighted sum of planar and sphere
def V_att_cyl(h,t,R,H):
    return (1/2)*V_att_planar(h,t,H)+(1/2)*V_att_sphere(h,R,H)
##http://lib.tkk.fi/Diss/2007/isbn9789512288861/isbn9789512288861.pdf
#attractive based on DLVO and potentials from model
def V_rep_sphere(h,R,P,kappa):
    return e**2



###############################################################
# behrens and grier j chem phys 2001
# sphere interaction energy of two charged spheres
#

def u(x,a,kappa,S,alpha):
    e=1.6
    lB=0.7
    sigma=alpha*e/S**2
    return (e/(lB*kappa))*(sigma**2*a**2/(1+a))*(np.exp(-x)/(2*a+x))#/(4.1*10**(-9)*10**(-12))



#Israelchivili ch. 13, ch 14
#p.254, 317
#V_rep (electrostatic) + V_att (vdw)
#Z coefficient
#def Z(lB,P):
#    return (16/lB)*(np.tanh(P/4))**2
#equal size sphere
#def E_sphere(h,a,R,H,P,lB):
#    return (R/2)*(Z(lB,P)*np.exp(-a*h/R)+(-H/(6*h)))
#parallel equal cylinders
#def E_cyl(h,kappa,a,R,H,P,lB):
#    return np.sqrt(R/2)*(np.sqrt(kappa/(2*np.pi))*Z(lB,P)*np.exp(-a*h/R)+(-H/(12*np.sqrt(2)*h**(3/2))))
#flat surfaces
#def E_planar(h,kappa,H,P,lB):
#    return (kappa/(2*np.pi))*Z(lB,P)*np.exp(-h)+(-H/(12*np.pi*h**2))


#########
#SPHERE
def psiS(x,a,P):
    return (2*a*(1 - np.exp(-3*a/10))*np.log((1 + np.exp(-x)*np.tanh(P/4))/(1 - np.exp(-x)*np.tanh(P/4))))/(a + x) + 2*np.exp(-3*a/10)*np.log((1 + np.exp(-x)*np.tanh((P*a)/((4*(a + x)))))/(1 - np.exp(-x)*np.tanh((P*a)/((4*(a + x))))))

def psiS_prime(x,a,P):
    return -2*a*(1 - np.exp(-3*a/10))*np.log((1 + np.exp(-x)*np.tanh(P/4))/(1 - np.exp(-x)*np.tanh(P/4)))/(a + x)**2 + (2*(1 - np.exp(-x)*np.tanh((P*a)/((4*(a + x)))))*(-(1 + np.exp(-x)*np.tanh((P*a)/((4*(a + x)))))*((a*P*(1/np.cosh((P*a)/(4*(a+x))))**2*np.exp(-x))/((4*(a + x)**2)) + np.exp(-x)*np.tanh((P*a)/((4*(a + x)))))/(1 - np.exp(-x)*np.tanh((P*a)/((4*(a + x)))))**2 + ((-a*P*(1/np.cosh((P*a)/(4*(a+x))))**2*np.exp(-x))/((4*(a + x)**2)) - np.exp(-x)*np.tanh((P*a)/((4*(a + x)))))/(1 - np.exp(-x)*np.tanh((P*a)/((4*(a + x))))))*np.exp(-3*a/10))/(1 + np.exp(-x)*np.tanh((P*a)/((4*(a + x))))) + (2*a*(1 - np.exp(-x)*np.tanh(P/4))*(1 - np.exp(-3*a/10))*(-np.exp(-x)*np.tanh(P/4)/(1 - np.exp(-x)*np.tanh(P/4)) - (1 + np.exp(-x)*np.tanh(P/4))*np.exp(-x)*np.tanh(P/4)/(1 - np.exp(-x)*np.tanh(P/4))**2))/(((1 + np.exp(-x)*np.tanh(P/4))*(a + x)))

def psiS_at0(a,P):
    #a=R/lambda_d
    #P=psi0
    return 2*(1 - np.exp(-3*a/10))*np.log((np.tanh(P/4) + 1)/(1 - np.tanh(P/4))) + 2*np.exp(-3*a/10)*np.log((np.tanh(P/4) + 1)/(1 - np.tanh(P/4)))

def psiS_prime_at0(a,P):
    return (2*(1 - np.exp(-3*a/10))*(1 - np.tanh(P/4))*(-(np.tanh(P/4) + 1)*np.tanh(P/4)/(1 - np.tanh(P/4))**2 - np.tanh(P/4)/(1 - np.tanh(P/4))))/(np.tanh(P/4) + 1) + (2*(1 - np.tanh(P/4))*(-(np.tanh(P/4) + 1)*(np.tanh(P/4) + (P*(1/(np.cosh(P/4))**2))/((4*a)))/(1 - np.tanh(P/4))**2 + (-np.tanh(P/4) - (P*((1/np.cosh(P/4))**2))/((4*a)))/(1 - np.tanh(P/4)))*np.exp(-3*a/10))/(np.tanh(P/4) + 1) - 2*(1 - np.exp(-3*a/10))*np.log((np.tanh(P/4) + 1)/(1 - np.tanh(P/4)))/a

##########
#CYLINDER
def psiC(x,a,P):
    return 2*np.exp(-3*a + x)*np.log((1 + np.exp(-x)*np.tanh((P*kn(0,(a + x)))/((4*kn(0,a))))/(1 - np.exp(-x)*np.tanh((P*kn(0,(a + x)))/((4*kn(0,a))))))) + (2*kn(0,a+x)*(1 - np.exp(-3*a))*np.exp(x)*np.log((1 + np.exp(-x)*np.tanh(P/4))/(1 - np.exp(-x)*np.tanh(P/4))))/((kn(0,a)))

def psiC_prime(x,a,P):
     return (2*np.log(-1 + (2*np.exp(x))/(np.exp(x) - np.tanh((P*kn(0,a+x))/((4*kn(0,a)))))) - 4*kn(0,a+x)*(np.exp(3*a) - 1)*np.exp(x)*np.tanh(P/4)/(kn(0,a)*((1/np.cosh(P/4))**2 + np.sinh(2*x) + np.cosh(2*x) - 1)) - (2*kn(0,a)*np.sinh((P*kn(0,a+x))/((2*kn(0,a)))) + P*kn(1,a+x))/(kn(0,a)*(np.sinh(x)*np.cosh((P*kn(0,a+x))/((2*kn(0,a)))) + np.cosh(x))) + (2*kn(0,a+x)*(np.exp(3*a) - 1)*np.log(1/(-1 + (2*np.exp(x))/(np.exp(x) + np.tanh(P/4)))))/((kn(0,a))) - 2*kn(1,a+x)*(np.exp(3*a) - 1)*np.log(1/(-1 + (2*np.exp(x))/(np.exp(x) + np.tanh(P/4))))/(kn(0,a)))*np.exp(-3*a + x)

def psiC_at0(a,P):
    #a=R/lambda_d
    #P=psi0
    return 2*(1 - np.exp(-3*a))*np.log((np.tanh(P/4) + 1)/(1 - np.tanh(P/4))) + 2*np.exp(-3*a)*np.log((np.tanh(P/4) + 1)/(1 - np.tanh(P/4)))

def psiC_prime_at0(a,P):
    return (2*(1 - np.exp(-3*a))*(1 - np.tanh(P/4))*(-(np.tanh(P/4) + 1)*np.tanh(P/4)/(1 - np.tanh(P/4))**2 - np.tanh(P/4)/(1 - np.tanh(P/4))))/(np.tanh(P/4) + 1) + (2*(1 - np.tanh(P/4))*(-(np.tanh(P/4) + 1)*(np.tanh(P/4) + (kn(1,a)*P*(1/np.cosh(P/4))**2)/((4*kn(0,a))))/(1 - np.tanh(P/4))**2 + (-np.tanh(P/4) - kn(1,a)*P*(1/np.cosh(P/4))**2/(4*kn(0,a)))/(1 - np.tanh(P/4)))*np.exp(-3*a))/(np.tanh(P/4) + 1) + 2*(1 - np.exp(-3*a))*np.log((np.tanh(P/4) + 1)/(1 - np.tanh(P/4))) + 2*np.exp(-3*a)*np.log((np.tanh(P/4) + 1)/(1 - np.tanh(P/4))) - 2*kn(1,a)*(1 - np.exp(-3*a))*np.log((np.tanh(P/4) + 1)/(1 - np.tanh(P/4)))/(kn(0,a))
##########
#PLANAR
def psiP(x,P):
    return 2*np.log((1+np.exp(-x)*np.tanh(P/4))/(1-np.exp(-x)*np.tanh(P/4)))

def psiP_prime(x,P):
    return -4*np.exp(x)*np.tanh(P/4)/(-1+np.exp(2*x)+(1/np.cosh(P/4))**2)

def psiP_at0(P):
    #P=psi0
    return 2*np.log((np.tanh(P/4) + 1)/(1 - np.tanh(P/4)))

def psiP_prime_at0(P):
    return (2*(1 - np.tanh(P/4))*(-(np.tanh(P/4) + 1)*np.tanh(P/4)/(1 - np.tanh(P/4))**2 - np.tanh(P/4)/(1 - np.tanh(P/4))))/(np.tanh(P/4) + 1)


#########
#NO ASSEMBLY ISOLATED MOLECULE
def psi_yukawa(x,P):
    return P*np.exp(-x)/x

def psi_yukawa_prime(x,P):
    return -np.exp(-x)*P*(1+x)/x**2

#Boundary condition
def BC_prefactor(lB,alpha,S):
    return 4*np.pi*lB*alpha/S

def f1(alpha):
    return alpha/(1-alpha)

def f2(pH,pKa):
    return 10**(pKa-pH)

#henderson hasselbach and hill for comparison
def hh(pH,pKa):
    return 1/(1+10**(pH-pKa))
    
def hill(pH,pKa,m):
    return 1/(1+10**(m*(pH-pKa)))



#electrostatic free energy of solvation
def mu_sol(kappa,lB):
    r_ion=0.2
    
    def f(k):
        return ((k-np.sqrt(kappa**2+k**2))/(k+np.sqrt(kappa**2+k**2)))*np.exp(-2*k*r_ion)

    return (lB/2)*integrate.quad(f,0,np.inf)[0]
    

#excess chemical potential of ions
#all functions scaled by beta *(i.e. mu_msa=beta*mu_msa)
#electrostatic contribution
#mean spherical approximation for chemical potential of ions
def mu_msa(kappa, conc_s):
    d=0.4 #ionic diameter
    conc_t=conc_s #molar concentrations of acid(base) and salt
    return (kappa*d*np.sqrt(1+2*kappa*d)-kappa*d-(kappa*d)**2)/(8*np.pi*conc_t*d**3)
#hard sphere
#carnahan-starling approximations
def mu_cs(conc_s):
    d=0.4 #ionic diameter
    conc_t=conc_s
    eta=np.pi*d**3*conc_t
    return (8*eta-9*eta**2+3*eta**3)/(1-eta)**3
def mu_ex(kappa,conc_s,A=1):
    return A*(mu_cs(conc_s)+mu_msa(kappa,conc_s))
#correction to proton activity at higher salt concentrations
def pH_activity(pH,kappa,conc_s):
     return np.log10(10**(pH)*np.exp(mu_ex(kappa,conc_s)))


#general
def calculate(pKa,pHrange,a,nacl,S0,h,shape='sphere',verbose=False):
    lB=0.7
    ld = 0.3/np.sqrt(nacl)
    kappa = (1/ld)
    R=a/kappa
    alpha0=0.99
    P0=18
    ent_scale=1 #scaling for entropy of ions term
    
    H=0.68 #hamaker (3/4)*((eps-epsw)/(eps+epsw))**2
    l0=0.18/ld#where the potential is evaluated
    #l0=0
    #l0=0
    #h=.5#height of headgroup for boundary condition, nm
    #S0=0.9#area per charge group, 1/nm^2

    sphere_sat=3.1
    cyl_sat=3.6
    planar_sat=4.0
    
    print(a,kappa,h)
    #def set_S0(S0,R):
    #    #different areas for appropriate boundary condition
    #    if shape=='sphere' or shape=="sphere_mu":
    #        S0=S0
    #    elif shape=='cyl' or shape=="cyl_mu":
    #        S0=S0/(1+(h/(h+2*R)))
    #    elif shape=='planar' or shape=='planar_mu':
    #        S0=S0/(1+(h/(h+2*R)))**2
    #    #elif shape=='iso':
    #    #    S0=lB**2
    #    else:
    #        print("Not a valid shape selection")
    #    return S0

    def set_S0(S0,R):
        #different areas for appropriate boundary condition
        if shape=='sphere' or shape=="sphere_mu":
            S0=S0
        elif shape=='cyl' or shape=="cyl_mu":
            S0=(R*S0/3)*2/radius_cyl
        elif shape=='planar' or shape=='planar_mu':
            S0=(R*S0/3)*4/thickness
        else:
            print("Not a valid shape selection")
        return S0

    
    potential_list=[]
    if shape=='sphere':
        print(shape)
        #S0=set_S0(S0,R)
        print(S0)
        result_list=[]
        #result_list2=[]
        for pH in pHrange:
            def func(variables):
                (alpha,P)=variables
                if P>sphere_sat:
                    eqn1=f1(alpha)-f2(pH,pKa)*np.exp(-psiS(l0,a,P))#+mu_ex(kappa,nacl,A=ent_scale))#+u(l0,a,kappa,S0,alpha))#-V_att_sphere(ld,R,H))+psiS(ld,a,P))
                    eqn2=BC_prefactor(lB,alpha,S0)+psiS_prime(0,a,P)
                else:
                    eqn1=f1(alpha)-f2(pH,pKa)*np.exp(-psiS(l0,a,P))#+mu_ex(kappa,nacl,A=ent_scale))#+u(l0,a,kappa,S0,alpha))#-V_att_sphere(ld,R,H))+psiS(ld,a,P))
                    eqn2=BC_prefactor(lB,alpha,S0)+psiS_prime(0,a,P)
                #eqn3=-BC_prefactor(lB,alpha,S0)+psiP_prime(x,P)
                return [eqn1,eqn2]#,eqn3]
            result = fsolve(func,(alpha0,P0))
            #if len(potential_list)==0 and result[1]<7:
            #print(result[1],psiS(l0,a,result[1]),l0)
            #print(result[1])
            if round(pH,4)==5.8:
                potential_list.append(result[1])
                print("pH={}".format(round(pH,4)))
                print("micelle potential={}".format(psiS(l0,a,result[1])))
            
                
            result_list.append(result[0])
            #result_list2.append(result[1])
        if verbose:
            print("{H,alpha}")
            for i in range(0,len(result_list)):
                H="{:.12f}".format(10**(-pHrange[i]))
                
                print("{"+str(H)+", "+ str(result_list[i])+"},")
        else:
            print("verbose=False")

        #xs=np.linspace(0,10*a,100)
        #ys=[psiS(x,a,potential) for x in xs]
        #plt.figure()
        #plt.scatter(xs,ys)
        #plt.show()
        return result_list
    elif shape=='cyl':
        print(shape)
        S0=set_S0(S0,R)
        print(S0)
        result_list=[]
        for pH in pHrange:
            def func(variables):
                (alpha,P)=variables
                if P>cyl_sat:
                    eqn1=f1(alpha)-f2(pH,pKa)*np.exp(-psiC(l0,a,P))#+mu_ex(kappa,nacl,A=ent_scale))#-V_att_cyl(ld,h,R,H))+psiC(ld,a,P))
                    eqn2=BC_prefactor(lB,alpha,S0)+psiC_prime(0,a,P)
                else:
                    eqn1=f1(alpha)-f2(pH,pKa)*np.exp(-psiC(l0,a,P))#+mu_ex(kappa,nacl,A=ent_scale))#-V_att_cyl(ld,h,R,H))+psiC(ld,a,P))
                    eqn2=BC_prefactor(lB,alpha,S0)+psiC_prime(0,a,P)
                #eqn3=-BC_prefactor(lB,alpha,S0)+psiP_prime(x,P)
                return [eqn1,eqn2]#,eqn3]
            result = fsolve(func,(alpha0,P0))
            #print(result[1])
            result_list.append(result[0])
        if verbose:
            print("{H,alpha}")
            for i in range(0,len(result_list)):
                H="{:.12f}".format(10**(-pHrange[i]))
                print("{"+str(H)+", "+ str(result_list[i])+"},")
        else:
            print("verbose=False")
        return result_list
    elif shape=='planar':
        print(shape)
        S0=set_S0(S0,R)
        print(S0)
        result_list=[]
        for pH in pHrange:
            def func(variables):
                (alpha,P)=variables
                if P>planar_sat:
                    eqn1=f1(alpha)-f2(pH,pKa)*np.exp(-psiP(l0,P))#+mu_ex(kappa,nacl,A=ent_scale))#-V_att_planar(ld,h,H))+psiP(ld,P))
                    eqn2=BC_prefactor(lB,alpha,S0)+psiP_prime(0,P)
                else:
                    eqn1=f1(alpha)-f2(pH,pKa)*np.exp(-psiP(l0,P))#+mu_ex(kappa,nacl,A=ent_scale))#-V_att_planar(ld,h,H))+psiP(ld,P))
                    eqn2=BC_prefactor(lB,alpha,S0)+psiP_prime(0,P)
                #eqn3=-BC_prefactor(lB,alpha,S0)+psiP_prime(x,P)
                return [eqn1,eqn2]#,eqn3]
            result = fsolve(func,(alpha0,P0))
            result_list.append(result[0])
        if verbose:
            print("{H,alpha}")
            for i in range(0,len(result_list)):
                H="{:.12f}".format(10**(-pHrange[i]))
                print("{"+str(H)+", "+ str(result_list[i])+"},")
        else:
            print("verbose=False")
        return result_list
    elif shape=='isolated':
        print(shape)
        #S0=set_S0(S0,R)
        result_list=[]
        for pH in pHrange:
            #def func(variables):
            #    (alpha,P)=variables
            #    eqn1=f1(alpha)-f2(pH,pKa)*np.exp(-psi_yukawa(lB/ld,P)+4*mu_ex(kappa,nacl))
            #    eqn2=BC_prefactor(lB,alpha,S0)+psi_yukawa_prime(lB/ld,P)
            #    return [eqn1,eqn2]#,eqn3]
            #result = fsolve(func,(alpha0,4))
            #print(result[1])
            #result_list.append(result[0])
            result_list.append(hh(pH,pKa))

        if verbose:
            print("{H,alpha}")
            for i in range(0,len(result_list)):
                H="{:.12f}".format(10**(-pHrange[i]))
                
                print("{"+str(H)+", "+ str(result_list[i])+"},")
        else:
            print("verbose=False")
        return result_list
    else:
        print("Not a valid shape selection")
        result_list=[]
        return result_list


def titration(Hs,alphas,label="N/A"):
    Ainit=0.008
    Kw=10**(-14)
    Va=5000
    Binit=0.1
    def Vb(H,alpha):
        return ((-H*(Ainit*(-1 + alpha) + H) + Kw)*Va)/(H*(Binit + H) - Kw)
    def Vb2(H,alpha):
        return -Va*(H**2 + H*Ainit*alpha - Kw)/(H*(H + Binit) - Kw) + (Ainit*Va)/Binit

    plt.scatter([Vb2(Hs[i],alphas[i]) for i in range(0,len(Hs))],[-np.log10(Hs[i]) for i in range(0,len(Hs))],label=label)

#plots titration back calculated along with a reference file
def titration2(Hs,alphas,fileName,Va,Ainit=0.008,Binit=0.1,Kw=10**(-14),equiv=250,label="N/A",color='#E41A1C',alpha=1):
    df=read_file(fileName,"0.1M NaOH")
    vol_lim=1.7
    equiv_index=equivalence_point(df)
    equiv_vol=df[df.columns[0]].values[equiv_index]
    #print(equiv_vol)
    df[df.columns[0]]= df[df.columns[0]]/equiv_vol*1.05#*0.95 #+/- 5% for the correct equiv volume
    df = df[df[df.columns[0]]<vol_lim]
    #plot_df(df)
    plt.scatter(df[df.columns[0]],df[df.columns[1]],color="#A65628")
    plt.xlabel(df.columns[0],fontsize=18)
    plt.ylabel(df.columns[1],fontsize=18)
    def Vb(H,alpha):
        return (((-H*(Ainit*(-1 + alpha) + H) + Kw)*Va)/(H*(Binit + H) - Kw))/equiv_vol
    def Vb2(H,alpha):
        return (-Va*(H**2 + H*Ainit*alpha - Kw)/(H*(H + Binit) - Kw) + (Ainit*Va)/Binit)/equiv_vol
    pHs_filt=[]
    Vb_filt=[]

    
    for i in range(1,len(alphas)):
        if all(Vb2(Hs[i],alphas[i])> Vb2(a,b) for a, b in zip(Hs[:i], alphas[:i])) and Vb2(Hs[i],alphas[i])>0:
            #make sure volume is greater than zero and check to make sure the next point is bigger than the one before (strictly increasig)
            pHs_filt.append(-np.log10(Hs[i]))
            Vb_filt.append(Vb2(Hs[i],alphas[i]))
        else:
            continue



        
    #print(all(i < j for i, j in zip(Vb_filt, Vb_filt[1:]))) # check if strictly increasing
    #plt.scatter([Vb2(Hs[i],alphas[i]) for i in range(0,len(Hs))],[-np.log10(Hs[i]) for i in range(0,len(Hs))],label=label)
    plt.scatter(Vb_filt,pHs_filt,label=label+"        ",c=color,alpha=alpha)
    plt.xlim(-0.1,vol_lim)
    
    #spline=us(Vb_filt,pHs_filt,s=0)
    #xs=np.linspace(Vb_filt[0],Vb_filt[-1],100)
    #plt.plot(xs,spline(xs))
    df_results=pd.DataFrame({'Vb_red':Vb_filt, 'pH':pHs_filt})
    return df_results
    
    
def run(pHrange,pKa,nacl,R,S0,h,shape='sphere',verbose=False):
    ld = 0.3/np.sqrt(nacl)
    kappa = (1/ld)
    a=R*kappa
    conc_a=0.008
    #if shape=="sphere" or shape=="sphere_mu":
    #    pHrange = np.arange(6.2,11,0.05)
    
    #if shape=="cyl" or shape=="cyl_mu":
    #    a=R*0.7*kappa
        
    #elif shape=="planar" or shape=="planar_mu":
    #    pHrange=np.arange(6.2,11,0.05)

    #adjust pKa of isolated lysine with change in salt concentration
    pKa=pH_activity(pKa,kappa,nacl)
    results=calculate(pKa,pHrange,a,nacl,S0,h,shape=shape,verbose=verbose)

    pH_corrected=[]
    results_filt=[]
    pHrange_filt=[]
    for i in range(0,len(results)):
        if results[i]<=1:
            results_filt.append(results[i])
            pHrange_filt.append(pHrange[i])
    for pH in pHrange_filt:
        pH_corrected.append(pH_activity(pH,kappa,nacl))
        
    #plt.scatter(pHrange_filt,results_filt,label=shape+", a="+str(np.round(a,2)))
    #plt.scatter(pHrange_filt,results_filt,label=shape)
    plt.scatter(pH_corrected,results_filt,label=shape)

        
    plt.xlabel("pH",fontsize=28)
    plt.ylabel("$\\alpha$, degree of ionization",fontsize=28)
    plt.title("[NaCl] = {}mM, a={}".format(np.round((nacl-conc_a)*1000,1),np.round(a,2)))
        
    #return [[10**(-pHrange[i]) for i in range(0,len(pHrange))],*results]
    return [[10**(-pHrange_filt[i]) for i in range(0,len(pHrange_filt))],*results_filt]
    #return [[10**(-pH_corrected[i]) for i in range(0,len(pH_corrected))],*results_filt]




main()

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
import pandas as pd
from scipy.optimize import curve_fit
from scipy.optimize import fsolve
from scipy.special import kv,kn
plt.style.use('veusz')



# constants
lB=0.7*10**(-9)
e=1.6*10**(-19)
eps0=8.85*10**(-12)
eps_r=80
kbT=4.1*10**(-21) #4.1 pN nM
conc_a=8 #mM



def psiS(x,a,P):
    return (2*a*(1 - np.exp(-3*a/10))*np.log((1 + np.exp(-x)*np.tanh(P/4))/(1 - np.exp(-x)*np.tanh(P/4))))/(a + x) + 2*np.exp(-3*a/10)*np.log((1 + np.exp(-x)*np.tanh((P*a)/((4*(a + x)))))/(1 - np.exp(-x)*np.tanh((P*a)/((4*(a + x))))))

def psiS_DH(x,a,P):
    return P*(a/(a+x))*np.exp(-x)




def psiS_prime(x,a,P):
    return -2*a*(1 - np.exp(-3*a/10))*np.log((1 + np.exp(-x)*np.tanh(P/4))/(1 - np.exp(-x)*np.tanh(P/4)))/(a + x)**2 + (2*(1 - np.exp(-x)*np.tanh((P*a)/((4*(a + x)))))*(-(1 + np.exp(-x)*np.tanh((P*a)/((4*(a + x)))))*((a*P*(1/np.cosh((P*a)/(4*(a+x))))**2*np.exp(-x))/((4*(a + x)**2)) + np.exp(-x)*np.tanh((P*a)/((4*(a + x)))))/(1 - np.exp(-x)*np.tanh((P*a)/((4*(a + x)))))**2 + ((-a*P*(1/np.cosh((P*a)/(4*(a+x))))**2*np.exp(-x))/((4*(a + x)**2)) - np.exp(-x)*np.tanh((P*a)/((4*(a + x)))))/(1 - np.exp(-x)*np.tanh((P*a)/((4*(a + x))))))*np.exp(-3*a/10))/(1 + np.exp(-x)*np.tanh((P*a)/((4*(a + x))))) + (2*a*(1 - np.exp(-x)*np.tanh(P/4))*(1 - np.exp(-3*a/10))*(-np.exp(-x)*np.tanh(P/4)/(1 - np.exp(-x)*np.tanh(P/4)) - (1 + np.exp(-x)*np.tanh(P/4))*np.exp(-x)*np.tanh(P/4)/(1 - np.exp(-x)*np.tanh(P/4))**2))/(((1 + np.exp(-x)*np.tanh(P/4))*(a + x)))

def psiC(x,a,P):
    return 2*np.exp(-3*a + x)*np.log((1 + np.exp(-x)*np.tanh((P*kn(0,(a + x)))/((4*kn(0,a))))/(1 - np.exp(-x)*np.tanh((P*kn(0,(a + x)))/((4*kn(0,a))))))) + (2*kn(0,a+x)*(1 - np.exp(-3*a))*np.exp(x)*np.log((1 + np.exp(-x)*np.tanh(P/4))/(1 - np.exp(-x)*np.tanh(P/4))))/((kn(0,a)))

def psiP(x,P):
    return 2*np.log((1+np.exp(-x)*np.tanh(P/4))/(1-np.exp(-x)*np.tanh(P/4)))

def sigma_sphere_fit(a,P,ld):
    return -e*psiS_prime(0,a,P)/(4*np.pi*lB*ld)

#USING GRAHAME EQUATION, ISRAELCHVILLI 1992, Ch 14
def sigma_grahame_israelchvilli(zeta,nacl,r):
    #zeta in mV
    zeta=zeta*1000
    return 0.117*np.sinh(zeta/51.4)*np.sqrt(nacl)

def sigma_grahame_curved(zeta,kappa,R):
    #having trouble with units here
    # zeta in V
    # applies for a>=0.5 (a=kappa*R)
    return (2*eps_r*eps0*kbT*kappa/e)*(np.sinh(zeta*e/(kbT*2))+(2/(kappa*R))*np.tanh((e*zeta)/(kbT*4)))

#USING GRAHAME EQUATION, ISRAELCHVILLI 1992, Ch 14
def sigma_grahame_flat(zeta,kappa):
    # zeta in V
    return (2*eps_r*eps0*kbT*kappa/e)*np.sinh(zeta*e/(kbT*2))


########################################################################

# ZETA POTENTIAL STUFF
#zetas=[42.74,57.33,62.1,59.84,54.97,40.98]
zetas=[42.74,57.33,62.1,59.84,54.97,40.98] #0 ,5, 10,25,50, 250 mM NaCl, measured by equipment
zetas=[42.74,57.33,62.1,59.84,54.97,40.98] #
#zetas=[41.8,52.7,62.3,59.8,55.0,48.5]

#for comparison
zetas_2=[53.44,52.70,62.61,60.6,53.4,29.5]#from single gaussian fits
zetas_2=[53.08,52.70,58.90,59.30,53.08,30.49]#from weighted 2 peak gaussian fits
zetas_2=[56.95,52.70,62.11,59.97,56.95,68.95]#from weighted 2 peak gaussian fits, weighted by peak area
zetas_2=[41.8,52.7,62.3,59.8,55.0,48.5]#from weighted 2 peak gaussian fits, weighted by peak area with bounds



zetas_2=[z/1000 for z in zetas_2]#V
zetas=[z/1000 for z in zetas]#V
conductivities=[0.347,1.054,1.179,2.584,4.755,19.79]#mS/cm
salts_zeta=[0,5,10,25,50,250]
salts_zeta=[s+conc_a for s in salts_zeta] #8 is for concentration of counter ions

lambdas_zeta=[0.3*10**(-9)/np.sqrt(s/1000) for s in salts_zeta]
kappas_zeta=[1/l for l in lambdas_zeta]


#FOR NACL CONC
rs_zeta=[0.15*(c)**(0.21)+1.62 for c in salts_zeta]
rs_zeta=[r*10**(-9) for r in rs_zeta]
hs_zeta=[-0.24*(c)**(0.14)+1.1 for c in salts_zeta]
hs_zeta=[h*10**(-9) for h in hs_zeta]


#FOR IONIC STRENGTH
rs_zeta=[0.52*(c)**(0.1)+1.24 for c in salts_zeta]
rs_zeta=[r*10**(-9) for r in rs_zeta]
hs_zeta=[-0.06*(c)**(0.28)+0.86 for c in salts_zeta]
hs_zeta=[h*10**(-9) for h in hs_zeta]


################################################################################################################
##############################################################################################################

#BASED ON ELECTROSTATIC MODEL
#salts=[0,5,50,250]
salts=list(np.linspace(0,18,18))
#salts=list(np.linspace(0,25,25))
salts.append(50)
salts.append(250)
salts=[s+conc_a for s in salts] #8 is for concentration of counter ions
rs=[0.15*(c)**(0.21)+1.6 for c in salts]
rs=[r*10**(-9) for r in rs]
lambdas=[0.3*10**(-9)/np.sqrt(s/1000) for s in salts]
kappas=[1/l for l in lambdas]
#salts for psi0s are 0,5,50,250 from electrostatic_model_noninteracting.py maximum value
#rhoh=450,rhotails=270,rhosolv=334
#psi0s=[5.05,5.33,5.88,6.18]
#apcs=[0.41234,0.40726, 0.3897327566176895, 0.36791198868260816]



#psi0s=[5.059472649478965,5.135941734488016,5.201928165329209,5.2580013830418775,5.2983653727402515,5.3419747305830985,5.380724498395619,5.415438692311508,5.4462546783138315,5.4752095540122125,5.493062032882116,5.521010540120405,5.545655062588104,5.556060116180537,5.574552946226468,5.591848110576301,5.616658883360142,5.632023680742685,5.6405137581782885,5.654154174069292,5.675693323419821,5.679351719476271,5.692243980705271,5.703379926871399,5.71401008693126]#, 5.88,6.18];

#20230404
#apcs=[0.42253298425399216, 0.4207685649580418, 0.4191910365272565, 0.41776417520962766, 0.41646140208965454,\
#      0.41526263043377104, 0.41415230575705025, 0.41311813465832536, 0.4121502307431008, 0.41124052336037287,\
#      0.4103823376374658, 0.40957008946348317, 0.40879905958946633, 0.4080652234102561, 0.4073651207163623, \
#      0.4066957546477167, 0.4060545123225869, 0.405439101786598, 0.40484750141034315, 0.4042779188955282, \
#      0.40372875777834405, 0.40319858984125057, 0.4026861322240351, 0.40219022830436757, 0.40170983162599233, \
#      0.39320136446296544, 0.3701844266957037]

apcs=[0.42253298425399216, 0.4207685649580418, 0.4191910365272565, 0.41776417520962766, 0.41646140208965454,\
      0.41526263043377104, 0.41415230575705025, 0.41311813465832536, 0.4121502307431008, 0.41124052336037287,\
      0.4103823376374658, 0.40957008946348317, 0.40879905958946633, 0.4080652234102561, 0.4073651207163623, \
      0.4066957546477167, 0.4060545123225869, 0.405439101786598, 0.39320136446296544, 0.3701844266957037]

pHs=[5.8,5.8,5.8,5.8,5.8,\
     5.8,5.8,5.8,5.8,5.8,\
     5.8,5.8,5.8,5.8,5.8,\
     5.8,5.8,5.8,5.8,5.8,\
     5.8,5.8,5.8,5.8,5.8,\
     5.8,5.8]
#AT SURFACE
psi0s=[4.973876021421224,
       5.08603670124291,
       5.175363383653271,
       5.248456628064572,
       5.309699095739395,
       5.361982487654276,
       5.407299283656255,
       5.4473027223589865,
       5.482600790587891,
       5.5141891396082245,
       5.542950820107754,
       5.568832727106386,
       5.592764850332888,
       5.614478209143536,
       5.634821424116785,
       5.653376902278928,
       5.670966878643306,
       5.687408897246772,
       5.868086284562173,
       6.191295410825807]

    
#pH=5.8
psi0s_dense=[4.9714024028416315,
             5.278383997884757,
       5.444991392011586,
       5.553641885579396,
       5.632412950406665,
       5.6925201726245955,
       5.740894321152391,
       5.780818426587852,
       5.814542520719794,
       5.843551333360791,
       5.872954356101171,
       5.895656874357402,
       5.915944124071621,
       5.934269856308285,
       5.950926516606541,
       5.966175816423687,
       5.98021306254818,
       5.993198853934921,
       6.005266753323928,
       6.016526403720161,
       6.027071638965548,
       6.036975850597383,
       6.046327413822146,
       6.0551251022515,
       6.063487428821352,
       6.071404456670089,
       6.078945988039045,
       6.08613445599701,
       6.092998447310869,
       6.099563612541925,
       6.105849395823427,
       6.111867242358255,
       6.11766497578128,
       6.12322444103964,
       6.128591363313143,
       6.13376559496252,
       6.138760469955284,
       6.143577189179929,
       6.148247308674593,
       6.152770097931885,
       6.157148575645126,
       6.16143010861322,
       6.165498788641396,
       6.169506568101272,
       6.173424578921408,
       6.17717855970513,
       6.1808378253266,
       6.184428180147557,
       6.187911783451833,
       6.191295410825807]


#apcs=[0.4123430267540477, 0.4111392047181438, 0.4100465393451609, 0.4090448401381566, 0.40811907564428984, 0.4072577362194517, 0.40645180226112193, 0.40569406672671976, 0.40497867477834615, 0.4043008017540605, 0.4036564222177759, 0.4030421407059458, 0.4024550653144417, 0.4018927116861335, 0.4013529289914734, 0.40083384209564726, 0.40033380582428457, 0.39985136839922664, 0.39938524191329994, 0.3989342782710571, 0.3984974494191184, 0.3980738309758091, 0.3976625885788436, 0.397262966424443, 0.39687427758696703]#,  0.3897327566176895, 0.36791198868260816]

apcs=[0.42253298425399216, 0.4153806696534619, 0.41055357860598596, 0.40689756170369956, 0.40395053365693984, 0.40147990205033735, 0.3993518241504846, 0.3974821138650452, 0.3958143300950349, 0.39430876737320253, 0.3929364217385807, 0.39167545479967814, 0.39050900859711724, 0.38942379541967104, 0.38840915432561846, 0.3874564003687997, 0.38655836383459874, 0.385709056530789, 0.38490342525480115, 0.38413716644457374, 0.3834065846392896, 0.38270848287525905, 0.3820400767371866, 0.381398926187203, 0.3807828809310566, 0.38019003621650954, 0.37961869675971477, 0.37906734706836404, 0.37853462684607386, 0.37801931046782533, 0.3775202897431875, 0.37703655935450014, 0.3765672044865142, 0.37611139026301615, 0.375668352682465, 0.3752373908042454, 0.37481785998389744, 0.3744091659926393, 0.37401075988590665, 0.373622133509192, 0.3732428155484403, 0.37287236804764246, 0.37251038332879605, 0.37215648125966794, 0.37181030682324, 0.3714715279497003, 0.3711398335776474, 0.3708149319160098, 0.370496548882234, 0.3701844266957037]



#MANNING CONDENSATION LAYER
manning=[np.sqrt(apc)/(lB*10**(9)) for apc in apcs]

a_list=[kappas[i]*rs[i] for i in range(0,len(rs))]



nacls_sat=np.linspace(0,250,50)
nacls_sat=[n+8 for n in nacls_sat]
lambdas_sat=[0.3*10**(-9)/np.sqrt(s/1000) for s in nacls_sat]
kappas_sat=[1/l for l in lambdas_sat]
#FOR IONIC STRENGTH
rs_sat=[0.52*(c)**(0.1)+1.24 for c in nacls_sat]
rs_sat=[r*10**(-9) for r in rs_sat]
a_list_sat=[kappas_sat[i]*rs_sat[i] for i in range(0,len(rs_sat))]


sigmas_model_dense=[sigma_sphere_fit(a_list_sat[i],psi0s_dense[i],lambdas_sat[i]) for i in range(0,len(a_list_sat))]



condensed_ion_layer=0.18
psi0s_sat=[psiS(condensed_ion_layer*10**(-9)/lambdas_sat[i],a_list_sat[i],psi0s_dense[i]) for i in range(0,len(psi0s_dense))]
print(psi0s_sat)


crossover=2
psi0s_mix=psi0s_dense[0:crossover]#ionic strength = 20mM
for i in psi0s_sat[crossover:]:
    psi0s_mix.append(i)
#print(psi0s_mix)

sigmas_model2=[sigma_sphere_fit(a_list_sat[i],psi0s_sat[i],lambdas_sat[i]) for i in range(0,len(a_list_sat))]

sigmas_model_mix=[sigma_sphere_fit(a_list_sat[i],psi0s_mix[i],lambdas_sat[i]) for i in range(0,len(a_list_sat))]



    
#to calculate location of diffuse layer
#diffuse layer is where \Psi=1 (i.e. e\psi=kT)
result_list=[]
gamma=1
for q in range(0,len(a_list)):
    def func(variables):
        (d)=variables
        eqn1= psiS(d,a_list[q],psi0s[q])-gamma
        return eqn1
    result = fsolve(func,1)
    result_list.append(result[0]*lambdas[q]*10**(9))
    #result_list.append(result[0])
#print("Diffuse layer location or lB_eff (nm): "+str(result_list))
#print("Diffuse layer location or lB_eff: "+str(result_list))



#PLOTTING POTENTIALS FOR DIFFERENT SALT CONCENTREATIONS
plt.figure()
xs=np.linspace(0,5,100)
for i in range(0, len(psi0s)):
    ys=[psiS(x,a_list[i],psi0s[i]) for x in xs]
    plt.plot([(x/kappas[i])/10**(-9) for x in xs],ys,label="{}".format(salts[i]-conc_a),linewidth=3)
plt.xlabel('r(nm)',fontsize=18)
plt.ylabel('$\Psi$',fontsize=18)
#plt.ylim(1)
plt.legend(title='[NaCl] (mM)',title_fontsize=18,fontsize=18)
plt.show()









######################################################################################################
######################################################################################################
#S0s (apls) used for model
#charge calculated as e*alpha/S
# for salts: 0,5,50,250
#S0s=[0.4123430267540477, 0.4074253126542586, 0.3897327566176895, 0.36791198868260816]
S0s=apcs
alpha=1
print(len(S0s))
sigmas_S0s=[e*alpha/(S0s[i]*10**(-18)) for i in range(0,len(S0s))]



#fig, ax1 = plt.subplots()
#ax1.scatter(salts,psi0s,s=100,color='k',marker='s')#,label='$\Psi$ (Left)')
#ax1.plot(salts,psi0s,linewidth=2,color='k',label='$\Psi_0$ (Left)')
##ax1.legend(loc=(0,1),fontsize=14)
#ax1.tick_params(axis='x',labelsize=18)
#ax1.tick_params(axis='y',labelsize=18)
##ax1.set_ylabel('Zeta Potential (mV)')
#ax2 = ax1.twinx()
#ax2.scatter(salts,[S*2 for S in S0s], s=100,color='k',marker='P')#,label='Area Per Lipid (Right)')
#ax2.plot(salts,[S*2 for S in S0s], linewidth=2,linestyle='dashed',color='k',label='Area Per Lipid (Right)')
##ax2.set_ylabel(None,size=30)
##ax2.legend(loc=(0.435,1),fontsize=14)
#ax2.tick_params(axis='y',labelsize=18)
#fig.legend(loc=(0.3,0.125),fontsize=14)
##fig.tight_layout()
#plt.show()


print(a_list_sat[crossover+1],psi0s_sat[crossover+1],nacls_sat[crossover+1])


fig, ax1 = plt.subplots()
#ax1.scatter(nacls_sat,psi0s_dense,s=100,color='k',marker='s')#,label='$\Psi$ (Left)')
ax1.scatter(nacls_sat[crossover+1],psi0s_sat[crossover+1],s=100,label='$\Psi(x=0.18/\lambda_D)_{peak}$ = (23.3, 4.16)')
ax1.plot(nacls_sat,psi0s_dense,linewidth=2,color='k',label='$\Psi_0$ (Left)',linestyle='--')
ax1.plot(nacls_sat,psi0s_sat,linewidth=2,color='k',label='$\Psi(x=0.18/\lambda_D)$ (Left)',linestyle='-.')
#ax1.legend(loc=(0,1),fontsize=14)
ax1.tick_params(axis='x',labelsize=18)
ax1.tick_params(axis='y',labelsize=18)
#ax1.set_ylabel('Zeta Potential (mV)')
ax2 = ax1.twinx()
#ax2.scatter(nacls_sat,[S*2 for S in S0s], s=100,color='k',marker='P')#,label='Area Per Lipid (Right)')
ax2.plot(nacls_sat,[S*2 for S in S0s], linewidth=2,color='k',label='Area Per Lipid (Right)')
#ax2.set_ylabel(None,size=30)
#ax2.legend(loc=(0.435,1),fontsize=14)
ax2.tick_params(axis='y',labelsize=18)
#fig.legend(loc=(0.3,0.125),fontsize=14)
fig.legend(loc=(0.3,0.525),fontsize=14)
#fig.tight_layout()
plt.show()


fig, ax1 = plt.subplots()
ax1.scatter(salts_zeta,[z*1000 for z in zetas],s=100,color='k',marker='s')#,label='Zeta Potential (Left)')
ax1.plot(salts_zeta,[z*1000 for z in zetas],linewidth=2,color='k',label='Zeta Potential (left)')
#ax1.legend(loc=(0,1),fontsize=14)
ax1.tick_params(axis='x',labelsize=18)
ax1.tick_params(axis='y',labelsize=18)
#ax1.set_ylabel('Zeta Potential (mV)')
ax2 = ax1.twinx()
ax2.scatter(salts_zeta,conductivities, s=100,color='k',marker='o')#,label='Conductivity (Right)')
ax2.plot(salts_zeta,conductivities, linewidth=2,linestyle='dashdot',color='k',label='Conductivity (right)')
#ax2.set_ylabel('Conductivity (mS/cm)')
ax2.tick_params(axis='y',labelsize=18)
#ax2.legend(loc=(0.532,1),fontsize=14)
#fig.tight_layout()
fig.legend(loc=(0.30,0.125),fontsize=14)
plt.show()





######################################################################################################
######################################################################################################


#CALCULATIONS

psi0s_zeta=[]
x=0.24*10**(-9) #slip plane
testxs=[]
testys=[]
for j in range(0,len(zetas)):
    def func(variables):
        (P)=variables
        eqn1= psiS((hs_zeta[j]+x)*kappas[j],rs_zeta[j]*kappas[j],P)-zetas[j]*e/(kbT)
        return eqn1
    result = fsolve(func,5)
    print((hs_zeta[j]+x),rs_zeta[j],result[0],zetas[j])
    testxs.append(hs_zeta[j]+x)
    print(psiS((hs_zeta[j]+x)*kappas[j],rs_zeta[j]*kappas[j],result[0])*kbT/(e))
    testys.append(psiS((hs_zeta[j]+x)*kappas[j],rs_zeta[j]*kappas[j],result[0])*kbT/(e))
    psi0s_zeta.append(result[0])
    #result_list.append(result[0])
    
psi0s_zeta_2=[]
for j in range(0,len(zetas)):
    def func(variables):
        (P)=variables
        eqn1= psiS((hs_zeta[j]+x)*kappas[j],rs_zeta[j]*kappas[j],P)-zetas_2[j]*e/(kbT)
        return eqn1
    result = fsolve(func,5)
    psi0s_zeta_2.append(result[0])
    #result_list.append(result[0])

    
    
#hs_flat=np.linspace(1.35,1.35,6)
    
psi0s_zeta_flat=[]
for j in range(0,len(zetas)):
    def func(variables):
        (P)=variables
        eqn1= psiP((hs_zeta[j]+x)*kappas[j],P)-zetas[j]*e/(kbT)
        #eqn1= psiP((hs_flat[j]*10**(-9)+x)*kappas[j],P)-zetas[j]*e/(kbT)
        return eqn1
    result = fsolve(func,5)
    psi0s_zeta_flat.append(result[0])


    
##salts_zeta=[0,5,10,25,50,250]
#diffuse_layer=[2.32,2.10,1.86,1.53,1.23,0.66]
#diffuse_layer=[d*10**(-9) for d in diffuse_layer]
#psi0s_diffuse=[]
#for j in range(0,len(zetas)):
#    def func(variables):
#        (P)=variables
#        eqn1= psiC((hs_zeta[j]+diffuse_layer[j])*kappas[j],rs_zeta[j]*kappas[j],P)-1
#        return eqn1
#    result = fsolve(func,4)
#    psi0s_diffuse.append(result[0])
#    print(result)


#psi0s_zeta_flat=[4*np.arctanh(np.exp(kappas[i]*(hs_zeta[i]+x))*np.tanh((zetas[i]*e/(kbT))/4)) for i in range(0,len(zetas))]
sigmas_zeta_curved=[sigma_grahame_curved((psi0s_zeta[i]*(kbT/e)),kappas_zeta[i],rs_zeta[i]) for i in range(0,len(rs_zeta))]
sigmas_zeta_curved_2=[sigma_grahame_curved((psi0s_zeta_2[i]*(kbT/e)),kappas_zeta[i],rs_zeta[i]) for i in range(0,len(rs_zeta))]
sigmas_zeta_flat=[sigma_grahame_flat((psi0s_zeta_flat[i]*(kbT/e)),kappas_zeta[i]) for i in range(0,len(rs_zeta))]
#grahame_diffuse=[sigma_grahame_curved(psi0s_diffuse[i]*(kbT/e),kappas_zeta[i],rs_zeta[i]) for i in range(0,len(psi0s_diffuse))]



sigmas_manning=[sigmas_S0s[i]*manning[i] for i in range(0,len(manning))]







plt.figure()
plt.plot(nacls_sat,sigmas_S0s,label="2$e$/APL",linewidth=4,linestyle="--",color='k')
plt.plot(nacls_sat,sigmas_manning,label="2$e\\alpha '$/APL",linewidth=4,linestyle="--",color='r')
print(sigmas_model_dense[-1],sigmas_model_dense[-2])
#plt.scatter(nacls_sat,sigmas_model_dense,label="model (bare charge)",linewidth=4,color='#72bf44',marker='P')
#plt.scatter(nacls_sat,sigmas_model2,label="model (eff charge)",linewidth=4,color='#72bf44',marker='x')
plt.scatter(salts_zeta,sigmas_zeta_flat,label="flat grahame",linewidth=4,color='#faa61a')
plt.scatter(salts_zeta,sigmas_zeta_curved,label="curved grahame",linewidth=4,color='#0066b3')
#plt.scatter(salts_zeta,sigmas_zeta_curved_2,linewidth=4,color='#0066b3')
#plt.plot(salts_zeta,grahame_diffuse,label="$\Psi$=1 (diffuse layer)",linewidth=4,color='#ed1c24')
#plt.scatter(salts_zeta,grahame_diffuse,label="$\zeta\equiv k_B T/e$",linewidth=4,color='#ed1c24')
plt.plot(nacls_sat,sigmas_model_dense,linewidth=4,color='#72bf44',linestyle='--')
plt.plot(salts_zeta,sigmas_zeta_flat,linewidth=4,color='#faa61a')
plt.plot(salts_zeta,sigmas_zeta_curved,linewidth=4,color='#0066b3')
plt.plot(nacls_sat,sigmas_model2,label="model (eff charge)",linewidth=4,color='#72bf44',linestyle='-')

#LETICIA's MCMD
plt.scatter([x+8 for x in [0,5,10,25,50,100,175,200]],[0.01,0.14,0.16,0.185,0.194,0.205,0.22,0.21],label="MC/MD",color='#72bf44')#,linewidth=4,color='#72bf44',linestyle='-.')
plt.plot([x+8 for x in [0,5,10,25,50,100,175,200]],[0.01,0.14,0.16,0.185,0.194,0.205,0.22,0.21],color='#72bf44',linewidth=4,linestyle='-.')
plt.xticks([0,50,100,150,200,250])
plt.yticks([0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7])


#plt.plot(nacls_sat,sigmas_model_mix,label="model (mix charge)",linewidth=5,color='k',alpha=0.3)
#print(rs,a_list,lB)
#plt.plot(salts,[(4.14*10**(-21)/e)*(rs[i]*(1+a_list[i])*1.5)/lB for i in range(0,len(a_list))],color='k')
#plt.plot(salts,[psi0s[i]*e/lB**2 for i in range(0,len(a_list))],color='k')
#plt.plot(salts_zeta,sigmas_zeta_curved_2,linewidth=2,color='#0066b3',linestyle="--",label="curved grahame Zeta peak")
#plt.plot(salts,grahame_diffuse,linewidth=4,color='#ed1c24')

plt.xlim(salts[0]-1,salts[-1]+1)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
#plt.plot(salts,[0.473 for s in salts], label="e/0.675 nm$^2$",linewidth=4,linestyle='--',color='k')
#plt.legend(fontsize=13)
plt.show()



plt.figure()
xs=np.linspace(0,5,100)
for i in range(0, len(psi0s)):
    ys=[psiS(x,a_list[i],psi0s[i]) for x in xs]
    plt.plot(xs,ys,label="a={}".format(a_list[i]))
plt.xlabel('x(kappa*r)')
plt.ylabel('$\Psi$')
#plt.ylim(1)
plt.legend()
plt.show()

#checking the different versions of grahame equation
# zeta in V, nacl in M to call functions
#print(sigma_grahame_flat(-0.075,0.15))
#print(sigma_grahame_israelchvilli(-0.075,0.15,0))
#print(sigma_grahame_curved(-0.075,0.15,2*10**(-9)))



plt.figure()
xs=np.linspace(0,5,100)
p0=5
a=0.5
#print(kappas[0])
ys1=[psiS(x,a,p0) for x in xs]
ys2=[psiC(x,a,p0) for x in xs]
ys3=[psiP(x,p0) for x in xs]
#plt.plot([(x/kappas[0])/10**(-9) for x in xs],ys1,label="Sphere",linewidth=4)
#plt.plot([(x/kappas[0])/10**(-9) for x in xs],ys2,label="Cylinder",linewidth=4)
#plt.plot([(x/kappas[0])/10**(-9) for x in xs],ys3,label="Planar",linewidth=4)
plt.plot(xs,ys1,label="Sphere",linewidth=4)
plt.plot(xs,ys2,label="Cylinder",linewidth=4)
plt.plot(xs,ys3,label="Planar",linewidth=4)
plt.xlabel('x ($\lambda_D$)',fontsize=24)
plt.ylabel('  (kT/e)',fontsize=24)
plt.legend(title='Assembly Shape',title_fontsize=24,fontsize=24)
plt.show()




################################################################################################3

print('5mM NaCl')
print(nacls_sat[1],lambdas_sat[1],a_list_sat[1],psi0s_dense[1],psi0s_sat[1])
plt.figure()
plt.title("[NaCl]=5mM")
p0=5.28
kappa=kappas_zeta[1]
r0=rs_zeta[1]
a=r0*kappa
xs2=np.linspace(0,0.98,100)
xs_nm=[x/kappa*10**(9) for x in xs2]
plt.plot(xs_nm,psiS(xs2,a,p0)*1000*kbT/e,label="Nonlinear",linewidth=4)
plt.plot(xs_nm,psiS_DH(xs2,a,p0)*1000*kbT/e,label="Linear (Debye Huckel)",linewidth=4)
#plt.plot(xs_nm,psiS(xs2,a,psi0s_zeta[1])*1000*kbT/e,label="Zeta")
#xzeta=0.977
#plt.plot([x for x in xs_nm if x <xzeta],[zetas[1]*1000 for x in xs_nm if x<xzeta],linestyle='--',color='k')
#plt.plot([xzeta-0.01 for x in xs_nm],np.linspace(0,zetas[1]*1000,len(xs_nm)),linestyle='--',color='k')
#plt.scatter(testxs[1]/10**(-9),testys[1]*1000,s=100,color='#4DAF4A')
#plt.plot([xzeta-0.01 for x in xs_nm],np.linspace(0,(p0+0.1)*25.4,len(xs_nm)),linestyle='--',color='k')
#xh=0.737/(1/kappa)
#plt.plot([x for x in xs_nm if x <xh*(1/kappa)],[psiS(xh,a,p0)*25.4 for x in xs_nm if x<xh*(1/kappa)],linestyle='--',color='k')
#plt.plot([xh*(1/kappa) for x in xs_nm],np.linspace(0,psiS(xh,a,p0)*25.4,len(xs_nm)),linestyle='--',color='k')
plt.plot([testxs[1]*10**(9) for x in xs_nm],np.linspace(0,(p0+0.1)*25.4,len(xs_nm)),linestyle='--',color='k')
plt.plot([testxs[1]*10**(9)-.24 for x in xs_nm],np.linspace(0,(p0+0.1)*25.4,len(xs_nm)),linestyle='--',color='k')
plt.legend(fontsize=14,title='PB Model')
plt.ylim(0,(p0+0.1)*25.4)
plt.ylim(0,(p0+0.1)*25.4)
plt.xlim(0,2.5)
#plt.xlabel('x(nm)')
#plt.ylabel('$\Psi$(kT)')
plt.show()

plt.figure()
print('250mM NaCl')
print(nacls_sat[-1],lambdas_sat[-1],a_list_sat[-1],psi0s_dense[-1],psi0s_sat[-1])
plt.title("[NaCl]=250mM")
p0=6.19
a=3.63
kappa=1.69
r0=2.146
xs2=np.linspace(0,5,50)
xs_nm=[x/kappa for x in xs2]
#xzeta=0.816
plt.plot(xs_nm,[psiS(x,a,p0)*25.4 for x in xs2],label="Nonlinear",color='#E41A1C',linestyle='--',linewidth=4)
plt.plot(xs_nm,[psiS_DH(x,a,p0)*25.4 for x in xs2],label="Linear (Debye Huckel)",color='#377EB8',linestyle='--',linewidth=4)
#plt.plot(xs2,[psiS(x,a,psi0s_zeta[-1])*25.4 for x in xs2],label="Zeta",color='#4DAF4A',linestyle='--')
#plt.scatter(testxs[-1]/10**(-9)*kappa,testys[-1]*1000*kbT/e,s=100,color='#4DAF4A')
plt.plot([testxs[-1]*10**(9) for x in xs_nm],np.linspace(0,(p0+0.1)*25.4,len(xs_nm)),linestyle='--',color='k')
plt.plot([testxs[-1]*10**(9)-.24 for x in xs_nm],np.linspace(0,(p0+0.1)*25.4,len(xs_nm)),linestyle='--',color='k')
plt.ylim(0,(p0+0.1)*25.4)
plt.ylim(0,(p0+0.1)*25.4)
plt.xlim(0,2.5)
#plt.xlabel('x(nm)')
#plt.ylabel('$\Psi$(kT)')
plt.show()



# DIFFERENT GEOMETRIES POTENTIAL AS A FUNCTION OF SALT
xs=np.linspace(0,5,100)
p0=5
a=[0.25,0.5,1]
kappas=[0.15,0.275,0.51]




plt.figure()
ys1=[psiS(x,a[0],p0) for x in xs]
ys2=[psiS(x,a[1],p0) for x in xs]
ys3=[psiS(x,a[2],p0) for x in xs]
plt.plot([x for x in xs],ys1,label="{}".format(a[0]),linewidth=4,color='#e41a1c')
plt.plot([x for x in xs],ys2,label="{}".format(a[1]),linewidth=4,color='#e41a1c',linestyle='dashed')
plt.plot([x for x in xs],ys3,label="{}".format(a[2]),linewidth=4,color='#e41a1c',linestyle='dotted')
plt.xlabel('x($\lambda_D$)',fontsize=24)
plt.ylabel('  (kT/e)',fontsize=24)
plt.legend(title='$a=R/\lambda_D$',title_fontsize=24,fontsize=24)
plt.title('Sphere')
plt.show()


ys1=[psiS(x,a[0],p0) for x in xs]
ys2=[psiS(x,a[1],p0) for x in xs]
ys3=[psiS(x,a[2],p0) for x in xs]
plt.plot([(x/kappas[0]) for x in xs],ys1,label="{}".format(a[0]),linewidth=4,color='#e41a1c')
plt.plot([(x/kappas[1]) for x in xs],ys2,label="{}".format(a[1]),linewidth=4,color='#e41a1c',linestyle='dashed')
plt.plot([(x/kappas[2]) for x in xs],ys3,label="{}".format(a[2]),linewidth=4,color='#e41a1c',linestyle='dotted')
plt.xlabel('r (nm)',fontsize=24)
plt.ylabel('  (kT/e)',fontsize=24)
plt.legend(title='$a=R/\lambda_D$',title_fontsize=24,fontsize=24)
plt.title('Sphere')
plt.show()

plt.figure()
ys1=[psiC(x,a[0],p0) for x in xs]
ys2=[psiC(x,a[1],p0) for x in xs]
ys3=[psiC(x,a[2],p0) for x in xs]
plt.plot([(x/kappas[0]) for x in xs],ys1,label="{}".format(a[0]),linewidth=4,color='#377eb8')
plt.plot([(x/kappas[1]) for x in xs],ys2,label="{}".format(a[1]),linewidth=4,color='#377eb8',linestyle='dashed')
plt.plot([(x/kappas[2]) for x in xs],ys3,label="{}".format(a[2]),linewidth=4,color='#377eb8',linestyle='dotted')
plt.xlabel('r (nm)',fontsize=24)
plt.ylabel('  (kT/e)',fontsize=24)
plt.legend(title='$a=R/\lambda_D$',title_fontsize=24,fontsize=24)
plt.title('Cylinder')
plt.show()

plt.figure()
ys1=[psiP(x,p0) for x in xs]
ys2=[psiP(x,p0) for x in xs]
ys3=[psiP(x,p0) for x in xs]
plt.plot([(x/kappas[0]) for x in xs],ys1,label="{}".format(kappas[0]),linewidth=4,color='#4daf4a')
plt.plot([(x/kappas[1]) for x in xs],ys2,label="{}".format(kappas[1]),linewidth=4,color='#4daf4a',linestyle='dashed')
plt.plot([(x/kappas[2]) for x in xs],ys3,label="{}".format(kappas[2]),linewidth=4,color='#4daf4a',linestyle='dotted')
plt.xlabel('r (nm)',fontsize=24)
plt.ylabel('  (kT/e)',fontsize=24)
plt.legend(title='$\kappa$',title_fontsize=24,fontsize=24)
plt.title('Planar')
plt.show()


plt.figure()
ax=plt.axes(projection="3d")
plt.xlim(salts[0]-1,salts[-1]+1)
x=salts
y=sigmas_model
z=psi0s
ax.plot3D(x,y,z,color='#72bf44',label="model",linewidth=4)
ax.tick_params(axis="x",labelsize=12)
ax.tick_params(axis="y",labelsize=12)
ax.tick_params(axis="z",labelsize=12)
ax.set_xlabel("Ionic Strength [mM]",fontsize=18)
ax.set_ylabel("$\sigma$",fontsize=18)
ax.set_zlabel("$\Psi_0$",fontsize=18)
ax.scatter3D(x, y, z, c=z, cmap='Blues',s=50)
plt.legend(fontsize=20)
plt.show()














########################################################



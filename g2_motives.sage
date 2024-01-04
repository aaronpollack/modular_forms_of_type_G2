#this file contains code used to compute Fourier coefficients of automorphic forms of type G_2, using results from the paper ``Exceptional theta functions and arithmeticity of modular forms on $G_2$"
#specifically, there are dual pairs G_2^c x Sp_6 in E_{7,3} and F_4^c x G_2 in E_{8,4}, and the code computes Fourier coefficients of level one exceptional theta lifts from algebraic modular forms on G_2^c and F_4^c to Sp_6 and G_2, respectively.
#See the accompanying document, ``Computation of Fourier coefficients of automorphic forms of type $G_2$" for more details about what this code does, and how it works
#in the case of G_2^c x Sp_6, the theta lifts are vector-valued holomorphic Siegel modular forms.  Fourier coefficients of vector valued holomorphic Siegel modular forms for Sp_6 are polynomials in 6 variables, v1, v2, v3, w1, w2, w3, one such polynomial for each half-integral symmetric matrix T.  The code can produce this polynomial when the upper left entry of T is 1.  (The code could be modified to handle general T, but this has not been implemented.)  The polynomials produced all have coefficents in Q(sqrt(-1)), so can be handled exactly by Sage
#in the case of F_4^c x G_2, the theta lifts are quaternionic modular forms.  Fourier coefficients of such forms are complex numbers, one for each integral binary cubic form au^3 + bu^2v + cuv^2 + dv^3.  The code can calculate this complex number, when the binary cubic has a=1.  (Again, the code could be modified to handle the general case, but this has not been implemented.). The Fourier coefficient produced always lives in Q(sqrt(-1)), so again this can be handled exactly by Sage

#HOW TO RUN THE G2^c-SP6 CODE:

#Tsplit_list=make_T_list_a1_vecs_split([1,b,c,d,e,f])
#a_Sp6_g_FC(Tsplit_list,k1,k2,[r1,r2,r3,r4,r5,r6]).expand()

#HOW TO RUN THE F4^c-G2 CODE FOR THE (J,I) THETA LIFT:

#the lists below list_1_2, list_1_3, list_2_1 are lists of numbers, at most of length 8.  they move Xoct and Yoct around by unipotent elements of F_4^c(K) to give a new singular pair
#X1,Y1=list_to_oct_pair_vec(list_1_2,list_1_3,list_2_1,Xoct,Yoct);

#the following computes the Fourier coefficients of the theta lifts Theta_I(m;X1,Y1) for m in the m_list.  The output is a dictionary, whose keys are tuples (b,c,d) corresponding to a monic binary cubic u^3+bu^2v+cuv^2+dv^3 with positive discriminant, b in {0,1} and b^2-2c <= n.
#my_coefficient_dict=G2_FC_dict_I(n, m_list, X1, Y1); 

#HOW TO RUN THE F4^c-G2 CODE FOR THE (J,E) THETA LIFT:

#pari.allocatemem(40737418240)
#n=5 #(you might run out of memory if n > 6)
#short_vecs=J_E_2.short_vector_list_up_to_length(n+1, up_to_sign_flag=False) #when n=5, this takes about 30 seconds on my laptop
#all_short_vecs=flatten(short_vecs,max_level=1)
#my_splitting_dictionary=vec_list_to_dict2(all_short_vecs) #this takes another 30 seconds on my laptop for n=5
#X1_E,Y1_E=list_to_oct_pair_E(list_1_2,list_1_3,list_2_1,Xoct,Yoct) #the list_1_2 etc must each have length at most 8
#G2_FC_dict_E(my_splitting_dictionary,X1_E,Y1_E,n,m_list) #this take about 7 minutes on my laptop for n=5

#if you want to know the dimension of the space of weight k >= 3, level one cuspidal quaternionic modular forms on G2, according to Rahul Dalal's explicit formula, call:
#Dalal_dim_k(k)

#THE FIELDS, POLYNOMIAL RINGS, AND QUATERNION ALGEBRAS USED TO DO AND OUTPUT COMPUTATIONS
H.<i,j,k>=QuaternionAlgebra(QQ,-1,-1);
K.<t>=QuadraticField(-1);
HK.<iK,jK,kK>=QuaternionAlgebra(K,-1,-1);

#create the polynomial ring which contains the Fourier coefficients of vector-valued Siegel modular forms on Sp6
BB.<v1,v2,v3,w1,w2,w3>=PolynomialRing(K)

#GENERAL CODE FOR OCTONIONS, BINARY CUBICS, AND THE EXCEPTIONAL CUBIC NORM STRUCTURE.

def disc_cubic(a,b,c,d):
    return -27*a^2*d^2+18*a*b*c*d+b^2*c^2-4*a*c^3-4*d*b^3

def oct_multiply(oct1,oct2,gamma):
    x1=oct1[0]
    y1=oct1[1]
    x2=oct2[0]
    y2=oct2[1]
    oct_prod=[x1*x2+gamma*y2.conjugate()*y1,y2*x1+y1*x2.conjugate()]
    return vector(oct_prod)

def oct_trace(oct1):
    return (oct1[0]).reduced_trace()

def oct_norm(oct1,gamma):
    return (oct1[0]).reduced_norm()-gamma*((oct1[1]).reduced_norm())

def oct_conjugate(oct1):
    x1=oct1[0]
    y1=oct1[1]
    return vector([x1.conjugate(),-y1])

def oct_pairing(oct1,oct2,gamma):
    oct3=oct_multiply(oct1,oct_conjugate(oct2),gamma)
    pairing=oct_trace(oct3)
    return pairing

def oct_trace_pairing(oct1,oct2,gamma):
    oct3=oct_multiply(oct1,oct2,gamma)
    pairing=oct_trace(oct3)
    return pairing

def oct_trilinear(oct1,oct2,oct3,gamma):
    prod=oct_multiply(oct_multiply(oct1,oct2,gamma),oct3,gamma)
    return oct_trace(prod)

def jordan_det(T,gamma):
    c1=T[0]
    c2=T[1]
    c3=T[2]
    x1=T[3]
    x2=T[4]
    x3=T[5]
    det=c1*c2*c3-c1*oct_norm(x1,gamma)-c2*oct_norm(x2,gamma)-c3*oct_norm(x3,gamma)+oct_trilinear(x1,x2,x3,gamma)
    return det

def jordan_sharp(T,gamma):
    c1=T[0]
    c2=T[1]
    c3=T[2]
    x1=T[3]
    x2=T[4]
    x3=T[5]
    b1=c2*c3-oct_norm(x1,gamma)
    b2=c3*c1-oct_norm(x2,gamma)
    b3=c1*c2-oct_norm(x3,gamma)
    y1=oct_multiply(oct_conjugate(x3),oct_conjugate(x2),gamma)-c1*x1
    y2=oct_multiply(oct_conjugate(x1),oct_conjugate(x3),gamma)-c2*x2
    y3=oct_multiply(oct_conjugate(x2),oct_conjugate(x1),gamma)-c3*x3
    Tsharp=[b1,b2,b3,y1,y2,y3]
    return Tsharp

def jordan_cross(T1,T2,gamma):
    Tsum=[T1[k]+T2[k] for k in range(6)]
    Tsum_sharp=jordan_sharp(Tsum,gamma)
    T1sharp=jordan_sharp(T1,gamma)
    T2sharp=jordan_sharp(T2,gamma)
    T1crossT2=[Tsum_sharp[k]-T1sharp[k]-T2sharp[k] for k in range(6)]
    return T1crossT2
    
def jordan_pairing(T1,T2,gamma):
    p1=T1[0]*T2[0]+T1[1]*T2[1]+T1[2]*T2[2]
    p2=oct_pairing(T1[3],T2[3],gamma)+oct_pairing(T1[4],T2[4],gamma)+oct_pairing(T1[5],T2[5],gamma)
    return p1+p2
    
def act_Phi(c,b,z,gamma):
    #-c x (b x z) + (c,z)b+(b,c)z
    V=jordan_cross(c,jordan_cross(b,z,gamma),gamma)
    p=jordan_pairing(c,z,gamma)
    q=jordan_pairing(b,c,gamma)
    T=[-V[k]+p*b[k]+q*z[k] for k in range(6)]
    return T
    
def act_wedge_Phi(c,b,z,gamma):
    V1=act_Phi(c,b,z,gamma)
    V2=act_Phi(b,c,z,gamma)
    V=[V1[k]-V2[k] for k in range(6)]
    return V

def oct_gram_basis(octonion_basis,gamma):
    basis_gram_matrix=matrix([[oct_pairing(octonion_basis[j],octonion_basis[k],gamma) for j in range(8)] for k in range(8)])
    return basis_gram_matrix
    
def oct_trace_pairing_matrix(octonion_basis,gamma):
    basis_gram_matrix=matrix([[oct_trace_pairing(octonion_basis[j],octonion_basis[k],gamma) for j in range(8)] for k in range(8)])
    return basis_gram_matrix
    
def vector_to_octonion(myVector,octonion_basis):
    octo=sum(myVector[j]*octonion_basis[j] for j in range(8))
    return octo

def octonion_to_vector(myOctonion,octonion_basis,gamma):
    basis_gram_matrix=oct_gram_basis(octonion_basis,gamma)
    inner_prods=vector([oct_pairing(myOctonion,octonion_basis[j],gamma) for j in range(8)])
    coords=basis_gram_matrix.inverse()*inner_prods
    return coords
    
    
#SPECIAL BASES AND MATRICES USED IN THE COMPUTATIONS
#create the Coxeter basis of the positive definite octonions
eETF=vector([H(0),H(1)]);hETF=1/2*vector([i+j+k,H(1)]);
coxETF0=oct_multiply(vector([j,H(0)]),hETF,-1);
coxETF1=eETF;
coxETF2=-hETF;
coxETF3=vector([j,H(0)]);
coxETF4=oct_multiply(vector([i,H(0)]),hETF,-1);
coxETF5=vector([H(1),H(0)]);
coxETF6=oct_multiply(eETF,hETF,-1);
coxETF7=oct_multiply(vector([k,H(0)]),eETF,-1);
coxETF=[coxETF0,coxETF1,coxETF2,coxETF3,coxETF4,coxETF5,coxETF6,coxETF7];
cox_Gram=matrix([[oct_pairing(coxETF[i],coxETF[j],-1) for i in range(8)] for j in range(8)])

oct1_vec=octonion_to_vector([H(1),H(0)],coxETF,-1)

# create the split basis of the positive definite octonions tensored up to K=Q(i)
e2ETF=1/2*(vector([HK(0),HK(1)])-t*vector([HK(0),iK]));
e3starETF=1/2*(vector([HK(0),jK])-t*vector([HK(0),kK]));
e3ETF=1/2*(vector([HK(0),-jK])-t*vector([HK(0),kK]));
e2starETF=1/2*(vector([HK(0),-HK(1)])-t*vector([HK(0),iK]));
ep1ETF=1/2*(vector([HK(1),HK(0)])-t*vector([iK,HK(0)]));
ep2ETF=1/2*(vector([HK(1),HK(0)])+t*vector([iK,HK(0)]));
e1ETF=1/2*(vector([jK,HK(0)])-t*vector([kK,HK(0)]));
e1starETF=1/2*(vector([-jK,HK(0)])-t*vector([kK,HK(0)]));
splitK_ETF=[ep1ETF,e1ETF,e2ETF,e3ETF,e1starETF,e2starETF,e3starETF,ep2ETF]

#rewrite the Coxeter basis in the octonion algebra tensored up to K
eETFK=vector([HK(0),HK(1)]);hETFK=1/2*vector([iK+jK+kK,HK(1)]);
coxETF0K=oct_multiply(vector([jK,HK(0)]),hETFK,-1);
coxETF1K=eETFK;
coxETF2K=-hETFK;
coxETF3K=vector([jK,HK(0)]);
coxETF4K=oct_multiply(vector([iK,HK(0)]),hETFK,-1);
coxETF5K=vector([HK(1),HK(0)]);
coxETF6K=oct_multiply(eETFK,hETFK,-1);
coxETF7K=oct_multiply(vector([kK,HK(0)]),eETFK,-1);
coxETFK=[coxETF0K,coxETF1K,coxETF2K,coxETF3K,coxETF4K,coxETF5K,coxETF6K,coxETF7K]

#a singular pair, used for the theta lifting from F4^c to G2
r1=1/2*vector([HK(0),1-t*iK]);
r2=1/2*vector([HK(0),1-t*iK]);
r3=-t*vector([iK,HK(0)]);
s1=-t*vector([HK(0),iK]);
s2=1/2*vector([HK(0),1+t*iK]);
s3=-1/2*vector([1+t*iK,HK(0)]);
R1=octonion_to_vector(r1,coxETFK,-1);
R2=octonion_to_vector(r2,coxETFK,-1);
R3=octonion_to_vector(r3,coxETFK,-1);
S1=octonion_to_vector(s1,coxETFK,-1);
S2=octonion_to_vector(s2,coxETFK,-1);
S3=octonion_to_vector(s3,coxETFK,-1);
Xoct=[1,-1,0,r1,r2,r3];
Yoct=[0,-1,1,s1,s2,s3];
XX=[1,-1,0,R1,R2,R3];
YY=[0,-1,1,S1,S2,S3];

#the Gram matrix for E8 root lattice, i.e., Q_E8_gram=oct_gram_basis(coxETF,-1)
Q_E8_gram=matrix([[2,0,0,0,0,-1,0,0],[0,2,-1,0,0,0,0,0],[0,-1,2,-1,0,0,0,0],[0,0,-1,2,-1,0,0,0],[0,0,0,-1,2,-1,0,0],[-1,0,0,0,-1,2,-1,0],[0,0,0,0,0,-1,2,-1],[0,0,0,0,0,0,-1,2]])

#the Gram matrix for the trace pairing of the Coxeter basis, i.e. Q_trace=oct_trace_pairing_matrix(coxETF,-1)
Q_trace=matrix([[-1,0,0,0,1,-1,1,0],[0,-2,1,0,0,0,0,0],[0,1,-2,1,0,0,0,0],[0,0,1,-2,1,0,0,0],[1,0,0,1,-1,-1,1,0], [-1,0,0,0,-1,2,-1,0],[1,0,0,0,1,-1,-1,1],[0,0,0,0,0,0,1,-2]])

#the transpose of the change of basis matrix from the coxeter basis to the split basis; one has R_cox_to_split=matrix([octonion_to_vector(coxETFK[j],splitK_ETF,-1) for j in range(8)])
R_cox_to_split=matrix([[1/2*t - 1/2, -1/2*t, 0, -1/2, -1/2*t, 0, 1/2, -1/2*t - 1/2],
 [0, 0, 1, 0, 0, -1, 0, 0],
 [-1/2*t, -1/2*t - 1/2, -1/2, 0, -1/2*t + 1/2, 1/2, 0, 1/2*t],
 [0, 1, 0, 0, -1, 0, 0, 0],
 [-1/2, 1/2*t - 1/2, 1/2*t, 0, 1/2*t + 1/2, 1/2*t, 0, -1/2],
 [1, 0, 0, 0, 0, 0, 0, 1],
 [-1/2, 0, -1/2*t, -1/2*t + 1/2, 0, -1/2*t, -1/2*t - 1/2, -1/2],
 [0, 0, 0, t, 0, 0, t, 0]])
 
#the Gram matrix for the norm pairing in the split basis, i.e., split_Gram=oct_gram_basis(splitK_ETF,-1)
split_Gram=matrix([[0,0,0,0,0,0,0,1],[0,0,0,0,-1,0,0,0],[0,0,0,0,0,-1,0,0],[0,0,0,0,0,0,-1,0],
                   [0,-1,0,0,0,0,0,0],[0,0,-1,0,0,0,0,0],[0,0,0,-1,0,0,0,0],[1,0,0,0,0,0,0,0]])
                   

#l1_mat is the matrix for the g2 nilpotent element E_{21} in the split basis
l1_mat=matrix([[0, 0, 0, 0, 0, 0, 0, 0],
 [0, 0, 0, 0, 0, 0, 0, 0],
 [0, 1, 0, 0, 0, 0, 0, 0],
 [0, 0, 0, 0, 0, 0, 0, 0],
 [0, 0, 0, 0, 0, -1, 0, 0],
 [0, 0, 0, 0, 0, 0, 0, 0],
 [0, 0, 0, 0, 0, 0, 0, 0],
 [0, 0, 0, 0, 0, 0, 0, 0]])
 
#l2_mat is the matrix for the g2 nilpotent element \delta_2 in the split basis
l2_mat=matrix([[0, 0, -1, 0, 0, 0, 0, 0],
 [0, 0, 0, 0, 0, 0, 1, 0],
 [0, 0, 0, 0, 0, 0, 0, 0],
 [0, 0, 0, 0, -1, 0, 0, 0],
 [0, 0, 0, 0, 0, 0, 0, 0],
 [1, 0, 0, 0, 0, 0, 0, -1],
 [0, 0, 0, 0, 0, 0, 0, 0],
 [0, 0, 1, 0, 0, 0, 0, 0]])
 
#here we create the other negative nilpotent elements:
l3_mat=l1_mat*l2_mat-l2_mat*l1_mat;
l4_mat=(l3_mat*l2_mat-l2_mat*l3_mat)/2;
l5_mat=(l4_mat*l2_mat-l2_mat*l4_mat)/3
l6_mat=l5_mat*l1_mat-l1_mat*l5_mat;


t3_cox_trilinear2=[matrix([[oct_trilinear(coxETF[j1],coxETF[j2],coxETF[j3],-1) for j3 in range(8)] for j2 in range(8)]) for j1 in range(8)]

def matrix_list_to_matrix(matrix_list,vector):
    L=len(vector)
    return sum(vector[j]*matrix_list[j] for j in range(L))
                   

#SOME HELPER FUNCTIONS FOR THE MAIN G2-SP6 FUNCTIONS
def mult_matrix(octo,oct_basis,gamma):
    mult_mat0=[]
    for j in range(8):
        vecj=octonion_to_vector(oct_multiply(octo,oct_basis[j],gamma),oct_basis,gamma)
        mult_mat0.append(vecj)
    mult_mat=matrix(mult_mat0).transpose()
    return mult_mat
    
list_of_mult_mats=[]
for k in range(8):
    Mk=mult_matrix(coxETF[k],coxETF,-1)
    list_of_mult_mats.append(Mk)
    
def mult_mat_vec_Cox(vec):
    M=sum([vec[j]*list_of_mult_mats[j] for j in range(8)])
    return M
    
def mult_cox_vecs(vec1,vec2):
    return mult_mat_vec_Cox(vec1)*vec2
    
def cox_vec_conjugate(vec):
    return (oct1_vec*Q_trace*vec)*oct1_vec-vec
    
#THE MAIN FUNCTIONS FOR G2-SP6

#this function takes in a list T0=[1,b,c,d,e,f] of integers, and finds all rank one T=[1,b,c,x1,x2,x3] in J_R whose projection is T0.  The output is a list of lists of size 3 [v1,v2,v3] where vj is the vector representation of xj in the ordered Coxeter basis

def make_T_list_a1_vecs_cox(T0):
    Tlist=[]
    Q_E8=QuadraticForm(ZZ,Q_E8_gram)
    a=T0[0]
    b=T0[1]
    c=T0[2]
    d=T0[3]
    e=T0[4]
    f=T0[5]
    if a!=1:
        raise ValueError("a must be 1 for this function")
    m=max(b,c)
    short_vecs=Q_E8.short_vector_list_up_to_length(m+1, up_to_sign_flag=False)
    vec_list_2=short_vecs[c]
    vec_list_3=short_vecs[b]
    for v2 in vec_list_2:
        if oct1_vec*Q_trace*v2==e:
            for v3 in vec_list_3:
                if oct1_vec*Q_trace*v3==f:
                    if v2*Q_trace*v3==d:
                        v1=cox_vec_conjugate(mult_cox_vecs(v2,v3))
                        Tlist.append([v1,v2,v3])
    return Tlist
    
#this function has the same input as make_T_list_a1_vecs_cox, but produces output in the ordered split basis instead of the ordered coxeter basis
    
def make_T_list_a1_vecs_split(T0):
    T_cox_list=make_T_list_a1_vecs_cox(T0)
    T_split_list=[]
    for T_vec in T_cox_list:
        v1_split=T_vec[0]*R_cox_to_split
        v2_split=T_vec[1]*R_cox_to_split
        v3_split=T_vec[2]*R_cox_to_split
        T_split_list.append([v1_split,v2_split,v3_split])
    return T_split_list
    

#this function makes the contribution to a Fourier coefficient on Sp6 of a single rank one T; here u,v are a null pair in the split basis, and u_pairing=u*split_Gram, v_pairing=v*split_Gram

def make_poly_Summand(Tvec,k1,k2,u_pairing,v_pairing):
    x1=Tvec[0]
    x2=Tvec[1]
    x3=Tvec[2]
    x1u=u_pairing*x1
    x2u=u_pairing*x2
    x3u=u_pairing*x3
    tvu=2*(x1u*v1+x2u*v2+x3u*v3)
    x1v=v_pairing*x1
    x2v=v_pairing*x2
    x3v=v_pairing*x3
    twuv=4*((x2u*x3v-x3u*x2v)*w1+(x3u*x1v-x1u*x3v)*w2+(x1u*x2v-x2u*x1v)*w3)
    return (tvu)^k1*(twuv)^k2
    
#this function computes the Fourier coefficient on Sp6 using the null pair u,v and a list of rank one T's in the split basis
    
def an_Sp6_FC(a_Tvec_list,k1,k2,u,v):
    u_pairing=u*split_Gram
    v_pairing=v*split_Gram
    total=0
    for Tvec in a_Tvec_list:
        total=total+make_poly_Summand(Tvec,k1,k2,u_pairing,v_pairing)
    return total
    
#this is our standard choice of null pair, corresponding to e3^* and e1, respectively
u_vec=vector([0,0,0,0,0,0,1,0]);v_vec=vector([0,1,0,0,0,0,0,0])


 
#this function creates a new null pair from a given one (U_vec, V_vec), using a list of 6 numbers r_vec
def new_pair(r_vec,U_vec,V_vec):
    nil=r_vec[0]*l1_mat+r_vec[1]*l2_mat+r_vec[2]*l3_mat+r_vec[3]*l4_mat+r_vec[4]*l5_mat+r_vec[5]*l6_mat
    g=exp(nil)
    return [g*U_vec,g*V_vec]
    
#this is the main function for the G2-Sp6 dual pair.  It starts with the standard null pair u_vec, v_vec, creates a new null pair using the list r_vec, then creates the Fourier coefficient on Sp6 from the new null pair and the list of rank one T's

def a_Sp6_g_FC(a_T_vec_list,k1,k2,r_vec):
    gpair=new_pair(r_vec,vector([0,0,0,0,0,0,1,0]),vector([0,1,0,0,0,0,0,0]))
    gFC=an_Sp6_FC(a_T_vec_list,k1,k2,gpair[0],gpair[1])
    return gFC
    
#HELPER FUNCTIONS FOR THE MAIN F4^c TO G2 FUNCTIONS
def Vnum(i,j,num):
    theta_basis=splitK_ETF
    if i==1:
        return [0,0,0,num*theta_basis[j],vector([HK(0),HK(0)]),vector([HK(0),HK(0)])]
    if i==2:
        return [0,0,0,vector([HK(0),HK(0)]),num*theta_basis[j],vector([HK(0),HK(0)])]
    if i==3:
        return [0,0,0,vector([HK(0),HK(0)]),vector([HK(0),HK(0)]),num*theta_basis[j]]
     
def list_to_UV_list(list_1_2,list_1_3,list_2_1):
    u_v_list=[]
    e11=[1,0,0,vector([HK(0),HK(0)]),vector([HK(0),HK(0)]),vector([HK(0),HK(0)])]
    e22=[0,1,0,vector([HK(0),HK(0)]),vector([HK(0),HK(0)]),vector([HK(0),HK(0)])]
    e33=[0,0,1,vector([HK(0),HK(0)]),vector([HK(0),HK(0)]),vector([HK(0),HK(0)])]
    L12=len(list_1_2)
    L13=len(list_1_3)
    L21=len(list_2_1)
    for j in range(L12):
        u_v_list.append([e11,Vnum(2,j,list_1_2[j])])
    for j in range(L13):
        u_v_list.append([e11,Vnum(3,j,list_1_3[j])])
    for j in range(L21):
        u_v_list.append([e22,Vnum(1,j,list_2_1[j])])
    return u_v_list
    
def Vj(j,my_oct):
    if j==1:
        return [0,0,0,my_oct,vector([HK(0),HK(0)]),vector([HK(0),HK(0)])]
    if j==2:
        return [0,0,0,vector([HK(0),HK(0)]),my_oct,vector([HK(0),HK(0)])]
    if j==3:
        return [0,0,0,vector([HK(0),HK(0)]),vector([HK(0),HK(0)]),my_oct]
        
def act_exp_Phi(c,b,z,gamma):
    #assuming Phi^3==0
    T0=z
    T1=act_Phi(c,b,T0,gamma)
    T2=act_Phi(c,b,T1,gamma)
    T=[T0[k]+T1[k]+T2[k]/2 for k in range(6)]
    return T    
    
def g_action2(u,v,Toct0,gamma):
    #assume Phi_{u wedge v}^3=0
    Toct1=act_wedge_Phi(u,v,Toct0,gamma)
    Toct2=act_wedge_Phi(u,v,Toct1,gamma)
    Toct_sum=Tsum=[Toct0[k]+Toct1[k]+Toct2[k]/2 for k in range(6)]
    return Toct_sum
    
def g_action_many(u_v_list,Toct,gamma):
    T=Toct
    L=len(u_v_list)
    for k in range(L):
        T=g_action2(u_v_list[k][0],u_v_list[k][1],T,gamma)
    return T
    
def Tvec_to_T(Tvec,theta_basis):
    T=[Tvec[0],Tvec[1],Tvec[2]]
    oct1=vector_to_octonion(Tvec[3],theta_basis)
    oct2=vector_to_octonion(Tvec[4],theta_basis)
    oct3=vector_to_octonion(Tvec[5],theta_basis)
    T.append(oct1)
    T.append(oct2)
    T.append(oct3)
    return T
    
def T_to_Tvec(T,theta_basis,gamma):
    v1=octonion_to_vector(T[3],theta_basis,gamma)
    v2=octonion_to_vector(T[4],theta_basis,gamma)
    v3=octonion_to_vector(T[5],theta_basis,gamma)
    Tvec=[T[0],T[1],T[2],v1,v2,v3]
    return Tvec
    
def list_to_oct_pair(list_1_2,list_1_3,list_2_1,X,Y):
    a_u_v_list=list_to_UV_list(list_1_2,list_1_3,list_2_1)
    gTX_oct=g_action_many(a_u_v_list,X,-1)
    gTY_oct=g_action_many(a_u_v_list,Y,-1)
    return [gTX_oct,gTY_oct]

    
    
#Suppose that f(u,v) = au^3+bu^2v+cuv^2+dv^3 is an integral binary cubic form with $a=1$.  Then there is GL_2(Z) translate of f(u,v) that has a=1 and b= 0 or b=1.  The following function finds all binary cubics f(u,v) with a=1, b in {0,1}, disc_cubic(f) > 0, and b^2 -2c at most the integer n.  (If lambda1, lambda2, lambda3 are the roots of f(u,1), then b^2-2c = lambda1^2 + lambda2^2 + lambda3^2.). It puts the tuples (b,c,d) into a dictionary as keys, with the corresponding value the 0 integer vector of length L.

def initialize_dict4(n,L):
    coeff_dict={}
    Blim=1
    Clim=n
    Dlim=floor((n/3)^(3/2))
    for b in range(0,Blim+1):
        for c in range(-Clim,Clim+1):
            for d in range(-Dlim,Dlim+1):
                if b^2-2*c <= n:
                    dis=disc_cubic(1,b,c,d)
                    if dis>0:
                        coeff_dict.update({(b,c,d):zero_vector(ZZ,L)})
    return coeff_dict
    
#SOME FUNCTIONS AND NOTATION FOR THE JORDAN PAIR (J,E):

beta_oct=1/2*vector([HK(-1)+iK+jK+kK,HK(1)+iK+jK+kK])
beta_oct_p_1=beta_oct+vector([HK(1),HK(0)])
oct1_oct=vector([HK(1),HK(0)])
E_oct=[2,2,2,beta_oct,beta_oct,beta_oct]
E_oct_sharp=[2,2,2,-(beta_oct_p_1),-(beta_oct_p_1),-(beta_oct_p_1)]

e11=[1,0,0,vector([HK(0),HK(0)]),vector([HK(0),HK(0)]),vector([HK(0),HK(0)])]
e22=[0,1,0,vector([HK(0),HK(0)]),vector([HK(0),HK(0)]),vector([HK(0),HK(0)])]
e33=[0,0,1,vector([HK(0),HK(0)]),vector([HK(0),HK(0)]),vector([HK(0),HK(0)])]
TcoxK=[e11,e22,e33]
for j in range(8):
    TcoxK.append(Vj(1,coxETFK[j]))
for j in range(8):
    TcoxK.append(Vj(2,coxETFK[j]))
for j in range(8):
    TcoxK.append(Vj(3,coxETFK[j]))

def delta_E(T_input):
    #implemented for quaternion algebra HK
    Temp1=act_exp_Phi(Vj(3,-beta_oct/2),e11,T_input,-1)
    Temp2=act_exp_Phi(e11,Vj(2,-(beta_oct_p_1)/2),Temp1,-1)
    Temp3=act_exp_Phi(e22,Vj(1,(-beta_oct)/2),Temp2,-1)
    Temp4=act_exp_Phi(Vj(1,oct1_oct),e22,Temp3,-1)
    Temp5=act_exp_Phi(e22,Vj(1,-oct1_oct),Temp4,-1)
    Temp6=act_exp_Phi(Vj(1,3/2*oct1_oct),e22,Temp5,-1)
    Temp7=act_exp_Phi(e22,Vj(1,-oct1_oct),Temp6,-1)
    return Temp7
    
def delta_E_inverse(T_input):
    #implemented for quaternion algebra HK
    Temp7=act_exp_Phi(e22,Vj(1,oct1_oct),T_input,-1)
    Temp6=act_exp_Phi(Vj(1,-3/2*oct1_oct),e22,Temp7,-1)
    Temp5=act_exp_Phi(e22,Vj(1,oct1_oct),Temp6,-1)
    Temp4=act_exp_Phi(Vj(1,-oct1_oct),e22,Temp5,-1)
    Temp3=act_exp_Phi(e22,Vj(1,(beta_oct)/2),Temp4,-1)
    Temp2=act_exp_Phi(e11,Vj(2,(beta_oct_p_1)/2),Temp3,-1)
    Temp1=act_exp_Phi(Vj(3,beta_oct/2),e11,Temp2,-1)
    return Temp1
    
def jordan_pairing_E(T1,T2,gamma):
    #(u,v)_E = (delta_E u, delta_E v)_I
    S1=delta_E(T1)
    S2=delta_E(T2)
    pairing=jordan_pairing(S1,S2,gamma)
    return pairing

def triple_norm(T1,T2,T3,gamma):
    T=jordan_cross(T1,T2,gamma)
    pair=jordan_pairing(T,T3,gamma)
    return pair
  
        
#THE MAIN FUNCTIONS FOR THE F4^c TO G2 LIFT

#this function produces an F_4^c(K)-translate of two elements X,Y in J_K.  If X,Y form a singular pair, then the output is acceptable to use as input into G2_FC_dict_I.  The input list_1_2, list_1_3, and list_2_1 must each be at most of length 8.

def list_to_oct_pair_vec(list_1_2,list_1_3,list_2_1,X,Y):
    gTX_oct,gTY_oct = list_to_oct_pair(list_1_2,list_1_3,list_2_1,X,Y)
    gTX_vec=T_to_Tvec(gTX_oct,coxETFK,-1)
    gTY_vec=T_to_Tvec(gTY_oct,coxETFK,-1)
    return [gTX_vec,gTY_vec]

#this function produces an F_4^c(K)-translate of two elements X,Y in J_K.  It then applies delta_E_inverse.  The output is acceptable to use at input into G2_FC_dict_E.  The input list_1_2, list_1_3, and list_2_1 must each be at most of length 8.

def list_to_oct_pair_E(list_1_2,list_1_3,list_2_1,X,Y):
    gTX_oct,gTY_oct=list_to_oct_pair(list_1_2,list_1_3,list_2_1,X,Y)
    gX_E=delta_E_inverse(gTX_oct)
    gY_E=delta_E_inverse(gTY_oct)
    return [gX_E,gY_E]

#for this function, x,y is a singular pair in the vector notation, such as x=XX and y=YY.  m_list is a list of non-negative integers. the function computes Fourier coefficients of Theta_I(m;x,y) for m in the m_list.  (This is a modular form on G2 of weight 4+m.)  It computes these Fourier coefficients for all positive discriminant monic integral binary cubics f(u,v) = u^3 + bu^2v+cuv^3+dv^3 with b in {0,1} and b^2-2c <= n.  The output is a dictionary with keys (b,c,d) and corresponding value a vector of the coefficients.  The j-th component of this vector corresponds to the j-th entry in m_list.  These coefficients are always elements of the quadratic field K.

def G2_FC_dict_I(n, m_list, x, y):
    L=len(m_list)
    coeff_dict = initialize_dict4(n,L)
    TVx3 = matrix_list_to_matrix(t3_cox_trilinear2, x[3])
    TVx4 = matrix_list_to_matrix(t3_cox_trilinear2, x[4])
    TVx5 = matrix_list_to_matrix(t3_cox_trilinear2, x[5])
    TVy3 = matrix_list_to_matrix(t3_cox_trilinear2, y[3])
    TVy4 = matrix_list_to_matrix(t3_cox_trilinear2, y[4])
    TVy5 = matrix_list_to_matrix(t3_cox_trilinear2, y[5])
    Qx3 = Q_E8_gram * x[3]
    Qx4 = Q_E8_gram * x[4]
    Qx5 = Q_E8_gram * x[5]
    Qy3 = Q_E8_gram * y[3]
    Qy4 = Q_E8_gram * y[4]
    Qy5 = Q_E8_gram * y[5]
    Q_E8 = QuadraticForm(ZZ, Q_E8_gram)
    short_vecs = Q_E8.short_vector_list_up_to_length(floor(n / 2) + 1, up_to_sign_flag=False)
    for n1 in range(floor(n / 2) + 1):
        for v1 in short_vecs[n1]:
            # do stuff that just depends on v1
            TVv1 = matrix_list_to_matrix(t3_cox_trilinear2, v1)
            v1Qx3 = v1 * Qx3
            v1Qy3 = v1 * Qy3
            TVy4_v1 = TVy4 * v1
            v1_TVy5 = v1 * TVy5
            TVx4_v1 = TVx4 * v1
            v1_TVx5 = v1 * TVx5
            for n2 in range(floor((n - 2 * n1) / 2) + 1):
                for v2 in short_vecs[n2]:
                    # do stuff that depends just on v1 and v2
                    TVv1v2 = v2 * TVv1
                    v2Qx4 = v2 * Qx4
                    v2Qy4 = v2 * Qy4
                    v2_TVy3 = v2 * TVy3  # slow
                    v1_TVy5_v2 = v1_TVy5 * v2
                    v2_TVx3 = v2 * TVx3  # slow
                    v1_TVx5_v2 = v1_TVx5 * v2
                    for n3 in range(floor((n - 2 * n1 - 2 * n2) / 2) + 1):
                        for v3 in short_vecs[n3]:
                            # do stuff that depends just on v1, v2, v3
                            TVv1v2_v3 = TVv1v2 * v3
                            v3Qx5 = v3 * Qx5
                            v3Qy5 = v3 * Qy5
                            TTy3 = v2_TVy3 * v3 + v3 * TVy4_v1 + v1_TVy5_v2
                            TTx3 = v2_TVx3 * v3 + v3 * TVx4_v1 + v1_TVx5_v2
                            temp1 = floor(sqrt(n - 2 * n1 - 2 * n2 - 2 * n3))
                            for u1 in range(-temp1, temp1 + 1):
                                temp2 = floor(sqrt(n - 2 * n1 - 2 * n2 - 2 * n3 - u1 ** 2))
                                for u2 in range(-temp2, temp2 + 1):
                                    temp3 = floor(sqrt(n - 2 * n1 - 2 * n2 - 2 * n3 - u1 ** 2 - u2 ** 2))
                                    for u3 in range(-temp3, temp3 + 1):
                                        d0 = u1 * u2 * u3 - u1 * n1 - u2 * n2 - u3 * n3
                                        TTx1 = (u2 * u3 - n1) * x[0] + (u3 * u1 - n2) * x[1] + (u1 * u2 - n3) * x[2]
                                        TTy1 = (u2 * u3 - n1) * y[0] + (u3 * u1 - n2) * y[1] + (u1 * u2 - n3) * y[2]
                                        Tx0 = u1 * x[0] + u2 * x[1] + u3 * x[2]
                                        Ty0 = u1 * y[0] + u2 * y[1] + u3 * y[2]
                                        b = u1 + u2 + u3
                                        c = u2 * u3 + u3 * u1 + u1 * u2 - (n1 + n2 + n3)
                                        d = d0 + TVv1v2_v3
                                        Tx = Tx0 + v1Qx3 + v2Qx4 + v3Qx5
                                        Ty = Ty0 + v1Qy3 + v2Qy4 + v3Qy5
                                        TTy2 = u1 * (v1Qy3) + u2 * (v2Qy4) + u3 * (v3Qy5)
                                        Tsharp_y = TTy1 - TTy2 + TTy3
                                        TTx2 = u1 * (v1Qx3) + u2 * (v2Qx4) + u3 * (v3Qx5)
                                        Tsharp_x = TTx1 - TTx2 + TTx3
                                        val = Tx * Tsharp_y - Ty * Tsharp_x
                                        val_vector=vector([val^(m_list[j]) for j in range(L)])
                                        if (b, c, d) in coeff_dict.keys():
                                            coeff_dict[(b, c, d)] = (coeff_dict[(b, c, d)] + val_vector)
    return coeff_dict
    
    
#QUADRATIC FORM (X,X)_E:

#J_E_gram=[[jordan_pairing_E(TcoxK[j],TcoxK[k],-1) for k in range(27)] for j in range(27)] #takes a while to compute, so this is hard-coded in 

J_E_gram=[[4,
  2,
  2,
  1,
  -1,
  2,
  -1,
  1,
  -3,
  3,
  -1,
  0,
  -2,
  4,
  -2,
  0,
  -2,
  4,
  -2,
  0,
  -2,
  4,
  -2,
  0,
  -2,
  4,
  -2],
 [2,
  4,
  2,
  0,
  -2,
  4,
  -2,
  0,
  -2,
  4,
  -2,
  1,
  -1,
  2,
  -1,
  1,
  -3,
  3,
  -1,
  0,
  -2,
  4,
  -2,
  0,
  -2,
  4,
  -2],
 [2,
  2,
  4,
  0,
  -2,
  4,
  -2,
  0,
  -2,
  4,
  -2,
  0,
  -2,
  4,
  -2,
  0,
  -2,
  4,
  -2,
  1,
  -1,
  2,
  -1,
  1,
  -3,
  3,
  -1],
 [1,
  0,
  0,
  4,
  0,
  0,
  0,
  0,
  -2,
  0,
  0,
  -1,
  -1,
  1,
  1,
  -1,
  0,
  1,
  0,
  -1,
  0,
  1,
  -2,
  1,
  0,
  1,
  -1],
 [-1,
  -2,
  -2,
  0,
  5,
  -4,
  1,
  0,
  1,
  -2,
  1,
  0,
  0,
  0,
  0,
  0,
  2,
  -4,
  2,
  -1,
  0,
  -3,
  2,
  -1,
  2,
  -1,
  0],
 [2,
  4,
  4,
  0,
  -4,
  8,
  -4,
  0,
  -2,
  4,
  -2,
  1,
  -3,
  3,
  -1,
  1,
  -4,
  5,
  -2,
  1,
  0,
  3,
  -2,
  1,
  -4,
  5,
  -2],
 [-1,
  -2,
  -2,
  0,
  1,
  -4,
  5,
  -2,
  1,
  -2,
  1,
  -2,
  2,
  -2,
  0,
  1,
  2,
  -2,
  0,
  1,
  0,
  -1,
  0,
  -1,
  2,
  -3,
  2],
 [1,
  0,
  0,
  0,
  0,
  0,
  -2,
  4,
  -2,
  0,
  0,
  1,
  -1,
  1,
  -1,
  -1,
  0,
  1,
  0,
  -1,
  0,
  1,
  1,
  -1,
  0,
  1,
  -1],
 [-3,
  -2,
  -2,
  -2,
  1,
  -2,
  1,
  -2,
  5,
  -4,
  1,
  0,
  2,
  -4,
  2,
  0,
  2,
  -4,
  2,
  0,
  2,
  -4,
  2,
  0,
  2,
  -4,
  2],
 [3,
  4,
  4,
  0,
  -2,
  4,
  -2,
  0,
  -4,
  8,
  -4,
  1,
  -1,
  5,
  -3,
  1,
  -4,
  5,
  -2,
  1,
  -4,
  5,
  -2,
  1,
  -4,
  5,
  -2],
 [-1,
  -2,
  -2,
  0,
  1,
  -2,
  1,
  0,
  1,
  -4,
  5,
  -1,
  0,
  -2,
  2,
  -1,
  2,
  -2,
  0,
  0,
  2,
  -2,
  0,
  0,
  2,
  -2,
  0],
 [0,
  1,
  0,
  -1,
  0,
  1,
  -2,
  1,
  0,
  1,
  -1,
  4,
  0,
  0,
  0,
  0,
  -2,
  0,
  0,
  -1,
  -1,
  1,
  1,
  -1,
  0,
  1,
  0],
 [-2,
  -1,
  -2,
  -1,
  0,
  -3,
  2,
  -1,
  2,
  -1,
  0,
  0,
  5,
  -4,
  1,
  0,
  1,
  -2,
  1,
  0,
  0,
  0,
  0,
  0,
  2,
  -4,
  2],
 [4,
  2,
  4,
  1,
  0,
  3,
  -2,
  1,
  -4,
  5,
  -2,
  0,
  -4,
  8,
  -4,
  0,
  -2,
  4,
  -2,
  1,
  -3,
  3,
  -1,
  1,
  -4,
  5,
  -2],
 [-2,
  -1,
  -2,
  1,
  0,
  -1,
  0,
  -1,
  2,
  -3,
  2,
  0,
  1,
  -4,
  5,
  -2,
  1,
  -2,
  1,
  -2,
  2,
  -2,
  0,
  1,
  2,
  -2,
  0],
 [0,
  1,
  0,
  -1,
  0,
  1,
  1,
  -1,
  0,
  1,
  -1,
  0,
  0,
  0,
  -2,
  4,
  -2,
  0,
  0,
  1,
  -1,
  1,
  -1,
  -1,
  0,
  1,
  0],
 [-2,
  -3,
  -2,
  0,
  2,
  -4,
  2,
  0,
  2,
  -4,
  2,
  -2,
  1,
  -2,
  1,
  -2,
  5,
  -4,
  1,
  0,
  2,
  -4,
  2,
  0,
  2,
  -4,
  2],
 [4,
  3,
  4,
  1,
  -4,
  5,
  -2,
  1,
  -4,
  5,
  -2,
  0,
  -2,
  4,
  -2,
  0,
  -4,
  8,
  -4,
  1,
  -1,
  5,
  -3,
  1,
  -4,
  5,
  -2],
 [-2,
  -1,
  -2,
  0,
  2,
  -2,
  0,
  0,
  2,
  -2,
  0,
  0,
  1,
  -2,
  1,
  0,
  1,
  -4,
  5,
  -1,
  0,
  -2,
  2,
  -1,
  2,
  -2,
  0],
 [0,
  0,
  1,
  -1,
  -1,
  1,
  1,
  -1,
  0,
  1,
  0,
  -1,
  0,
  1,
  -2,
  1,
  0,
  1,
  -1,
  4,
  0,
  0,
  0,
  0,
  -2,
  0,
  0],
 [-2,
  -2,
  -1,
  0,
  0,
  0,
  0,
  0,
  2,
  -4,
  2,
  -1,
  0,
  -3,
  2,
  -1,
  2,
  -1,
  0,
  0,
  5,
  -4,
  1,
  0,
  1,
  -2,
  1],
 [4,
  4,
  2,
  1,
  -3,
  3,
  -1,
  1,
  -4,
  5,
  -2,
  1,
  0,
  3,
  -2,
  1,
  -4,
  5,
  -2,
  0,
  -4,
  8,
  -4,
  0,
  -2,
  4,
  -2],
 [-2,
  -2,
  -1,
  -2,
  2,
  -2,
  0,
  1,
  2,
  -2,
  0,
  1,
  0,
  -1,
  0,
  -1,
  2,
  -3,
  2,
  0,
  1,
  -4,
  5,
  -2,
  1,
  -2,
  1],
 [0,
  0,
  1,
  1,
  -1,
  1,
  -1,
  -1,
  0,
  1,
  0,
  -1,
  0,
  1,
  1,
  -1,
  0,
  1,
  -1,
  0,
  0,
  0,
  -2,
  4,
  -2,
  0,
  0],
 [-2,
  -2,
  -3,
  0,
  2,
  -4,
  2,
  0,
  2,
  -4,
  2,
  0,
  2,
  -4,
  2,
  0,
  2,
  -4,
  2,
  -2,
  1,
  -2,
  1,
  -2,
  5,
  -4,
  1],
 [4,
  4,
  3,
  1,
  -1,
  5,
  -3,
  1,
  -4,
  5,
  -2,
  1,
  -4,
  5,
  -2,
  1,
  -4,
  5,
  -2,
  0,
  -2,
  4,
  -2,
  0,
  -4,
  8,
  -4],
 [-2,
  -2,
  -1,
  -1,
  0,
  -2,
  2,
  -1,
  2,
  -2,
  0,
  0,
  2,
  -2,
  0,
  0,
  2,
  -2,
  0,
  0,
  1,
  -2,
  1,
  0,
  1,
  -4,
  5]]

J_E_gram_mat=matrix(J_E_gram)
J_E_gram_mat2=2*J_E_gram_mat
J_E_2=QuadraticForm(ZZ,J_E_gram_mat2)
Q_E8=QuadraticForm(ZZ,Q_E8_gram)

pair_E_sharp=vector([jordan_pairing(E_oct_sharp,TcoxK[j],-1) for j in range(27)])
    
def vec_list_to_dict2(list_of_vecs):
    splitting_dict={}
    for vec in list_of_vecs:
        B=pair_E_sharp*vec
        if B>-1 and B<2:
            vec1=vec[0:11]
            vec1.set_immutable()
            vec2=vec[11:27]
            if vec1 in splitting_dict.keys():
                splitting_dict[vec1].append(vec2)
            else:
                splitting_dict.update({vec1:[vec[0:11],vec2]})
    return splitting_dict

    
#for this function, X,Y is a singular pair in the octonion notation, for the pointed cubic norm structure (J,E).  For example, one can take X=delta_E_inverse(Xoct) and Y=delta_E_inverse(Y_oct).  m_list is a list of non-negative integers. the function computes Fourier coefficients of Theta_E(m;x,y) for m in the m_list.  (This is a modular form on G2 of weight 4+m.)  It computes these Fourier coefficients for all positive discriminant monic integral binary cubics f(u,v) = u^3 + bu^2v+cuv^3+dv^3 with b in {0,1} and b^2-2c <= n.  The output is a dictionary with keys (b,c,d) and corresponding value a vector of the coefficients.  The j-th component of this vector corresponds to the j-th entry in m_list.  These coefficients are always elements of the quadratic field K.

def G2_FC_dict_E(a_split_dict,X,Y,n,m_list):
    L_m_list=len(m_list)
    coeff_dict=initialize_dict4(n,L_m_list)
    vec_X_1=octonion_to_vector(X[3],coxETFK,-1);vec_X1_Q=vec_X_1*Q_E8_gram
    vec_X_2=octonion_to_vector(X[4],coxETFK,-1);vec_X2_Q=vec_X_2*Q_E8_gram
    vec_X_3=octonion_to_vector(X[5],coxETFK,-1);vec_X3_Q=vec_X_3*Q_E8_gram
    vec_Y_1=octonion_to_vector(Y[3],coxETFK,-1);vec_Y1_Q=vec_Y_1*Q_E8_gram
    vec_Y_2=octonion_to_vector(Y[4],coxETFK,-1);vec_Y2_Q=vec_Y_2*Q_E8_gram
    vec_Y_3=octonion_to_vector(Y[5],coxETFK,-1);vec_Y3_Q=vec_Y_3*Q_E8_gram
    pair_E_X=vector([jordan_pairing_E(X,TcoxK[j],-1) for j in range(27)])
    pair_E_X0=pair_E_X[0:11];pair_E_X1=pair_E_X[11:27]
    pair_I_X=vector([jordan_pairing(X,TcoxK[j],-1) for j in range(27)])
    pair_I_X0=pair_I_X[0:11];pair_I_X1=pair_I_X[11:27]
    pair_E_Y=vector([jordan_pairing_E(Y,TcoxK[j],-1) for j in range(27)])
    pair_E_Y0=pair_E_Y[0:11];pair_E_Y1=pair_E_Y[11:27]
    pair_I_Y=vector([jordan_pairing(Y,TcoxK[j],-1) for j in range(27)])
    pair_I_Y0=pair_I_Y[0:11];pair_I_Y1=pair_I_Y[11:27]
    pair_I_E=vector([jordan_pairing(E_oct,TcoxK[j],-1) for j in range(27)])
    pair_E_sharp=vector([jordan_pairing(E_oct_sharp,TcoxK[j],-1) for j in range(27)])
    beta_vec=octonion_to_vector(beta_oct,coxETFK,-1)
    beta_star_vec=-octonion_to_vector(beta_oct_p_1,coxETFK,-1)
    beta_pair=beta_vec*Q_E8_gram
    beta_star_pair=beta_star_vec*Q_E8_gram
    beta_trilinear=matrix_list_to_matrix(t3_cox_trilinear2,beta_vec)
    X1_trilinear=matrix_list_to_matrix(t3_cox_trilinear2,vec_X_1)
    X2_trilinear=matrix_list_to_matrix(t3_cox_trilinear2,vec_X_2)
    X3_trilinear=matrix_list_to_matrix(t3_cox_trilinear2,vec_X_3)
    Y1_trilinear=matrix_list_to_matrix(t3_cox_trilinear2,vec_Y_1)
    Y2_trilinear=matrix_list_to_matrix(t3_cox_trilinear2,vec_Y_2)
    Y3_trilinear=matrix_list_to_matrix(t3_cox_trilinear2,vec_Y_3)
    for key in a_split_dict.keys():
        vec_list=a_split_dict[key]
        vec1_long=vec_list[0]
        L=len(vec_list)
        c1,c2,c3,x1=vec1_long[0],vec1_long[1],vec1_long[2],vec1_long[3:11]
        nx1=1/2*x1*Q_E8_gram*x1
        B0=pair_E_sharp[0:11]*vec1_long
        C0=2*(c1*c2+c2*c3+c3*c1-nx1)-c1*beta_pair*x1
        D0=c1*(c2*c3-nx1)
        X_E_T=pair_E_X0*vec1_long
        X_T_sharp=X[0]*(c2*c3-nx1)+X[1]*c1*c3+X[2]*c1*c2-c1*vec_X1_Q*x1
        Y_E_T=pair_E_Y0*vec1_long
        Y_T_sharp=Y[0]*(c2*c3-nx1)+Y[1]*c1*c3+Y[2]*c1*c2-c1*vec_Y1_Q*x1
        beta_trilinear_vec1=beta_trilinear*x1;vec1_beta_trilinear=x1*beta_trilinear
        vec1_trilinear=matrix_list_to_matrix(t3_cox_trilinear2,x1)
        for p in range(1,L):
            vecp=vec_list[p]
            vec2=vecp[0:8];vec3=vecp[8:16]
            nx2=1/2*vec2*Q_E8_gram*vec2;nx3=1/2*vec3*Q_E8_gram*vec3
            X_E_T_p=pair_E_X1*vecp;X_E_T_sum=X_E_T+X_E_T_p
            Y_E_T_p=pair_E_Y1*vecp;Y_E_T_sum=Y_E_T+Y_E_T_p
            B1=beta_star_pair*(vec2+vec3);B=B0+B1
            C2=-2*(nx2+nx3)+vec2*beta_trilinear*vec3 
            C1=beta_pair*(-c2*vec2-c3*vec3)+vec3*beta_trilinear_vec1+vec1_beta_trilinear*vec2
            C=C0+C1+C2
            D1=-c2*nx2-c3*nx3+vec2*vec1_trilinear*vec3;D=D0+D1
            X_T_p_sharp=-X[1]*nx2-X[2]*nx3+vec2*X1_trilinear*vec3
            Y_T_p_sharp=-Y[1]*nx2-Y[2]*nx3+vec2*Y1_trilinear*vec3
            X_T_T_p=-c2*vec_X2_Q*vec2-c3*vec_X3_Q*vec3+vec3*X2_trilinear*x1+x1*X3_trilinear*vec2
            Y_T_T_p=-c2*vec_Y2_Q*vec2-c3*vec_Y3_Q*vec3+vec3*Y2_trilinear*x1+x1*Y3_trilinear*vec2
            X_T_sharp_sum=X_T_sharp+X_T_T_p+X_T_p_sharp
            Y_T_sharp_sum=Y_T_sharp+Y_T_T_p+Y_T_p_sharp
            val=X_E_T_sum*Y_T_sharp_sum-Y_E_T_sum*X_T_sharp_sum
            val_vector=vector([val^(m_list[j]) for j in range(L_m_list)])
            if (B,C,D) in coeff_dict.keys():
                coeff_dict[(B,C,D)]=(coeff_dict[(B,C,D)]+val_vector)
    return coeff_dict

#Rahul Dalal's explicit formula for the dimension of the space of level one cuspidal quaternionic modular forms on G2 of weight k at least 3

def Dalal_dim_k(k):
    n=k-2
    term1=1/(12096*120)*(n+1)*(3*n+4)*(n+2)*(3*n+5)*(2*n+3)
    term2=1/(216*6)*(n+1)*(n+2)*(2*n+3)
    if n%2==0:
        term3=5/192*1/8*(n+2)*(3*n+4)
    if n%2==1:
        term3=5/192*1/8*(-1)*(n+1)*(3*n+5)
    if n%3==0:
        term4=1/18*(2*n/3+1)
    if (n%3)!=0:
        term4=1/18*(-1)*(floor(n/3)+1)
    if n%4==0:
        term5=1/32*(3*n/2+10)
    if n%4==1:
        term5=1/32*(6*floor(n/4)-4)
    if n%4==2 or n%4==3:
        term5=1/32*(-1)*(2*floor(n/4)+2)
    if n%6==0 or n%6==1:
        term6=1/24*(3*floor(n/6)+5)
    if n%6==2 or n%6==3:
        term6=1/24*(3*floor(n/6)-2)
    if n%6==4 or n%6==5:
        term6=1/24*(3*floor(n/6)+3)
    if n%7==0:
        term7=1/7*1
    if n%7==4:
        term7=-1/7
    if n%7!=0 and n%7!=4:
        term7=0
    if n%8==0:
        term8=1/4
    if n%8==5:
        term8=-1/4
    if n%8!=0 and n%8!=5:
        term8=0
    if n%12==0:
        term9=floor((n+2)/4)*(floor((n+2)/12)-1)
    if n%12==2 or n%12==4 or n%12==6 or n%12==8 or n%12==10:
        term9=floor((n+2)/4)*(floor((n+2)/12))
    if n%12==11:
        term9=-(floor((3*n+5)/12)-1)*(floor((n+3)/12)-1)
    if n%12==3 or n%12==7:
        term9=-(floor((3*n+5)/12)-1)*(floor((n+3)/12))
    if n%12==1 or n%12==5 or n%12==9:
        term9=-(floor((3*n+5)/12))*(floor((n+3)/12))
    return term1+term2+term3+term4+term5+term6+term7+term8+term9
    

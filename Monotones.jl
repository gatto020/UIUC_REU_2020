#Calculating J of Triangle
using LinearAlgebra
using QBase
#|0>
k0 = [1,0]
#k1
k1 = [0,1]
#|00>
k00 = [1,0,0,0]
#|11>
k11 = [0,0,0,1]
#|01>
k01 = [0,1,0,0]
#|10>
k10 = [0,0,1,0]
#|001>
k001 = kron(k00,k1)
#|000>
k000 = [1,0,0,0,0,0,0,0]
#|011>
k011 = [0,0,0,1,0,0,0,0]
#|101>
k101 = [0,0,0,0,0,1,0,0]
#|k110>
k110 = [0,0,0,0,0,0,1,0];
#|010>
k010 = [0,0,1,0,0,0,0,0]
#|111>
k111 = [0,0,0,0,0,0,0,1];
k000000 = kron(k000,k000)
k000101 = kron(k000,k101)
k010010 = kron(k010,k010)
k010111 = kron(k010,k111)
k101000 = kron(k101,k000)
k101101 = kron(k101,k101)
k111010 = kron(k111,k010)
k111111 = kron(k111,k111)
k000110 = kron(k000,k110)
k101110 = kron(k101,k110)
k010001 = kron(k010,k001)
k111001 = kron(k111,k001)

Triangle_ket = k000000 + k000110 + k101000 + k101110 + k010001 + k010111 + k111001 + k111111

Triangle = (1/8)*Triangle_ket*Triangle_ket'

eigenvalues_Triangle = eigvals(Triangle)

for i = 1:64
    if(eigenvalues_Triangle[i] > 0.001)
        print(i," ",eigenvalues_Triangle[i],"\n")
    end
end
#So this is really the only nonzero eigenvalue.This does not seem right. S(AA'BB'CC') = 0.
#How would I find S(AA',BB',CC')? Are they the same? I would say so.
#IS THIS RIGHT????

#Now trace out AA' to get S(BB'CC').
#Tracing out A first:
Triangle_A_traced_out = QBase.partial_trace(Triangle,[2,2,2,2,2,2],1)
#find a way to check if this worked

#Tracing out A' now:
Triangle_AA_traced_out = QBase.partial_trace(Triangle_A_traced_out,[2,2,2,2,2],1)
#find a way to check if this worked
#This is now the subsystem BB'CC'.

#Finding entropy 
QBase.Information.von_neumann_entropy(QBase.States.DensityMatrix(Triangle))#confere
S_BBCC = QBase.Information.von_neumann_entropy(QBase.States.DensityMatrix(Triangle_AA_traced_out))
#By symmetry, this is also S(AA'BB') and S(AA'CC'). 

#Sanity check of the eigenvalues of Triangle_AA'_traced_out
eigvals(Triangle_AA_traced_out) #checks out!



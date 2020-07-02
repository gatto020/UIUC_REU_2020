#Triangle bilocal network. Alice, Bob and Charlie share Bell states
#Each one of them makes parity measurements.
#Let's just try to test if the entropies were right first.
using QBase

#|0>
k0 = [1,0]
#k1
k1 = [0,1]
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

rho_clas_ABC = (1/4)*(k000*k000' + k011*k011' + k101*k101' + k110*k110');

rho_ABC = QBase.States.DensityMatrix(rho_clas_ABC)

S_ABC = QBase.Information.von_neumann_entropy(rho_ABC)

QBase.States.is_density_matrix(rho_clas_ABC)

using QBase: QMath
rho_BC = QMath.partial_trace(rho_clas_ABC,[2,2,2],1)#tracing out A

rho_C0 = QMath.partial_trace(rho_BC,[2,2],1)#tracing out B
rho_C = QBase.States.DensityMatrix(rho_C0)

S_C = QBase.Information.von_neumann_entropy(rho_C)

#Parity measurements
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

#Even Parity:
P0 = k00*k00' + k11*k11'
#Odd Parity
P1 = k10*k10' + k01*k01'

#Density matrix
k000000 = kron(k000,k000)
k000101 = kron(k000,k101)
k010010 = kron(k010,k010)
k010111 = kron(k010,k111)
k101000 = kron(k101,k000)
k101101 = kron(k101,k101)
k111010 = kron(k111,k010)
k111111 = kron(k111,k111)

g_state = (k000000 + k000101 + k010010 + k010111 + k101000 + k101101 + k111010 + k111111)

g_den_matrix0 = (1/8)*g_state*g_state'

#eigvals(g_den_matrix0) #if this is a pure state, should have given only 1 and 0 for eigvals

QBase.States.is_density_matrix(g_den_matrix0)

g_den_matrix = QBase.States.DensityMatrix(g_den_matrix0)

SqABC = QBase.Information.von_neumann_entropy(g_den_matrix)

#Checking if it is a pure state:
using LinearAlgebra
tr(g_den_matrix*g_den_matrix)

#Finding p(0,0,0)
P000 = kron(kron(P0,P0),P0)#has expected dimensions

p000 = tr(P000*g_den_matrix) #matches expectations!

P111 = kron(kron(P1,P1),P1)
p111 = tr(P111*g_den_matrix) #also matches expectations!

#Ok, so it works. What now? We see if we get the same results with the classical density matrix.
#Getting only diagonal terms:
class_den_matrix0 = ones(64,64)
n = 0
for i = 1:64
    for j = 1:64
        if(i==j)
            class_den_matrix0[i,j] = g_den_matrix[i,j];
            if(g_den_matrix[i,j] != 0)
                n = n+1; #n is the number of nonzero diagonal terms. we have 8.
            end            
        else
            class_den_matrix0[i,j] = 0;
        end
    end
end

class_den_matrix = QBase.States.DensityMatrix(class_den_matrix0)

p000_class = tr(P000*class_den_matrix)

p111_class = tr(P111*class_den_matrix)

P110 = kron(kron(P1,P1),P0)
P101 = kron(kron(P1,P0),P1)
P011 = kron(kron(P0,P1),P1)
P001 = kron(kron(P0,P0),P1);

p110_class = tr(P110*class_den_matrix)

p011_class = tr(P011*class_den_matrix)

p101_class = tr(P101*class_den_matrix)

p001_class = tr(P001*class_den_matrix)

dens_beta0 = (1/2)*(k00*k00' + k00*k11' + k11*k00' + k11*k11')
dens_beta1 = (1/2)*((k01 + k10)*(k01 + k10)')
dens_beta2 = (1/2)*((k00 - k11)*(k00 - k11)')
dens_beta3 = (1/2)*((k01 - k10)*(k01 - k10)')

Bell000 = kron(kron(dens_beta0,dens_beta0),dens_beta0)
Bell001 = kron(kron(dens_beta0,dens_beta0),dens_beta1)

p_Bell000 = tr(Bell000*g_den_matrix) #1/16

p_Bell001 = tr(Bell001*g_den_matrix)

outcomes = [dens_beta0,dens_beta1,dens_beta2,dens_beta3]

#Measurement matrix
function measurement(a,b,c)
    kron(kron(a,b),c)
end

#Testing if it works as I think it should.
testmeasurement = measurement(outcomes[1],outcomes[1],outcomes[1])
testp000B = tr(testmeasurement*g_den_matrix)
#It does.

#Calculates probabilities, given outcomes a,b,c (dens matrix form) and a density matrix
function prob(a,b,c,dens_matrix)
    tr(measurement(a,b,c)*dens_matrix)
end
prob(outcomes[1],outcomes[1],outcomes[1],g_den_matrix)#works just fine

#Calculating all probs!!!
bell_probs = ones(4,4,4)
for i = 1:4
    for j = 1:4
        for k = 1:4
            bell_probs[i,j,k] = prob(outcomes[i],outcomes[j],outcomes[k],g_den_matrix);
            print("p(",i,",",j,",",k,") = ",bell_probs[i,j,k],"\n")
        end
    end
end

           

#Calculating all probabilities for classical state
bell_probs_class = ones(4,4,4)
nonzero = 0;
for i = 1:4
    for j = 1:4
        for k = 1:4
            bell_probs_class[i,j,k] = prob(outcomes[i],outcomes[j],outcomes[k],class_den_matrix);
            print("p(",i,",",j,",",k,") = ",bell_probs_class[i,j,k],"\n")
            if(bell_probs_class[i,j,k] != 0)
                nonzero = nonzero+1;
            end
        end
    end
end


QBase.bell_kets

#Placeholder for angles. Will probably change in a loop or sth.
theta = 0; 
phi = 0;

#Generates probability distributions. It would be more useful if the dimensions were flexible,
#but right now this will do.
#Prints density matrices
function print_prob_dist(outcomes,dens_matrix) #outcomes will be density matrices
    probs_matrix = ones(4,4,4); #how to make dimensions flexible?
    for i = 1:4
        for j = 1:4
            for k = 1:4
                probs_matrix[i,j,k] = prob(outcomes[i],outcomes[j],outcomes[k],dens_matrix);
                print("p(",i,",",j,",",k,") = ",probs_matrix[i,j,k],"\n")
            end
        end
    end
end

#Generates probability distributions. Does not print.
function prob_dist(outcomes,dens_matrix) #outcomes need to be dens matrices
    probs_matrix = ones(4,4,4); 
    for i = 1:4
        for j = 1:4
            for k = 1:4
                probs_matrix[i,j,k] = prob(outcomes[i],outcomes[j],outcomes[k],dens_matrix);
            end
        end
    end
    return probs_matrix;
end

#Testing prob_dist function
print_prob_dist(outcomes,class_den_matrix)#seems to be working.

#first, let's see if it will work with our usual bell kets
theta = pi/4;
phi = pi/4;
#General Bell states
gen_bell_ket1 = cos(theta)*k00 + sin(theta)*k11
gen_bell_ket2 = -sin(theta)*k00 + cos(theta)*k11
gen_bell_ket3 = cos(phi)*k10 + sin(phi)*k01
gen_bell_ket4 = -sin(phi)*k10 + cos(phi)*k01

gen_bell_dens1 = gen_bell_ket1*gen_bell_ket1'
gen_bell_dens2 = gen_bell_ket2*gen_bell_ket2'
gen_bell_dens3 = gen_bell_ket3*gen_bell_ket3'
gen_bell_dens4 = gen_bell_ket4*gen_bell_ket4'
gen_bell_outcomes = [gen_bell_dens1,gen_bell_dens2,gen_bell_dens3,gen_bell_dens4]

print_prob_dist(gen_bell_outcomes,class_den_matrix) 
#Same distribution as before. I just changed the order of the bell basis.

#Returns probability distribution of making generalized Bell measurements on a given state. Prints.
function bell_distr_angles(theta,phi,dens_matrix)#input angles and state. 
    #General Bell states
    gen_bell_ket1 = cos(theta)*k00 + sin(theta)*k11
    gen_bell_ket2 = -sin(theta)*k00 + cos(theta)*k11
    gen_bell_ket3 = cos(phi)*k10 + sin(phi)*k01
    gen_bell_ket4 = -sin(phi)*k10 + cos(phi)*k01
    #Density matrices
    gen_bell_dens1 = gen_bell_ket1*gen_bell_ket1'
    gen_bell_dens2 = gen_bell_ket2*gen_bell_ket2'
    gen_bell_dens3 = gen_bell_ket3*gen_bell_ket3'
    gen_bell_dens4 = gen_bell_ket4*gen_bell_ket4'
    gen_bell_outcomes = [gen_bell_dens1,gen_bell_dens2,gen_bell_dens3,gen_bell_dens4]
    #Printing prob distr
    print_prob_dist(gen_bell_outcomes,dens_matrix) 
    return prob_dist(gen_bell_outcomes,dens_matrix);
end


function no_print_bell_distr_angles(theta,phi,dens_matrix) #Same as before, but does not print.
    #General Bell states
    gen_bell_ket1 = cos(theta)*k00 + sin(theta)*k11
    gen_bell_ket2 = -sin(theta)*k00 + cos(theta)*k11
    gen_bell_ket3 = cos(phi)*k10 + sin(phi)*k01
    gen_bell_ket4 = -sin(phi)*k10 + cos(phi)*k01
    #Density matrices
    gen_bell_dens1 = gen_bell_ket1*gen_bell_ket1'
    gen_bell_dens2 = gen_bell_ket2*gen_bell_ket2'
    gen_bell_dens3 = gen_bell_ket3*gen_bell_ket3'
    gen_bell_dens4 = gen_bell_ket4*gen_bell_ket4'
    gen_bell_outcomes = [gen_bell_dens1,gen_bell_dens2,gen_bell_dens3,gen_bell_dens4]
    
    return prob_dist(gen_bell_outcomes,dens_matrix);
end

testing = bell_distr_angles(pi/4,pi/4,class_den_matrix);

#I know there was already a function to calculate the entropy, but I wanted to do it myself for practice.
function total_entropy(probs_matrix)#WARNING: works only for these particular cases of classical distributions. NOT general entropy
    SUM = 0;
    for i = 1:4
        for j = 1:4
            for k = 1:4
                if (probs_matrix[i,j,k] != 0)
                    if(log2(probs_matrix[i,j,k]) != NaN)
                        SUM = SUM - probs_matrix[i,j,k]*log2(probs_matrix[i,j,k])
                    end
                end
            end
        end
    end
    return SUM;
end

total_entropy(testing)

function p1_outcome_distr_coeffs(probs_matrix)# NONZERO coefficients in outcome probability distribution
    SUM = 0*ones(4);
    for i = 1:4
        for j = 1:4
            for k = 1:4
                if (probs_matrix[i,j,k] != 0)
                    SUM[i] = SUM[i] + probs_matrix[i,j,k]
                end
            end
        end
    end
    return SUM; #ONLY NONZERO COEFFS
end

testing_coeffs = p1_outcome_distr_coeffs(testing)

function partial_entropy(coeffs)#gives the entropy of CLASSICAL cases where the coefficients are the probs and the eigenvalues.
    S = 0;
    for i =  1:length(coeffs)
        S = S - coeffs[i]*log2(coeffs[i])
    end
    return S;
end

partial_entropy(testing_coeffs)

#This gives the coefficients of the probability density matrix of the two-party subsystem
function p2_outcome_distr_coeffs(probs_matrix)# NONZERO coefficients in outcome probability distribution
    SUM = 0*ones(16);
    n = 1;
    for i = 1:4
        for j = 1:4
            for k = 1:4
                if (probs_matrix[i,j,k] != 0)
                    SUM[n] = SUM[n] + probs_matrix[i,j,k]
                end
            end
            n = n+1;
        end
    end
    return SUM; #ONLY NONZERO COEFFS
end

partial_entropy(p2_outcome_distr_coeffs(testing))

#3-party MMI 
function I3(probs_matrix)# WARNING: FOR TOTALLY SYMMETRIC BILOCAL NETWORKS ONLY, WITH CLASSICAL DISTRIBUTIONS
    return total_entropy(probs_matrix)-3*partial_entropy(p2_outcome_distr_coeffs(probs_matrix)) + 3*partial_entropy(p1_outcome_distr_coeffs(probs_matrix));
end
    

I3(testing) #checks out!

function print_gen_bell_I3(theta,phi,den_matrix)
    testing = bell_distr_angles(theta,phi,den_matrix);
    return I3(testing);
end

print_gen_bell_I3(pi/4,pi/4,class_den_matrix)

#print_gen_bell_I3(pi/4,pi/4,g_den_matrix) #também funciona!!

function gen_bell_I3(theta,phi,den_matrix)
    testing = no_print_bell_distr_angles(theta,phi,den_matrix);
    return I3(testing);
end

#gen_bell_I3(pi/4,pi/4,g_den_matrix) #THIS SHOULD NOT BE DONE. NOT CLASSICAL, SO IM NOT SURE THIS WORKS.

gen_bell_I3(pi/3,pi/6,class_den_matrix)

theta_axis = ones(100)#did not find a way to simply convert LinRange to an array
phi_axis = LinRange(0,pi/2,100)
for i = 1:100
    theta_axis[i] = i*pi/600+pi/6;
end


I3_axis = ones(99)
for i=2:100
    I3_axis[i-1] = gen_bell_I3(pi/4,theta_axis[i],g_den_matrix)
end
#can't do this because of the logs. Will need to find some representative angles instead.

print_gen_bell_I3(pi/3,pi/6,g_den_matrix) #DONT TRUST I3, but prob distr is right.

pi3pi6 = no_print_bell_distr_angles(pi/3,pi/6,g_den_matrix)

rhoA_pi3pi6 = p1_outcome_distr_coeffs(pi3pi6)#finding \rhoA

function rhoB_outcome_distr_coeffs(probs_matrix)# NONZERO coefficients in outcome probability distribution
    SUM = 0*ones(4);
    for j = 1:4
        for i = 1:4
            for k = 1:4
                if (probs_matrix[i,j,k] != 0)
                    SUM[j] = SUM[j] + probs_matrix[i,j,k]
                end
            end
        end
    end
    return SUM; #ONLY NONZERO COEFFS
end

rhoB_pi3pi6 = rhoB_outcome_distr_coeffs(pi3pi6)

function rhoC_outcome_distr_coeffs(probs_matrix)# NONZERO coefficients in outcome probability distribution
    SUM = 0*ones(4);
    for k = 1:4
        for i = 1:4
            for j = 1:4
                if (probs_matrix[i,j,k] != 0)
                    SUM[k] = SUM[k] + probs_matrix[i,j,k]
                end
            end
        end
    end
    return SUM; #ONLY NONZERO COEFFS
end

rhoC_pi3pi6 = rhoC_outcome_distr_coeffs(pi3pi6)

rhoAB_pi3pi6 = p2_outcome_distr_coeffs(pi3pi6)

function rhoBC_outcome_distr_coeffs(probs_matrix)# NONZERO coefficients in outcome probability distribution
    SUM = 0*ones(16);
    n = 1;
    for j = 1:4
        for k = 1:4
            for i = 1:4
                if (probs_matrix[i,j,k] != 0)
                    SUM[n] = SUM[n] + probs_matrix[i,j,k]
                end
            end
            n = n+1;
        end
    end
    return SUM; #ONLY NONZERO COEFFS
end

rhoBC_pi3pi6 = rhoBC_outcome_distr_coeffs(pi3pi6)

print_gen_bell_I3(pi/3,pi/6,class_den_matrix)

print_gen_bell_I3(pi/3,pi/6,g_den_matrix)

print_gen_bell_I3(pi/4,pi/6,g_den_matrix)

print_gen_bell_I3(pi/4,pi/6,class_den_matrix)

print_gen_bell_I3(pi/4,0,g_den_matrix)

print_gen_bell_I3(pi/4,0,class_den_matrix)

print_gen_bell_I3(0,0,g_den_matrix)

print_gen_bell_I3(0,0,class_den_matrix)

print_gen_bell_I3(0,pi/4,g_den_matrix)



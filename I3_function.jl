using QBase

#Calculates 3-party Multivariate Mutual Information for a given density matrix of a 3-party systemerror
#I3(A:B:C) = S(ABC)-S(AB)-S(BC)-S(AC)+S(A)+S(B)+S(C)
#Only works if subsystems share qubits (otherwise partial_trace would not work)
function I3(density_matrix,n) #n is the total number of qubits
    if (n == 3)
        marginal_AB = QBase.States.DensityMatrix(QBase.partial_trace(density_matrix,[2,2,2],3))
        marginal_BC = QBase.States.DensityMatrix(QBase.partial_trace(density_matrix,[2,2,2],1))
        marginal_AC = QBase.States.DensityMatrix(QBase.partial_trace(density_matrix,[2,2,2],2))
        marginal_A = QBase.States.DensityMatrix(QBase.partial_trace(marginal_AB,[2,2],2))
        marginal_B = QBase.States.DensityMatrix(QBase.partial_trace(marginal_BC,[2,2],2))
        marginal_C = QBase.States.DensityMatrix(QBase.partial_trace(marginal_AC,[2,2],1))
        i3 = QBase.Information.von_neumann_entropy(density_matrix) - QBase.Information.von_neumann_entropy(marginal_AB) - QBase.Information.von_neumann_entropy(marginal_BC)-QBase.Information.von_neumann_entropy(marginal_AC) + QBase.Information.von_neumann_entropy(marginal_A)+QBase.Information.von_neumann_entropy(marginal_B)+QBase.Information.von_neumann_entropy(marginal_C)
    elseif(n == 6)
        marginal_AABB = QBase.States.DensityMatrix(QBase.partial_trace(QBase.partial_trace(density_matrix,[2,2,2,2,2,2],6),[2,2,2,2,2],5))
        marginal_BBCC = QBase.States.DensityMatrix(QBase.partial_trace(QBase.partial_trace(density_matrix,[2,2,2,2,2,2],1),[2,2,2,2,2],1))
        marginal_AACC = QBase.States.DensityMatrix(QBase.partial_trace(QBase.partial_trace(density_matrix,[2,2,2,2,2,2],3),[2,2,2,2,2],3))
        marginal_AA = QBase.States.DensityMatrix(QBase.partial_trace(QBase.partial_trace(marginal_AABB,[2,2,2,2],4),[2,2,2],3))
        marginal_BB = QBase.States.DensityMatrix(QBase.partial_trace(QBase.partial_trace(marginal_BBCC,[2,2,2,2],4),[2,2,2],3))
        marginal_CC = QBase.States.DensityMatrix(QBase.partial_trace(QBase.partial_trace(marginal_AACC,[2,2,2,2],1),[2,2,2],1))
        i3 = QBase.Information.von_neumann_entropy(density_matrix) - QBase.Information.von_neumann_entropy(marginal_AABB) - QBase.Information.von_neumann_entropy(marginal_BBCC)-QBase.Information.von_neumann_entropy(marginal_AACC) + QBase.Information.von_neumann_entropy(marginal_AA)+QBase.Information.von_neumann_entropy(marginal_BB)+QBase.Information.von_neumann_entropy(marginal_CC)
    else
        print("Invalid number of qubits")
    end
    return i3 
end



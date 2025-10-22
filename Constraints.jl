module Constraints_Conserveds

include("Parameters.jl")
using .Parameters
using DelimitedFiles

export constraints_checker
function constraints_checker(Z,dZdt)

    P = dZdt

    #Constraints of CQC
    constraint_1 = conj.(P)*transpose(P) - P*P'
    constraint_2 = conj.(Z)*transpose(Z) - Z*Z'
    constraint_3 = im.*(Z*P' - conj.(Z)*transpose(P))

    #Check if the constraints satisfied
    c1 = all(abs(constraint_1[J,K]) < 1e-8 for J in 1:N2 , K in 1:N2)
    c2 = all(abs(constraint_2[J,K]) < 1e-8 for J in 1:N2 , K in 1:N2)
    c3 = all( J==K ? (abs(real(constraint_3[J,K])-1) < 1e-8  && abs(imag(constraint_3[J,K])) < 1e-8) 
    : abs(constraint_3[J,K]) < 1e-8 for J in 1:N2 , K in 1:N2)

    #Check if the constraints satisfied
    c1 = zero_checker(constraint_1)
    c2 = zero_checker(constraint_2)
    c3 = identity_checker(constraint_3)

    #Reporting
    if (c1==false || c2==false || c3==false)
        error("Constraint not satified! \n--> c1=$c1, c2=$c2 , c3=$c3 <-- Stopping the run." )
    end
       
    
    # #!Tester
    # open("data/constraints.dat","r") do io
    #     data = [readline(io) for _ in 1:4]
    #     constraint_3 .= [ parse(ComplexF64, split(data[j],"\t")[k]) for j in 1:4, k in 1:4 ]
    #     # open("data/anani.dat","w") do io2
    #     #     writedlm(io2,constraint_3)
    #     # end
    #     #Constraint-2   
    #     # c2 = all(abs(constraint_2[J,K]) < 1e-8 for J in 1:N2 , K in 1:N2)
    #     # println("Constraint-2: --> ",c2," <---")
    # end

    # open("data/constraints.dat","w") do io
    #     # write(io,"Constraint-1: \n")
    #     writedlm(io,constraint_3)
    # end

end




#Checks the conserved quantities of CQC.
#By CQC:  J=conserved_1=1 // J_bar=conserved_2=0 
export conserved_checker
function conserved_checker(Z,dZdt)

    P = dZdt

    #Conversed quantities of CQC
    conserved_1 = im.*(P'*Z - Z'*P)
    conserved_2 = im.*(P'*conj.(Z) - Z'*conj.(P))

    #Check their values 
    j1 = identity_checker(conserved_1)
    j2 = zero_checker(conserved_2)

    #Reporting
    if (j1==false || j2==false)
        error("Conserved quantity is not conserved! \n--> j1=$j1, j2=$j2  <-- Stopping the run." )
    end

end

#Checks if the given matrix is zero with a numerical tolerance of 1e-8
#Matrix M can be a complex, doesn't have to be real.
function zero_checker(M)
    return all(abs(M[J,K]) < 1e-8 for J in 1:(size(M,1)) , K in 1:size(M,2))
end

#Checks if the given matrix is an identity matrix with a numerical tolerance of 1e-8
#Matrix M can be complex, doesn't have to be real.
function identity_checker(M)
    return all( J==K ? (abs(real(M[J,K])-1) < 1e-8  && abs(imag(M[J,K])) < 1e-8) 
                : abs(M[J,K]) < 1e-8 for J in 1:size(M,1) , K in 1:size(M,2))
end

end
module FindVortex

include("Parameters.jl")
using .Parameters



export vortex_finder
function vortex_finder(ϕ)

    vortex_count = 0
    anti_vortex_count = 0
    #Go over each of the plaquettes:
    for j=lx:rx-1
        for k=ly:ry-1
           
            #extract the phase
            θ = zeros(4)
            θ[1] = angle(ϕ[j,k])
            θ[2] = angle(ϕ[j+1,k])
            θ[3] = angle(ϕ[j+1,k+1])
            θ[4] = angle(ϕ[j,k+1])
            
            #test the plaquette for winding number
            wind = winding_number(θ) #give 4 phases and get the winding number

            tolerance = 0.01
            #Check if the winding is a multiple of 2π
            if 1-tolerance  < wind < 1+tolerance
                println("Found a vortex!")
                vortex_count += 1
            elseif -1-tolerance  < wind < -1+tolerance
                println("Found an anti-vortex!")
                anti_vortex_count += 1
            # else
                # println("No vortex :(")
            end

        end
    end

    if vortex_count > 0  || anti_vortex_count > 0 
        println("Number of vortices: ",vortex_count)
        println("Number of anti-vortices: ", anti_vortex_count)
    else
        println("No vortex :(")
    end


end



#Calculates the winding number of a given plaquette.
#This routine can only find the windings of ±1 -> thus, this can only be used from bottom-up approach for finding vortices.
export winding_number #!remove after testing done.
function winding_number(θ)

    #Do the integral for a single plaquette
    integral = 0.
    for i=1:4

        #Calculate Δθ
        nn = circular_bc(i) #nn= next neighbour
        Δθ = θ[nn]-θ[i]

        #Check and correct for Δθ's not in the range of [-π,π]
        if Δθ < -π
            Δθ = Δθ + 2π
        elseif Δθ > π
            Δθ = Δθ - 2π
        end

        integral += Δθ/(2π)

    end

return integral
end




function circular_bc(i)
    if i==4
        next_element = 1
    else
        next_element = i+1
    end
return next_element
end




end #module
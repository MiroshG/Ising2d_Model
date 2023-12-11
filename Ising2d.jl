using Plots, Random, Statistics, LaTeXStrings, GLM, DataFrames, DelimitedFiles



# Creation of the lattice of the network and the initial condition ################################################

function Initial_Lattice(L)

    s=zeros(L^2,1)
    for i in 1:length(s)
        if rand()<0.5
            s[i]=1
        else 
            s[i]=-1
        end
    end
    return s
end


# function for the list of neighbours ##################################################################################

function neighbours(
    i, # central spin 
    L # dimensions of the square lattice
    )
    
    n_array=collect(0:L-1)
    up=0
    down=0
    left=0
    right=0

    # top and bottom boundaries

    if L < i <= L^2-L # the selected spin is not in the top or botton boundary

        up=i+L
        down=i-L

    elseif  i <= L # the selected spin is in the botton boundary

        up=i+L
        down=i+L^2-L

    elseif i > L^2-L             # the selected spin is in the top boundary
       
        up=i-L^2+L
        down=i-L

    end

    

    # left and right boundaries

    for n in n_array
        if n*L+1 < i < (n+1)*L # the selected spin is not in the left or right boundary
            
            right=i+1
            left=i-1
        
        elseif i==n*L+1           # the selected spin in in the left boundary

            right=i+1
            left=i+L-1
            
        elseif i==(n+1)*L             # te selected spin is in the right boundary 

            right=i-L+1
            left=i-1

        end
    end
    return(right, left, up, down)
end


# Applyin of the metropolis method to obtain the evolution of the system until a time t #########################################

function Metropolis_MC_Method(
    T, # temperature of the system (in units of J/k_B) [0,5] ΔT=0.1
    L, # dimensions of the square lattice L=(4, 8, 16, 32, 64)
    s  #lattice to use 
    )

    n=L^2

    for _ in 1:n
        Diff_E=0
        i=rand(1:n)                                      # Selecting a random spin of the Lattice
        right, left, up, down=neighbours(i, L )          # calling of the neighbours
        Diff_E=2*s[i]*(s[right]+s[left]+s[up]+s[down])
        if Diff_E<=0          # Energy is to high and need to get lower
            s[i]*=-1
        elseif rand()<min(1,exp(-Diff_E/T))
            s[i]*=-1
        end
    end

end



function System_Energy(s,L)
    H=0
    for i in 1:length(s)
        right, left, up, down=neighbours(i, L )
        H += -s[i]*(s[right]+s[left]+s[up]+s[down])
    end
    return H/(4)
end

function Magnetization(s,L)
    m=0
    m=sum(s)
    return(abs(m))
end

# [8,16,32,64]

L_array=[8,16,32,64,128]

termalization_time=1000 # checked


# TEMPERATURES TO USE #########################################################################################################

# Low temperatures

low_temperature_array_1= [x for x in 0.1:0.1:2]
#low_temperature_array_2= [x for x in 1:0.1:2]

# High temperatures

high_temperature_array= [x for x in 2.7:0.1:5]

# Medium temperature
medium_temperature_array= [x for x in 2.1:0.01:2.6]
#medium_temperature_array_less_than_Tc=[x for x in 2.1:0.01:2.25]

# Temperatures to use

temperature_array=append!(low_temperature_array_1, medium_temperature_array, high_temperature_array)
temperature_array_reverse=reverse(temperature_array)

#temperature_array_less_than_Tc=append!(low_temperature_array_2, medium_temperature_array_less_than_Tc)

temperature=[x for x in 0.1:0.1:5]
temperature_reverse=reverse(temperature)
temperature_2=[4,3,2,1]


# CALCULATION OF THE DIFFERENT PHYSICAL QUANTITIES #####################################################################################

function Quantities(L_array, temperature)
    
    # Final dimensions of the matrix: 
    # -columns:length L_array
    # -rows: length temperature_array

    Energy_matrix=zeros(length(temperature),0)
    Magnetization_matrix=zeros(length(temperature),0) 
    Specific_Heat_matrix=zeros(length(temperature),0) 
    Susceptibility_matrix=zeros(length(temperature),0)


    for L in L_array

        time_compute=@elapsed begin
            s=Initial_Lattice(L)    # Initialize de lattice for every L

            Energy_mean_array=[]                # array with the mean of H for each T
            Magnetization_mean_array=[]         # array with the mean of M for each T
            Specific_Heat_array=[]              # array with the mean of C for each T
            Susceptibility_array=[]             # array with the mean of X for each T

            for T in temperature #temperature_array_reverse   # Use the initialized lattice to evaluate in the higher T     

                beta= 1/T
                Energy_array=[]                 # array for H at temperature T
                Magnetization_array=[]          # array for m at temperature T
                Energy_square_array=[]          # array for H^2 at temperature T
                Magnetization_square_array=[]   # array for m^2 at temperature T

                for _ in 1:termalization_time
                    Metropolis_MC_Method(T,L,s)      # just update the system until thermalization
                end

                for _ in 1:10*termalization_time
                        
                    Metropolis_MC_Method(T,L,s)      # we update for a montecarlo step  

                    
                    H=System_Energy(s,L)        # with the configuration at the correlation time steps we obtain the energy H
                    m=Magnetization(s,L)          #  and magnetization m 
                    
                    push!(Energy_array, H)                      # save the values of H in each correlation time step
                    push!(Magnetization_array, m)               # save the values of m in each correlation time step
                    push!(Energy_square_array, H^2)             # save the values of H^2 in each correlation time step
                    push!(Magnetization_square_array, m^2)      # save the values of m^2 in each correlation time step

                end
                H_mean=mean(Energy_array)
                M_mean=mean(Magnetization_array)
                H_2_mean=mean(Energy_square_array)
                M_2_mean=mean(Magnetization_square_array)

                C=beta^2*( H_2_mean - H_mean^2)/L^2             # save the mean value of C for each T
                X=beta*( M_2_mean - M_mean^2)/L^2               # save the mean value of X for each T

                push!(Energy_mean_array, H_mean/L^2)                # save the mean value of H for each T
                push!(Magnetization_mean_array, M_mean/L^2)         # save the mean value of m for each T


                push!(Specific_Heat_array, C)
                push!(Susceptibility_array, X)
            end
            
            Energy_matrix=hcat(Energy_matrix, Energy_mean_array)                      # save the array with the mean values of H(T) for each L
            Magnetization_matrix=hcat(Magnetization_matrix, Magnetization_mean_array)        # save the array with the mean values of m(T) for each L
            Susceptibility_matrix=hcat(Susceptibility_matrix, Susceptibility_array)           # save the array with the mean values of X(T) for each L
            Specific_Heat_matrix=hcat(Specific_Heat_matrix, Specific_Heat_array)             # save the array with the mean values of C(T) for each L

        end
        println("L: ", L)
        println("time: ", time_compute)
    
    end
    return(Energy_matrix, Magnetization_matrix, Susceptibility_matrix, Specific_Heat_matrix)
end

Energy_matrix, Magnetization_matrix, Susceptibility_matrix, Specific_Heat_matrix=Quantities(L_array, temperature_array_reverse)

writedlm("Energy_matrix.txt", Energy_matrix)

writedlm("Magnetization_matrix.txt", Magnetization_matrix)

writedlm("Susceptibility_matrix.txt", Susceptibility_matrix)

writedlm("Specific_Heat_matrix.txt", Specific_Heat_matrix)



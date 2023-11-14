using Plots, Random



# Creation of the lattice of the network and the initial condition

function Lattice(L)

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


# function for the list of neighbours

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

neighbours(12,4)


# Applyin of the metropolis method to obtain the evolution of the system until a time t

function Metropolis_MC_Method(
    T, # temperature of the system (in units of J/k_B) [0,5] ΔT=0.1
    t, # units of time
    L # dimensions of the square lattice L=(4, 8, 16, 32, 64)
    )

    # initial parameters
    n=L^2
    MC_steps=10

    # calling of the function Lattice

    s=Lattice(L)

    # array for plotting the time vs the Energy

    H_array=[]
    time_array=[]


    # Updating of the system 

    for _ in 1:t

        for _ in 1:MC_steps*n
            s_i=0
            B_i=0
            i=rand(1:n)                                      # Selecting a random spin of the Lattice
            s_i=s[i]
            right, left, up, down=neighbours(i, L )          # calling of the neighbours
            B_i=s[i]*(s[right]+s[left]+s[up]+s[down])
            if B_i<0
                s[i]=-s[i]
                s_i=s[i]
            elseif (rand()<min(1, exp(-2*B_i/T)))
                s[i]=-s[i]
                s_i=s[i]
            end

        end

        # Energy of the system 

        H=s[1]*last(s)
        for i in 1:n-1
            H+=s[i]*s[i+1]
        end

        push!(H_array,H)
        push!(time_array, t)

    end

    plot(time_array, H_array, xlabel="t", ylabel="H", label="")

    savefig("Thermalization for T=$T, L=$L.png")
    
end

Metropolis_MC_Method(5, 1000, 4)
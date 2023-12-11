using Plots, Random, Statistics, LaTeXStrings, GLM, DataFrames, StatsBase

Random.seed!(0)

# Creation of the lattice of the network and the initial condition

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



# Applyin of the metropolis method to obtain the evolution of the system until a time t

function Metropolis_MC_Method(
    T, # temperature of the system (in units of J/k_B) [0,5] ΔT=0.1
    L, # dimensions of the square lattice L=(4, 8, 16, 32, 64)
    s  #lattice to use 
    )

    n=L^2

    for _ in n
        B_i=0
        i=rand(1:n)                                      # Selecting a random spin of the Lattice
        right, left, up, down=neighbours(i, L )          # calling of the neighbours
        B_i=2*s[i]*(s[right]+s[left]+s[up]+s[down])
        if B_i<=0          # Energy is to high and need to get lower
            s[i]*=-1
        elseif rand()<min(1,exp(-B_i/T))
            s[i]*=-1
        end
    end

end

function Magnetization(s,L)
    m=sum(s)/L^2
    return(abs(m))
end

function Termalization(T,L, t)

    #s=Initial_Lattice(L)
    s=ones(L^2)
    m_array=[]
    t_array=[]
    
    for i in 1:t

        Metropolis_MC_Method(T,L,s)

        m=Magnetization(s,L)
        push!(m_array, m)
        push!(t_array, i)

    end

    plot(t_array, m_array, xlabel="t", ylabel="m", ylim=(0, 1), label=false, grid=false, color=:black)
    savefig("Termalization.png")
    
    
end

Termalization(4, 16, 5000)

termalization_time=1000

# CALCULATION OF THE DIFFERENT PHYSICAL QUANTITIES ##################################################

function Quantities(L, T)
    
    # Final dimensions of the matrix: 
    # -columns:length L_array
    # -rows: length temperature_array

    Magnetization_array=[]

    s=Initial_Lattice(L)    # Initialize de lattice for every L 

    for _ in 1:termalization_time
        Metropolis_MC_Method(T,L,s)      # just update the system until thermalization
    end

    for _ in 1:10*termalization_time
            
        Metropolis_MC_Method(T,L,s)      # we update for a montecarlo step  

        m=Magnetization(s,L)          #  and magnetization m 
        
        push!(Magnetization_array, m)               # save the values of m in each correlation time step
    end
    return(Magnetization_array)
end

Magnetization_array=Quantities(16,3)

# Autocorrelation time for T=3, L=16 is 2500

step=[1,2,3,4,5,10,15,20,30,40,50,100,150,200,300]
function Correlation_time(Magnetization_array, step,L)
    mean_M=mean(Magnetization_array)
    corr_array = autocor(Magnetization_array.- mean_M^2, step ;demean=true)

    k=50
    new_corr_array = corr_array[1:11]
    new_step=step[1:11]
    
    model = lm(@formula(y ~ x), DataFrame(x=new_step , y=log.(new_corr_array)))

    fit_params = coef(model)
    slope, intercept = fit_params[2], fit_params[1]

    # Obtain the errors 
    errors = stderror(model)
    slope_error, intercept_error = errors[2], errors[1]

    r_squared=r2(model)

    println("correlation time: ", -1/(slope), "±", slope_error/(slope^2))
    println("r2= " , r_squared)
    println(" ")

    y=slope .*step .+ intercept

    plot(step, y, label="1/tc= $slope ± $slope_error")
    scatter!(step, log.(corr_array), xlabel=L"\tau", ylabel=L"\ln(\left\langle m \right\rangle)", grid=false, label="")
    savefig("Correlation_time.png")
end

Correlation_time(Magnetization_array, step, 16)

# CORRELATION LENGTH #####################################################################################################

function other_neighbours(
    i, # central spin 
    L, # dimensions of the square lattice
    m, # order of the neighbours
    )
    
    
    up=0
    down=0
    left=0
    right=0

    # top and bottom boundaries

    if L*m < i <= (L^2-L)-(m-1)*L # the selected spin is not in the m top or m botton rows

        up=i+m*L
        down=i-m*L


    elseif  i <= L*m # the selected spin is in the  m botton boundary

        up=i+m*L
        down=i+(L^2-L)-(m-1)*L

    elseif i > (L^2-L)-(m-1)*L             # the selected spin is in the top boundary
       
        up=i-(L^2-L)+(m-1)*L
        down=i-m*L

    end

    

    # left and right boundaries

    n_array=collect(0:L-1) # collect(m:L-m)

    for n in n_array
        if right!=0 && left!=0 
            return(right, left, up, down)
        else
            if n*L+m < i < (n+1)*L-(m-1) # the selected spin is not in the m left or m right columns
                
                right=i+m
                left=i-m
            
            elseif n*L< i <=n*L+m           # the selected spin in in the m left columns

                right=i+m
                left=i+L-m
                
            elseif (n+1)*L-(m-1)<= i <L*(n+1)+1           # te selected spin is in the m right columns 

                right=i-L+m
                left=i-m

            end
        end
    end
    
end

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

function System_after_quantities(L, temperature)

    s=Initial_Lattice(L)    # Initialize de lattice for every L
  

    for _ in 1:11*termalization_time
        Metropolis_MC_Method(T,L,s)      # just update the system until thermalization
    end



end

function Correlation_length(L, points, temperature)
    n=L^2

    for j in neighbours_distance
        for i in 1:points # loop over the number of points  I want to use to compute the correlation length


        end
    end

    
end


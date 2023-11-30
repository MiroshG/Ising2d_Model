using Plots, Random, Statistics, LaTeXStrings, GLM, DataFrames



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

a=neighbours(56,8)



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
    return H/(2)
end

function Magnetization(s,L)
    m=0
    m=sum(s)
    return(abs(m))
end

# [8,16,32,64]

L_array=[8,16,32,64,20,30,40]

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
    Fourth_order_cumulant_matrix=zeros(length(temperature),0)




    for L in L_array

        time_compute=@elapsed begin
            s=Initial_Lattice(L)    # Initialize de lattice for every L

            Energy_mean_array=[]                # array with the mean of H for each T
            Magnetization_mean_array=[]         # array with the mean of M for each T
            Specific_Heat_array=[]              # array with the mean of C for each T
            Susceptibility_array=[]             # array with the mean of X for each T
            Fourth_order_cumulant_mean_array=[]

            for T in temperature #temperature_array_reverse   # Use the initialized lattice to evaluate in the higher T     

                beta= 1/T
                Energy_array=[]                 # array for H at temperature T
                Magnetization_array=[]          # array for m at temperature T
                Energy_square_array=[]          # array for H^2 at temperature T
                Magnetization_square_array=[]   # array for m^2 at temperature T
                Magnetization_4=[]

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
                    push!(Magnetization_4, m^4)

                end
                H_mean=mean(Energy_array)
                M_mean=mean(Magnetization_array)
                H_2_mean=mean(Energy_square_array)
                M_2_mean=mean(Magnetization_square_array)

                C=beta^2*( H_2_mean - H_mean^2)/L^2             # save the mean value of C for each T
                X=beta*( M_2_mean - M_mean^2)/L^2               # save the mean value of X for each T

                push!(Energy_mean_array, H_mean/L^2)                # save the mean value of H for each T
                push!(Magnetization_mean_array, M_mean/L^2)         # save the mean value of m for each T

                Fourth_order=mean(Magnetization_4)/mean(Magnetization_square_array)^2

                push!(Specific_Heat_array, C)
                push!(Susceptibility_array, X)
                push!(Fourth_order_cumulant_mean_array, Fourth_order)
            end
            
            Energy_matrix=hcat(Energy_matrix, Energy_mean_array)                      # save the array with the mean values of H(T) for each L
            Magnetization_matrix=hcat(Magnetization_matrix, Magnetization_mean_array)        # save the array with the mean values of m(T) for each L
            Susceptibility_matrix=hcat(Susceptibility_matrix, Susceptibility_array)           # save the array with the mean values of X(T) for each L
            Specific_Heat_matrix=hcat(Specific_Heat_matrix, Specific_Heat_array)             # save the array with the mean values of C(T) for each L
            Fourth_order_cumulant_matrix=hcat(Fourth_order_cumulant_matrix, Fourth_order_cumulant_mean_array)
        end
        println("L: ", length(s))
        println("time: ", time_compute)
    
    end
    return(Energy_matrix, Magnetization_matrix, Susceptibility_matrix, Specific_Heat_matrix, Fourth_order_cumulant_matrix)
end

Energy_matrix, Magnetization_matrix, Susceptibility_matrix, Specific_Heat_matrix, Fourth_order_cumulant_matrix=Quantities(L_array, temperature_array_reverse)

# CRITICAL TEMPERATURE ########################################################################################################

function Critical_temperature(Specific_Heat_matrix, Susceptibility_matrix, L_array, temperature_array_reverse)
    Inverse_L=1 ./ L_array
    Tc_L_C=[]
    Tc_L_X=[]

    for i in 1:length(L_array)
        j=argmax(Specific_Heat_matrix[:,i])
        push!(Tc_L_C, temperature_array_reverse[j])
    end
    for i in 1:length(L_array)
        j=argmax(Susceptibility_matrix[:,i])
        push!(Tc_L_X, temperature_array_reverse[j])
    end

    Tc_L_C=convert(Vector{Float64}, Tc_L_C)
    Tc_L_X=convert(Vector{Float64}, Tc_L_X)

    # Linear regression to obtain the value of Tc because Tc(L)=Tc+C*L^(-1)

    model_C = lm(@formula(y ~ x), DataFrame(x=Inverse_L , y=Tc_L_C))
    model_X = lm(@formula(y ~ x), DataFrame(x=Inverse_L , y=Tc_L_X))
   
    fit_params_C = coef(model_C)
    fit_params_X = coef(model_X)
    slope_C, intercept_C = fit_params_C[2], fit_params_C[1]
    slope_X, intercept_X = fit_params_X[2], fit_params_X[1]

    # Obtain the errors for C
    errors_C = stderror(model_C)
    errors_X = stderror(model_X)
    slope_error_C, intercept_error_C = errors_C[2], errors_C[1]
    slope_error_X, intercept_error_X = errors_X[2], errors_X[1]



    r_squared_C=r2(model_C)
    r_squared_X=r2(model_X)

    println("REGRESSION RESULTS")
    println(" ")
    println("Tc_C= ", intercept_C, " ± ", intercept_error_C )
    println("Tc_X= ", intercept_X, " ± ", intercept_error_X )
    println(" ")
    println("r2_C= " , r_squared_C)
    println("r2_X= " , r_squared_X)
    println(" ")

    y_C=slope_C .*Inverse_L .+ intercept_C
    y_X=slope_X .*Inverse_L .+ intercept_X

    plot(Inverse_L, y_C, label="Tc(C)= $intercept_C ± $intercept_error_C")
    plot(Inverse_L, y_X, label="Tc(χ)= $intercept_X ± $intercept_error_X")
    scatter!(Inverse_L, Tc_L_C, xlabel="1/L", ylabel=L"Tc", grid=false, label=L"Tc(C_V)")
    scatter!(Inverse_L, Tc_L_X, xlabel="1/L", ylabel=L"Tc", grid=false, label=L"Tc(\chi)")
    savefig("Tc-Inverse_L.png")
    return(intercept_C, intercept_X)
end

Tc_2_C,Tc_2_X=Critical_temperature(Specific_Heat_matrix, Susceptibility_matrix, L_array, temperature_array_reverse)



# PLOTING THE DIFFERENT QUANTITIES ##########################################################################################


# Fourth order cumulant 
#plot(markershape = :circle, ms=3)
#for i in 1:length(L_array)
#   L=L_array[i]
#   plot!(temperature_array_reverse, Fourth_order_cumulant_matrix[:, i], label="L=$L")
#end
#savefig("Fourth_order_cumulant.png")

Tc=2.269


# Energy vs T

x_H=collect(range(1, stop=4, length=1000))
mean_H= -coth.(2 ./x_H)
plot(xlabel="T", ylabel="H", grid=false, label="Theoretical result")
vline!([Tc], label="Tc", color=:black, linestyle=:dash)
ylims!(-2.3, -.05)
for i in 1:length(L_array)
    L=L_array[i]
    scatter!(temperature_array_reverse, Energy_matrix[:, i], label="L=$L", ms=3)
end
savefig("H-T.png")

# Magnetization vs T

x_m=collect(range(1, stop=2.269, length=1000))
mean_m=(1 .-(sinh.(2 ./x_m)).^(-4)).^(1/8)
plot(x_m, mean_m, xlabel="T", ylabel="m", grid=false)
vline!([Tc], label="Tc", color=:black, linestyle=:solid)
ylims!(0, 1.1)
for i in 1:length(L_array)
    L=L_array[i]
    scatter!(temperature_array_reverse, Magnetization_matrix[:, i], label="L=$L", ms=3)
end
savefig("m-T.png")

# Susceptibility vs T

plot( xlabel="T", ylabel=L"\chi_T", grid=false)
vline!([Tc], label="Tc", color=:black, linestyle=:solid)
ylims!(0, 40)
for i in 1:length(L_array)
    L=L_array[i]
    scatter!(temperature_array_reverse, Susceptibility_matrix[:, i], label="L=$L", ms=3)
end
savefig("X-T.png")

# Specific Heat vs T 

x_C_1=collect(range(1.5, stop=2.25, length=1000))
x_C_2=collect(range(2.28, stop=3, length=1000))
x_C=append!(x_C_1,x_C_2)
mean_C= (2/pi).*(2 ./x_C).^2 .*(-1 .*log.(abs.(1 .- x_C ./Tc)) .+ log(Tc/2) .- (1+pi/4))
plot(xlabel="T", ylabel=L"C_V", grid=false)
vline!([Tc], label="Tc", color=:black, linestyle=:dash)
ylims!(0, 3)
for i in 1:length(L_array)
    L=L_array[i]
    scatter!(temperature_array_reverse, Specific_Heat_matrix[:, i], label="L=$L", ms=3)
end
savefig("C-T.png")





# CRITICAL EXPONENTS ###################################################################################################


function nu_value(Specific_Heat_matrix, L_array, Tc, temperature_array_reverse)

    Tc_L_C=[]
    
    for i in 1:length(L_array)
        j=argmax(Specific_Heat_matrix[:,i])
        push!(Tc_L_C, temperature_array_reverse[j])
    end

    Tc_L_C=convert(Vector{Float64}, Tc_L_C)
    
    model = lm(@formula(y ~ x), DataFrame(x=log.(L_array) , y=log.(abs.(Tc_L_C .-Tc))))

    fit_params = coef(model)
    slope, intercept = fit_params[2], fit_params[1]

    # Obtain the errors 
    errors = stderror(model)
    slope_error, intercept_error = errors[2], errors[1]

    r_squared=r2(model)

    println("REGRESSION RESULTS FOR NU")
    println(" ")
    println("nu= ", slope, " ± ", slope_error )
    println(" ")
    println("r2_nu= " , r_squared)
    println(" ")

    y=slope .*log.(L_array) .+ intercept

    plot(log.(L_array), y, label="ν= $slope ± $slope_error")
    scatter!(log.(L_array), log.(abs.(Tc_L_C .-Tc)), xlabel="ln(L)", ylabel="ln(Tc(L)-Tc)", grid=false, label="")
    savefig("nu.png")
    return(slope)
end

println("NU FROM SPECIFIC HEAT")
println(" ")
nu_C=nu_value(Specific_Heat_matrix, L_array, Tc_2_C, temperature_array_reverse)
println("NU FROM SUSCEPTIBILITY")
println(" ")
nu_X=nu_value(Susceptibility_matrix, L_array, Tc_2_X, temperature_array_reverse)

function beta_value(Magnetization_matrix, Susceptibility_matrix, L_array)

    m=[]

    for i in 1:length(L_array)
        j=argmax(Susceptibility_matrix[:,i])
        push!(m, Magnetization_matrix[j,i])
    end

    m=convert(Vector{Float64}, m)


    model = lm(@formula(y ~ x), DataFrame(x=log.(L_array) , y=log.(m)))

    fit_params = coef(model)
    slope, intercept = fit_params[2], fit_params[1]

    # Obtain the errors 
    errors = stderror(model)
    slope_error, intercept_error = errors[2], errors[1]

    r_squared=r2(model)

    println("REGRESSION RESULTS FOR BETA")
    println(" ")
    println("beta= ", slope, " ± ", slope_error )
    println(" ")
    println("r2_beta= " , r_squared)
    println(" ")

    y=slope .*log.(L_array) .+ intercept

    plot(log.(L_array), y, label="β= $slope ± $slope_error")
    scatter!(log.(L_array), log.(m), xlabel="ln(L)", ylabel="ln(m)", grid=false, label="")
    savefig("beta.png")


    

end

beta_value(Magnetization_matrix, Susceptibility_matrix, L_array)

function alpha_value(Specific_Heat_matrix, L_array)

    C=[]

    for i in 1:length(L_array)
        j=argmax(Specific_Heat_matrix[:,i])
        push!(C, Specific_Heat_matrix[j,i])
    end

    C=convert(Vector{Float64}, C)


    model = lm(@formula(y ~ x), DataFrame(x=log.(L_array) , y=log.(C)))

    fit_params = coef(model)
    slope, intercept = fit_params[2], fit_params[1]

    # Obtain the errors 
    errors = stderror(model)
    slope_error, intercept_error = errors[2], errors[1]

    r_squared=r2(model)

    println("REGRESSION RESULTS FOR ALPHA")
    println(" ")
    println("alpha= ", slope, " ± ", slope_error )
    println(" ")
    println("r2_alpha= " , r_squared)
    println(" ")

    y=slope .*log.(L_array) .+ intercept

    plot(log.(L_array), y, label="α= $slope ± $slope_error")
    scatter!(log.(L_array), log.(C), xlabel="ln(L)", ylabel=L"ln(C_V)", grid=false, label="")
    savefig("alpha.png")
end

alpha_value(Specific_Heat_matrix, L_array)

function gamma_value(Susceptibility_matrix, L_array)

    X=[]

    for i in 1:length(L_array)
        j=argmax(Susceptibility_matrix[:,i])
        push!(X, Susceptibility_matrix[j,i])
    end

    X=convert(Vector{Float64}, X)


    model = lm(@formula(y ~ x), DataFrame(x=log.(L_array) , y=log.(X)))

    fit_params = coef(model)
    slope, intercept = fit_params[2], fit_params[1]

    # Obtain the errors 
    errors = stderror(model)
    slope_error, intercept_error = errors[2], errors[1]

    r_squared=r2(model)

    println("REGRESSION RESULTS FOR GAMMA")
    println(" ")
    println("gamma= ", slope, " ± ", slope_error )
    println(" ")
    println("r2_alpha= " , r_squared)
    println(" ")

    y=slope .*log.(L_array) .+ intercept

    plot(log.(L_array), y, label="γ= $slope ± $slope_error")
    scatter!(log.(L_array), log.(X), xlabel="ln(L)", ylabel=L"ln(χ)", grid=false, label="")
    savefig("gamma.png")
    return(slope)
end

gamma=gamma_value(Susceptibility_matrix, L_array)


# FINITE SIZE SCALING ################################################################################################

function Scaling_Susceptibility(Susceptibility_matrix, gamma, L_array, temperature, nu_X, Tc)
    New_Susceptibility_matrix=zeros(length(temperature),0)
    for i in 1:length(L_array)
        XL=Susceptibility_matrix[:,i].*L_array[i].^(-gamma/nu_X)
        New_Susceptibility_matrix=hcat(New_Susceptibility_matrix, XL)
    end 

    plot( xlabel="T", ylabel=L"\chi_T", grid=false)
    vline!([Tc], label="Tc", color=:black, linestyle=:solid)
    for i in 1:length(L_array)
        L=L_array[i]
        scatter!(temperature, New_Susceptibility_matrix[:, i], label="L=$L", ms=3)
    end
    savefig("Finite_Size_Scaling.png")
    
end

Scaling_Susceptibility(Susceptibility_matrix, gamma, L_array, temperature_array_reverse, nu_X, Tc_2_X)


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
        return(right, left, up, down)
    end
    
end

e=other_neighbours(49, 8, 3)



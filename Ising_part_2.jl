using DelimitedFiles, Plots, Random, Statistics, LaTeXStrings, GLM, DataFrames


L_array=[8,16,32,64,128]

Energy_matrix=readdlm("Energy_matrix.txt")

Magnetization_matrix=readdlm("Magnetization_matrix.txt")

Susceptibility_matrix=readdlm("Susceptibility_matrix.txt")

Specific_Heat_matrix=readdlm("Specific_Heat_matrix.txt")

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


# CRITICAL TEMPERATURE ########################################################################################################

function Critical_temperature(Specific_Heat_matrix, Susceptibility_matrix, L_array, temperature)
    Inverse_L=1 ./ L_array
    Tc_L_C=[]
    Tc_L_X=[]

    for i in 1:length(L_array)
        j=argmax(Specific_Heat_matrix[:,i])
        push!(Tc_L_C, temperature[j])
    end
    for i in 1:length(L_array)
        j=argmax(Susceptibility_matrix[:,i])
        push!(Tc_L_X, temperature[j])
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

    plot(Inverse_L, y_C, label="Tc(C)= 2.2716 ± 0.0043", color=:red)
    plot!(Inverse_L, y_X, label="Tc(χ)= 2.2463 ± 0.0076", color=:blue)
    scatter!(Inverse_L, Tc_L_C, xlabel="1/L", ylabel=L"Tc", grid=false, label="Tc(C)", color=:red)
    scatter!(Inverse_L, Tc_L_X, xlabel="1/L", ylabel=L"Tc", grid=false, label=L"Tc(\chi)", color=:blue)
    savefig("Tc-Inverse_L.png")
    return(intercept_C, intercept_X)
end

Tc_2_C,Tc_2_X=Critical_temperature(Specific_Heat_matrix, Susceptibility_matrix, L_array, temperature_array_reverse)

# PLOTING THE DIFFERENT QUANTITIES ##########################################################################################


Tc=2.269


# Energy vs T

x_H=collect(range(1, stop=4, length=1000))
mean_H= -coth.(2 ./x_H)
plot(xlabel="T", ylabel=L"\left\langle H \right\rangle", grid=false, label="Theoretical result")
vline!([Tc], label="Tc", color=:black, linestyle=:dash)
ylims!(-1.1, -.06)
for i in 1:length(L_array)
    L=L_array[i]
    scatter!(temperature_array_reverse, Energy_matrix[:, i], label="L=$L", ms=3)
end
savefig("H-T.png")

# Magnetization vs T

x_m=collect(range(0, stop=2.269, length=1000))
mean_m=(1 .-(sinh.(2 ./x_m)).^(-4)).^(1/8)
plot(x_m, mean_m, xlabel="T", ylabel=L"\left\langle M \right \rangle", grid=false, label= "Onsager", color=:black, linestyle=:solid)
vline!([Tc], label="Tc", color=:black, linestyle=:dash)
ylims!(0, 1.1)
for i in 1:length(L_array)
    L=L_array[i]
    scatter!(temperature_array_reverse, Magnetization_matrix[:, i], label="L=$L", ms=3)
end
savefig("m-T.png")

# Susceptibility vs T

plot( xlabel="T", ylabel=L"\chi", grid=false)
vline!([Tc], label="Tc", color=:black, linestyle=:dash)
ylims!(0, 43)
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
ylims!(0, 1)
for i in 1:length(L_array)
    L=L_array[i]
    scatter!(temperature_array_reverse, Specific_Heat_matrix[:, i], label="L=$L", ms=3)
end
savefig("C-T.png")


# CRITICAL EXPONENTS ###################################################################################################


function nu_value_Specific_Heat(Specific_Heat_matrix, L_array, Tc, temperature)

    Tc_L_C=[]
    
    for i in 1:length(L_array)
        j=argmax(Specific_Heat_matrix[:,i])
        push!(Tc_L_C, temperature[j])
    end

    Tc_L_C=convert(Vector{Float64}, Tc_L_C)
    
    model = lm(@formula(y ~ x), DataFrame(x=log.(L_array) , y=log.(abs.(Tc_L_C .-Tc))))

    fit_params = coef(model)
    slope, intercept = fit_params[2], fit_params[1]

    # Obtain the errors 
    errors = stderror(model)
    slope_error, intercept_error = errors[2], errors[1]

    r_squared=r2(model)

    println("nu_C= ", 1/slope, " ± ", slope_error/slope^2 )
    println(" ")
    println("r2_nu= " , r_squared)
    println(" ")

    return(abs(slope))
end

println("NU FROM SPECIFIC HEAT")
println(" ")
nu_C=nu_value_Specific_Heat(Specific_Heat_matrix, L_array, Tc_2_C, temperature_array_reverse)

function nu_value_Suscepitibility(Specific_Heat_matrix, L_array, Tc, temperature)

    Tc_L_C=[]
    
    for i in 1:length(L_array)
        j=argmax(Specific_Heat_matrix[:,i])
        push!(Tc_L_C, temperature[j])
    end

    Tc_L_C=convert(Vector{Float64}, Tc_L_C)
    
    model = lm(@formula(y ~ x), DataFrame(x=log.(L_array) , y=log.(abs.(Tc_L_C .-Tc))))

    fit_params = coef(model)
    slope, intercept = fit_params[2], fit_params[1]

    # Obtain the errors 
    errors = stderror(model)
    slope_error, intercept_error = errors[2], errors[1]

    r_squared=r2(model)

    println("nu_X= ", 1/slope, " ± ", slope_error/slope^2 )
    println(" ")
    println("r2_nu= " , r_squared)
    println(" ")

    return(abs(slope))
end
println("NU FROM SUSCEPTIBILITY")
println(" ")
nu_X=nu_value_Suscepitibility(Susceptibility_matrix, L_array, Tc_2_X, temperature_array_reverse)

function nu_value(Specific_Heat_matrix,Susceptibility_matrix, L_array, Tc_C, Tc_X, temperature)

    Tc_L_C=[]
    Tc_L_X=[]
    
    for i in 1:length(L_array)
        j=argmax(Specific_Heat_matrix[:,i])
        k=argmax(Susceptibility_matrix[:,i])
        push!(Tc_L_C, temperature[j])
        push!(Tc_L_X, temperature[k])
    end

    Tc_L_C=convert(Vector{Float64}, Tc_L_C)
    Tc_L_X=convert(Vector{Float64}, Tc_L_X)
    
    model_C = lm(@formula(y ~ x), DataFrame(x=log.(L_array) , y=log.(abs.(Tc_L_C .-Tc_C))))
    model_X = lm(@formula(y ~ x), DataFrame(x=log.(L_array) , y=log.(abs.(Tc_L_X .-Tc_X))))

    fit_params_C = coef(model_C)
    slope_C, intercept_C = fit_params_C[2], fit_params_C[1]
    fit_params_X = coef(model_X)
    slope_X, intercept_X = fit_params_X[2], fit_params_X[1]


    y_C=slope_C .*log.(L_array) .+ intercept_C
    y_X=slope_X .*log.(L_array) .+ intercept_X

    plot(log.(L_array), y_C, label="ν(C)= 0.82 ± 0.20", color=:red)
    plot!(log.(L_array), y_X, label="ν(χ)= 1.16 ± 0.16", color=:blue)

    scatter!(log.(L_array), log.(abs.(Tc_L_C .-Tc_C)), xlabel="ln(L)", ylabel="ln(Tc(L)-Tc)", grid=false, label="Specific heat", color=:red)
    scatter!(log.(L_array), log.(abs.(Tc_L_X .-Tc_X)), xlabel="ln(L)", ylabel="ln(Tc(L)-Tc)", grid=false, label="Susceptibility", color=:blue)
    savefig("nu.png")
end

nu_value(Specific_Heat_matrix,Susceptibility_matrix, L_array, Tc_2_C, Tc_2_X, temperature_array_reverse)

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

    plot(log.(L_array), y, label="β= -0.115 ± 0.024")
    scatter!(log.(L_array), log.(m), xlabel="ln(L)", ylabel="ln(m)", grid=false, label="")
    savefig("beta.png")
    return(abs(slope))

end

beta=beta_value(Magnetization_matrix, Susceptibility_matrix, L_array)

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

    plot(log.(L_array), y, label="α= 0.288 ± 0.025")
    scatter!(log.(L_array), log.(C), xlabel="ln(L)", ylabel=L"ln(C_V)", grid=false, label="")
    savefig("alpha.png")
    return(abs(slope))
end

alpha=alpha_value(Specific_Heat_matrix, L_array)

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
    println("r2_gamma= " , r_squared)
    println(" ")

    y=slope .*log.(L_array) .+ intercept

    plot(log.(L_array), y, label="γ= 1.715 ± 0.055")
    scatter!(log.(L_array), log.(X), xlabel="ln(L)", ylabel=L"ln(χ)", grid=false, label="")
    savefig("gamma.png")
    return(abs(slope))
end

gamma=gamma_value(Susceptibility_matrix, L_array)


# FINITE SIZE SCALING ################################################################################################

function Scaling_Susceptibility(Susceptibility_matrix, gamma, L_array, temperature, nu_X, Tc)

    New_Susceptibility_matrix=zeros(length(temperature),0)

    for i in 1:5
        XL=Susceptibility_matrix[:,i].*L_array[i]^(-gamma/nu_X)
        New_Susceptibility_matrix=hcat(New_Susceptibility_matrix, XL)
    end 

    New_temperature=zeros(length(temperature),0)

    for i in 1:5
        T=(temperature ./Tc.-1) .* L_array[i] .^ (1/nu_X)
        New_temperature=hcat(New_temperature, T)
    end 

    plot( xlabel=L"tL^{1/\nu}", ylabel=L"\chi L^{-\gamma/\nu}", grid=false)
    vline!([0], label="Tc", color=:black, linestyle=:dash)
    ylims!(-0.005, 0.05)

    for i in 1:5
        L=L_array[i]
        scatter!(New_temperature[:,i], New_Susceptibility_matrix[:, i], label="L=$L", ms=3)
    end
    savefig("Finite_Size_Scaling_Susceptibility.png")
    
end

Scaling_Susceptibility(Susceptibility_matrix, gamma, L_array, temperature_array_reverse, nu_X, Tc_2_X)

function Scaling_Specific_Heat(Specific_Heat_matrix, alpha, L_array, temperature, nu_C, Tc)

    New_Specific_Heat_matrix=zeros(length(temperature),0)

    for i in 1:5
        XL=Specific_Heat_matrix[:,i].*L_array[i]^(-alpha/nu_C)
        New_Specific_Heat_matrix=hcat(New_Specific_Heat_matrix, XL)
    end 

    New_temperature=zeros(length(temperature),0)

    for i in 1:5
        T=(temperature ./Tc.-1) .* L_array[i] .^ (1/nu_C)
        New_temperature=hcat(New_temperature, T)
    end 


    plot( xlabel=L"tL^{1/\nu}", ylabel=L"C_VL^{-\alpha/\nu}", grid=false)
    vline!([0], label="Tc", color=:black, linestyle=:dash)
    ylims!(-0.005, 0.25)

    for i in 1:5
        L=L_array[i]
        scatter!(New_temperature[:,i], New_Specific_Heat_matrix[:, i], label="L=$L", ms=3)
    end
    savefig("Finite_Size_Scaling_Specific_Heat.png")
    
end

Scaling_Specific_Heat(Specific_Heat_matrix, alpha, L_array, temperature_array_reverse, nu_C, Tc_2_C)

function Scaling_Magnetization(Magnetization_matrix, beta, L_array, temperature, nu_X, Tc)

    New_Magnetization_matrix=zeros(length(temperature),0)

    for i in 1:5
        XL=Magnetization_matrix[:,i].*L_array[i]^(beta/nu_X)
        New_Magnetization_matrix=hcat(New_Magnetization_matrix, XL)
    end 

    New_temperature=zeros(length(temperature),0)

    for i in 1:5
        T=(temperature ./Tc .-1) .* L_array[i] .^ (1/nu_X)
        New_temperature=hcat(New_temperature, T)
    end 


    plot( xlabel=L"tL^{1/\nu}", ylabel=L"\left\langle M \right\rangle L^{\beta/\nu}", grid=false)
    vline!([0], label="Tc", color=:black, linestyle=:dash)

    for i in 1:5
        L=L_array[i]
        scatter!(New_temperature[:,i], New_Magnetization_matrix[:, i], label="L=$L", ms=3)
    end
    savefig("Finite_Size_Scaling_Magnetization.png")
    
end

Scaling_Magnetization(Magnetization_matrix, beta, L_array, temperature_array_reverse, nu_X, Tc_2_X)



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
using Plots
gr(size=(600,600))
default(fmt = :png)

function rungekutta4(tf, y0, f, t0; n=100)
    h = (tf - t0) / n
    t = linspace(t0, tf, n+1)
    y = zeros(4, n + 1)
    y[:,1] = y0
    for i = 1:n
        k1 = f(t[i], y[:,i])
        k2 = f(t[i] + 0.5 * h, y[:,i] + 0.5 * h * k1)
        k3 = f(t[i] + 0.5 * h, y[:,i] + 0.5 * h * k2)
        k4 = f(t[i] + h, y[:,i] + h * k3)
        Y[:,i+1] = y[:,i] + ((1/6) * h) * (k1 + 2 * k2 + 2 * k3 + k4)
    end
    return Y
end

F(t, Y) = [Y[3] ; Y[4];
        (-g*(2*m₁+m₂)*sin(θ₁))-(m₂*g*sin(θ₁-2*θ₂))-(2*sin(θ₁-θ₂)*m₂*(Y[4]^2)*l₂+(Y[3]^2)*l₁*cos(θ₁-θ₂)) / (l₁*(2*m₁+m₂-m₂*cos(2*θ₁-2*θ₂))) ;
        2*sin(θ₁-θ₂)*((Y[3]^2)*l₁*(m₁+m₂)+g*(m₁+m₂)*cos(θ₁)+((Y[4]^2)*l₂*m₂*cos(2*θ₁-2*θ₂))) / (l₂*(2*m₁+m₂-m₂*cos(2*θ₁-2*θ₂)))]

#θ₁' = ω₁
#θ₂' = ω₂
#ω₁' = .......ω₂^2 e ω₁^2
#ω₂' = .......ω₁^2 e ω₂^2
m₁ = m₂ = 1.0
g = 9.81
l₁ = l₂ = 3.0
θ₁ = 0.75 * π
θ₂ = π
ω₁ = ω₂ = 0.0
ω₀ = 0.0
θ₁ = π/2
θ₂ = π/5
θ₀ = 0.75 * π
t0, tf = 0.0, 60.0

n = 10000
Y0 = [θ₀; θ₀; ω₀; ω₀]
Y = zeros(4, 10001)
Y = rungekutta4(tf, Y0, F, t0, n=n)
θ₁ = Y[1,:](?)
θ₂ = Y[2,:](?)
ω₁ = Y[3,:](?)
ω₂ = Y[4,:](?)
#plot
plot!(linspace(t0, tf, n+1), (?))

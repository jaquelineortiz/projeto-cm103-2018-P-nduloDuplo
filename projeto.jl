using Plots
gr(size=(400,300))
default(fmt = :png)

function rungekutta4(tf, Y0, f, t0; n=100)
    h = (tf - t0) / n
    t = linspace(t0, tf, n+1)
    Y = zeros(length(Y0), n + 1)
    Y[:,1] = Y0
    for i = 1:n
        k1 = f(t[i], Y[:,i])
        k2 = f(t[i] + 0.5 * h, Y[:,i] + 0.5 * h * k1)
        k3 = f(t[i] + 0.5 * h, Y[:,i] + 0.5 * h * k2)
        k4 = f(t[i] + h, Y[:,i] + h * k3)
        Y[:,i+1] = Y[:,i] + ((1/6) * h) * (k1 + 2 * k2 + 2 * k3 + k4)
    end
    return Y
end


m₁ = m₂ = 1.0
g = 1.62
l₁ = l₂ = 3.0
F(t, Y) = [Y[3];
           Y[4];
           (-m₂ * l₁ * Y[4]^2 * sin(Y[1] - Y[2]) * cos(Y[1] - Y[2]) + g * m₂ * sin(Y[2]) * cos(Y[1] - Y[2]) - m₂ * l₂ * Y[4]^2 * sin(Y[1] - Y[2]) - (m₁ + m₂) * g * sin(Y[1])) / (l₁ * (m₁ + m₂) - m₂ * l₁ * cos(Y[1] - Y[2])^2);
           (m₂ * l₂ * Y[4]^2 * sin(Y[1] - Y[2]) * cos(Y[1] - Y[2]) + g * sin(Y[1]) * cos(Y[1] - Y[2]) * (m₁ + m₂) + l₁ * Y[4]^2 * sin(Y[1] - Y[2]) * (m₁ + m₂) - g * sin(Y[2]) * (m₁ + m₂)) / (l₂ * (m₁ + m₂) - m₂ * l₂ * cos(Y[1] - Y[2])^2)]

#θ₁' = ω₁
#θ₂' = ω₂
#ω₁' = .......ω₂^2 e ω₁^2
#ω₂' = .......ω₁^2 e ω₂^2


ω₀ = 0.0
#θ₀ = 0.5 * π
t0, tf = 0.0, 30.0

n = 4000
Y0 = [pi/2; pi/2; θ₀; ω₀; ω₀]
Y = zeros(4, n+1)
Y = rungekutta4(tf, Y0, F, t0, n=n)
θ₁ = Y[1,:]
θ₂ = Y[2,:]
ω₁ = Y[3,:]
ω₂ = Y[4,:]
#plot(linspace(t0, tf, n+1),)

gr(size=(400,300))
anim = Animation()

x₁ = l₁* sin.(θ₁)
y₁ = -l₁ * cos.(θ₁)
x₂ = x₁ + l₂ * sin.(θ₂)
y₂ = y₁ - l₂ * cos.(θ₂)
l = l₁ + l₂

skip = div(n,200)
for i = 1:skip:n
    scatter([x₁[i], x₂[i]], [y₁[i], y₂[i]], m=(stroke(1,:black),:black,10))
    plot!([0; x₁[i]; x₂[i]], [0; y₁[i]; y₂[i]], c=:black)
    xlims!(-1.2 * l, 1.2 * l)
    ylims!(-1.2 * l, 1.2 * l)

    frame(anim)
end
gif(anim, "pendulo.gif",fps=24)

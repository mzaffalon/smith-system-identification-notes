using LinearAlgebra: rank
using ToeplitzMatrices: Toeplitz
using ControlSystemsBase: tf, lsim

Ts, T = 1, 100
τmax = 5
G = tf([1], [1, 0.2], Ts)
G = tf([0.8, 2.28, -6.8, 6.03, 25.5], [1, 0, 0, 0, 0, 0], Ts)
u(t) = sin(2*pi*0.1*t)
r = lsim(G, (x,t)->u(t), T)

Φu = Toeplitz(u.(0:Ts:T), zeros(τmax+1))
@assert rank(Φu) == 6 # It is a step function

# Bug in Toeplitz
gest = Float64.(Φu) \ vec(r.y)

Φu_per = Φu[10:end,:]
@assert rank(Φu_per) == 2
y_per = vec(r.y[:,10:end])
gest_per = Float64.(Φu_per) \ y_per

gest = Float64.(Φu[3:end,:]) \ vec(r.y[:,3:end])

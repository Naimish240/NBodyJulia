using Random
using Plots
using Printf

function getAcc(pos, mass, G)
	# positions r = [x,y] for all particles
	x = pos[:,1]
	y = pos[:,2]

	# all pairwise particle separations
	Δx = x' .- x
	Δy = y' .- y

	# 1/r^3 for all particle pairwise particle separations
	inv_r3 = (Δx.^2 + Δy.^2)
	inv_r3[inv_r3.>0] = inv_r3[inv_r3.>0].^(-1.5)

	ax = G .* (Δx .* inv_r3) * mass
	ay = G .* (Δy .* inv_r3) * mass

	# return acceleration components
	a = [ax ay]
	return a;
end


function main(N, tEnd, dt, G, show=false, save=false)
	# Initialize time
	t = 0

	# Generate Initial Conditions
	seed = MersenneTwister(42)
	mass = randn(seed, N, 1)
	pos  = randn(seed, N, 2)
	vel  = randn(seed, N, 2)

	# calculate initial gravitational accelerations
	acc = getAcc(pos, mass, G);

	# number of timesteps
	Nt = convert(Int64,ceil(tEnd/dt));

	# Simulation Loop
	for i in 1:Nt
		vel += acc * dt/2
		pos += vel * dt
		acc = getAcc(pos, mass, G)
		vel += acc * dt/2
		t += dt;

		x = pos[:, 1]
		y = pos[:, 2]

		if show == true
			display(scatter(x, y, xlims=(-5, 5), ylims=(-5,5)))
		end
		if save == true
			savefig(scatter(x, y, xlims=(-5, 5), ylims=(-5,5)), "frames/$i.png")
		end
	end
end

main(25, 1, 0.01, 1, true, true)

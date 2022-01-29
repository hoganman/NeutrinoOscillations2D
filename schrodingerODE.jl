using LinearAlgebra
using Plots

function get_mass2_mat(nu1_mass::Float16, nu2_mass::Float16)::Array
    """Create a 2D mass-squared diagonal matrix"""
    ret_array::Array{Float16} = Diagonal([nu1_mass^2, nu2_mass^2])
    return ret_array
end

function get_mixing_matrix(mixing_angle::Float16, isdeg::Bool = true)::Array
    """Create a 2D mixing/rotation matrix with a given angle"""
    if isdeg
        mixing_angle = deg2rad(mixing_angle)
    end
    mixing_mat::Array{Float16} = [
        +cos(mixing_angle) -sin(mixing_angle)
        +sin(mixing_angle) +cos(mixing_angle)
    ]
    return mixing_mat
end

function time_propagator_mat(
        energy_nu::Float16,
        time_step::Float16,
        delta_mat::Array{Float16}
)::Array{ComplexF16}
    """Calculate the time propagator. Due to a bug in Julia, complex 
    matrix arguments for the exponential function are not supported. 
    So instead, I use the Euler relation instead.
    """
    exp_arg::Array{ComplexF16} = time_step / (2.0 * energy_nu) * delta_mat
    exp_mat::Array{ComplexF16} = cos(exp_arg) + (0 + 1im) * sin(exp_arg)
    return exp_mat
end


function schrodingerODE()
    """Run a 2D toy example of neutrino oscillations"""

    # Define scalars
    nu1_mass::Float16 = 1.0
    nu2_mass::Float16 = 2.0
    mixing_angle_deg::Float16 = 40
    time_steps::Int = 100
    time_step_size::Float16 = 0.01
    energy_nu::Float16 = 10.0

    # Define physics matrices
    mass2_mat::Array{Float16} = get_mass2_mat(nu1_mass, nu2_mass)
    mixing_mat::Array{Float16} = get_mixing_matrix(mixing_angle_deg, true)
    delta_mat::Array{Float16} = transpose(conj(mixing_mat)) * mass2_mat * mixing_mat

    # Store results in arrays
    # initial state is nu1 = 100%, nu2 = 0%
    initial_state_vector::Vector{ComplexF16} = Vector([1 + 0im, 0 + 0im])
    state_vectors::Array{ComplexF16} = zeros(ComplexF16, (2, time_steps))
    state_fractions::Array{Float16} = zeros(Float16, (2, time_steps))
    state_vectors[:, 1] = initial_state_vector
    state_fractions[:, 1] = Array([1.0; 0.0])

    # Evaluate the wave-function over time and save the state over each step
    for step_index::Int = 2:time_steps
        previous_step_index::Int = step_index - 1
        previous_state_vector::Array{ComplexF16} = state_vectors[:, previous_step_index]
        time_step::Float16 = (step_index - 1) * time_step_size
        time_propagator::Array{ComplexF16} =
            time_propagator_mat(energy_nu, time_step, delta_mat)
        new_state::Vector{ComplexF16} = time_propagator * previous_state_vector
        new_state_fraction::Vector{Float16} =
            Array([abs2(new_state[1]); abs2(new_state[2])])
        state_vectors[:, step_index] = new_state
        state_fractions[:, step_index] = new_state_fraction
    end

    # Plot the results
    Plots.pyplot()
    fig = plot(1:time_steps, state_fractions[1, :], label="\$\\nu_1\$")
    plot!(1:time_steps, state_fractions[2, :], label="\$\\nu_2\$")
    xlabel!("Time step")
    ylabel!("Wave Function Component")
    title!("Flavor Oscillations")
    @show fig
end

schrodingerODE()
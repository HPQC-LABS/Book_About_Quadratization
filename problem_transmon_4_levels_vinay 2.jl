using LinearAlgebra
using QuantumAnnealingTools, OrdinaryDiffEq


const PHI_0 = 2.067833831e-15
const H_PLANCK = 6.62607015e-34;
const E_CHARGE = 1.602176634e-19;


################# Fucntions #############################
function transmon_H(e_j, e_c, n_max, ϕx, d)
    n2 = Diagonal([(i)^2 for i = -n_max:n_max])
    tot_n = 2 * n_max + 1
    ϕ0 = atan(d * tan(ϕx / 2))

    d1 = Bidiagonal(zeros(Int64, tot_n), [1 for i = 1:tot_n-1], :L)
    cosphi = 1 / 2 * (exp(1im * ϕ0) * d1 + exp(-1im * ϕ0) * Transpose(d1))

    H = (
        4 * e_c * n2 -
        e_j * cos(ϕx / 2) * sqrt(1 + d^2 * tan(ϕx / 2)^2) * cosphi
    )     #### here I am taking ng = 0 i.e. we are always at the sweet spot
    return H  ## in GHz
end


function A_current(e_j_eff, d, phi_x, n_max)
    tot_n = 2 * n_max + 1
    d1 = Bidiagonal(zeros(Int64, tot_n), [1 for i = 1:tot_n-1], :L)
    cosphi = 1 / 2 * (d1 + Transpose(d1))
    sinphi = 1 / 2im * (d1 - Transpose(d1))


    A =
        H_PLANCK *
        1e9 *
        (
            e_j_eff * pi / PHI_0 *
            (sin(phi_x / 2) * cosphi - d * cos(phi_x / 2) * sinphi)
        )
    return A
end


function make_drive_H(n_max, ωq, η)
    ## This function returns the H_drive matrix already in the diagonal basis of H_sys
    tot_n = 2 * n_max + 1
    a_dagger = Bidiagonal(
        zeros(Float64, tot_n),
        [sqrt(i) * (1 - (i - 1) / 2 * η / ωq) for i = 1:tot_n-1],
        :L,
    )
    a = Transpose(a_dagger)
    H_x = (a + a_dagger) ./ 2
    H_y = 1im * (-a + a_dagger) ./ 2
    N_hat = Diagonal([i - 1 for i = 1:tot_n])
    return H_x, H_y, N_hat#,a_dagger*a
end
########################################################

function build_hamiltonian_and_coupling(
    E_J,
    E_C,
    n,
    phi_x_value,
    d_value;
    n_trunc = 4,
)
    H_trans = transmon_H(E_J, E_C, n, phi_x_value, d_value)
    e_list = eigvals(H_trans)
    e_list_new = e_list - e_list[1] * ones(2 * n + 1)
    ω10 = e_list_new[2] - e_list_new[1]
    ω21 = e_list_new[3] - e_list_new[2]
    η_1 = -ω21 + ω10
    A = A_current(E_J, d_value, phi_x_value, n)

    U_eigen = eigvecs(H_trans)
    H_trans_diag = U_eigen' * H_trans * U_eigen
    A_eigen = U_eigen' * A * U_eigen
    H_drive_x, H_drive_y, num_op = make_drive_H(n, ω10, η_1)

    H_trans_trunc = H_trans_diag[1:n_trunc, 1:n_trunc]
    A_eigen_trunc = A_eigen[1:n_trunc, 1:n_trunc]
    Ay = A_eigen_trunc .- Diagonal(diag(A_eigen_trunc))
    Az = zeros(4, 4) .+ Diagonal(diag(A_eigen_trunc))

    H_trans_trunc =
        H_trans_trunc - H_trans_trunc[1, 1] * Diagonal([1 for i = 1:n_trunc])

    H_trans_trunc, A_eigen_trunc, Ay, Az
end

## setting the parameters
### We want for qubit 2, wq = 5.23817555736479 GHz  and Anharmoicity = 331.9 MHz
E_J = 13.4## in GHz
E_C = 0.287  ## in GHz
n = 20      ## number of levels
d_value = 0.03
phi_x_value = 0.04 * pi

const num_op_trunc = Diagonal([i - 1 for i = 1:4])
const Hm, Am, Ay, Az =
    build_hamiltonian_and_coupling(E_J, E_C, n, phi_x_value, d_value)
const ωd = real(Hm[2, 2] - Hm[1, 1] + 0.004)
const H_rotated = Hm - ωd * num_op_trunc

build_unitary(tf) = (s) -> exp(-2im * π * s * tf * Diagonal(H_rotated))

calc_unitary(s, tf) = Diagonal(exp.(-1im * ωd * s * tf * 2 * π * (0:3)))

function build_rotate_coupling(tf, mat)
    coupling = function (s)
        U = calc_unitary(s, tf)
        2 * π * g * PHI_0 / H_PLANCK * 1e-9 * U' * mat * U
    end
    CustomCouplings([coupling])
end

const v0 = [1.0, 0im, 0, 0]
const v1 = [0im, 1.0, 0, 0]
const vp = (v0 + v1) / sqrt(2)
const v2 = [0.0, 0, 1.0, 0]
const v3 = [0.0, 0, 1.0, 0]
const vm = (v0 - v1) / sqrt(2)
const vip = (v0 + 1im * v1) / sqrt(2)
const vim = (v0 - 1im * v1) / sqrt(2)

## parameters for both the  bath
const η = 1e-4
const T = 20
const g = 0.03

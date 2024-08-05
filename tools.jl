using BitBasis
# using Combinatorics
using Plots
using LinearAlgebra

### note: indices are backwards
# e.g. 0000000111 is nonzero at indices 1, 2, 3, equiv to int 7
# basis for single d atom: ml: 2, 1, 0, -1, -2, 2, 1, 0, -1, -2,
#                          ms: u, u, u,  u,  u, d, d, d,  d,  d
#                         idx: 10, 9, 8, 7,  6, 5, 4, 3,  2,  1

include("tables.jl")
include("params.jl")
k_table = allowed_k_table()
gaunt_table = gaunt_coeff_table()
gaunt_heads = gaunt_coeff_headings()
V_CEF = V_CFE_d_cubic();

"""
This function creates a single-shell Hilbert space for a system with n_orb orbitals and n_elec electrons.
Each state in the Hilbert space is an integer, which when unpacked in binary is in (ml, ms) basis as described below:
basis for single d atom: ml: 2, 1, 0, -1, -2, 2, 1, 0, -1, -2,
                         ms: u, u, u,  u,  u, d, d, d,  d,  d
                        idx: 10, 9, 8, 7,  6, 5, 4, 3,  2,  1
This function calls make_Hilbert_helper() which builds the Hilbert space recursively
params:
    n_elec: number of electrons in the shell (e.g. 4 for a d4 atom)
    n_orb: number of orbitals in the shell, counting spins separately (e.g. 10 for a d shell)
returns:
    h_space: a vector of ints, each representing one valid state in the Hilbert space
"""
function make_Hilbert_space(n_elec, n_orb)
    if n_orb < 16
        h_space = Vector{Int16}([])
    else
        h_space = Vector{Int32}([])
    end
    make_Hilbert_helper(n_elec, n_orb, h_space, 0)

    return h_space
end

"""
Helper function for make_Hilbert_space(). Uses recursive algorithm to create basis states
"""
function make_Hilbert_helper(n_elec, n_orb, h_space, state)
    # base case: all electrons added to state
    if n_elec == 0
        append!(h_space, state)
        return
    end

    if n_elec == n_orb
        for i in 1:n_orb
            state += 2^(i-1)
        end
        append!(h_space, state)
        return
    end
    
    # add electron to last orbital:
    state_tmp = state + 2^(n_orb-1)
    make_Hilbert_helper(n_elec-1, n_orb-1, h_space, state_tmp)

    # don't add electron to last orbital
    make_Hilbert_helper(n_elec, n_orb-1, h_space, state)
end

"""
Get the total values for Lz and Sz for a given state
params:
    state: an integer representing a basis state
    n_orb: number of orbitals in the state (e.g. 10 for a d orbital state)
returns:
    Lz: sum of ml for all component electrons in the state
    Sz: sum of ms for all component electrons
"""
function get_Lz_Sz(state, n_orb)
    l = (n_orb/2 - 1)/2
    Lz = 0
    Sz = 0.0
    elecs = baddrs(state)
    for e in elecs
        Lz += rem(e-1, n_orb/2) - l 
        Sz += fld(e-1, n_orb/2) - 0.5
    end
    return Lz, Sz
end

"""
This functino returns the Gaunt coefficient for electrons with quantum numbers (l1, m1) and (l2, m2). 
Values are stored in "tables.jl", but this function includes any sign changes as well
params:
    l1, m1, l2, m2: quantum numbers l and m_l for two electrons (ints)
returns:
    c_k: Gaunt coefficient (float)
"""
function get_ck(l1,m1,l2,m2)
    # Here we get ck(l1,m1,l2,m2), aray of gaunt coefficients (non-zero)
    if l1 == l2
        if abs(m1) == abs(m2) || abs(m1) > abs(m2)
            gaunt_12 = gaunt_table[l1][l2][m1][m2]./gaunt_heads[l1][l2]
        else#if abs(m1) < abs(m2)
            gaunt_12 = ((-1)^(m1-m2)*gaunt_table[l1][l2][m2][m1])./gaunt_heads[l1][l2]
        end

    elseif l1 < l2
        gaunt_12 = gaunt_table[l1][l2][m1][m2]./gaunt_heads[l1][l2]

    else#if l1 > l2
        gaunt_12 = ((-1)^(m1-m2)*gaunt_table[l2][l1][m2][m1])./gaunt_heads[l2][l1]
    end
        
    gaunt_12
end

"""
This function calculates the Coulomb interaction matrix element between two states 
params:
    state1, state2: basis states as in the format created by make_Hilbert_space()
    norb: number of orbitals represented by the basis states (d -> 10, p -> 6)
    F, G: Direct and Exchange radial integrals (Slater-Condon parameters), defined in "params.jl"
returns:    
    H_elem: matrix element of coloumb interaction Hamiltonian between the two states (float)
"""
function h_int_elem(state1, state2, norb, F, G)
    if get_Lz_Sz(state1, norb) != get_Lz_Sz(state2, norb)
        return 0
    end

    l::Int8 = (norb/2 - 1) / 2
    # k_arr = k_table[l, l]

    # Diagonal matrix element
    if state1 == state2
        H_elem = 0.0
        indices = baddrs(state1)
        for i in indices
            ml1::Int8 = rem(i-1, norb/2) - l
            ms1 = fld(i-1, norb/2) - 0.5
            for j in indices
                ml2::Int8 = rem(j-1, norb/2) - l
                ms2 = fld(j-1, norb/2) - 0.5
                H_elem += sum(get_ck(l, ml1, l, ml1) .* get_ck(l, ml2, l, ml2) .* F[1:l+1]) # Direct
                if ms1 == ms2
                    H_elem -= sum(get_ck(l, ml1, l, ml2) .* get_ck(l, ml1, l, ml2) .* G[1:l+1])    # Exchange
                end
            end
        end
        return H_elem / 2
    end

    # off-diagonal element
    H_elem = 0.0
    elecs1 = baddrs(state1)
    elecs2 = baddrs(state2)
    in_both = state1 & state2   # integer -- intersection of both states
    only_in_1 = baddrs(xor(state1, in_both))
    only_in_2 = baddrs(xor(state2, in_both))

    num_ex = 0
    for e in only_in_1
        num_ex += length(elecs2[elecs2 .< e])
    end
    for e in only_in_2
        num_ex += length(elecs1[elecs1 .< e])
    end
    sign = (-1)^num_ex
    # print(sign)

    

    if (length(only_in_1) != 2) && (length(only_in_2) != 2)
        return 0
    end

    m1a::Int8, s1a = rem(only_in_1[1] - 1, norb/2) - l, fld(only_in_1[1]-1, norb/2) - 0.5
    m1b::Int8, s1b = rem(only_in_1[2] - 1, norb/2) - l, fld(only_in_1[2]-1, norb/2) - 0.5
    m2a::Int8, s2a = rem(only_in_2[1] - 1, norb/2) - l, fld(only_in_2[1]-1, norb/2) - 0.5
    m2b::Int8, s2b = rem(only_in_2[2] - 1, norb/2) - l, fld(only_in_2[2]-1, norb/2) - 0.5

    if s1a == s2a
        H_elem += sign .* sum(get_ck(l, m1a, l, m2a) .* get_ck(l, m2b, l, m1b) .* F[1:l+1])
        # H_elem += sign .* sum(get_ck(l, m1a, l, m2a) .* get_ck(l, m1b, l, m2b) .* F)
    end

    if s1a == s2b
        H_elem -= sign .* sum(get_ck(l, m1a, l, m2b) .* get_ck(l, m2a, l, m1b) .* G[1:l+1])
        # H_elem -= sign .* sum(get_ck(l, m1a, l, m2b) .* get_ck(l, m1b, l, m2a) .* G)
    end

    return H_elem

end

"""
Returns the matrix element for the crystal field Hamiltonian between two D states. 
Reads V_CEF from "tables.jl", which is currently only defined for a d atom in a octahedral environment
params:
    state1, state2: basis states (ints)
    norb: number of orbitals represented by state1 and state2 (int)
    Dq10: value of 10Dq, which represents the energy difference between the eg and t2g manifolds
returns:
    (float) value of matrix element
"""
function crystal_field_elem(state1, state2, norb, Dq10)
    # Diagonal element:
    l::Int8 = ((norb/2) - 1)/2
    if state1 == state2
        H_elem = 0
        for e in baddrs(state1)
            ml::Int8 = rem(e - 1, norb/2) - l
            H_elem += V_CEF[ml, ml] * Dq10 / 10
        end
        if hole_lang
            H_elem *= -1
        end
        return H_elem
    end

    # Off diagonal element
    elecs1 = baddrs(state1)
    elecs2 = baddrs(state2)
    both = state1 & state2
    only1 = baddrs(xor(state1, both))
    only2 = baddrs(xor(state2, both))
    if (length(only1) != 1) || (length(only2) != 1)
        return 0    # should only move one electron
    end
    
    if fld(only1[1] - 1, norb/2) != fld(only2[1] -1, norb/2)
        return 0    # no spin flips
    end

    num_ex = length(elecs1[elecs1 .< only2[1]]) + length(elecs2[elecs2 .< only1[1]])
    sign = (-1)^num_ex

    if hole_lang
        sign *= -1
    end

    ml1::Int8 = rem(only1[1] - 1, norb/2) - l
    ml2::Int8 = rem(only2[1] - 1, norb/2) - l

    return sign*V_CEF[ml1, ml2] * Dq10 / 10
end

"""
In some texts, the radial integrals are defined with Racah parameters instead of Slater-Condon parameters (F).
This helper function converts between the two
"""
function F_to_Racah(F)
    A = F[1] - 49 * F[3] / 441
    B = F[2] / 49 - 5 * F[3] / 441
    C = 35 * F[3] / 441
    return A, B, C
end

"""
This function builds a Hamiltonian containing elements from Coulomb interaction and crystal fields
    H = H_int + H_cef 
    Useful for building Tunabe Sugano diagrams
params:
    n_elec: number of electrons in the system
    norb: number of orbitals
    F, G: Slater-Condon parameters
    Dq10: 10Dq (energy difference btween eg and t2g)
returns:
    real-valued Hamiltonian, a NxN matrix, where N is the size of the Hilbert space
"""
function coulomb_and_cef(nelec, norb, F, G, Dq10)
    hil = make_Hilbert_space(nelec, norb)
    num = length(hil)
    ham = zeros(num, num)
    for i in 1:num
        for j in 1:num
            ham[i, j] += h_int_elem(hil[i], hil[j], norb, F, G)
            ham[i, j] += crystal_field_elem(hil[i], hil[j], norb, Dq10)
        end
    end
    return real(eigvals(ham))
end

"""
This function generates data and plots Tunabe-Sugano diagrams, which show how the eigenenergies change relative to one another as 10Dq increases
This function generates plots energies relative to the ground state vs. 10Dq in units of Racah parameter B 
params:
    nelec, norb: number of electrons and orbitals in the system
    F, G: Slater-Condon parameters
"""
function tunabe_sugano_data(nelec, norb, F, G)
    hil = make_Hilbert_space(nelec, norb)
    num = length(hil)
    A, B, C = F_to_Racah(F)
    Dq10 = LinRange(0, 35*B, 30)
    eigs_to_plot = zeros(length(Dq10), num)
    for cef in 1:length(Dq10)
        ham = zeros(num, num)
        for i in 1:num
            for j in 1:num
                ham[i, j] += h_int_elem(hil[i], hil[j], norb, F, G)
                ham[i, j] += crystal_field_elem(hil[i], hil[j], norb, Dq10[cef])
            end
        eigs_to_plot[cef, :] = real(eigvals(ham))
            
        end
    end

    ground_states = []
    for i in 1:length(Dq10)
        append!(ground_states, minimum(eigs_to_plot[i, :]))
    end

    println(ground_states)

    p = plot(legend=false)
    for i in 1:num
        plot!(p, Dq10./B, (eigs_to_plot[:, i] .- ground_states) ./ B)
    end
    display(p)
    return Dq10, eigs_to_plot
end

"""
Generates Hamiltonian matrix element for spin-orbit coupling
params:
    state1, state2: basis states (ints)
    norb: number of orbitals
    l_soc: strength of spin-orbit coupling
returns:
    (float) Hamiltonian matrix element 
"""
function spin_orbit_elem(state1, state2, norb, l_soc)
    # Diagonal element:
    l::Int8 = ((norb/2) - 1)/2
    if state1 == state2
        H_elem = 0
        for e in baddrs(state1)
            ml::Int8 = rem(e - 1, norb/2) - l
            sigma = fld(e-1, norb/2) * 2 - 1   # +1 if spin up, -1 if spin down
            H_elem += sigma*ml
        end
        if hole_lang
            H_elem *= -1
        end
        return H_elem * l_soc / 2
    end

    # Off diagonal element
    elecs1 = baddrs(state1)
    elecs2 = baddrs(state2)
    both = state1 & state2
    only1 = baddrs(xor(state1, both))
    only2 = baddrs(xor(state2, both))
    if (length(only1) != 1) || (length(only2) != 1)
        return 0    # should only move one electron
    end
    

    ml1::Int8 = rem(only1[1] - 1, norb/2) - l
    ms1::Int8 = fld(only1[1] - 1, norb/2)   # 0 if spin down, 1 if spin up
    ml2::Int8 = rem(only2[1] - 1, norb/2) - l
    ms2::Int8 = fld(only2[1] - 1, norb/2)   # 0 if spin down, 1 if spin up
    if ml1 + ms1 != ml2 + ms2
        return 0    # no spin flips
    end

    num_ex = length(elecs1[elecs1 .< only2[1]]) + length(elecs2[elecs2 .< only1[1]])
    sign = (-1)^num_ex

    if hole_lang
        sign *= -1
    end

    m = min(ml1, ml2)

    return sign*l_soc*sqrt((l - m) * (l + m + 1)) / 2
end

"""
Generic function for building a single-shell Hamiltonian. Can turn on/off what you want to include
params:
    nelec: number of electrons
    norb: number of orbitals
    F, G: Slater-Condon params
    Dq10: value of 10Dq
    l_soc: strength of spin-orbit coupling
    coulomb: turns on H_int (Coulomb interaction)(defaults to true)
    cef: turns on H_cef (crystal field) (defaults to false)
    soc: turns on H_soc (spin-orbit coupling) (defaults to false)
returns: 
    ham: NxN matrix of floats, where N is the size of the Hilbert space. Diagonalize to get eigenvalues/eigenvectors
"""
function build_hamiltonian(nelec, norb, F, G; Dq10=0, l_soc=0, coulomb=true, cef=false, soc=false)
    hil = make_Hilbert_space(nelec, norb)
    num = length(hil)
    ham = zeros(num, num)
    for i in 1:num
        for j in 1:num
            # these three do alot of the same calculations -- can probably optimize later
            if coulomb
                ham[i, j] += h_int_elem(hil[i], hil[j], norb, F, G)
            end
            if cef
                ham[i, j] += crystal_field_elem(hil[i], hil[j], norb, Dq10)
            end
            if soc
                ham[i, j] += spin_orbit_elem(hil[i], hil[j], norb, l_soc)
            end
        end
    end
    return ham
end

# only supports 2 shells (core and valence), need recursive to include ligands
"""
Creates a multi-shell Hilbert space, e.g. d-electron valence and p-electron core.
currently only 2 shells is implemented
params:
    num_elecs: (vector of ints) number of electrons in each shell in order
    l_vals: value of l for both shells. XAS functions assume the default order = [2, 1]
returns
    hilbert: A vector of length-2 tuples of ints representing the Hilbert space
"""
function multi_shell_hilbert(num_elecs, l_vals=[2, 1])
    @assert length(num_elecs) == length(l_vals)
    num_orbs = 2 .* (2 .* l_vals .+ 1)
    hilbert::Vector{Tuple{Int16, Int16}} = []
    subhil1 = make_Hilbert_space(num_elecs[1], num_orbs[1])
    subhil2 = make_Hilbert_space(num_elecs[2], num_orbs[2])
    for i in subhil1
        for j in subhil2
            push!(hilbert, ((i, j)))
        end
    end

    return hilbert
end


### hil1, hil2 consist of tuples, e.g. (Int describing d shell occupation, Int describing p shell occupation)
# polarization in [ex, ey, ez]
# only supports two shells (core and valence) -- do we need more?
"""
Creates a complex-values matrix for dipole transitions between basis states of two multi-shelled Hilbert spaces
params:
    hil1, hil2: two multi-shell hilbert spaces (e.g. for ground state and excited state)
    polarization: photon polarization as 1d vector [ex, ey, ez]
    l_vals: l values for corresponding shells in hil1 and hil2, Use [2, 1] for d valence electrons and p core electrons
returns:
    table: a N1xN2 complex-values matrix describing dipole transitions betweens states in hil1 and hil2 
"""
function dipole_in_basis(hil1, hil2, polarization, l_vals=[2, 1])
    n_vals = Dict{Int8, Float64}(-1 => (polarization[1] + 1im*polarization[2])/sqrt(2),
                                0 => polarization[3],
                                1 => (-1*polarization[1] + 1im*polarization[2])/sqrt(2))
    @assert length(hil1[1]) == length(l_vals) && length(hil2[1]) == length(l_vals) "mismatch in number of shells"
    num_shells = length(l_vals)
    n_elec1 = [length(baddrs(hil1[1][x])) for x in 1:num_shells]
    n_elec2 = [length(baddrs(hil2[1][x])) for x in 1:num_shells]
    table::Matrix{ComplexF64} = zeros(length(hil1), length(hil2))
    for i in 1:length(hil1)
        for j in 1:length(hil2)
            ml::Vector{Int8} = [0, 0]
            l::Vector{Int8} = [0, 0]
            ms = [0, 0]
            check1, check2 = false, false
            num_ex::Int64 = 0
            for k in 1:num_shells
                both = hil1[i][k] & hil2[j][k]
                all1 = baddrs(hil1[i][k])
                all2 = baddrs(hil2[j][k])
                only1 = baddrs(xor(both, hil1[i][k]))
                only2 = baddrs(xor(both, hil2[j][k]))
                if length(only1) == 1 && length(only2) == 0
                    l[1] = l_vals[k]
                    ml[1] = rem(only1[1] - 1, 2*l_vals[k] + 1) - l_vals[k]
                    ms[1] = fld(only1[1] -1, 2*l_vals[k] + 1)       # 1 if spin up, 0 if down
                    check1 = true
                    num_ex += length(all2[all2 .< only1[1]])
                    if k == 2
                        num_ex += n_elec2[1]
                    end
                end
                if length(only2) == 1 && length(only1) == 0
                    l[2] = l_vals[k]
                    ml[2] = rem(only2[1] - 1, 2*l_vals[k] + 1) - l_vals[k]
                    ms[2] = fld(only2[1] - 1, 2*l_vals[k] + 1)       # 1 if spin up, 0 if down
                    check2 = true
                    num_ex += length(all1[all1 .< only2[1]])
                    if k == 2
                        num_ex += n_elec1[1]
                    end
                end
            end

            sign = (-1)^num_ex

            if check1 && check2 && (ms[2] == ms[1]) && (abs(ml[2] - ml[1]) < 2)
                if hole_lang
                    table[i, j] = sign * (-1)^(ml[1] - ml[2] + 1) * get_ck(l[2], ml[2], l[1], ml[1])[1] * n_vals[ml[1] - ml[2]]
                else 
                    table[i, j] = sign * get_ck(l[1], ml[1], l[2], ml[2])[1] * n_vals[ml[1] - ml[2]]
                end
            end
        end
    end
    return table
end

                
"""
Calculates a Coulomb interactions matrix element between core holes and valence electrons
params:
    state1, state2: tuples of ints representing electron/hole occupation in two shells
    Fcv, Gcv: core-valence Slater-Condon parameters (defined in "params.jl")
returns:
    H_elem: a matrix element (float)
"""
function core_valence_int_elem(state1, state2, Fcv, Gcv; l_vals=[2, 1])
    n_orbs = 4 .* l_vals .+ 2
    # must conserve Lz and Sz
    if get_Lz_Sz(state1[1], n_orbs[1]) .+ get_Lz_Sz(state1[2], n_orbs[2]) != 
        get_Lz_Sz(state2[1], n_orbs[1]) .+ get_Lz_Sz(state2[2], n_orbs[2])
        return 0
    end

    num_k = fld(sum(l_vals), 2) + 1

    # diagonal elements
    if state1 == state2
        H_elem = 0.0
        indices_d = baddrs(state1[1])
        indices_p = baddrs(state1[2])
        for i in indices_d
            ml1::Int8 = rem(i-1, n_orbs[1]/2) - l_vals[1]
            ms1 = fld(i-1, n_orbs[1]/2) - 0.5
            for j in indices_p
                ml2::Int8 = rem(j-1, n_orbs[2]/2) - l_vals[2]
                ms2 = fld(j-1, n_orbs[2]/2) - 0.5
                H_elem += sum(get_ck(l_vals[1], ml1, l_vals[1], ml1)[1:num_k] .* get_ck(l_vals[2], ml2, l_vals[2], ml2)[1:num_k] .* Fcv)
                if ms1 == ms2
                    H_elem -= sum(get_ck(l_vals[1], ml1, l_vals[2], ml2) .* get_ck(l_vals[1], ml1, l_vals[2], ml2) .* Gcv) 
                end
            end
        end
        return H_elem
    end

    # off-diagonal element
    H_elem = 0.0
    elecs1 = [baddrs(x) for x in state1]
    elecs2 = [baddrs(x) for x in state2]
    in_both = [state1[x] & state2[x] for x in 1:2]
    only_in_1 = [baddrs(xor(state1[x], in_both[x])) for x in 1:2]
    only_in_2 = [baddrs(xor(state2[x], in_both[x])) for x in 1:2]
    
    if !(length(only_in_1[1]) == length(only_in_1[2]) == length(only_in_2[1]) == length(only_in_2[2]))
        return 0    # both states need to move only 1 core and 1 valence elec/hole 
    end
    num_ex = 0
    for i in 1:length(l_vals)
        num_ex += length(elecs2[i][elecs2[i] .< only_in_1[i][1]])
        num_ex += length(elecs1[i][elecs1[i] .< only_in_2[i][1]])
    end
    sign = (-1)^num_ex

    m1a::Int8, s1a = rem(only_in_1[1][1] - 1, n_orbs[1]/2) - l_vals[1], fld(only_in_1[1][1] - 1, n_orbs[1] / 2) - 0.5 # d
    m1b::Int8, s1b = rem(only_in_1[2][1] - 1, n_orbs[2]/2) - l_vals[2], fld(only_in_1[2][1] - 1, n_orbs[2] / 2) - 0.5 # p
    m2a::Int8, s2a = rem(only_in_2[1][1] - 1, n_orbs[1]/2) - l_vals[1], fld(only_in_2[1][1] - 1, n_orbs[1] / 2) - 0.5
    m2b::Int8, s2b = rem(only_in_2[2][1] - 1, n_orbs[2]/2) - l_vals[2], fld(only_in_2[2][1] - 1, n_orbs[2] / 2) - 0.5 # p

    

    if s1a == s2a
        H_elem += sign .* sum(get_ck(l_vals[1], m1a, l_vals[1], m2a)[1:num_k] .* get_ck(l_vals[2], m2b, l_vals[2], m1b)[1:num_k] .* Fcv)
    end

    if s1a == s2b
        H_elem -= sign .* sum(get_ck(l_vals[1], m1a, l_vals[2], m2b) .* get_ck(l_vals[1], m2a, l_vals[2], m1b) .* Gcv)
    end

    return H_elem

end

### hilbert should be a vector of 2D tuples, i.e. (valencestate, corestate)
# SOC only applies to p orbitals(core)
# CEF only applies to d orbitalas (valence)
"""
This function builds a Hamiltonian matrix that includes a core shell. Basis is the multi-shell Hilbert space
params:
    hil: vector of tuples of ints, e.g. output of multi_shell_hilbert()
    Fd, Gd: Slater-Condon parameters for d-shell Coulomb interactions (vectors of floats)
    Dq10: value of 10Dq (float)
    l_soc: SOC strength (float) 
    coulombd: turns on d-electron coulomb interactions
    cef: turns on d-electron crystal field splitting
    soc: turns on p-orbital spin orbit coupling
    cv_int: turns on core-valence interactions
return:
    ham: a NxN real-valued matrix, N=size of Hilbert space 
"""
function build_hamiltonian_with_core(hil, Fd, Gd; Dq10=0, l_soc=0, Fcv=[0., 0.], Gcv=[0., 0.], coulombd=true, cef=false, soc=false, cv_int=false)
    num = length(hil)
    ham = zeros(num, num)
    for i in 1:num
        for j in 1:num
            same_core = (hil[i][2] == hil[j][2])
            same_valence = (hil[i][1] == hil[j][1])
            # these three do alot of the same calculations -- can probably optimize later
            if coulombd && same_core   # same core states
                ham[i, j] += h_int_elem(hil[i][1], hil[j][1], 10, Fd, Gd)
            end
            # if coulombp && same_valence     # same valence state
            #     ham[i, j] += h_int_elem(hil[i][2], hil[j][2], 6, Fp, Gp)
            # end
            if cef && same_core
                ham[i, j] += crystal_field_elem(hil[i][1], hil[j][1], 10, Dq10)
            end
            if soc && same_valence
                ham[i, j] += spin_orbit_elem(hil[i][2], hil[j][2], 6, l_soc)
            end
            if cv_int
                ham[i, j] += core_valence_int_elem(hil[i], hil[j], Fcv, Gcv)
            end
        end
    end
    return ham
end


"""
This function calcualtes the L-edge XAS spectrum for a single atom in a crystal field. 
params:
    n_elec: number of valence (d-shell) electrons (or holes if hole_lang==true)
    polarization: the incident light polarization as [ex, ey, ez]
    Fd, Gd: Slater-Condon parameters for d-shell Coulomb interactions (vectors of floats)
    Dq10: value of 10Dq (float)
    l_soc: SOC strength (float) 
    coulombd: turns on d-electron coulomb interactions
    cef: turns on d-electron crystal field splitting
    soc: turns on p-orbital spin orbit coupling
    cv_int: turns on core-valence interactions
returns: 
    omegas: vector of floats containing all the energy differences captured by XAS => E_excited - E_gs 
    intensities: strength of the signal for each omega. Vector of floats with the same length as omegas
    Outputs should be plotted with some Gaussian or Lorentzian smearing
"""
function XAS_L_edge(n_elec, polarization, Fd, Gd; Dq10=0, l_soc=0, Fcv=[0, 0], Gcv=[0, 0], coulombd=true, cef=false, soc=false, cv_int=false)
    if hole_lang
        gs_hil = multi_shell_hilbert([n_elec, 0])
        es_hil = multi_shell_hilbert([n_elec-1, 1])
    else
        gs_hil = multi_shell_hilbert([n_elec, 6])
        es_hil = multi_shell_hilbert([n_elec+1, 5])
    end
        dipole_table = dipole_in_basis(gs_hil, es_hil, polarization)

    eig_gs = eigen(build_hamiltonian_with_core(gs_hil, Fd, Gd, Dq10=Dq10, l_soc=l_soc, coulombd=coulombd, cef=cef, soc=soc, cv_int=false))
    eig_es = eigen(build_hamiltonian_with_core(es_hil, Fd, Gd, Dq10=Dq10, l_soc=l_soc, Fcv=Fcv, Gcv=Gcv, coulombd=coulombd, cef=cef, soc=soc, cv_int=cv_int))
    # println(x for x in eig_gs.values)
    E0 = minimum(eig_gs.values)
    println(E0)

    omegas::Vector{Float64} = []
    intensities::Vector{Float64} = []

    for i in 1:length(eig_gs.values)
        if abs(eig_gs.values[i] - E0) > 1E-6 
            continue
        end
        for j in 1:length(eig_es.values)
            dE = eig_es.values[j] - (eig_gs.values[i])
            vec1 = eig_gs.vectors[:, i]
            vec2 = eig_es.vectors[:, j]
            dip::ComplexF64 = 0
            for k in 1:length(vec1)
                for l in 1:length(vec2)
                    dip += vec1[k] * vec2[l] * dipole_table[k, l]
                end
            end
            push!(omegas, dE)
            push!(intensities, dip * conj(dip))
        end
    end

    return omegas, intensities
end










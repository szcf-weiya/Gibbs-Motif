using StatsBase
using Distributions
using SpecialFunctions
# read file into string
s = open("hw2_part2_data.txt") do file
    read(file, String)
end
# split string into substring (DNA) via `>`
R = split(s, ">", keepempty = false)
# number of DNA
N = length(R) # 30
# extract the DNA sequence name
DNAnames = String[]
L = ones(Int64, N)
for i = 1:N
    r = split(R[i], "\r\n", keepempty = false)
    push!(DNAnames, r[1])
    R[i] = join(r[2:end])
    L[i] = length(R[i])
end

# prior
W = 18
ALPHA = ones(4) / 4
BETA = ones(W, 4)


# count the number of AGCT
function countAGCT(s::T) where {T <: AbstractString}
    numA = 0
    numG = 0
    numC = 0
    numT = 0
    for c in collect(s)
        if c == 'A' || c == 'a'
            numA += 1
        elseif c == 'G' || c == 'g'
            numG += 1
        elseif c == 'C' || c == 'c'
            numC += 1
        elseif c == 'T' || c == 't'
            numT += 1
        end
    end
    return [numA, numG, numC, numT]
end

function countAGCT(s::Array)
    countAGCT(join(s))
end


function get_site(A::Array, k::Int64) # k = 0 represents no except
    site = String[]
    for i = 1:N
        if i == k
            push!(site, "")
        else
            push!(site, R[i][A[i]:A[i]+W-1])
        end
    end
    return site
end

# sample theta0
function spl_theta0(A::Array)
    # extract site region
    site = get_site(A, 0)
    # count Ac
    freq_all = countAGCT(R)
    freq_site = countAGCT(site)
    freq_nonsite = freq_all - freq_site
    param = freq_nonsite + ALPHA
    return rand(Dirichlet(param))
end

# sample A = (a_1, \ldots, a_N)
function spl_A!(A::Array{Int64}, theta0)
    for k = 1:N
        #spl_ak2!(k, A, theta0)
        spl_ak!(k, A, theta0)
    end
end

function spl_ak!(k::Int64, A::Array{Int64}, theta0)
    # extract current motif region except k
    site_nok = get_site(A, k)
    pos = join(site_nok)
    # calculate probs at position i
    probs = ones(L[k]-W+1)
    for i = 1:L[k]-W+1
        sitek = R[k][i:i+W-1]
        freq_sitek = countAGCT(sitek)
        for j = 1:W
            freq = countAGCT(pos[(0:(N-2))*W .+ j])
            theta = freq .+ BETA[j,:]
            if sitek[j] == 'A' || sitek[j] == 'a'
                probs[i] *= theta[1] / theta0[1]
            elseif sitek[j] == 'G' || sitek[j] == 'g'
                probs[i] *= theta[2] / theta0[2]
            elseif sitek[j] == 'C' || sitek[j] == 'c'
                probs[i] *= theta[3] / theta0[3]
            elseif sitek[j] == 'T' || sitek[j] == 't'
                probs[i] *= theta[4] / theta0[4]
            end
        end
    end
    # sample ak according to the probs.
    A[k] = sample(1:L[k]-W+1, pweights(probs))
end

function compare_site(site1, site2; delta = 1)
    pos1 = join(site1)
    pos2 = join(site2)
    if delta == 1
        freq1 = countAGCT(pos1[(0:(N-1))*W .+ 1])
        freq2 = countAGCT(pos2[(0:(N-1))*W .+ W])
    else
        freq1 = countAGCT(pos1[(0:(N-1))*W .+ W])
        freq2 = countAGCT(pos2[(0:(N-1))*W .+ 1])
    end                
    lg2 = sum(loggamma.(filter(!iszero, freq2)))
    lg1 = sum(loggamma.(filter(!iszero, freq1)))
    if lg2 > lg1
        p = 1
    else
        p = exp(lg2 - lg1)
    end
    return p
end

# shift
function shift!(A)
    site = get_site(A, 0)
    if rand() < 1/2  # delta = 1
        while maximum(A) < L[1] - W + 1 # here L are the same.
            newsite = get_site(A .+ 1, 0)
            p = compare_site(site, newsite)
            if rand() < p
                A .= A .+ 1
            else
                break
            end
        end
    else # delta = -1
        while minimum(A) > 1
            newsite = get_site(A .- 1, 0)
            p = compare_site(site, newsite, delta = -1)
            if rand() < p
                A .= A .- 1
            else
                break
            end
        end
    end
end

# gibbs sampler
function gibbs(M = 3000, Burnin = 1000)
    # initial
    A = ones(Int64, N)
    resA = ones(Int64, M, N)
    restheta0 = ones(M, 4)
#    theta0 = [0.208, 0.302, 0.319, 0.171]
    for t = 1:M + Burnin
        println("t = $t...")
        # sample theta0
        theta0 = spl_theta0(A)
        # sample A
        spl_A!(A, theta0)
        # shifting
        shift!(A)
        if t > Burnin
            resA[t - Burnin,:] .= A
            restheta0[t - Burnin,:] .= theta0
        end
    end
    return resA, restheta0
end

# estimate A
function estA(resA)
    A = ones(Int, N)
    for i = 1:N
        A[i] = mode(resA[:,i])
    end
    return A
end

# calculate theta_jk
function cal_theta(A)
    site = get_site(A, 0)
    pos = join(site)
    theta = ones(W, 4)
    for j = 1:W
        freq = countAGCT(pos[(0:(N-1))*W .+ j])
        theta[j, :] = freq ./ sum(freq)
    end
    return theta
end

# run
M = 1000
B = 1000
resA, restheta0 = gibbs(M, B)
using DelimitedFiles
using Dates
timestamp = Dates.now()
# save results
# writedlm("resA$(timestamp).txt", resA)
# writedlm("restheta0$(timestamp).txt", restheta0)

# estimate A
Ahat = estA(resA)
motif = get_site(Ahat,0)
# estimate theta0
theta0 = mean(restheta0, dims=1)
# eatimate theta
theta = cal_theta(Ahat)
# save motif
writedlm("resmotif$(timestamp).txt", motif)

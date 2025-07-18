using Combinatorics
using Hecke
using SaferIntegers
using GAP


######################################################
# THRESHOLD/RESONANCE ARRANGEMENT
######################################################


function threshold_hyperplanes(n, ring=nothing)
    M=[]
    for i in 0:2^n-1
        M=vcat(M,vcat([1],digits(i,base=2,pad=n)))
    end
    M = convert(Array{Int64,2},reshape(M,n+1,2^n))
    if ring != nothing
        M = [ring(i) for i in M]
    else
        M = convert(Array{SafeInt64,2}, M)
    end
    return M
end

function resonance_hyperplanes(n, ring=nothing)
    M=[]
    for i in 1:2^n-1
        M=vcat(M,digits(2^n-i,base=2,pad=n))
    end
    M = convert(Array{Int64,2},reshape(M,n,2^n-1))
    if ring != nothing
        M = [ring(i) for i in M]
    else
        M = convert(Array{SafeInt64,2}, M)
    end
    return M
end

#This function will define the permutations that multiply a row by -1
function reflections(step::Int64, n::Int64)
    incr = 2^step
    pairs = [[1, 1+incr]]
    current = 2
    while size(pairs, 1) < 2^(n-1)
        if current in vcat(pairs...)
            current = current+1
            continue
        end
        push!(pairs, [current, current+incr])
        current = current+1
    end
    return pairs
end

#This will define the permutation that shift the coordintates cyclically
function index_shift_perm(k::Int64, n::Int64)
    k1 = k-1
    if 2*k1 >=2^n
        return 2*k1-2^n+2
    else
        return 2*k1+1
    end
end

#This will define the permutation that swaps two coordinates
function index_swap(n::Int64)
    pairs = []
    current = 1
    for k in 1:2^(n-2)
        push!(pairs, [current, current+3])
        current = current+4
    end
    return pairs
end

function extra_symmetry(n::Int64)
    pairs = []
    if n == 1
        return pairs
    end

    pair_sum = 2^(n-1)+1
    for k in 1:2^(n-2)
        push!(pairs, [k, pair_sum-k])
    end
    return pairs
end

function symmetry_threshold(n::Int64)
    perms = []
    p_cyc = GAP.Globals.MappingPermListList(GapObj(2:(2^n-1)), GapObj([index_shift_perm(k, n) for k in 2:(2^n-1)]))
    push!(perms, p_cyc)
    for step in 0:(n-1)
        push!(perms, cycle_list_to_perm(reflections(step, n)))
    end
    push!(perms, cycle_list_to_perm(index_swap(n)))
    push!(perms, cycle_list_to_perm(extra_symmetry(n)))
    return GAP.Globals.Group(GapObj(perms))
end


function symmetry_resonance(n::Int64)
    return GAP.Globals.Stabilizer(symmetry_threshold(n), GapObj(1:(2^n)-1), GAP.Globals.OnSets)
end

######################################################
# PERMUTOHEDRON
######################################################

function permutohedron_hyperplanes(n::Int64, ring=nothing)
    M = [vcat([1],v) for v in permutations(Array(1:n))]
    M = [v[i] for i in 1:n+1, v in M]
    if ring != nothing
        M = [ring(i) for i in M]
    else
        M = convert(Array{SafeInt128,2}, M)
    end
    return M
end

function symmetry_permutohedron(n::Int64)
    matrix = permutohedron_hyperplanes(n)
    permuted_row_permutations = GapObj(find_all_permutations_from_row_permutations(matrix,2))

    return GAP.Globals.Group(permuted_row_permutations)
end

######################################################
# CROSSPOLYTOPES
######################################################

function crosspolytope_hyperplanes(n::Int64, ring=nothing)
    M = []
    for i in 0:n-1
        M=vcat(M, vcat([1], digits(2^i, base = 2, pad = n)))
        M=vcat(M, vcat([1], digits(-(2^i), base = 2, pad = n)))
    end
    M = convert(Array{Int64,2},reshape(M,n+1,2*n))
    if ring != nothing
        M = [ring(i) for i in M]
    else
        M = convert(Array{SafeInt64,2}, M)
    end
    return M
end

function symmetry_crosspolytope(n::Int64)
    return GAP.Globals.WreathProduct(GAP.Globals.SymmetricGroup(2),GAP.Globals.SymmetricGroup(n))
end

######################################################
# DEMICUBES
######################################################

# Vertices: The vertices of the n-cube with an odd number of 1's
function indices_demicube(n::Int64)::Array{Int64,1}
    dc_indices = Array{Int64,1}()
    for i in 1:2^n-1
        if isodd(sum(digits(i,base=2)))
            push!(dc_indices, i)
        end
    end
    return dc_indices
end

function demicube_hyperplanes(n::Int64, ring=nothing)
    M = []
    for i in indices_demicube(n)
        M = vcat(M,vcat([1],digits(i,base=2,pad=n)))
    end
    M = convert(Array{Int64,2},reshape(M,n+1,2^(n-1)))
    if ring != nothing
        M = [ring(i) for i in M]
    else
        M = convert(Array{SafeInt64,2}, M)
    end
    return M
end

function symmetry_demicube(n)
    G = symmetry_threshold(n)
    ind = indices_demicube(n)
    ind_shifted = [i+1 for i in ind]
    AutDC = GAP.Globals.Stabilizer(G, GapObj(ind_shifted), GAP.Globals.OnSets)
    perms = []
    for j in 1:GAP.Globals.Length(GAP.Globals.GeneratorsOfGroup(AutDC))
        perm = GAP.Globals.GeneratorsOfGroup(AutDC)[j]
        permuted_inds = [GAP.Globals.ListPerm(perm, 2^n)[ind_shifted[i]] for i in 1:2^(n - 1)]
        permuted_inds = GapObj([findfirst(isequal(i), ind_shifted) for i in permuted_inds])
        new_perm = GAP.Globals.PermList(permuted_inds)
        push!(perms, new_perm)
    end
    return GAP.Globals.Group(GapObj(perms))
end

######################################################
# REGULAR POLYTOPES
######################################################

function pos_neg_combinations(l, symmetry_group)
    M = []
    group_size = GAP.Globals.Size(symmetry_group)
    len = length(l)
    type = typeof(l[1])

    for negative_indices in powerset(1:len)
        current_l=copy(l)
        for ind in negative_indices
            current_l[ind] = -1*current_l[ind]
        end
        for i in 1:group_size
             push!(M,Vector{type}(GAP.Globals.Permuted(GapObj(current_l),GAP.Globals.List(symmetry_group)[i])))
        end
    end
    return unique(M)
end

function dodecahedron_hyperplanes()
    Qx, x = PolynomialRing(FlintQQ, "x");
    K, phi = NumberField(x^2 -x -1, "phi");

    M = []
    S = GAP.Globals.SymmetricGroup(3)
    A = GAP.Globals.AlternatingGroup(3)

    l = nf_elem[K(1), K(1), K(1)]
    M = vcat(M,pos_neg_combinations(l, S))

    l = nf_elem[K(phi), K(phi+1), K(0)]
    M = vcat(M,pos_neg_combinations(l, A))

    h20 = [vcat(nf_elem[K(1)],v) for v in M]
    h20 = [K(j) for k in 1:20 for j in h20[k]]
    h20 = reshape(h20, 4,20)
    return h20
end

function symmetry_dodecahedron()
    M = dodecahedron_hyperplanes()
    return find_all_permutations(M, 2)
end

function icosahedron_hyperplanes()
    Qx, x = PolynomialRing(FlintQQ, "x");
    K, phi = NumberField(x^2 -x -1, "phi");

    M = []
    A = GAP.Globals.AlternatingGroup(3)

    l = nf_elem[K(1), K(phi), K(0)]
    M = vcat(M,pos_neg_combinations(l, A))

    h12 = [vcat(nf_elem[K(1)],v) for v in M]
    h12 = [K(j) for k in 1:12 for j in h12[k]]
    h12 = reshape(h12, 4,12)
    return h12
end

function symmetry_icosahedron()
    M = icosahedron_hyperplanes()
    return find_all_permutations(M, 2)
end

# Vertices: all permuations of  (±1,±1,0,0)
function regular_24_cell_hyperplanes()
    S = GAP.Globals.SymmetricGroup(4)
    M = pos_neg_combinations([1,1,0,0], S)
    M = [vcat([1],v) for v in M]
    hyperplanes_24 = [v[i] for i in 1:5, v in M]
    return hyperplanes_24
end

function symmetry_24_cell()
    M = regular_24_cell_hyperplanes()
    return find_all_permutations(M, 2)
end

function regular_600_cell_hyperplanes()
    Qx, x = PolynomialRing(FlintQQ, "x");
    K, phi = NumberField(x^2 -x -1, "phi");

    M = []
    S = GAP.Globals.SymmetricGroup(4)
    A = GAP.Globals.AlternatingGroup(4)

    l = nf_elem[K(2), K(0), K(0), K(0)]
    M = vcat(M,pos_neg_combinations(l, S))

    l = nf_elem[K(1), K(1), K(1), K(1)]
    M = vcat(M,pos_neg_combinations(l, S))

    l = nf_elem[phi, K(1), phi^-1,K(0)]
    M = vcat(M,pos_neg_combinations(l, A))

    h600 = [vcat(nf_elem[K(1)],v) for v in M]
    h600 = [K(j) for k in 1:120 for j in h600[k]]
    h600 = reshape(h600, 5,120)
    return h600
end

function symmetry_600_cell()
    M = regular_600_cell_hyperplanes()
    return find_all_permutations(M, 2)
end

function regular_120_cell_hyperplanes()
    Qx, x = PolynomialRing(FlintQQ, "x");
    K, phi = NumberField(x^2 -x -1, "phi");
    r5 = 2*phi-1

    M = []
    S = GAP.Globals.SymmetricGroup(4)
    A = GAP.Globals.AlternatingGroup(4)
    
    #all permutations of (±2,±2,0,0)
    l = nf_elem[K(2),K(2),K(0),K(0)]
    M = vcat(M, pos_neg_combinations(l, S))

    #all permutations of (±1,±1,±1,±sqrt(5))
    l = nf_elem[K(1),K(1),K(1),r5]
    M = vcat(M, pos_neg_combinations(l, S))

    #all permutations of (±phi^-2,±phi,±phi,±phi)
    l = nf_elem[phi^(-2),phi,phi,phi]
    M = vcat(M, pos_neg_combinations(l, S))

    #all permutations of (±phi^-1,±phi^-1,±phi^-1,±phi^2)
    l = nf_elem[phi^(-1),phi^(-1),phi^(-1),phi^2]
    M = vcat(M, pos_neg_combinations(l, S))

    #all even permutations of (0,±phi^-2,±1,±phi^2)
    l = nf_elem[K(0),phi^(-2),K(1),phi^2]
    M = vcat(M, pos_neg_combinations(l, A))

    #all even permutations of (0,±phi^-1,±phi,±sqrt(5))
    l = nf_elem[K(0),phi^(-1),phi,r5]
    M = vcat(M, pos_neg_combinations(l, A))

    #all even permutations of (±phi^-1,±1,±phi,±2)
    l = nf_elem[phi^(-1),K(1),phi,K(2)]
    M = vcat(M, pos_neg_combinations(l, A))

    h120 = [vcat(nf_elem[K(1)],v) for v in M]
    h120 = [K(j) for k in 1:600 for j in h120[k]]
    h120 = reshape(h120, 5,600)
    return h120
end

function symmetry_120_cell()
    M = regular_120_cell_hyperplanes()
    return find_all_permutations(M, 2)
end


######################################################
# DISCRIMINANT/SOFT DISCRIMINANT ARRANGEMENT
######################################################


function discriminantal_hyperplanes(MyPoints::Array{QQFieldElem,2})
    Npoints=size(MyPoints,1)
    PointDim=size(MyPoints,2)
    MyHomPoints=matrix(QQ,hcat(Array{Int,1}(ones(Npoints)),MyPoints))
    HypIterator=powerset(collect(1:Npoints),PointDim,PointDim)
    Hyperplanes=[QQFieldElem(i) for i in hcat([kernel(MyHomPoints[K,:])[2] for K in HypIterator])]
    ConstantTerms=-Hyperplanes[PointDim+1,1:size(Hyperplanes,2)]
    Hyperplanes=Hyperplanes[1:PointDim,1:size(Hyperplanes,2)]
    return [Hyperplanes,ConstantTerms]
end


function symmetry_discriminantal(d::Int64,n::Int64)
    Npoints=n
    PointDim=d
    HypIterator=powerset(collect(1:Npoints),PointDim,PointDim)
    C=collect(HypIterator)
    function myperm(i,j)
        if i==Npoints && j==i
            return(1)
        elseif i==Npoints && j==1
            return(Npoints)
        end
        if j==i
            return(i+1)
        elseif j==i+1
            return(i)
        else
            return(j)
        end
    end
    permBucket=[]
    Id=GapObj([i for i in 1:length(C)+1])
    for i in 1:Npoints
        Cmoved=[[myperm(i,j) for j in c] for c in C]
        actingOnPlanes=[findfirst(x->sort(x)==c,Cmoved) for c in C]
        push!(actingOnPlanes,length(actingOnPlanes)+1)
        GapPerm=GAP.Globals.MappingPermListList(Id,GapObj(actingOnPlanes))
        push!(permBucket,GapPerm)
    end

    DiscriminantSymmetry=GAP.Globals.Group(GapObj(permBucket))
    return DiscriminantSymmetry
end


function discriminantal_hyperplanes(d::Int64,n::Int64)
    MyPoints=[QQFieldElem(rand(-1000:1000)//rand(-1000:1000)) for i in zeros(n,d)]
    MyFloatPoints=[float(i*1.0) for i in MyPoints]
    D=discriminantal_hyperplanes(MyPoints)
    return D
end

function soft_discriminantal_hyperplanes(d::Int64,n::Int64)
    d=d-1
    n=n-1
    MyPoints=[QQFieldElem(rand(-1000:1000)//rand(-1000:1000)) for i in zeros(n-d,d)]
    I=zeros(d,d)
    for i in 1:d
        I[i,i]=1
    end
    MyID=vcat([QQFieldElem(Int(i)) for i in zeros(d)'],[QQFieldElem(i) for i in Array{Int,2}(I)])
    MyID=MyID'
    MyHomPoints=matrix(QQ,vcat(MyID,hcat(Array{Int,1}(ones(n-d)),MyPoints)))
    Npoints=n
    HypIterator=reverse(collect(powerset(collect(1:Npoints),d,d)))
    Hyperplanes=[QQFieldElem(i) for i in hcat([kernel(MyHomPoints[K,:])[2] for K in HypIterator])]
    ConstantTerms=-Hyperplanes[1,1:size(Hyperplanes,2)-1]
    Hyperplanes=Hyperplanes[2:size(Hyperplanes,1),1:size(Hyperplanes,2)-1]
    return [Hyperplanes,ConstantTerms]
end


function symmetry_soft_discriminantal(d::Int64,n::Int64)
    d=d-1
    n=n-1
    MyPoints=[QQFieldElem(rand(-1000:1000)//rand(-1000:1000)) for i in zeros(n-d,d)]
    I=zeros(d,d)
    for i in 1:d
        I[i,i]=1
    end
    MyID=vcat([QQFieldElem(Int(i)) for i in zeros(d)'],[QQFieldElem(i) for i in Array{Int,2}(I)])
    MyID=MyID'
    MyHomPoints=matrix(QQ,vcat(MyID,hcat(Array{Int,1}(ones(n-d)),MyPoints)))
    Npoints=n
    HypIterator=powerset(collect(1:Npoints),d,d)
    C=reverse(collect(HypIterator))
    function myperm(i,j)
        if i==Npoints && j==i
            return(1)
        elseif i==Npoints && j==1
            return(Npoints)
        end
        if j==i
            return(i+1)
        elseif j==i+1
            return(i)
        else
            return(j)
        end
    end
    permBucket=[]
    Id=GapObj([i for i in 1:length(C)])
    for i in 1:Npoints
        Cmoved=[[myperm(i,j) for j in c] for c in C]
        actingOnPlanes=[findfirst(x->sort(x)==c,Cmoved) for c in C]
        GapPerm=GAP.Globals.MappingPermListList(Id,GapObj(actingOnPlanes))
        push!(permBucket,GapPerm)
    end

    G=GAP.Globals.Group(GapObj(permBucket))

    G = GAP.Globals.Stabilizer(G, GapObj(length(C)))
    return G
end


######################################################
# AUXILIARY FUNCTIONS TO COMPUTE THE SYMMETRY GROUP
######################################################

function cycle_list_to_perm(cycles)
    current_p = GAP.Globals.CycleFromList(GapObj([]))
    for p in cycles
        current_p = current_p*GAP.Globals.CycleFromList(GapObj(p))
    end
    return current_p
end

## Given two matrices we are trying to find a permutation of the columns identifying these two matrices
function find_permutation(matrix, new_matrix)
    perm = []

    for i in 1:size(matrix,2)
        found_matching_column = false
        for j in 1:size(matrix,2)
            if view(matrix,:,i) == view(new_matrix,:,j)
                push!(perm, j)
                found_matching_column = true
                break;
            end
        end
        if found_matching_column == false
            return nothing
        end
    end
    if length(perm) < size(matrix,2)
        return nothing
    end
    return GAP.Globals.PermList(GapObj(perm))
end

## We are trying to find a permutation that recovers the matrix after scaling some row by some scaling factor
function find_permutation_from_scaling(matrix, scaling_factor, row_index::Int64)
    new_matrix = matrix[:,:]
    new_matrix[row_index,:] = matrix[row_index,:]*scaling_factor
    return find_permutation(matrix, new_matrix)
end

function find_all_permutations_from_scaling(matrix, scaling_factors)
    perms = []
    for i in scaling_factors
        for row_index in 1:size(matrix,1)
            new_matrix = matrix[:,:]
            new_matrix[row_index,:] = matrix[row_index,:]*i
            perm = find_permutation(matrix, new_matrix)
            if perm != nothing
                push!(perms, perm)
            end
        end
    end
    return perms
end

## We are trying to find a permutation that recovers the matrix after permuting its rows by row_perm
## We assume row_perm to be an array of length equal to the number of rows of the matrix
function find_permutation_from_permuteted_rows(matrix, row_perm::Array{Int64,1})
    @assert length(unique(row_perm)) == size(matrix,1)
    new_matrix = matrix[row_perm,:]
    return find_permutation(matrix, new_matrix)
end

function find_all_permutations_from_row_permutations(matrix,first_relevant_row=1)
    perms = []
    for p in permutations(Array(first_relevant_row:size(matrix,1)))
        padded_p = p
        if first_relevant_row != 1
            padded_p = vcat(1:first_relevant_row-1,padded_p)
        end
        perm = find_permutation_from_permuteted_rows(matrix, padded_p)
        if perm != nothing
            push!(perms, perm)
        end
    end
    return perms
end

function find_all_permutations(matrix,first_relevant_row=1)
    elts = unique(matrix)
    elts = [i for i in elts if i != 0 && i != 1]
    scaling_permutations = GapObj(find_all_permutations_from_scaling(matrix, [-1]))
    permuted_row_permutations = GapObj(find_all_permutations_from_row_permutations(matrix, first_relevant_row))

    return GAP.Globals.Group(GAP.Globals.Concatenation(scaling_permutations,permuted_row_permutations))
end

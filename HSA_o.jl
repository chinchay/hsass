
# using TimerOutputs # for `to`
## Create a TimerOutput, this is the main type that keeps track of everything.
# const to = TimerOutput()

# using DelimitedFiles # for writedlm()

# include("./customreadlammps.jl")
# using .readLammpsModule
# using PyPlot
#
# using StatsBase
# using Statistics # for `std`
### \notit for `not in` : https://stackoverflow.com/questions/59978282/is-there-an-elegant-way-to-do-not-in-in-julia

# using Random # for `randperm`
using Test # for @test

include("./common.jl")
include("./utils.jl")

function printListAndMean(
    n::Int,
    list::Array{Float64,1},
    )
    for i in 1:n
        print(list[i]); print(", ")
    end
    print(mean(list)); print(", ")
    println(std(list))
end

#-------------------------------------------------------------------------------

# found = is_r_in_2Refill(r, refill, L)
function is_r_in_2Refill(
    r::Int,
    refill::Array{Int64,2},
    L::Int
    )
    #
    found = false
    for k in 1:L
        if r == refill[2, k]
            found = true
        end
    end
    return found
end

# found = is_r_in_il(r, il, L)
function is_r_in_il(
    r::Int,
    il::Array{Int64,1},
    L::Int
    )
    #
    found = false
    for k in 1:L
        if r == il[k]
            found = true
        end
    end
    return found
end

# e0 = getRefillTotalEnergy(L, energyBase, r, U, UionAionBlist)
function getRefillTotalEnergy(
    L::Int,
    energyBase::Float64,
    r::Array{Int64,2},
    U::Array{Float64, 2},
    UionAionBlist::Array{Float64,4}
    )
    #
    Ut   = getEnergy_ion_atoms(L, r, U)
    Uion = getEnergy_ion_ion(L, r, UionAionBlist)
    return energyBase + Ut + Uion
end

# Ut = get_ilEnergy_ion_atoms(L, il, U, ion1, ion2)
function get_ilEnergy_ion_atoms(
    L::Int,
    il::Array{Int64,1},
    U::Array{Float64,2},
    ion1::Int,
    ion2::Int
    )
    L1plus1 = L1 + 1
    s = 0.0
    for i in 1:L1
        ia = il[i]
        if ia != 0
            s += U[ia, ion1]
        end
    end
    for i in L1plus1:L
        ia = il[i]
        if ia != 0
            s += U[ia, ion2]
        end
    end
    return s
end



# Uion = get_ilEnergy_ion_ion(L1, L, il, UionAionBlist, ion1, ion2)
function get_ilEnergy_ion_ion(
    L1::Int,
    L::Int,
    il::Array{Int64,1},
    UionAionBlist::Array{Float64,4},
    ion1::Int,
    ion2::Int
    )
    L1plus1 = L1 + 1
    Uion = 0.0
    for i in 1:L1
        ia = il[i]
        if ia != 0
            for j in i+1:L
                ib = il[j]
                if ib != 0
                    if j <= L1
                        Uion += UionAionBlist[ia, ib, ion1, ion1]
                    else
                        Uion += UionAionBlist[ia, ib, ion1, ion2]
                    end
                end
            end
        end
    end
    for i in L1plus1:L
        ia = il[i]
        if ia != 0
            for j in i+1:L
                ib = il[j]
                if ib != 0
                    if j <= L1
                        Uion += UionAionBlist[ia, ib, ion2, ion1]
                    else
                        Uion += UionAionBlist[ia, ib, ion2, ion2]
                    end
                end
            end
        end
    end
    #
    return Uion
end


# function get_energy_from_il(
#     L1::Int,
#     L::Int,
#     il::Array{Int64,1},
#     ion1::Int,
#     ion2::Int,
#     )
#     Ut   = get_ilEnergy_ion_atoms(L, il, U, ion1, ion2)
#     Uion = get_ilEnergy_ion_ion(L1, L, il, UionAionBlist, ion1, ion2)

# end

# mutate_notRefill!(refill, notRefill, L, Nv)
function mutate_notRefill!(
    refill::Array{Int64,2},
    notRefill::Array{Int64,1},
    L::Int,
    Nv::Int
    )
    #
    j = 0
    for r in 1:Nv
        if !is_r_in_2Refill(r, refill, L)
            j += 1
            notRefill[j] = r
        end
    end
end

# mutate_refill!(refill, L1, L, Nv, ion1, ion2, removedSites)
function mutate_refill!(
        refill::Array{Int64,2},
        L1::Int,
        L::Int,
        Nv::Int,
        ion1::Int,
        ion2::Int,
        removedSites::Array{Int64,1}
    )
    # create refill
    il = sample( 1:Nv, L, replace=false ) # `using StatsBase` is needed
    for (l, r) in enumerate(il)
        refill[ 1, l] = l
        refill[ 2, l] = r
        refill[ 3, l] = removedSites[r]
        refill[ 4, l] = l <= L1 ? ion1 : ion2
    end
end

# copyRefillIntoMem!(L, refill, M_r, memRow)
function copyRefillIntoMem!(
    L::Int,
    refill::Array{Int64,2},
    M_r::Array{Int64,3},
    memRow::Int
    )
    #
    for i in 1:L
        for j in 1:4
            M_r[j, i, memRow] = refill[j, i]
        end
    end
end

# copyNotRefillIntoMem!(Ne, notRefill, M_n)
function copyNotRefillIntoMem!(
    Ne::Int,
    notRefill::Array{Int64,1},
    M_n::Array{Int64,2},
    memRow::Int
    )
    # copy notRefill into M_n and M_n_opt
    for i in 1:Ne
        M_n[i, memRow] = notRefill[i]
    end
end

# # M_r, M_n, M_r_opt, M_n_opt, M_e0, M_dEsum, M_dEsum_opt, M_E, M_acc, worstMemIndex, worstEnergy = getInitMemory(L1, L, Nv, Ne, ion1, ion2, Nmem, maximumMoves, UionAionBlist, removedSites, energyBase)
# function getInitMemory(
#     L1::Int,
#     L::Int,
#     Nv::Int,
#     Ne::Int,
#     ion1::Int,
#     ion2::Int,
#     Nmem::Int,
#     maximumMoves::Int,
#     UionAionBlist::Array{Float64,4},
#     removedSites::Array{Int64,1},
#     energyBase::Float64
#     )
#     #
#     M_r     = zeros(Int64, 4, L, Nmem) # 4×L×nWalkers ::Array{Float64,3}:
#     M_n     = zeros(Int64, Ne, Nmem)
#     M_r_opt = zeros(Int64, 4, L, Nmem)
#     M_n_opt = zeros(Int64, Ne, Nmem)
#     M_e0        = zeros(Nmem)
#     M_dEsum     = zeros(Nmem)
#     M_dEsum_opt = zeros(Nmem)
#     M_E     = zeros(maximumMoves, Nmem) # 0×1 Array{Float64,2}
#     M_acc   = zeros(Bool, maximumMoves, Nmem) # maximumMoves×2 Array{Bool,2}
#     #
#     refill    = zeros(Int, 4, L)
#     notRefill = zeros(Int, Ne)

#     worstMemIndex = -1 # int that can take a value from 1 to Nmem
#     worstEnergy = 0.0

#     #
#     for memRow in 1:Nmem
#         # get new refill
#         mutate_refill!(refill, L1, L, Nv, ion1, ion2, removedSites)

#         #
#         copyRefillIntoMem!(L, refill, M_r, memRow)
#         copyRefillIntoMem!(L, refill, M_r_opt, memRow)

#         # get energy of refill
#         M_e0[m] = getRefillTotalEnergy(L, energyBase, r, U, UionAionBlist)

#         # mutate notRefill from refill
#         mutate_notRefill!(refill, notRefill, L, Nv)

#         # copy notRefill into M_n and M_n_opt
#         copyNotRefillIntoMem!(Ne, notRefill, M_n, memRow)
#         copyNotRefillIntoMem!(Ne, notRefill, M_n_opt, memRow)

#         if isWorstEnergy(memRow, e, worstEnergy)
#             worstMemIndex = memRow
#             worstEnergy   = M_e0[m]
#         end

#     end
#     #
#     return M_r, M_n, M_r_opt, M_n_opt, M_e0, M_dEsum, M_dEsum_opt, M_E, M_acc, worstMemIndex, worstEnergy



# end


# hM_il, hM_U, worstMemIndex, worstEnergy = getInitMemory(L, Nmem, Nv_list, il, ion1, ion2, U, UionAionBlist)
function getInitMemory(
    L::Int,
    Nmem::Int,
    Nv_list::Array{Int64,1},
    il::Array{Int64,1},
    ion1::Int,
    ion2::Int,
    U::Array{Float64,2},
    UionAionBlist::Array{Float64,4}
    )
    # harmony memory
    hM_il = zeros(Int64, Nmem, L) # 4×L×nWalkers ::Array{Float64,3}:
    hM_ev = zeros(Int64, Nmem, Ne) # 4×L×nWalkers ::Array{Float64,3}:
    hM_U  = zeros(Float64, Nmem)
    worstMemIndex = -1 # int that can take a value from 1 to Nmem
    worstEnergy = 0.0

    bestMemIndex = -1 # int that can take a value from 1 to Nmem
    bestEnergy  = 0.0
    #
    for memRow in 1:Nmem
        get_il_proposal!(L, Nv_list, il)
        copy_il_to_hM!(L, il, hM_il, memRow)
        #
        # determine energy of il
        u = get_il_Ut_plus_Uion(L1, L, il, U, UionAionBlist, ion1, ion2, energyBase)

        hM_U[memRow] = u
        if isWorstEnergy(memRow, u, worstEnergy)
            worstMemIndex = memRow
            worstEnergy   = hM_U[memRow]
        end

        if isBestEnergy(memRow, u, bestEnergy)
            bestMemIndex = memRow
        end
    end
    #
    return hM_il, hM_U, worstMemIndex, worstEnergy, bestMemIndex
end


function check_if_taken()

end

# # zerar!(L, vec)
# function zerar!(L::Int, vec::Array{Int,1})
#     for i in 1:L
#         vec[i] = 0
#     end
# end


# get_il_proposal!(L, Nv_list, il)
function get_il_proposal!(
    L::Int,
    Nv_list::Array{Int64,1},
    il::Array{Int64,1}
    )
    shuffle!(Nv_list)
    for i in 1:L
        il[i] = Nv_list[i]
    end
end



# notRefill = get_notRefillFromIl(il, L, Ne, Nv)
function get_notRefill(refill::Array{Int64,2}, L::Int, Ne::Int, Nv::Int)
    notRefill = zeros(Int, Ne)
    j = 0
    # found = false
    for r in 1:Nv
        if !is_r_in_il(r, il, L)
            j += 1
            notRefill[j] = r
        end
    end
    return notRefill
end


# copy_il_to_hM!(L, il, hM_il, memRow)
function copy_il_to_hM!(
    L::Int,
    il::Array{Int64,1},
    hM_il::Array{Int64,2},
    memRow::Int
    )
    for i in 1:L
        hM_il[memRow, i] = il[i]
    end
end


# b = isBestEnergy(memRow, e, bestEnergy)
function isBestEnergy(
    memRow::Int,
    e::Float64,
    bestEnergy::Float64
    )
    return (memRow > 1) ? (e < bestEnergy) : true
end


# b = isWorstEnergy(memRow, e, worstEnergy)
function isWorstEnergy(
    memRow::Int,
    e::Float64,
    worstEnergy::Float64
    )
    return (memRow > 1) ? (e > worstEnergy) : true
    # if memRow > 1
    #     return (e > worstEnergy)
    # else # memRow=1
    #     return true
    # end
end


# ia, randRow = getRandElementFromMemory!(Nmem, col, list, hM_il) ##mutates temporalList
function getRandElementFromMemory!(
    Nmem::Int,
    col::Int,
    list::Array{Int64,1},
    hM_il::Array{Int64,2}
    )
    for row in 1:Nmem
        list[row] = hM_il[row, col]
    end

    randRow = rand(1:Nmem)
    return list[randRow], randRow
end

# x_wrapped = wrap(x, Nv)
function wrap( x::Int, Nv::Int )
    if x > Nv
        return x - Nv
    elseif x < 1
        return x + Nv
    else
        return x
    end
end

# a1, a2 = getPitchedValues(ia, delta, Nv)
function getPitchedValues(
    ia::Int,
    delta::Int,
    Nv::Int
    )
    #
    if rand() < 0.5
        ia1 = wrap( ia - delta, Nv)
        ia2 = wrap( ia + delta, Nv)
        # println("---", il[col])
    else
        ia1 = wrap( ia + delta, Nv)
        ia2 = wrap( ia - delta, Nv)
        # println("+++", il[col])
    end
    return ia1, ia2
end

# ia1, ia2 = pitchAdjust(ia, par, delta)
function pitchAdjust(
    ia::Int,
    par::Float64,
    delta::Int
    )
    #
    # # probability: HMCR * (1 - PAR)
    # ia, randRow = getRandElementFromMemory!(Nmem, col, lTemp, hM_il)


    ia1 = ia
    ia2 = ia
    # STAGE 2
    if rand() <= par # probability: HMCR * PAR
        a1, a2 = getPitchedValues(ia, delta, Nv)
    end
    # println(ia, "  ", ia1, "  ", ia2)
    # return ia
    return ia1, ia2
end

# ia = pitchAdjust_ensureUnrepeatedElement(ia, Nv, par, il)
function pitchAdjust_ensureUnrepeatedElement(
    ia::Int,
    Nv::Int,
    par::Float64,
    il::Array{Int,1}
    )
    #
    repeat = true
    delta = 0

    ia1, ia2 = pitchAdjust(ia, par, delta)

    while repeat && (delta < Nv)
        delta += 1
        ia1, ia2 = pitchAdjust(ia, par, delta)

        ia = ia1
        # println("delta, ia, il[col]: ", delta, "  ", ia, "  ", il[col], "  ", l[col])
        repeat = ia in il
        if repeat
            if ia2 != ia1
                ia = ia2
                # println("delta, ia, il[col]: ", delta, "  ", ia, "  ", il[col], "  ", l[col])
                repeat = ia in il
            end
        end
    end
    return ia
end


# getShuffled_MemColumn!(Nmem, col, memList, hM_il) # mutates memList
function getShuffled_MemColumn!(
    Nmem::Int,
    col::Int,
    memList::Array{Int,1},
    hM_il::Array{Int,2},
    )
    #
    for row in 1:Nmem
        memList[row] = hM_il[row, col]
    end
    shuffle!(memList)
end

# get_MemColumn!(Nmem, col, memList, hM_il) # mutates memList
function get_MemColumn!(
    Nmem::Int,
    col::Int,
    memList::Array{Int,1},
    hM_il::Array{Int,2},
    )
    #
    for row in 1:Nmem
        memList[row] = hM_il[row, col]
    end
end

# newHarmony!(L, Nmem, hcmr, par, il, Nv_list, memoryList, hM_il) # mutates il
function newHarmony!(
    L::Int,
    Nv::Int,
    Nmem::Int,
    hcmr::Float64,
    par::Float64,
    il::Array{Int,1},
    Nv_list::Array{Int,1},
    memoryList::Array{Int,1},
    hM_il::Array{Int,2},
    )
    #
    # l = zeros(Int, 52)
    # for i in 1:L
    #     l[i] = il_[i]
    # end

    # STAGE 1
    # get_il_proposal!(L, Nv_list, il_) # probability: (1 - HMCR), if doing the following stages:
    # println(il_)

    # for i in 1:L
    #     il_[i] = hM_il[1, i]
    # end

    # # initialize lTemp2
    # for i in 1:L
    #     lTemp2[i] = -1
    # end


    boolNotRefilled = ones(Bool, Nv) # = [true, true, true, ...]
    shuffle!(Nv_list)


    for col in 1:L
        a_final = -1
        if rand() <= hcmr
            # memList will be a column of hM_il, shuffled
            getShuffled_MemColumn!(Nmem, col, memoryList, hM_il) # mutates memList

            pitch = (rand() <= par)

            if !pitch
                repeat = true
                count = 0
                while repeat && (count < Nmem)
                    count += 1
                    a = memoryList[count]
                    if boolNotRefilled[a]
                        a_final = a
                        repeat = false
                    end
                end
                pitch = repeat
            end
            #
            if pitch
                a = memoryList[1]
                repeat = true
                delta = 0
                while repeat && (delta < Nv)
                    delta += 1
                    a1, a2 = getPitchedValues(a, delta, Nv)
                    if boolNotRefilled[a1]
                        a_final = a1
                        repeat = false
                    elseif boolNotRefilled[a2]
                        a_final = a2
                        repeat = false
                    end
                end
            end
        else
            count = 0
            repeat = true
            while repeat && (count < Nv)
                count += 1
                a = Nv_list[count]
                if boolNotRefilled[a]
                    a_final = a
                    repeat = false
                end
            end
        end
        il[col] = a_final
        boolNotRefilled[a_final] = false
    end






    # col = rand(1:L)

    # # for col in 1:L
    #     # STAGE 1
    #     if rand() <= hcmr

    #         count = 0
    #         while repeat && (count < Nv)
    #             count += 1

    #             # probability: HMCR * (1 - PAR)
    #             ia, randRow = getRandElementFromMemory!(Nmem, col, lTemp, hM_il)

    #             # make sure the element from memory is not present in the current list
    #             repeat = ia in lTemp2
    #             if !repeat
    #                 lTemp2[col] = ia
    #         end
    #         if repeat

    #         # # pitchAdjust
    #         # ia = pitchAdjust_ensureUnrepeatedElement(ia, Nv, par, il_)

    #         #
    #         il_[col] = ia

    #     end
    # # end
    # # println(il_)

        # STAGE 1
        ## the code above was simplified from:
        # if rand() > hcmr
        #     # "the decision variable of the new harmony is randomly generated"
        #     # already done when shuffling: `get_il_proposal!(L, Nv_list, il)`
        #     # saved into il from get_il_proposal!(L, Nv_list, il)
        #     # do nothing else
        # else
        #     # "one of the harmonies storred in HM is randomly selected (see eq 5)"
        #     # now il[j] is going to be changed:
        #     ia = getRandElementFromMemory!(Nmem, col, temporalList, hM_il)

        #     # # if ia is already present in `il`, then just grab the
        #     # # proposal il[i], i.e., leave it as it was, do nothing.

        #     # # otherwise il[i] = ia
        #     # il[i] = ia

        #     # Actually, don't worry if `ia` is already present in `il`, because
        #     # the proposal is going to have a very high energy (modify
        #     # interaction energy if required) so that the configuration will be
        #     # replaced with another configuration eventually.
        #     il[i] = ia
        # end


end

# u = getInteractionWithProposalSite(L1, L, siteProposal, ionTypeProposal, il, UionAionBlist)
function getInteractionWithProposalSite(
	L1::Int,
	L::Int,
    siteProposal::Int,
    ionTypeProposal::Int,
	il::Array{Int,1},
    UionAionBlist::Array{Float64,4}
    )
    #
    u = 0.0
	for l in 1:L
		siteX  = il[l]
		ionX = getIonType(l, L1)
		if siteX != siteProposal
			if siteX != 0
				u += UionAionBlist[siteProposal, siteX, ionTypeProposal, ionX]
			end
		end
	end
    return u
end

# du = get_dE_byFillingAnEmptySite(L1, L, col, ia, il, U, UionAionBlist)
function get_dE_byFillingAnEmptySite(
	L1::Int,
	L::Int,
	i_empty::Int,
	ia_site::Int,
	il::Array{Int,1},
    U::Array{Float64,2},
    UionAionBlist::Array{Float64,4}
	)
	#
	ionA = (i_empty <= L1) ? ion1 : ion2
	dUt = U[ia_site, ionA]
	#

	UiB = 0.0
	for l in 1:L
		ilk  = il[l]
		ionX = (l <= L1) ? ion1 : ion2
		if ilk != ia_site
			if ilk != 0
				UiB += UionAionBlist[ia_site, ilk, ionA, ionX]
			end
		end
	end
	dUion = UiB
	dUdistr = dUt + dUion
	return dUdistr
end

# ionType = getIonType(i, L1)
function getIonType(
    i::Int,
    L1::Int,
    )
    #
    return (i <= L1) ? 1 : 2 #ion1=1, ion2=2
end


# repair_and_add!(l_init, l_end, ionType, L1, L, il, vacs, du2List, U, UionAionBlist) # mutates il
function repair_and_add!(
    l_init::Int,
    l_end::Int,
    ionType::Int,
    L1::Int,
    L::Int,
    il::Array{Int,1},
    vacs::Array{Bool,1},
    du2List::Array{Float64,1},
    U::Array{Float64,2},
    UionAionBlist::Array{Float64,4},
    )
    #
    # initialize lists
    for site in 1:Nv
        vacs[site] = true
        du2List[site] = 0.0
    end

    # `vacs` is a list of sites that have not been used yet (trues) to fill `il`
    for site in il
        if site != 0
            vacs[site] = false #it's occupied
        end
    end

    # initializing variables to be modified later
    m = 0 # a counter
    site0 = 0 # a site
    # ionType0 = 0 # an ionType=ion1 or ion2
    #
    for l in l_init:l_end
        if il[l] == 0
            # updating variables
            m += 1
            # ionType = getIonType(l, L1)

            # initializing variables to be modifiend later
            site = 0
            cont = 0
            siteMin = 0
            duMin = 0

            for bool in vacs
                site += 1
                if bool
                    cont += 1

                    # initialize siteMin
                    if cont == 1
                        siteMin = site
                    end

                    # use previous calculations to avoid redundant work
                    if m > 1
                        # du2 = du2List[site] + UionAionBlist[site, site0, ionType, ionType0]
                        du2 = du2List[site] + UionAionBlist[site, site0, ionType, ionType]
                    else
                        du2 = getInteractionWithProposalSite(L1, L, site, ionType, il, UionAionBlist)
                    end

                    # save the the previous calculation for future use
                    du2List[site] = du2

                    # get the total energy change if `site` was introduced at `il`
                    du = U[site, ionType] + du2

                    # save the best site that minimizes the energy the most
                    if du < duMin
                        duMin = du
                        siteMin = site
                    end
                end
            end
            il[l] = siteMin
            vacs[siteMin] = false
            site0 = siteMin
            # ionType0 = ionType

        end
    end
    #

end

# printIl(L, il)
function printIl(L::Int, il::Array{Int,1})
    println("")
    for i in 1:L
        print(il[i])
        print(" ")
    end
    println("")
end

# getBestHarmony!(Nmem, bestRow, bestHarmony, hM_il) # mutates bestHarmony
function get_listAvoidBest!(
    Nmem::Int,
    col::Int,
    bestRow::Int,
    listAvoidBest::Array{Int,1}, #listAvoidBest has length = L-1
    hM_il::Array{Int,2},
    )
    #
    cont = 0
    for row in 1:Nmem
        if row != bestRow
            cont += 1
            listAvoidBest[cont] = hM_il[row, col]
        end
    end
end

# initializeEmpties!( emptiesDict ) # mutates emptiesDict
function initializeEmpties!(
    emptiesDict::Dict,
    )
    #
    for (key, value) in emptiesDict
        emptiesDict[key] = true
    end
end

# site = getEmptySite(Nv_list, emptiesDict)
function getEmptySite(
    Nv_list::Array{Int,1},
    emptiesDict::Dict
    )
    #
    for site in Nv_list
        if emptiesDict[site]
            return site
        end
    end
end

# getBestHarmony!(L, bestRow, bestHarmony, hM_il) # mutates bestHarmony 
function getBestHarmony!(
    L::Int,
    bestRow::Int,
    bestHarmony::Array{Int,1},
    hM_il::Array{Int,2},
    ) 
    #
    for col in 1:L
        bestHarmony[col] = hM_il[bestRow, col]
    end
end


# It will mutate `il`
# newHarmony_algo3_dghsa!(L, Nmem, hcmr, par, il, hM_il, bestRow, bestHarmony, listAvoidBest, emptiesDict)
function newHarmony_algo3_dghsa!(
    L::Int,
    Nmem::Int,
    hcmr::Float64,
    par::Float64,
    il::Array{Int,1},
    hM_il::Array{Int,2},
    bestRow::Int,
    bestHarmony::Array{Int,1},
    listAvoidBest::Array{Int,1},
    emptiesDict::Dict,
    Nv_list::Array{Int,1},
    )
    #
    # pick the best harmony
    getBestHarmony!(L, bestRow, bestHarmony, hM_il) # mutates bestHarmony

    initializeEmpties!( emptiesDict )
    
    # Nv_list will be used as random keys for emptiesDict
    shuffle!(Nv_list)

    for col in 1:L
        a = -1
        if rand() <= hcmr
            # memory consideration
            a = bestHarmony[col]            
            if a in il
                a = 0
            end
        else
            # Generate a random integer number 𝑎 ∈ [1..HMS] ∧ 𝑎 ≠ bestRow
            # listAvoidBest will be a column of hM_il, avoiding best harmony
            get_listAvoidBest!(Nmem, col, bestRow, listAvoidBest, hM_il) # mutates listAvoidBest
            a = rand(listAvoidBest)

            # pitch adjustment for the pitch chosen randomly
            if rand() <= par
                a = getEmptySite(Nv_list, emptiesDict)
            end

            if a in il
                a = 0
            end

        end

        # If `a=0`, we would need to repair.
        # Instead of repairing ( by using repair_and_add() ) we would just take
        # an available site:
        if a == 0
            a = getEmptySite(Nv_list, emptiesDict)
        end

        # now insert into `il`:
        il[col] = a

        # update empty sites dictionary
        emptiesDict[a] = false
    end



end

# newHarmony_dghsa!( L, Nv, Nmem, hcmr, par, il, memList, hM_il, U, UionAionBlist, ion1, ion2, energyBase, vacs, du2List )
function newHarmony_dghsa!(
    L1::Int,
    L::Int,
    Nv::Int,
    Nmem::Int,
    hcmr::Float64,
    par::Float64,
    il::Array{Int,1},
    memList::Array{Int,1},
    hM_il::Array{Int,2},
    U::Array{Float64,2},
    UionAionBlist::Array{Float64,4},
    ion1::Int,
    ion2::Int,
    energyBase::Float64,
    vacs::Array{Bool,1},
    du2List::Array{Float64,1}
    )
    #
    for col in 1:L
        a = -1
        if rand() <= hcmr
            # memList will be a column of hM_il, shuffled
            get_MemColumn!(Nmem, col, memList, hM_il) # mutates memList
            a = rand(memList)

            if rand() <= par
                if rand() < 0.5
                # if rand(Bool)
                    a -= 1
                else
                    a += 1
                end
                a = wrap(a, Nv)
            end
            if a in il
                a = 0
            end
        else
            a = rand(1:Nv)
            if a in il
                a = 0
            end
        end
        il[col] = a
    end
    #
    # # printIl(L, il)
    # il0 = copy(il)

    # for col in 1:L
    #     if il[col] == 0
    #         du_min = 0
    #         ia_min = 0
    #         cont = 0
    #         for ia in 1:Nv
    #             if ia ∉ il
    #                 cont += 1
    #                 if cont == 1
    #                     ia_min  = ia
    #                 end
    #                 # il[col] = ia
    #                 # uf = get_il_Ut_plus_Uion(L1, L, il, U, UionAionBlist, ion1, ion2, energyBase)
    #                 # du = uf - u
    #                 du = get_dE_byFillingAnEmptySite(L1, L, col, ia, il, U, UionAionBlist)
    #                 if du < du_min
    #                     du_min = du
    #                     ia_min = ia
    #                 end
    #             end
    #         end

    #         # now put ia_min in col
    #         il[col] = ia_min
    #     end
    # end

    # # printIl(L, il)
    # for i in 1:L
    #     il[i] = il0[i]
    # end

    ############################################################################
    ionType = ion1
    l_init = 1
    l_end = L1
    repair_and_add!(l_init, l_end, ionType, L1, L, il, vacs, du2List, U, UionAionBlist) # mutates il

    ############################################################################
    ionType = ion2
    l_init = L1 + 1
    l_end = L
    repair_and_add!(l_init, l_end, ionType, L1, L, il, vacs, du2List, U, UionAionBlist) # mutates il

    ############################################################################

    # printIl(L, il)
    # println("............................")

end

# e = get_il_Ut_plus_Uion(L1, L, il, U, UionAionBlist, ion1, ion2)
function get_il_Ut_plus_Uion(
    L1::Int,
    L::Int,
    il::Array{Int64,1},
    U::Array{Float64,2},
    UionAionBlist::Array{Float64,4},
    ion1::Int,
    ion2::Int,
    energyBase::Float64
    )
    # determine energy of il
    Ut   = get_ilEnergy_ion_atoms(L, il, U, ion1, ion2)
    Uion = get_ilEnergy_ion_ion(L1, L, il, UionAionBlist, ion1, ion2)
    return Ut + Uion + energyBase
end

# loop!(L, Nmem, Nv_list, il, ion1, ion2, U, UionAionBlist, hcmr, par, temporalList) # mutates bestEnergyMem, worstEnergyMem
function loop!(
    L::Int,
    Nmem::Int,
    Nv_list::Array{Int,1},
    il::Array{Int,1},
    ion1::Int,
    ion2::Int,
    U::Array{Float64,2},
    UionAionBlist::Array{Float64,4},
    hcmr::Float64,
    par::Float64,
    tempListNmem::Array{Int,1},
    tempListL::Array{Int,1},
    nIterations::Int,
    energyBase::Float64,
    bestEnergyMem::Array{Float64,1},
    worstEnergyMem::Array{Float64,1},
    )
    #

    hM_il, hM_U, worstRow, worstEnergy = getInitMemory(L, Nmem, Nv_list, il, ion1, ion2, U, UionAionBlist)
    # M_r, M_n, M_r_opt, M_n_opt, M_e0, M_dEsum, M_dEsum_opt, M_E, M_acc, worstMemIndex, worstEnergy = getInitMemory(L1, L, Nv, Ne, ion1, ion2, Nmem, maximumMoves, UionAionBlist, removedSites, energyBase)

    for i in 1:nIterations

        # `il` will be modified:
        # newHarmony!(L, il, hM_il, hcmr, par, Nmem, Nv_list, tempListNmem)
        newHarmony!(L, Nmem, hcmr, par, il, Nv_list, tempListNmem, hM_il) # mutates il
        # println(hM_il)

        # for jj in 1:L
        #     tempListL[jj] = hM_il[1, jj]
        # end
        # print(tempListL)
        # println("....")
        # println(il)

        # determine energy of il
        u = get_il_Ut_plus_Uion(L1, L, il, U, UionAionBlist, ion1, ion2, energyBase)

        # # println(hM_U)
        # println("i, worstRow, worstEnergy: ", i, "  ", worstRow, "  ", worstEnergy, " ", u)
        # println("")


        # if il's energy is best than the worst row in memory, replace it
        if u < worstEnergy
            copy_il_to_hM!(L, il, hM_il, worstRow)
            hM_U[worstRow] = u
        end

        # determine again the new worstRow and worstEnergy
        worstRow = argmax(hM_U)
        worstEnergy = hM_U[worstRow]

        bestEnergyMem[i] = hM_U[ argmin(hM_U) ]
        worstEnergyMem[i] = worstEnergy

    end

    # # println(hM_il)
    # for i in 1:L
    #     il[i] = hM_il[1,i]
    # end
    # u = get_il_Ut_plus_Uion(L1, L, il, U, UionAionBlist, ion1, ion2, energyBase)
    # println(u)

    println(worstEnergyMem[end])
    println(bestEnergyMem[end])
end

# loop_par!(L, Nmem, Nv_list, il, ion1, ion2, U, UionAionBlist, hcmr, tempListNmem, energyBase, bestEnergyMem, worstEnergyMem, lpar) # mutates bestEnergyMem, worstEnergyMem
function loop_par!(
    L1::Int,
    L::Int,
    Nv::Int,
    Nmem::Int,
    Nv_list::Array{Int,1},
    il::Array{Int,1},
    ion1::Int,
    ion2::Int,
    U::Array{Float64,2},
    UionAionBlist::Array{Float64,4},
    tempListNmem::Array{Int,1},
    nIterations::Int,
    energyBase::Float64,
    bestEnergyMem::Array{Float64,1},
    worstEnergyMem::Array{Float64,1},
    acceptance_hist::Array{Bool,1},
    vacs::Array{Bool,1},
    du2List::Array{Float64,1},
    bestHarmony::Array{Int,1},
    listAvoidBest::Array{Int,1},
    emptiesDict::Dict,
    parValues::Array{Float64,1},
    hcmrValues::Array{Float64,1},
    )
    #

    hM_il, hM_U, worstRow, worstEnergy, bestRow = getInitMemory(L, Nmem, Nv_list, il, ion1, ion2, U, UionAionBlist)
    # M_r, M_n, M_r_opt, M_n_opt, M_e0, M_dEsum, M_dEsum_opt, M_E, M_acc, worstMemIndex, worstEnergy = getInitMemory(L1, L, Nv, Ne, ion1, ion2, Nmem, maximumMoves, UionAionBlist, removedSites, energyBase)

    for i in 1:nIterations
        par  = parValues[i]
        hcmr = hcmrValues[i]

        # `il` will be modified:
        # newHarmony!(L, il, hM_il, hcmr, par, Nmem, Nv_list, tempListNmem)

        # I commented this line to replace it with newHarmony_dghsa!()
        #newHarmony!(L, Nv, Nmem, hcmr, par, il, Nv_list, tempListNmem, hM_il) # mutates il

        # propose new harmony with a "greedy" consideration
        # newHarmony_dghsa!( L1, L, Nv, Nmem, hcmr, par, il, tempListNmem, hM_il, U, UionAionBlist, ion1, ion2, energyBase )
        # newHarmony_dghsa!( L1, L, Nv, Nmem, hcmr, par, il, tempListNmem, hM_il, U, UionAionBlist, ion1, ion2, energyBase, vacs, du2List )

        newHarmony_algo3_dghsa!(L, Nmem, hcmr, par, il, hM_il, bestRow, bestHarmony, listAvoidBest, emptiesDict, Nv_list)

        # println(hM_il)

        # for jj in 1:L
        #     tempListL[jj] = hM_il[1, jj]
        # end
        # print(tempListL)
        # println("....")
        # println(il)

        # determine energy of il
        u = get_il_Ut_plus_Uion(L1, L, il, U, UionAionBlist, ion1, ion2, energyBase)


        # u2 = getEnergyLammps(L1, L, il)
        # if abs(u - u2) > 0.1
        #     throw(error())
        #     println(u, " /// ", u2)
        # # else
        # #     println(u, " ... ", u2)
        # end

        # # println(hM_U)
        # println("i, worstRow, worstEnergy: ", i, "  ", worstRow, "  ", worstEnergy, " ", u)
        # println("")


        # if il's energy is best than the worst row in memory, replace it
        if u < worstEnergy
            copy_il_to_hM!(L, il, hM_il, worstRow)
            hM_U[worstRow] = u
            acceptance_hist[i] = true
        else
            acceptance_hist[i] = false
        end

        # determine again the new worstRow and worstEnergy
        worstRow = argmax(hM_U)
        worstEnergy = hM_U[worstRow]

        bestEnergyMem[i] = hM_U[ argmin(hM_U) ]
        worstEnergyMem[i] = worstEnergy


        bestRow = argmin(hM_U) #  newHarmony_algo3_dghsa! will use bestRow

    end

    # # println(hM_il)
    # for i in 1:L
    #     il[i] = hM_il[1,i]
    # end
    # u = get_il_Ut_plus_Uion(L1, L, il, U, UionAionBlist, ion1, ion2, energyBase)
    # println(u)

    println(worstEnergyMem[end])
    println(bestEnergyMem[end])
    
    # println(hM_il)
    # println("")
    # println(hM_U)
end

# plot1(namePlot1, nIterations, bestE_hist, worstE_hist, par_hist)
function plot1(
    namePlot1::String,
    nIterations::Int,
    bestE_hist::Array{Float64,1},
    worstE_hist::Array{Float64,1},
    par_hist::Array{Float64,1}
    )
    #
    x = collect(1:nIterations)
    #
    PyPlot.close() # it gave me errors saying that other plots were open.
    myColors = ["r", "k", "g", "m", "y"]
    fig, ax1 = PyPlot.subplots()
    ax1.scatter(x, bestE_hist, s=0.5, color=myColors[1])
    ax1.scatter(x, worstE_hist, s=0.5, color=myColors[2])
    # xlabel = string("move x ", factNumPoints)
    ax1.set_xlabel("iteration")
    ax1.set_ylabel("energy (eV)", color="k")
    # ax1.set_ylim([-5600, -5150])
    ax1.set_ylim([-43000, -40000])
    #
    ax2 = ax1.twinx()
    ax2.plot(x, par_hist, "--b", linewidth=1.4)
    # plt.errorbar(x, yTimeMean, yerr=yTimeDev, color="b")
    # ax2.scatter(x, yTimeMean, s=0.5, color="b")
    ax2.set_ylabel("PAR", color="b")
    ax2.set_ylim([0.0, 1.05])
    ax2.tick_params(axis="y", colors="blue")
    #
    # fig.savefig("hsa_")
    fig.savefig(namePlot1)
    PyPlot.close(fig)
end

# plot2(x_vals, y1mean, y1dev, y2mean, y2dev, y3mean, y3dev, y3max, xlabel)
function plot2(
    x_vals::Array{Float64,1},
    y1mean::Array{Float64,1},
    y1dev::Array{Float64,1},
    y2mean::Array{Float64,1},
    y2dev::Array{Float64,1},
    y3mean::Array{Float64,1},
    y3dev::Array{Float64,1},
    y1max::Float64,
    y3max::Float64,
    xlabel::String,
    )
    #
    PyPlot.close() # it gave me errors saying that other plots were open.
    myColors = ["r", "k", "g", "m", "y"]
    fig, ax1 = PyPlot.subplots()
    # ax1.scatter(x, yWorstMean, s=0.5, color=myColors[1])
    # ax1.plot(x, yWorstMean, "--o", color=myColors[1])
    # ax1.errorbar(x, yWorstMean, yerr=yWorstDev, "--o", color=myColors[1])

    # ax1.scatter(x, yBestMean, s=0.5, color=myColors[2])
    # ax1.plot(x, yBestMean, "--o", color=myColors[2])
    plt.errorbar(x_vals, y2mean, yerr=y2dev, color=myColors[2])
    plt.errorbar(x_vals, y1mean, yerr=y1dev, color=myColors[1])
    # xlabel = string("move x ", factNumPoints)

    # ax1.set_xlabel("nIterations")
    ax1.set_xlabel(xlabel)

    ax1.set_ylabel("mean energy (eV)", color="k")
    ax1.set_ylim([0, y1max])
    fig.savefig("hsa_evolution")
    #
    ax2 = ax1.twinx()
    # ax2.plot(x_vals, y3mean, "--ob", linewidth=0.4)
    plt.errorbar(x_vals, y3mean, yerr=y3dev, color="b")
    # ax2.scatter(x, yTimeMean, s=0.5, color="b")
    ax2.set_ylabel("time (seconds)", color="b")
    ax2.set_ylim([0.0, y3max])
    ax2.tick_params(axis="y", colors="blue")


    fig.savefig("hsa_evolution")
    #
    PyPlot.close(fig)

end


# experiment10(nruns, L, Nmem, Nv_list, il, ion1, ion2, U, UionAionBlist, hcmr, tempListNmem, nIterations, energyBase, x_vals, yWorstMean, yWorstDev, yBestMean, yBestDev, yTimeMean, yTimeDev, bestE_hist, worstE_hist, par_hist)
function experiment10!(
    nruns::Int,
    L::Int,
    Nmem::Int,
    Nv_list::Array{Int,1},
    il::Array{Int,1},
    ion1::Int,
    ion2::Int,
    U::Array{Float64,2},
    UionAionBlist::Array{Float64,4},
    hcmr::Float64,
    tempListNmem::Array{Int,1},
    nIterations::Int,
    energyBase::Float64,
    x_vals::Array{Float64,1},
    yBestMean::Array{Float64,1},
    yBestDev::Array{Float64,1},
    yWorstMean::Array{Float64,1},
    yWorstDev::Array{Float64,1},
    yTimeMean::Array{Float64,1},
    yTimeDev::Array{Float64,1},
    bestE_hist::Array{Float64,1},
    worstE_hist::Array{Float64,1},
    par_hist::Array{Float64,1},
    )
    #
    # nruns = 20
    #
    if nruns > 1
        for (k, par) in enumerate(x_vals)
            println("k, par: ", k, " ", par)
            #
            results = zeros(5, nruns )
            for i in 1:nruns
                # println("i: ", i)
                # time = @elapsed loop_par!(L, Nmem, Nv_list, il, ion1, ion2, U, UionAionBlist, hcmr, tempListNmem, nIterations, energyBase, bestE_hist, worstE_hist, par_hist)
                time = @elapsed loop!(L, Nmem, Nv_list, il, ion1, ion2, U, UionAionBlist, hcmr, par, tempListNmem, tempListL, nIterations, energyBase, bestE_hist, worstE_hist)
                results[1, i] = bestE_hist[end]
                results[2, i] = worstE_hist[end]
                results[3, i] = results[1, i] + 5591.3
                results[4, i] = results[2, i] + 5591.3
                results[5, i] = time
            end
            #
            yBestMean[k]  = mean( results[3, 1:5] )
            yBestDev[k]   = std(  results[3, 1:5] )

            yWorstMean[k] = mean( results[4, 1:5] )
            yWorstDev[k]  = std(  results[4, 1:5] )

            yTimeMean[k]  = mean( results[5, 1:5] )
            yTimeDev[k]   = std(  results[5, 1:5] )

            println("k, par, WorstMean[k]: ", k, " ", par, " ", yWorstMean[k], " ", yWorstDev[k])
        end
    else

        time = @elapsed loop_par!(L, Nmem, Nv_list, il, ion1, ion2, U, UionAionBlist, hcmr, tempListNmem, nIterations, energyBase, bestE_hist, worstE_hist, par_hist)
    end
end


# mean_lastNenergies_relative, std_lastNenergies_relative, mean_lastNacceptances, std_lastNacceptances = getData(lastN, listE, listAccepted)
function getData(
    lastN::Int,
    listE::Array{Float64,1},
    listAccepted::Array{Bool,1}
    )
    #
    # lastNaccepted_energies = getLastAcceptedNpoints(lastN, listE, listAccepted)

    lastNenergies = listE[end - lastN + 1: end]
    mean_lastNenergies = mean(lastNenergies)
    std_lastNenergies  = std(lastNenergies)

    # differences between the ground state energy and each of 10 last energies
    lastNenergies .-= -5591.3
    mean_lastNenergies_relative = mean(lastNenergies)
    std_lastNenergies_relative  = std(lastNenergies)

    # same thing for acceptances
    lastNacceptances = listAccepted[end - lastN + 1: end]
    mean_lastNacceptances = std(lastNacceptances)
    std_lastNacceptances = mean(lastNacceptances)

    return mean_lastNenergies, std_lastNenergies, mean_lastNenergies_relative, std_lastNenergies_relative, mean_lastNacceptances, std_lastNacceptances
end

# vals = getDynamicVals(Xconstant, is_X_dynamic, Xmin, Xmax, nIterations)
function getDynamicVals(
    Xconstant::Float64,
    is_X_dynamic::Bool,
    Xmin::Float64,
    Xmax::Float64,
    nIterations::Int
    )
    #
    vals = zeros(nIterations)
    
    if is_X_dynamic
        # Xmin = 0.3
        # Xmax = 0.0 #1.0
        dX = (Xmax - Xmin) / nIterations

        for i in 1:nIterations
            X = Xmin + (i * dX)
            vals[i] = X
        end
    else
        for i in 1:nIterations
            vals[i] = Xconstant
        end
    end
    return vals
end


# experiment11!( nruns, L, Nv, Nv_list, il, ion1, ion2, U, UionAionBlist, nIterations, energyBase, x_vals, yBestMean, yBestDev, yWorstMean, yWorstDev, yTimeMean, yTimeDev, bestE_hist, worstE_hist, par_hist, bestHarmony, emptiesDict, par, is_PAR_dynamic, parMin, parMax, hcmr, is_HCMR_dynamic, hcmrMin, hcmrMax )
function experiment11!(
    nruns::Int,
    L::Int,
    Nv::Int,
    Nv_list::Array{Int,1},
    il::Array{Int,1},
    ion1::Int,
    ion2::Int,
    U::Array{Float64,2},
    UionAionBlist::Array{Float64,4},
    nIterations::Int,
    energyBase::Float64,
    x_vals::Array{Int,1},
    yBestMean::Array{Float64,1},
    yBestDev::Array{Float64,1},
    yWorstMean::Array{Float64,1},
    yWorstDev::Array{Float64,1},
    yTimeMean::Array{Float64,1},
    yTimeDev::Array{Float64,1},
    bestE_hist::Array{Float64,1},
    worstE_hist::Array{Float64,1},
    par_hist::Array{Float64,1},
    bestHarmony::Array{Int,1},
    emptiesDict::Dict,
    par::Float64,
    is_PAR_dynamic::Bool,
    parMin::Float64,
    parMax::Float64,
    hcmr::Float64,
    is_HCMR_dynamic::Bool,
    hcmrMin::Float64,
    hcmrMax::Float64,
    )
    #
    # nruns = 20
    #

    println(x_vals)

    lastN = 10
    acceptance_hist = zeros(Bool, nIterations)
    l_mean_lastNenergies = zeros(nruns)
    l_std_lastNenergies = zeros(nruns)
    l_mean_lastNenergies_relative = zeros(nruns)
    l_std_lastNenergies_relative = zeros(nruns)
    l_mean_lastNacceptances = zeros(nruns)
    l_std_lastNacceptances = zeros(nruns)
    l_time = zeros(nruns)
    l_lastE = zeros(nruns)
    l_lastAcc = zeros(nruns)

    vacs = ones(Bool, Nv)
    du2List = zeros(Nv)

    parValues  = getDynamicVals(par,  is_PAR_dynamic,  parMin,  parMax,  nIterations)
    hcmrValues = getDynamicVals(hcmr, is_HCMR_dynamic, hcmrMin, hcmrMax, nIterations)

    for (k, Nmem) in enumerate(x_vals)
        println("k, Nmem: ", k, " ", Nmem)
        #
        tempListNmem = zeros(Int, Nmem)
        results = zeros(5, nruns )

        listAvoidBest = zeros(Int, Nmem - 1)
        println("length(listAvoidBest): ", length(listAvoidBest))

        if nruns > 1
            for i in 1:nruns
                # println("i: ", i)
                # time = @elapsed loop_par!(L1, L, Nv, Nmem, Nv_list, il, ion1, ion2, U, UionAionBlist, hcmr, tempListNmem, nIterations, energyBase, bestE_hist, worstE_hist, par_hist, acceptance_hist)
                # time = @elapsed loop_par!(L1, L, Nv, Nmem, Nv_list, il, ion1, ion2, U, UionAionBlist, hcmr, tempListNmem, nIterations, energyBase, bestE_hist, worstE_hist, par_hist, acceptance_hist, vacs, du2List )
                # time = @elapsed loop_par!(L1, L, Nv, Nmem, Nv_list, il, ion1, ion2, U, UionAionBlist, hcmr, tempListNmem, nIterations, energyBase, bestE_hist, worstE_hist, par_hist, acceptance_hist, vacs, du2List, bestHarmony, listAvoidBest, emptiesDict)
                time = @elapsed loop_par!(L1, L, Nv, Nmem, Nv_list, il, ion1, ion2, U, UionAionBlist, tempListNmem, nIterations, energyBase, bestE_hist, worstE_hist, acceptance_hist, vacs, du2List, bestHarmony, listAvoidBest, emptiesDict, parValues, hcmrValues)
                println(il)
                # time = @elapsed loop!(L, Nmem, Nv_list, il, ion1, ion2, U, UionAionBlist, hcmr, par, tempListNmem, tempListL, nIterations, energyBase, bestE_hist, worstE_hist)
                results[1, i] = bestE_hist[end]
                results[2, i] = worstE_hist[end]
                results[3, i] = results[1, i] + 5591.3
                results[4, i] = results[2, i] + 5591.3
                results[5, i] = time
                #
                ################################################################
                mean_lastNenergies, std_lastNenergies, mean_lastNenergies_relative, std_lastNenergies_relative, mean_lastNacceptances, std_lastNacceptances = getData(lastN, bestE_hist, acceptance_hist)
                l_mean_lastNenergies[i] = mean_lastNenergies
                l_std_lastNenergies[i]  = std_lastNenergies
                l_mean_lastNenergies_relative[i] = mean_lastNenergies_relative
                l_std_lastNenergies_relative[i]  = std_lastNenergies_relative
                l_mean_lastNacceptances[i] = mean_lastNacceptances
                l_std_lastNacceptances[i]  = std_lastNacceptances
                l_time[i] = time
                #
                lastE = bestE_hist[end]
                lastAcc = acceptance_hist[end]
                l_lastE[i] = lastE
                l_lastAcc[i] = lastAcc
                namePlot1 = string("hsa___", i)
                plot1(namePlot1, nIterations, bestE_hist, worstE_hist, par_hist)
                ################################################################            
            end
            #
            yBestMean[k]  = mean( results[3, 1:5] )
            yBestDev[k]   = std(  results[3, 1:5] )

            yWorstMean[k] = mean( results[4, 1:5] )
            yWorstDev[k]  = std(  results[4, 1:5] )

            yTimeMean[k]  = mean( results[5, 1:5] )
            yTimeDev[k]   = std(  results[5, 1:5] )

            println("k, Nmem, WorstMean[k]: ", k, " ", Nmem, " ", yWorstMean[k], " ", yWorstDev[k])

            #
            ################################################################
            println("*******************************************************")
            println("length(l_mean_lastNenergies_relative), lastN: ", length(l_mean_lastNenergies_relative), " ", lastN)
            printListAndMean(nruns, l_mean_lastNenergies)
            printListAndMean(nruns, l_std_lastNenergies)
            printListAndMean(nruns, l_mean_lastNenergies_relative)
            printListAndMean(nruns, l_std_lastNenergies_relative)
            printListAndMean(nruns, l_mean_lastNacceptances)
            printListAndMean(nruns, l_std_lastNacceptances)
            printListAndMean(nruns, l_time)
            printListAndMean(nruns, l_lastE)
            printListAndMean(nruns, l_lastAcc)
            println("*******************************************************")
            ################################################################
        else
            tempListNmem = zeros(Int, Nmem)
            # time = @elapsed loop_par!(L1, L, Nv, Nmem, Nv_list, il, ion1, ion2, U, UionAionBlist, hcmr, tempListNmem, nIterations, energyBase, bestE_hist, worstE_hist, par_hist, acceptance_hist)
            # time = @elapsed loop_par!(L1, L, Nv, Nmem, Nv_list, il, ion1, ion2, U, UionAionBlist, hcmr, tempListNmem, nIterations, energyBase, bestE_hist, worstE_hist, par_hist, acceptance_hist, vacs, du2List )
            # time = @elapsed loop_par!(L1, L, Nv, Nmem, Nv_list, il, ion1, ion2, U, UionAionBlist, hcmr, tempListNmem, nIterations, energyBase, bestE_hist, worstE_hist, par_hist, acceptance_hist, vacs, du2List, bestHarmony, listAvoidBest, emptiesDict )
            time = @elapsed loop_par!(L1, L, Nv, Nmem, Nv_list, il, ion1, ion2, U, UionAionBlist, tempListNmem, nIterations, energyBase, bestE_hist, worstE_hist, acceptance_hist, vacs, du2List, bestHarmony, listAvoidBest, emptiesDict, parValues, hcmrValues)
            println(il)    
        end
    end


end

##


##
# str = getContentInputLammps(L1, L, il)
function getContentInputLammps(
    L1::Int,
    L::Int,
    il::Array{Int,1}
    )
    #
    str = """
    PRIM to lammps
    
    188 atoms
    6 atom types
    
       0.0000000000    12.9827000000 xlo xhi
       0.0000000000    12.9827000000 ylo yhi
       0.0000000000    12.9827000000 zlo zhi
    
    Masses
    
    1 138.90547   # La
    2 91.224      #Zr	
    3 15.99977    # O
    4 6.941       # Li_Td
    5 6.941       #Li_Oh
    6 69.723
    
    Atoms # full
    
       1  0  1   3.000000     6.4913500000     3.2456750000    11.3598625000
       2  0  1   3.000000    11.3598625000     6.4913500000     3.2456750000
       3  0  1   3.000000     6.4913500000     9.7370250000     8.1141875000
       4  0  1   3.000000     1.6228375000     0.0000000000     3.2456750000
       5  0  1   3.000000     0.0000000000     9.7370250000     4.8685125000
       6  0  1   3.000000     0.0000000000     3.2456750000     1.6228375000
       7  0  1   3.000000     3.2456750000    11.3598625000     6.4913500000
       8  0  1   3.000000     9.7370250000     4.8685125000     0.0000000000
       9  0  1   3.000000     4.8685125000     0.0000000000     9.7370250000
      10  0  1   3.000000     9.7370250000     8.1141875000     6.4913500000
      11  0  1   3.000000     3.2456750000     1.6228375000     0.0000000000
      12  0  1   3.000000     8.1141875000     6.4913500000     9.7370250000
      13  0  1   3.000000     0.0000000000     3.2456750000     8.1141875000
      14  0  1   3.000000     6.4913500000     9.7370250000     1.6228375000
      15  0  1   3.000000     0.0000000000     9.7370250000    11.3598625000
      16  0  1   3.000000     6.4913500000     3.2456750000     4.8685125000
      17  0  1   3.000000     3.2456750000     8.1141875000     0.0000000000
      18  0  1   3.000000     9.7370250000     1.6228375000     6.4913500000
      19  0  1   3.000000     9.7370250000    11.3598625000     0.0000000000
      20  0  1   3.000000     3.2456750000     4.8685125000     6.4913500000
      21  0  1   3.000000     8.1141875000     0.0000000000     3.2456750000
      22  0  1   3.000000     1.6228375000     6.4913500000     9.7370250000
      23  0  1   3.000000    11.3598625000     0.0000000000     9.7370250000
      24  0  1   3.000000     4.8685125000     6.4913500000     3.2456750000
      25  0  2   4.000000     9.7370250000     3.2456750000     9.7370250000
      26  0  2   4.000000     0.0000000000     6.4913500000     6.4913500000
      27  0  2   4.000000     3.2456750000     9.7370250000     9.7370250000
      28  0  2   4.000000     0.0000000000     0.0000000000     0.0000000000
      29  0  2   4.000000     9.7370250000     9.7370250000     3.2456750000
      30  0  2   4.000000     3.2456750000     3.2456750000     3.2456750000
      31  0  2   4.000000     6.4913500000     6.4913500000     0.0000000000
      32  0  2   4.000000     6.4913500000     0.0000000000     6.4913500000
      33  0  2   4.000000     3.2456750000     9.7370250000     3.2456750000
      34  0  2   4.000000     3.2456750000     3.2456750000     9.7370250000
      35  0  2   4.000000     9.7370250000     3.2456750000     3.2456750000
      36  0  2   4.000000     9.7370250000     9.7370250000     9.7370250000
      37  0  2   4.000000     0.0000000000     0.0000000000     6.4913500000
      38  0  2   4.000000     0.0000000000     6.4913500000     0.0000000000
      39  0  2   4.000000     6.4913500000     0.0000000000     0.0000000000
      40  0  2   4.000000     6.4913500000     6.4913500000     6.4913500000
      41  0  3  -2.000000     7.7883217300     3.9441442600     9.3263821990
      42  0  3  -2.000000     0.4106428010     7.1898192600     4.5426467300
      43  0  3  -2.000000     5.1943782700    10.4354942600    10.1476678010
      44  0  3  -2.000000    12.5720571990     0.6984692600     1.9487032700
      45  0  3  -2.000000    11.6857282700     9.0385557400     2.8350321990
      46  0  3  -2.000000     1.2969717300     2.5472057400     3.6563178010
      47  0  3  -2.000000     3.9441442600     9.3263821990     7.7883217300
      48  0  3  -2.000000     9.0385557400     2.8350321990    11.6857282700
      49  0  3  -2.000000     6.9019928010    12.2842307400     8.4400532700
      50  0  3  -2.000000    10.4354942600    10.1476678010     5.1943782700
      51  0  3  -2.000000     2.5472057400     3.6563178010     1.2969717300
      52  0  3  -2.000000     9.3263821990     7.7883217300     3.9441442600
      53  0  3  -2.000000     6.0807071990     5.7928807400    11.0339967300
      54  0  3  -2.000000     2.8350321990    11.6857282700     9.0385557400
      55  0  3  -2.000000    10.1476678010     5.1943782700    10.4354942600
      56  0  3  -2.000000     3.6563178010     1.2969717300     2.5472057400
      57  0  3  -2.000000    11.6857282700     2.5472057400    10.1476678010
      58  0  3  -2.000000     6.9019928010     0.6984692600    11.0339967300
      59  0  3  -2.000000     7.7883217300     2.5472057400     2.8350321990
      60  0  3  -2.000000     7.7883217300    10.4354942600     3.6563178010
      61  0  3  -2.000000     6.0807071990     7.1898192600     8.4400532700
      62  0  3  -2.000000     1.2969717300     3.9441442600    10.1476678010
      63  0  3  -2.000000     1.2969717300     9.0385557400     9.3263821990
      64  0  3  -2.000000     5.1943782700     3.9441442600     2.8350321990
      65  0  3  -2.000000     2.5472057400    10.1476678010    11.6857282700
      66  0  3  -2.000000    10.4354942600     3.6563178010     7.7883217300
      67  0  3  -2.000000     0.4106428010     5.7928807400     1.9487032700
      68  0  3  -2.000000     2.5472057400     2.8350321990     7.7883217300
      69  0  3  -2.000000     9.0385557400     9.3263821990     1.2969717300
      70  0  3  -2.000000     3.9441442600     2.8350321990     5.1943782700
      71  0  3  -2.000000    10.1476678010    11.6857282700     2.5472057400
      72  0  3  -2.000000    12.5720571990    12.2842307400     4.5426467300
      73  0  3  -2.000000    10.1476678010     1.2969717300     3.9441442600
      74  0  3  -2.000000     3.6563178010     7.7883217300    10.4354942600
      75  0  3  -2.000000     9.3263821990     1.2969717300     9.0385557400
      76  0  3  -2.000000     2.8350321990     5.1943782700     3.9441442600
      77  0  3  -2.000000     4.5426467300     6.0807071990     0.6984692600
      78  0  3  -2.000000     5.7928807400    11.0339967300     6.0807071990
      79  0  3  -2.000000     0.4106428010    12.2842307400    11.0339967300
      80  0  3  -2.000000     8.4400532700     6.9019928010    12.2842307400
      81  0  3  -2.000000     7.1898192600     1.9487032700     6.9019928010
      82  0  3  -2.000000     8.4400532700    12.5720571990     5.7928807400
      83  0  3  -2.000000    12.2842307400     8.4400532700     6.9019928010
      84  0  3  -2.000000     1.9487032700     6.9019928010     7.1898192600
      85  0  3  -2.000000     7.1898192600     4.5426467300     0.4106428010
      86  0  3  -2.000000    11.0339967300     0.4106428010    12.2842307400
      87  0  3  -2.000000     0.6984692600     1.9487032700    12.5720571990
      88  0  3  -2.000000     0.6984692600     4.5426467300     6.0807071990
      89  0  3  -2.000000     4.5426467300     0.4106428010     7.1898192600
      90  0  3  -2.000000     5.7928807400     8.4400532700    12.5720571990
      91  0  3  -2.000000    11.0339967300     6.0807071990     5.7928807400
      92  0  3  -2.000000    12.2842307400    11.0339967300     0.4106428010
      93  0  3  -2.000000     1.9487032700    12.5720571990     0.6984692600
      94  0  3  -2.000000     6.0807071990     0.6984692600     4.5426467300
      95  0  3  -2.000000    12.5720571990     5.7928807400     8.4400532700
      96  0  3  -2.000000     6.9019928010     7.1898192600     1.9487032700
      97  0  3  -2.000000     1.2969717300    10.4354942600     2.8350321990
      98  0  3  -2.000000     5.1943782700     2.5472057400     9.3263821990
      99  0  3  -2.000000    11.6857282700     3.9441442600     3.6563178010
     100  0  3  -2.000000     7.7883217300     9.0385557400    10.1476678010
     101  0  3  -2.000000    10.4354942600     2.8350321990     1.2969717300
     102  0  3  -2.000000     2.5472057400     9.3263821990     5.1943782700
     103  0  3  -2.000000     3.9441442600     3.6563178010    11.6857282700
     104  0  3  -2.000000     9.0385557400    10.1476678010     7.7883217300
     105  0  3  -2.000000     2.8350321990     1.2969717300    10.4354942600
     106  0  3  -2.000000     9.3263821990     5.1943782700     2.5472057400
     107  0  3  -2.000000     3.6563178010    11.6857282700     3.9441442600
     108  0  3  -2.000000    10.1476678010     7.7883217300     9.0385557400
     109  0  3  -2.000000     5.1943782700     9.0385557400     3.6563178010
     110  0  3  -2.000000    11.6857282700    10.4354942600     9.3263821990
     111  0  3  -2.000000     9.0385557400     3.6563178010     5.1943782700
     112  0  3  -2.000000     3.9441442600    10.1476678010     1.2969717300
     113  0  3  -2.000000    10.4354942600     9.3263821990    11.6857282700
     114  0  3  -2.000000     3.6563178010     5.1943782700     9.0385557400
     115  0  3  -2.000000     2.8350321990     7.7883217300     2.5472057400
     116  0  3  -2.000000     9.3263821990    11.6857282700    10.4354942600
     117  0  3  -2.000000    11.0339967300    12.5720571990     7.1898192600
     118  0  3  -2.000000     1.9487032700     6.0807071990    12.2842307400
     119  0  3  -2.000000     8.4400532700     0.4106428010     0.6984692600
     120  0  3  -2.000000     4.5426467300     6.9019928010     5.7928807400
     121  0  3  -2.000000     7.1898192600    11.0339967300    12.5720571990
     122  0  3  -2.000000    12.2842307400     1.9487032700     6.0807071990
     123  0  3  -2.000000     0.6984692600     8.4400532700     0.4106428010
     124  0  3  -2.000000     5.7928807400     4.5426467300     6.9019928010
     125  0  3  -2.000000    12.5720571990     7.1898192600    11.0339967300
     126  0  3  -2.000000     6.0807071990    12.2842307400     1.9487032700
     127  0  3  -2.000000     0.4106428010     0.6984692600     8.4400532700
     128  0  3  -2.000000     6.9019928010     5.7928807400     4.5426467300
     129  0  3  -2.000000     1.9487032700     0.4106428010     5.7928807400
     130  0  3  -2.000000    11.0339967300     6.9019928010     0.6984692600
     131  0  3  -2.000000     4.5426467300    12.5720571990    12.2842307400
     132  0  3  -2.000000     8.4400532700     6.0807071990     7.1898192600
     133  0  3  -2.000000     5.7928807400     1.9487032700     0.4106428010
     134  0  3  -2.000000     0.6984692600    11.0339967300     6.9019928010
     135  0  3  -2.000000    12.2842307400     4.5426467300    12.5720571990
     136  0  3  -2.000000     7.1898192600     8.4400532700     6.0807071990
    """
    indx = 136

    for i in 1:L
        if i <= L1
            ion = 6
            Qion = 3.0
        else
            ion = 4
            Qion = 1.0
        end
        indx += 1

        site = il[i] + 136

        x = positions[site, 1]
        y = positions[site, 2]
        z = positions[site, 3]
        str *= " " * string(indx) * "  0  " * string(ion) * "  " * string(Qion) * "  " * string(x) * " " * string(y) * " " * string(z) * "\n"
    end
    return str
end

# getEnergyLammps(L1, L, il)
function getEnergyLammps(
    L1::Int,
    L::Int,
    il::Array{Int,1}
    )
    #
    str = getContentInputLammps(L1, L, il)
    # str = "TEST"
    io = open("data.lammps", "w")
    println(io, str)
    close(io)

    run(`./myscript.sh`)

    # open("energy.dat") do io
    #     flg = 0
    #     while !eof(io)
    #         line = readline(io)
    #         if "Step PotEng E_coul E_long E_vdwl E_tail" == line
    #             flg = 1
    #             continue
    #         end
    #         if flg == 1
    #             if "Loop time" in line
    #                 break
    #             end
    #             data = split(line, " ")
    #             return float(data[2])  
    #         end
    #     end
    # end
    # return -1.0
    io = open("energy.dat")
    data = split(readline(io))
    close(io)
    run(`rm energy.dat`)
    # println(data)
    return parse( Float64, data[1] )
end


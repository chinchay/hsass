#
include("./HSA.jl")
#
# Base algorithm:
L1 = 10 
L2 = 418
Nv = 576
# maximumMoves = 10
cutOff = 10.0 #4.3 #10.0 #4.3 #6.0
# file = "testVac/data.lammps"
file = "testVac/data.lammps_222"
L  = L1 + L2
ion1=1
ion2=2
Ne = Nv - L # Number of empty vacancies
# #
distSite2Col, neighbors_of, charges, ion1, ion2, removedSites, U, UionAionBlist = initialize(L1, L2, Nv, cutOff, file)
#
energyBase = getEnergyBase(removedSites, charges, distSite2Col, neighbors_of)
println("energyBase: ", energyBase)
#
# ##
# il = zeros(Int, L)
# Nv_list = collect(1:Nv)
# tempListL = zeros(Int,L)
# tempListNmem = zeros(Int, Nmem)
# hcmr = 1.0    #1.0
# par = 1.0   #0.1
# nIterations = 100_000 #200_000
# Nmem = 10
# # bestEnergyMem = zeros(nIterations)
# # worstEnergyMem = zeros(nIterations)
# # for i in 1:5
# #     @time loop!(L, Nmem, Nv_list, il, ion1, ion2, U, UionAionBlist, hcmr, par, tempListNmem, tempListL, nIterations, energyBase, bestEnergyMem, worstEnergyMem)
# # end
#
# ##
# hcmr = 0.9   #1.0
# # nIterations = 10_000 #200_000
# nIterations = 200_000
# nruns = 5
# # x_vals = collect( 10:10:100 ) # values of Nmem
# x_vals = collect( 10:10 ) # values of Nmem
# n = length(x_vals)
# yBestMean  = zeros(n)
# yBestDev   = zeros(n)
# yWorstMean = zeros(n)
# yWorstDev  = zeros(n)
# yTimeMean  = zeros(n)
# yTimeDev   = zeros(n)
# # for EACH simulation, we need to use:
# bestE_hist  = zeros(nIterations)
# worstE_hist = zeros(nIterations)
# par_hist    = zeros(nIterations) # used if nruns == 1
# #
# # experiment10!( nruns, L, Nmem, Nv_list, il, ion1, ion2, U, UionAionBlist, hcmr, tempListNmem, nIterations, energyBase, x_vals, yBestMean, yBestDev, yWorstMean, yWorstDev, yTimeMean, yTimeDev, bestE_hist, worstE_hist, par_hist)
# experiment11!(nruns, L, Nv, Nv_list, il, ion1, ion2, U, UionAionBlist, hcmr, nIterations, energyBase, x_vals, yBestMean, yBestDev, yWorstMean, yWorstDev, yTimeMean, yTimeDev, bestE_hist, worstE_hist, par_hist)
# # answer = string("")
# # run(`sayProgramFinished $answer`)
#
# ##
# if nruns > 1
#     y1max = 20.0 #140
#     y3max = 3.0 #10.0
#     xlabel = "HMS"
#     plot2(1.0*x_vals, yBestMean, yBestDev, yWorstMean, yWorstDev, yTimeMean, yTimeDev, y1max, y3max, xlabel)
#     plot1("hsa_", nIterations, bestE_hist, worstE_hist, par_hist)
# else
#     plot1("hsa_", nIterations, bestE_hist, worstE_hist, par_hist)
# end


##


##
##
il = zeros(Int, L)
Nv_list = collect(1:Nv)
tempListL = zeros(Int,L)
# tempListNmem = zeros(Int, Nmem)
# hcmr = 1.0    #1.0
# par = 1.0   #0.1
# nIterations = 100_000 #200_000
# Nmem = 10
# bestEnergyMem = zeros(nIterations)
# worstEnergyMem = zeros(nIterations)
# for i in 1:5
#     @time loop!(L, Nmem, Nv_list, il, ion1, ion2, U, UionAionBlist, hcmr, par, tempListNmem, tempListL, nIterations, energyBase, bestEnergyMem, worstEnergyMem)
# end



##
# if `is_PAR_dynamic=false` then `parMin` and `parMax` are not used
# if `is_PAR_dynamic=true` then parMin < parMax belongs to [0.0, 1.0]
par = 1.0 #
is_PAR_dynamic = false # options: false, true
parMin = 0.0 # 0.0
parMax = 0.0 # 0.0

# if `is_HCMR_dynamic=false` then `hcmrMin` and `hcmrMax` are not used
# if `is_HCMR_dynamic=true` then hcmrMin < hcmrMax belongs to [0.0, 1.0]
hcmr = 1.0 #
is_HCMR_dynamic = false
hcmrMin = 0.0 # 0.0
hcmrMax = 0.0 # 0.0

nIterations = 100_000 # 200_000
nruns = 1 #5
memory_vals = collect( 10:10:100 ) # values of Nmem
# memory_vals = collect( 10:10 ) # values of Nmem
n = length(memory_vals)
yBestMean  = zeros(n)
yBestDev   = zeros(n)
yWorstMean = zeros(n)
yWorstDev  = zeros(n)
yTimeMean  = zeros(n)
yTimeDev   = zeros(n)
# for EACH simulation, we need to use:
bestE_hist  = zeros(nIterations)
worstE_hist = zeros(nIterations)
par_hist    = zeros(nIterations) # used if nruns == 1


bestHarmony = zeros(Int, L)
emptiesDict = Dict{Int,Bool}( site => true for site in 1:Nv )



using Profile
# @profiler experiment11!(nruns, L, Nv, Nv_list, il, ion1, ion2, U, UionAionBlist, hcmr, nIterations, energyBase, memory_vals, yBestMean, yBestDev, yWorstMean, yWorstDev, yTimeMean, yTimeDev, bestE_hist, worstE_hist, par_hist)
experiment11!(nruns, L, Nv, Nv_list, il, ion1, ion2, U, UionAionBlist, nIterations, energyBase, memory_vals, yBestMean, yBestDev, yWorstMean, yWorstDev, yTimeMean, yTimeDev, bestE_hist, worstE_hist, par_hist, bestHarmony, emptiesDict, par, is_PAR_dynamic, parMin, parMax, hcmr, is_HCMR_dynamic, hcmrMin, hcmrMax )

# u1 = get_il_Ut_plus_Uion(L1, L, il, U, UionAionBlist, ion1, ion2, energyBase)
# u2 = getEnergyLammps(L1, L, il)
# if abs(u1 - u2) > 0.1
#     throw(error())
#     println(u1, " /// ", u2)
# # else
# #     println(u1, " ... ", u2)
# end

##



if nruns > 1
    y1max = 15.0
    y3max = 10.0
    xlabel = "HMS"
    plot2(1.0*memory_vals, yBestMean, yBestDev, yWorstMean, yWorstDev, yTimeMean, yTimeDev, y1max, y3max, xlabel)
else
    plot1("hsa_", nIterations, bestE_hist, worstE_hist, par_hist)
end

# answer = string("")
# run(`sayProgramFinished $answer`)
##
# if nruns > 1
#     y1max = 140.0
#     y3max = 6.0
#     xlabel = "PAR"
#     plot2(memory_vals, yBestMean, yBestDev, yWorstMean, yWorstDev, yTimeMean, yTimeDev, y1max, y3max, xlabel)
# else
#     plot1(nIterations, bestE_hist, worstE_hist, par_hist)
# end




#
#
# ##
# bestEnergyMem = zeros(nIterations)
# worstEnergyMem = zeros(nIterations)
# lpar = zeros(nIterations)
# # @time loop!(L, Nmem, Nv_list, il, ion1, ion2, U, UionAionBlist, hcmr, par, tempListNmem, tempListL, nIterations, energyBase, bestEnergyMem, worstEnergyMem)
# # @time loop_experiment1!(L, Nmem, Nv_list, il, ion1, ion2, U, UionAionBlist, hcmr, tempListNmem, nIterations, energyBase, bestEnergyMem, worstEnergyMem, lpar)
# ##
# # plot1(nIterations, bestEnergyMem, worstEnergyMem, lpar)
# ##
#
#
# ##
# # vals_nIterations = 1_000:20_000:210_000
# # temporal = collect(vals_nIterations)
#
# vals_Nmem = 1:2:100
# temporal = collect(vals_Nmem)
#
# x = zeros(length(temporal))
# yWorstMean = zeros(length(temporal))
# yWorstDev  = zeros(length(temporal))
# yBestMean  = zeros(length(temporal))
# yBestDev   = zeros(length(temporal))
# yTimeMean  = zeros(length(temporal))
# yTimeDev   = zeros(length(temporal))
#
# # for (k, nIterations) in enumerate(vals_nIterations)
# # println("k, nIterations: ", k, " ", nIterations)
#
# for (k, Nmem) in enumerate(vals_Nmem)
#     tempListNmem = zeros(Int, Nmem)
#     tempListL = zeros(Int,L)
#     println("k, Nmem: ", k, " ", Nmem)
#
#     bestEnergyMem = zeros(nIterations)
#     worstEnergyMem = zeros(nIterations)
#     nruns = 5
#     # nruns = 100
#     results = zeros(5, nruns )
#     for i in 1:nruns
#         println("i: ", i)
#         # @time loop!(L, Nmem, Nv_list, il, ion1, ion2, U, UionAionBlist, hcmr, par, tempListNmem, tempListL, nIterations, energyBase, bestEnergyMem, worstEnergyMem)
#
#         time = @elapsed loop!(L, Nmem, Nv_list, il, ion1, ion2, U, UionAionBlist, hcmr, par, tempListNmem, tempListL, nIterations, energyBase, bestEnergyMem, worstEnergyMem)
#         results[1, i] = worstEnergyMem[end]
#         results[2, i] = bestEnergyMem[end]
#         results[3, i] = results[1, i] + 5591.3
#         results[4, i] = results[2, i] + 5591.3
#         results[5, i] = time
#     end
#
#     # x[k] = nIterations
#     x[k] = Nmem
#
#     yWorstMean[k] = mean( results[3, 1:5] )
#     yWorstDev[k]  = std(  results[3, 1:5] )
#
#     yBestMean[k]  = mean( results[4, 1:5] )
#     yBestDev[k]   = std(  results[4, 1:5] )
#
#     yTimeMean[k]  = mean( results[5, 1:5] )
#     yTimeDev[k]   = std(  results[5, 1:5] )
# end
# ##
#
# yWorstMean
# PyPlot.close() # it gave me errors saying that other plots were open.
# myColors = ["k", "r", "g", "m", "y"]
# fig, ax1 = PyPlot.subplots()
# # ax1.scatter(x, yWorstMean, s=0.5, color=myColors[1])
# # ax1.plot(x, yWorstMean, "--o", color=myColors[1])
# # ax1.errorbar(x, yWorstMean, yerr=yWorstDev, "--o", color=myColors[1])
# plt.errorbar(x, yWorstMean, yerr=yWorstDev, color=myColors[1])
# # ax1.scatter(x, yBestMean, s=0.5, color=myColors[2])
# # ax1.plot(x, yBestMean, "--o", color=myColors[2])
# plt.errorbar(x, yBestMean, yerr=yBestDev, color=myColors[2])
# # xlabel = string("move x ", factNumPoints)
#
# # ax1.set_xlabel("nIterations")
# ax1.set_xlabel("HMS")
#
# ax1.set_ylabel("mean energy (eV)", color="k")
# ax1.set_ylim([0, 140])
# fig.savefig("hsa_evolution")
# #
# ax2 = ax1.twinx()
# # ax2.plot(x, yTimeMean, "--ob", linewidth=0.4)
# plt.errorbar(x, yTimeMean, yerr=yTimeDev, color="b")
# # ax2.scatter(x, yTimeMean, s=0.5, color="b")
# ax2.set_ylabel("time (seconds)", color="b")
# ax2.set_ylim([0.0, 7.0])
# ax2.tick_params(axis="y", colors="blue")
#
#
# fig.savefig("hsa_evolution")
# #
# PyPlot.close(fig)
#
#
#
# # println()
# # for row in 1:5
# #     avg = mean( results[row, 1:5] )
# #     std = std( results[row, 1:5] )
# #     println(avg, ", ", std)
# #     # for col in 1:nruns
# #     #     # print(results[row, col], ",")
# #     # end
#
# #     println()
# # end
# # println()
# # println()
#
#
#
# ##
# PyPlot.close() # it gave me errors saying that other plots were open.
# myColors = ["k", "r", "g", "m", "y"]
# fig, ax1 = PyPlot.subplots()
# ax1.scatter(collect(1:nIterations), bestEnergyMem, s=0.5, color=myColors[1])
# ax1.scatter(collect(1:nIterations), worstEnergyMem, s=0.5, color=myColors[2])
# # xlabel = string("move x ", factNumPoints)
# ax1.set_xlabel("iteration")
# ax1.set_ylabel("energy (eV)", color="k")
# ax1.set_ylim([-5600, -5150])
# fig.savefig("hsa_")
# PyPlot.close(fig)
#
#
# #
# # aa = [2 17 13 5 10 54 52 23 27 4 16 22 3 38 11 48 14 30 12 9 56 44 19 62 29 6 36 47 43 7 60 68 8 55 20 49 31 59 51 21 15 72 18 40 1 25 57 45 35 69 24 65]
# # bb = zeros(Int,52)
# # for i in 1:52
# #     bb[i] = aa[i]
# # end
# #
#
# alist = rand(collect(1:72), 100)
# blist = rand(100, 72, 2, 2)
# function fill_(
#     n::Int,
#     alist::Array{Int,1},
#     blist::Array{Float64,4}
#     )
#     for j in 1:n
#         sum_a = 0
#         sum_b = 0.0
#         for i in 1:100
#             sum_a += alist[i]
#             sum_b += blist[i, 50, 1, 2]
#         end
#     end
#     nothing
# end
# # using ProfileView
# # @profview fill_(1, alist, blist)  # run once to trigger compilation (ignore this one)
# # @profview fill_(10, alist, blist)
#
# using Profile
# Profile.init()
# @profiler fill_(100, alist, blist)

##
include("./customreadlammps.jl")
file = "testVac/data.lammps"
structure, charges = getStructure(file)
positions = structure.cart_coords #position of each site

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



##
aa = [28, 31, 21, 10, 45, 36, 71, 32, 48, 40, 63, 49, 39, 12, 19, 30, 37, 51, 17, 6, 67, 18, 22, 64, 70, 7, 27, 47, 26, 57, 38, 5, 16, 55, 14, 61, 60, 46, 8, 1, 66, 13, 42, 53, 43, 58, 15, 59, 34, 2, 29, 72]
# getEnergyLammps(L1, L, il)
getEnergyLammps(L1, L, aa)


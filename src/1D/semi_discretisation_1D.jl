


# @views av(A) = 0.5*(A[1:end-1].+A[2:end])
# # macros to avoid array allocation
# macro qx(D,C,ix, Δx_ad) esc(:((0.5*($D[$ix]+$D[$ix+1])) * ($C[$ix+1]-$C[$ix])/Δx_ad )) end
# macro sum_D(CMg, CFe, CMn, D0, ix) esc(:($D0[1] * $CMg[$ix] + $D0[2] * $CFe[$ix] + $D0[3] * $CMn[$ix] +
#                                        $D0[4] * (1 - $CMg[$ix] - $CFe[$ix] - $CMn[$ix]))) end



# function stencil_diffusion_1D!(dtCMg, dtCFe, dtCMn, CMg, CFe ,CMn, DMgMg, DMgFe, DMgMn, DFeMg, DFeFe, DFeMn, DMnMg, DMnFe, DMnMn, Δx_ad)


#     for I in eachindex((dtCMg))
#         if I<size(dtCMg,1) && I> 1
#             dtCMg[ix] = (@qx(DMgMg,CMg,ix,Δx_ad) - @qx(DMgMg,CMg,ix-1,Δx_ad)) / Δx_ad +
#                         (@qx(DMgFe,CFe,ix,Δx_ad) - @qx(DMgFe,CFe,ix-1,Δx_ad)) / Δx_ad +
#                         (@qx(DMgMn,CMn,ix,Δx_ad) - @qx(DMgMn,CMn,ix-1,Δx_ad)) / Δx_ad
#             dtCFe[ix] = (@qx(DFeMg,CMg,ix,Δx_ad) - @qx(DFeMg,CMg,ix-1,Δx_ad)) / Δx_ad +
#                         (@qx(DFeFe,CFe,ix,Δx_ad) - @qx(DFeFe,CFe,ix-1,Δx_ad)) / Δx_ad +
#                         (@qx(DFeMn,CMn,ix,Δx_ad) - @qx(DFeMn,CMn,ix-1,Δx_ad)) / Δx_ad
#             dtCMn[ix] = (@qx(DMnMg,CMg,ix,Δx_ad) - @qx(DMnMg,CMg,ix-1,Δx_ad)) / Δx_ad +
#                         (@qx(DMnFe,CFe,ix,Δx_ad) - @qx(DMnFe,CFe,ix-1,Δx_ad)) / Δx_ad +
#                         (@qx(DMnMn,CMn,ix,Δx_ad) - @qx(DMnMn,CMn,ix-1,Δx_ad)) / Δx_ad
#         else
#             dtCMg[ix] = 0
#             dtCFe[ix] = 0
#             dtCMn[ix] = 0
#         end
#     end
# end


# test = ones(100,100)

# dtCMg = ones(100)
# dtCMn = ones(100)
# dtCFe = ones(100)


# for I in eachindex((dtCMg))
#     if I<size(dtCMg,1) && I> 1
#         dtCMg[ix] = (@qx(DMgMg,CMg,ix,Δx_ad) - @qx(DMgMg,CMg,ix-1,Δx_ad)) / Δx_ad +
#                     (@qx(DMgFe,CFe,ix,Δx_ad) - @qx(DMgFe,CFe,ix-1,Δx_ad)) / Δx_ad +
#                     (@qx(DMgMn,CMn,ix,Δx_ad) - @qx(DMgMn,CMn,ix-1,Δx_ad)) / Δx_ad
#         dtCFe[ix] = (@qx(DFeMg,CMg,ix,Δx_ad) - @qx(DFeMg,CMg,ix-1,Δx_ad)) / Δx_ad +
#                     (@qx(DFeFe,CFe,ix,Δx_ad) - @qx(DFeFe,CFe,ix-1,Δx_ad)) / Δx_ad +
#                     (@qx(DFeMn,CMn,ix,Δx_ad) - @qx(DFeMn,CMn,ix-1,Δx_ad)) / Δx_ad
#         dtCMn[ix] = (@qx(DMnMg,CMg,ix,Δx_ad) - @qx(DMnMg,CMg,ix-1,Δx_ad)) / Δx_ad +
#                     (@qx(DMnFe,CFe,ix,Δx_ad) - @qx(DMnFe,CFe,ix-1,Δx_ad)) / Δx_ad +
#                     (@qx(DMnMn,CMn,ix,Δx_ad) - @qx(DMnMn,CMn,ix-1,Δx_ad)) / Δx_ad
#     else
#         dtCMg[ix] = 0
#         dtCFe[ix] = 0
#         dtCMn[ix] = 0
#     end
# end


# for i in eachindex(test)

#     @show i

# end
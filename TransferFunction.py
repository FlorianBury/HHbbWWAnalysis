import os.path
import sys
from bamboo.analysismodules import NanoAODModule, NanoAODHistoModule, DataDrivenBackgroundHistogramsModule
from bamboo import treefunctions as op
from bamboo import treeoperations as to
from bamboo.analysisutils import makeMultiPrimaryDatasetTriggerSelection
from bamboo.plots import Plot, EquidistantBinning, SummedPlot, Selection, Skim

sys.path.append(os.path.dirname(os.path.abspath(__file__))) # Add scripts in this directory

from BaseHHtobbWW import BaseNanoHHtobbWW


lambda_hadronFlavour = lambda subj: op.multiSwitch(
    (subj.nBHadrons>0, op.c_int(5)),
    (subj.nCHadrons>0, op.c_int(4)),
    op.c_int(0)
)
lambda_partonFlavour = lambda gen: op.multiSwitch(
    (op.abs(gen.pdgId)==5, op.c_int(5)),
    (op.abs(gen.pdgId)==4, op.c_int(4)),
    op.c_int(0)
)



def plotNumber(contName,cont,sel,Nmax,xTitle):
    return [Plot.make1D(contName,
                        op.rng_len(cont),
                        sel,
                        EquidistantBinning(Nmax,0.,Nmax),
                        xTitle=xTitle)]

def plotMatching(contName,cont,sel,type,Egen_binning=(3000,0,3000),dE_binning=(8000,-1000,1000)):
    plots = []
    assert type in ['lepton','jet','tau']
    # [0] is the last copy gen object = post FSR
    # [1] is the gen object found in the recursion = pre FSR
    # [2] is the reco object


    plots.append(Plot.make1D(contName+"_preFSR_postFSR_dR",
                             op.map(cont,lambda m : op.deltaR(m[0][0].p4,m[0][1].p4)),
                             sel,
                             EquidistantBinning(100,0.,1.0)))
    plots.append(Plot.make1D(contName+"_postFSR_reco_dR",
                             op.map(cont,lambda m : op.deltaR(m[0][1].p4,m[1].p4)),
                             sel,
                             EquidistantBinning(100,0.,1.0)))
    plots.append(Plot.make1D(contName+"_preFSR_reco_dR",
                             op.map(cont,lambda m : op.deltaR(m[0][0].p4,m[1].p4)),
                             sel,
                             EquidistantBinning(100,0.,1.0)))

    plots.append(Plot.make1D(contName+"_preFSR_postFSR_dE",
                             op.map(cont,lambda m : 2*(m[0][0].p4.E()-m[0][1].p4.E())/(m[0][0].p4.E()+m[0][1].p4.E())),
                             sel,
                             EquidistantBinning(200,-5.,5.)))
    plots.append(Plot.make1D(contName+"_postFSR_reco_dE",
                             op.map(cont,lambda m : 2*(m[0][1].p4.E()-m[1].p4.E())/(m[0][1].p4.E()+m[1].p4.E())),
                             sel,
                             EquidistantBinning(200,-5.,5.)))
    plots.append(Plot.make1D(contName+"_preFSR_reco_dE",
                             op.map(cont,lambda m : 2*(m[0][0].p4.E()-m[1].p4.E())/(m[0][0].p4.E()+m[1].p4.E())),
                             sel,
                             EquidistantBinning(200,-5.,5.)))

    if type in ['lepton','tau']:
        plots.append(Plot.make2D(contName+"_preFSR_postFSR_flavour",
                                 [op.map(cont,lambda m : m[0][0].pdgId),op.map(cont,lambda m : m[0][1].pdgId)],
                                 sel,
                                 [EquidistantBinning(25,0,25),EquidistantBinning(25,0,25)]))
        if type == 'lepton':
            plots.append(Plot.make2D(contName+"_postFSR_reco_flavour",
                                     [op.map(cont,lambda m : m[0][1].pdgId),op.map(cont,lambda m : m[1].genPartFlav)],
                                     sel,
                                     [EquidistantBinning(25,0,25),EquidistantBinning(25,0,25)]))
    if type in ['jet','tau']:
        plots.append(Plot.make1D(contName+"_genjet_reco_dR",
                                 op.map(cont,lambda m : op.deltaR(m[1].genJet.p4,m[1].p4)),
                                 sel,
                                 EquidistantBinning(100,0.,1.0)))
        plots.append(Plot.make1D(contName+"_postFSR_genJet_dR",
                                 op.map(cont,lambda m : op.deltaR(m[0][1].p4,m[1].genJet.p4)),
                                 sel,
                                 EquidistantBinning(100,0.,1.0)))
        plots.append(Plot.make1D(contName+"_preFSR_genJet_dR",
                                 op.map(cont,lambda m : op.deltaR(m[0][0].p4,m[1].genJet.p4)),
                                 sel,
                                 EquidistantBinning(100,0.,1.0)))
        plots.append(Plot.make1D(contName+"_genjet_reco_dE",
                                 op.map(cont,lambda m : 2*(m[1].genJet.p4.E()-m[1].p4.E())/(m[1].genJet.p4.E()+m[1].p4.E())),
                                 sel,
                                 EquidistantBinning(200,-5.,5.)))
        if type == 'jet':
            plots.append(Plot.make2D(contName+"_flavour",
                                     [op.map(cont,lambda m : lambda_partonFlavour(m[0][1])),op.map(cont,lambda m : m[1].hadronFlavour)],
                                     sel,
                                     [EquidistantBinning(6,0,6),EquidistantBinning(6,0,6)]))

    plots.append(Plot.make2D(contName+"_TF",
                             (op.map(cont,lambda m : m[0][0].p4.E()),op.map(cont,lambda m : m[1].p4.E()-m[0][0].p4.E())),
                             sel,
                             [EquidistantBinning(*Egen_binning),EquidistantBinning(*dE_binning)]))

    return plots

def plotFatjetMatching(cont,sel,sub1,sub2,Egen_binning=(3000,0,3000),dE_binning=(8000,-1000,1000)):
    plots = []

    # Subjet 1 #
    plots.append(Plot.make1D(f"fatjet{sub1}{sub2}_postFSR_sub1_dR",
                             op.map(cont,lambda m : op.deltaR(m[2].subJet1.p4,m[0][1].p4)),
                             sel,
                             EquidistantBinning(100,0.,1.0)))
    plots.append(Plot.make1D(f"fatjet{sub1}{sub2}_preFSR_sub1_dR",
                             op.map(cont,lambda m : op.deltaR(m[2].subJet1.p4,m[0][0].p4)),
                             sel,
                             EquidistantBinning(100,0.,1.0)))
    plots.append(Plot.make1D(f"fatjet{sub1}{sub2}_postFSR_sub1_dE",
                             op.map(cont,lambda m : 2*(m[0][1].p4.E()-m[2].subJet1.p4.E())/(m[0][1].p4.E()+m[2].subJet1.p4.E())),
                             sel,
                             EquidistantBinning(200,-5.,5.0)))
    plots.append(Plot.make1D(f"fatjet{sub1}{sub2}_preFSR_sub1_dE",
                             op.map(cont,lambda m : 2*(m[0][0].p4.E()-m[2].subJet1.p4.E())/(m[0][0].p4.E()+m[2].subJet1.p4.E())),
                             sel,
                             EquidistantBinning(200,-5.,5.0)))
    plots.append(Plot.make2D(f"fatjet{sub1}{sub2}_sub1_flavour",
                             [op.map(cont, lambda m: lambda_partonFlavour(m[0][1])), op.map(cont, lambda m: lambda_hadronFlavour(m[2].subJet1))],
                             sel,
                             [EquidistantBinning(6,0.,6.),EquidistantBinning(6,0.,6.)]))
    plots.append(Plot.make2D(f"fatjet{sub1}{sub2}_sub1_TF",
                             [op.map(cont,lambda m : m[0][0].p4.E()),op.map(cont,lambda m : m[2].subJet1.p4.E()-m[0][0].p4.E())],
                             sel,
                             [EquidistantBinning(*Egen_binning),EquidistantBinning(*dE_binning)],
                             xTitle="E^{{parton}}({sub1})",
                             yTitle=f"#Delta E = E^{{reco}}({sub1})-E^{{parton}}({sub1})"))
    # Subjet 2 #
    plots.append(Plot.make1D(f"fatjet{sub2}{sub2}_postFSR_sub2_dR",
                             op.map(cont,lambda m : op.deltaR(m[2].subJet2.p4,m[1][1].p4)),
                             sel,
                             EquidistantBinning(100,0.,1.0)))
    plots.append(Plot.make1D(f"fatjet{sub2}{sub2}_preFSR_sub2_dR",
                             op.map(cont,lambda m : op.deltaR(m[2].subJet2.p4,m[1][0].p4)),
                             sel,
                             EquidistantBinning(100,0.,1.0)))
    plots.append(Plot.make1D(f"fatjet{sub2}{sub2}_postFSR_sub2_dE",
                             op.map(cont,lambda m : 2*(m[1][1].p4.E()-m[2].subJet2.p4.E())/(m[1][1].p4.E()+m[2].subJet2.p4.E())),
                             sel,
                             EquidistantBinning(200,-5.,5.0)))
    plots.append(Plot.make1D(f"fatjet{sub2}{sub2}_preFSR_sub2_dE",
                             op.map(cont,lambda m : 2*(m[1][0].p4.E()-m[2].subJet2.p4.E())/(m[1][0].p4.E()+m[2].subJet2.p4.E())),
                             sel,
                             EquidistantBinning(200,-5.,5.0)))
    plots.append(Plot.make2D(f"fatjet{sub2}{sub2}_sub2_flavour",
                             [op.map(cont, lambda m: lambda_partonFlavour(m[1][1])), op.map(cont, lambda m: lambda_hadronFlavour(m[2].subJet2))],
                             sel,
                             [EquidistantBinning(6,0.,6.),EquidistantBinning(6,0.,6.)]))
    plots.append(Plot.make2D(f"fatjet{sub2}{sub2}_sub2_TF",
                             [op.map(cont,lambda m : m[1][0].p4.E()),op.map(cont,lambda m : m[2].subJet2.p4.E()-m[1][0].p4.E())],
                             sel,
                             [EquidistantBinning(*Egen_binning),EquidistantBinning(*dE_binning)],
                             xTitle="E^{{parton}}({sub2})",
                             yTitle=f"#Delta E = E^{{reco}}({sub2})-E^{{parton}}({sub2})"))


    # Compare subjets #
    plots.append(Plot.make1D(f"fatjet{sub1}{sub2}_dR",
                             op.map(cont,lambda m : op.deltaR(m[2].subJet1.p4,m[2].subJet2.p4)),
                             sel,
                             EquidistantBinning(100,0.,1.0)))
    plots.append(Plot.make2D(f"fatjet{sub1}{sub2}_E2D",
                             [op.map(cont,lambda m : m[2].subJet1.p4.E()),op.map(cont,lambda m : m[2].subJet2.p4.E())],
                             sel,
                             [EquidistantBinning(300,0.,3000.),EquidistantBinning(300,0.,3000.)]))
    plots.append(Plot.make2D(f"fatjet{sub1}{sub2}_dE2D",
                             [op.map(cont,lambda m : m[2].subJet1.p4.E()-m[0][0].p4.E()),op.map(cont,lambda m : m[2].subJet2.p4.E()-m[1][0].p4.E())],
                             sel,
                             [EquidistantBinning(1000,-1000.,1000.),EquidistantBinning(1000,-1000.,1000.)],
                             xTitle=f"subJet1 {sub1} E - gen {sub1} E",
                             yTitle=f"subJet2 {sub2} E - gen {sub2} E"))
    plots.append(Plot.make2D(f"fatjet{sub1}{sub2}_flavour2D",
                             [op.map(cont,lambda m : lambda_hadronFlavour(m[2].subJet1)),op.map(cont,lambda m : lambda_hadronFlavour(m[2].subJet2))],
                             sel,
                             [EquidistantBinning(6,0.,6.),EquidistantBinning(6,0.,6.)]))


    return plots


class TransferFunction(BaseNanoHHtobbWW,DataDrivenBackgroundHistogramsModule):
    def definePlots(self, t, noSel, sample=None, sampleCfg=None):
        noSel = super(TransferFunction,self).prepareObjects(t, noSel, sample, sampleCfg, 'DL')
        noSel = self.beforeJetselection(noSel)
        era = sampleCfg['era']
        plots = []

        # protection against data #
        if not self.is_MC:
            return plots

        # Hack the selections DF to make all weights at 1 #
        for selName,sel in noSel._fbe.selDFs.items():
            if 'nominal' in sel.weight.keys():
                noSel._fbe.selDFs[selName].weight['nominal'] = to.Const('float',1.)

        ##########################################################
        #                      Gen Leptons                       #
        ##########################################################
        #----- Gen particles -----#
        gen_e = op.select(
            t.GenPart,
            lambda p : op.AND(
                op.abs(p.pdgId) == 11,
                # eg some pair creation radiated form e or mu
                # in ttbar sample, less than 0.1% reco e have pt<~2
                op.rng_any(
                    p.ancestors,
                    lambda a: op.AND(
                        a.parent.idx >= 0, # Not initial state
                        op.OR(
                            op.abs(a.pdgId) == 23, # Z
                            op.abs(a.pdgId) == 24, # W+/W-
                        )
                    )
                ),
                p.statusFlags & ( 0x1 << 13 ), # isLastCopy
                op.OR(
                    p.statusFlags & ( 0x1 << 0), # isPrompt
                    p.statusFlags & ( 0x1 << 5), # isDirectPromptTauDecayProduct
                        # Prompt : tau from W decay
                        # Direct : we do not want a hadronic tau with D/else -> lepton
                ),
            )
        )
        gen_m = op.select(
            t.GenPart,
            lambda p : op.AND(
                op.abs(p.pdgId) == 13,
                # in ttbar sample reco muons have pt > 3
                op.rng_any(
                    p.ancestors,
                    lambda a: op.AND(
                        a.parent.idx >= 0, # Not initial state
                        op.OR(
                            op.abs(a.pdgId) == 23, # Z
                            op.abs(a.pdgId) == 24, # W+/W-
                        )
                    )
                ),
                p.statusFlags & ( 0x1 << 13 ), # isLastCopy
                op.OR(
                    p.statusFlags & ( 0x1 << 0), # isPrompt
                    p.statusFlags & ( 0x1 << 5), # isDirectPromptTauDecayProduct
                        # Prompt : tau from W decay
                        # Direct : we do not want a hadronic tau with D/else -> lepton
                ),
            )
        )
        gen_t = op.select(
            t.GenPart,
            lambda p : op.AND(
                op.abs(p.pdgId) == 15,
                op.rng_any(
                    p.ancestors,
                    lambda a: op.AND(
                        a.parent.idx >= 0, # Not initial state
                        op.OR(
                            op.abs(a.pdgId) == 23, # Z
                            op.abs(a.pdgId) == 24, # W+/W-
                        )
                    )
                ),
                p.statusFlags & ( 0x1 << 13 ), # isLastCopy
                p.statusFlags & ( 0x1 << 0), # isPrompt
                # We are only interested in hadronic decaying taus
                # All decays contain at least neural and charged pions
                # Check all the gen particles that have the tau has parent
                op.rng_any(
                    t.GenPart,
                    lambda s: op.AND(
                        s.idx > p.idx,          # must be successor, so after in idx
                        s.parent.idx == p.idx,  # parent is our tau candidate
                        op.OR(
                            op.abs(s.pdgId) == 111, # is a neutral pion
                            op.abs(s.pdgId) == 211, # is a charged pion
                        ),
                    ),
                ),
            )
        )

        ##########################################################
        #                      Gen Partons                       #
        ##########################################################
        gen_q = op.select(
            t.GenPart,
            lambda p : op.AND(
                op.OR(
                    op.abs(p.pdgId) == 1,
                    op.abs(p.pdgId) == 2,
                    op.abs(p.pdgId) == 3,
                    op.abs(p.pdgId) == 4,
                    op.abs(p.pdgId) == 5,
                ),
                p.statusFlags & ( 0x1 << 13 ), # isLastCopy
                p.statusFlags & ( 0x1 << 0),   # isPrompt
                op.OR(
                    # ISR cases #
                    op.AND(
                        # List potential origins of quark #
                        op.OR(
                            # Sometimes, a jet can emerge with mother idx = -1 (assume it's ISR), but not being the initial state (eg, the two first particles)
                            # But can also be the initial state itself (in which case, pt~=0, phi=0 and eta=~+/-20000), must not consider them
                            op.AND(
                                p.parent.idx == -1,
                                op.abs(p.eta) < 10,
                            ),
                            # Any parton (including b) or gluon (most likely) in initial state (idx=-1)
                            # Can give rise to a b quark that leads to a bjet
                            # Make sure to find the ancestor at idx=-1
                            op.AND(
                                p.parent.idx >= 0, # No initial state
                                # Find the initial state parton in the chain
                                op.rng_any(
                                    p.ancestors,
                                    lambda a: op.AND(
                                        a.parent.idx == -1,            # initial state
                                        op.OR(
                                            op.abs(a.pdgId) == 1,   # d quark
                                            op.abs(a.pdgId) == 2,   # u quark
                                            op.abs(a.pdgId) == 3,   # s quark
                                            op.abs(a.pdgId) == 4,   # c quark
                                            op.abs(a.pdgId) == 5,   # b quark
                                            op.abs(a.pdgId) == 21,  # gluon
                                        )
                                    )
                                ),
                            )
                        ),
                        # But also that the entire chain is just quarks or gluons
                        # Optimally should only be the same q for q(idx=-1)->q->...->b->...->b(last copy)
                        # But would need to check flavour by flavour
                        # Since one can have q -> q q' q'(bar) from pair creation (still in the q cone, so fine to dismiss q')
                        # But we want to keep q (beam) -> q' (because that is VBF typical)
                        # For gluons we only want cases where g(beam)->q->...−>q (last copy)
                        # cases such as g(idx=-1)(-> g(idx>=0) ->)...−>q (last copy)
                        # should be considered as the gluon only with FSR quarks (hence around the same cone)
                        # Make sure nothing else is in the ancestor chain
                        op.NOT(
                            op.rng_any(
                                p.ancestors,
                                lambda a: op.OR(
                                    op.AND(
                                        a.pdgId != p.pdgId,     # different particle than last copy candidate
                                        op.abs(a.eta) < 10,     # not in beam
                                    ),
                                )
                            )
                        )
                    ),
                    # top, Z, W and H decay cases #
                    op.AND(
                        # check if ancestor is top, Z, W and H
                        op.rng_any(
                            p.ancestors,
                            lambda a: op.AND(
                                a.parent.idx >= 0, # is not initial state
                                op.OR(
                                    op.abs(a.pdgId) == 6,  # top
                                    op.abs(a.pdgId) == 23, # Z
                                    op.abs(a.pdgId) == 24, # W
                                    op.abs(a.pdgId) == 25, # H
                                )
                            )
                        ),
                        # Need to remove pair-creations (same as above)
                        # Keep only cases with x -> t/Z/W/H -> q -> ... -> q (with no q' in the chain)
                        op.NOT(
                            op.rng_any(
                                p.ancestors,
                                lambda a: op.AND(
                                    # different quark than last copy candidate
                                    a.pdgId != p.pdgId,
                                    # Not any of the resonance we look at
                                    op.abs(a.pdgId) != 6,  # top
                                    op.abs(a.pdgId) != 23, # Z
                                    op.abs(a.pdgId) != 24, # W
                                    op.abs(a.pdgId) != 25, # H
                                    # Do not look at anything before t,Z,W,H (eg beam or previous resonances)
                                    # equivalent to : need to have either of those as ancestor (downstream)
                                    op.rng_any(
                                        a.ancestors,
                                        lambda aa: op.OR(
                                            op.abs(aa.pdgId) == 6,  # top
                                            op.abs(aa.pdgId) == 23, # Z
                                            op.abs(aa.pdgId) == 24, # W
                                            op.abs(aa.pdgId) == 25, # H
                                        ),
                                    ),
                                ),
                            )
                        )
                    ),
                ),
            )
        )
        # Select different flavour of quarks (for filling gen tree)
        gen_d = op.select(gen_q,lambda q: op.abs(q.pdgId)==1)
        gen_u = op.select(gen_q,lambda q: op.abs(q.pdgId)==2)
        gen_s = op.select(gen_q,lambda q: op.abs(q.pdgId)==3)
        gen_c = op.select(gen_q,lambda q: op.abs(q.pdgId)==4)
        gen_b = op.select(gen_q,lambda q: op.abs(q.pdgId)==5)

        gen_g = op.select(
            # Only ISR gluons that propagate or from an initial state quark (more rare)
            # FSR gluons (ie emitted from a quark) are usually colinear and end up in the same cone
            # -> Not a problem because we go back to the pre FSR parton to make TF
            # But for ISR gluons, they produce a big bunch of gluons in the final state
            # - some hard and some soft gluons
            # - More or less in a broad cone
            t.GenPart,
            lambda p : op.AND(
                p.pdgId == 21,
                p.statusFlags & (0x1 << 13), # isLastCopy
                p.statusFlags & (0x1 << 0),   # isPrompt
                # In ttbar sample, minimum reco jet energy is 15 GeV
                op.OR(
                    # Either a gluon form ISR (not initial state)
                    # Should be last copy because of above
                    op.AND(
                        p.parent.idx == -1,
                        op.abs(p.eta) < 10000,
                    ),
                    op.AND(
                        # Checking ancestor chain
                        # Make sure the initial state gluon is in the ancestor chain
                        p.parent.idx >= 0, # is not initial state
                        op.rng_any(
                            p.ancestors,
                            lambda a: op.AND(
                                a.parent.idx < 0,     # initial state
                                op.OR(
                                    op.abs(a.pdgId) == 1,   # b quark
                                    op.abs(a.pdgId) == 2,   # u quark
                                    op.abs(a.pdgId) == 3,   # s quark
                                    op.abs(a.pdgId) == 4,   # c quark
                                    op.abs(a.pdgId) == 5,   # b quark
                                    op.abs(a.pdgId) == 21,  # gluon
                                )
                            )
                        ),
                        # Also make sure the entire ancestor chain is only made up of gluons
                        # Except for initial state
                        op.NOT(
                            op.rng_any(
                                p.ancestors,
                                lambda a: op.AND(
                                    a.parent.idx >= 0, # Not initial state
                                    a.pdgId != 21,     # Is not gluon
                                )
                            )
                        )

                    )
                )
            )
        )
        ##########################################################
        #                      Overlapping                       #
        ##########################################################
        # Remove the objects in same tight DR cone (keeping the highest energy)
        # For example, e+e- creation off an e or mu, so only keep the highest energy one
        # Since we still recurse through the ancestor,
        #   - if it is off a muon, then the electron might be at reco level
        #   - if it is off an electron, then it does not matter since this electron will be the ancestor
        def remove_gen_overlap(cont):
            cont = op.select(
                cont,
                lambda o1: op.NOT(
                    op.rng_any(
                        cont,
                        lambda o2: op.AND(
                            o1.idx != o2.idx,              # Cannot be the same element
                            o1.pdgId == o2.pdgId,          # Has to be the same pdgid
                            op.deltaR(o1.p4,o2.p4) < 0.01, # Objects in same cone
                            o1.p4.E() < o2.p4.E(),         # Not taken if there is highest energy in the cone
                        )
                    )
                )
            )
            return cont


        gen_e = remove_gen_overlap(gen_e)
        gen_m = remove_gen_overlap(gen_m)
        gen_t = remove_gen_overlap(gen_t)

        gen_g = op.select(
            # Only keep a single final-state gluon from ISR within a DR=0.4 cone
            # Note : can happen that two gluons with parent idx = -1 (ISR gluons) are in the same cone
            # -> need to ask that they have a common ancestor
            # Select the highest pt of them all
            gen_g,
            lambda g1: op.NOT(
                op.rng_any(
                    gen_g,
                    lambda g2: op.AND(
                        g1.idx != g2.idx,               # Not the same object
                        op.deltaR(g1.p4,g2.p4) < 0.4,   # Objects in same cone
                        g1.pt < g2.pt,                  # Only consider highest pt gluon
                        op.rng_any(                     # Must have at least one common ancestor
                            g1.ancestors,
                            lambda a1: op.rng_any(
                                g2.ancestors,
                                lambda a2 : a1.idx == a2.idx,
                            )
                        )
                    )
                )
            )
        )

        # Clear overlaps between leptons and partons #
        # This is because electrons can be identified as jets
        # and any overlap will result in overestimation of the object energy
        def return_lambda_clean_overlap(containers,DR):
            return lambda t: op.NOT(
                op.OR(
                    *[
                        op.rng_any(
                            container,
                            lambda c: op.AND(
                                c.idx != t.idx,
                                op.deltaR(c.p4,t.p4) < DR,
                            )
                        )
                        for container in containers
                    ]
                )
            )
        # To avoid some segfaults (?) and precompute such values, use map
        # Still useful to put in the gen tree for debugging
        lambda_clean_e = return_lambda_clean_overlap([gen_e,gen_m,gen_t,gen_q,gen_g],0.4)
        # We ask the electron to be isolated from partons because a jet can leave ECAL deposits
        lambda_clean_m = lambda m: op.c_bool(True)
        # We always consider muons as the cleanest object
        # Since they only rely on tracking for energy measurement
        # no mismeasurement if coincides with a jet
        lambda_clean_t = return_lambda_clean_overlap([gen_e,gen_m,gen_t,gen_q,gen_g],0.4)
        lambda_clean_parton = return_lambda_clean_overlap([gen_e,gen_m,gen_t],0.4)
        # We do not clean between quarks and gluons (yet) because they might end up in AK8 jets
        # lambda is not used for gen used in the matching, because we still want to check "uniqueness" (see later)

        map_clean = op.map(
            t.GenPart,
            lambda p : op.multiSwitch(
                (
                    op.abs(p.pdgId) == 11,
                    lambda_clean_e(p),
                ),
                (
                    op.abs(p.pdgId) == 13,
                    lambda_clean_m(p),
                ),
                (
                    op.abs(p.pdgId) == 15,
                    lambda_clean_t(p),
                ),
                (
                    op.OR(
                        op.AND(
                            op.abs(p.pdgId) >= 1,
                            op.abs(p.pdgId) <= 5,
                        ),
                        op.abs(p.pdgId) == 21,
                    ),
                    lambda_clean_parton(p),
                ),
                # Default is False (should never be called for these particles anyway)
                op.c_bool(False),
            )
        )

        # Select partons (quarks and gluons) that are not within the same DR cone
        # Needs to be compared to gen_q/g before cleaning, otherwise if two quarks are in the same pair
        # but one is not clean, the other quark is considered resolved
        # Those will be treated later as AK8 jets
        lambda_resolved_parton = lambda p1 : op.NOT(
            op.OR(
                op.rng_any(
                    gen_q,
                    lambda p2: op.AND(
                        p1.idx != p2.idx,               # Different particle
                        op.deltaR(p1.p4,p2.p4) < 0.4,   # Same DR=4 cone
                    )
                ),
                op.rng_any(
                    gen_g,
                    lambda p2: op.AND(
                        p1.idx != p2.idx,               # Different particle
                        op.deltaR(p1.p4,p2.p4) < 0.4,   # Same DR=4 cone
                    )
                ),
            )
        )
        map_parton_resolved = op.map(
            t.GenPart,
            lambda p: op.multiSwitch(
                (
                    op.OR(
                        op.abs(p.pdgId) == 11,
                        op.abs(p.pdgId) == 13,
                        op.abs(p.pdgId) == 15,
                    ),
                    op.c_bool(True),
                ),
                (
                    op.OR(
                        op.AND(
                            op.abs(p.pdgId) >= 1,
                            op.abs(p.pdgId) <= 5,
                        ),
                        op.abs(p.pdgId) == 21,
                    ),
                    lambda_resolved_parton(p),
                ),
                op.c_bool(False),
            )
        )

        ##########################################################
        #                       Gen trees                        #
        ##########################################################

        # Define variable dicts #
        flavours = ['e','m','t','g','d','u','s','c','b']
        ranges   = [4,4,3,8,5,5,5,5,5]
        gens     = [gen_e,gen_m,gen_t,gen_g,gen_d,gen_u,gen_s,gen_c,gen_b]
        varDicts = [
            {
                'event': None,
                f'n_{flav}_gen': op.static_cast("UInt_t",op.rng_len(gen)),
            }
            for flav, gen in zip(flavours,gens)
        ]

        # Loop to fill trees #
        for flav, gen, varDict, rng in zip(flavours,gens,varDicts,ranges):
            for i in range(rng):
                varDict[f"{flav}_{i}_clean"] = op.switch(op.rng_len(gen)>i,op.static_cast("Int_t",map_clean[gen[i].idx]),op.c_int(-9999))
                varDict[f"{flav}_{i}_resolved"] = op.switch(op.rng_len(gen)>i,op.static_cast("Int_t",map_parton_resolved[gen[i].idx]),op.c_int(-9999))
                varDict[f"{flav}_{i}_idx"]   = op.switch(op.rng_len(gen)>i,op.static_cast("Int_t",gen[i].idx),op.c_int(-9999))
                varDict[f"{flav}_{i}_pdgId"] = op.switch(op.rng_len(gen)>i,gen[i].pdgId,op.c_int(-9999))
                varDict[f"{flav}_{i}_E"]     = op.switch(op.rng_len(gen)>i,gen[i].p4.E(),op.c_float(-9999))
                varDict[f"{flav}_{i}_pt"]    = op.switch(op.rng_len(gen)>i,gen[i].pt,op.c_float(-9999))
                varDict[f"{flav}_{i}_eta"]   = op.switch(op.rng_len(gen)>i,gen[i].eta,op.c_float(-9999))
                varDict[f"{flav}_{i}_phi"]   = op.switch(op.rng_len(gen)>i,gen[i].phi,op.c_float(-9999))
                varDict[f"{flav}_{i}_M"]     = op.switch(op.rng_len(gen)>i,gen[i].mass,op.c_float(-9999))
            # Add skim to plots #
            plots.append(Skim(f'gen_{flav}',varDict,noSel))

        ##########################################################
        #                       Matching                         #
        ##########################################################
        # Specific note on FSR
        # In rare cases a pre-FSR (taken from the ancestor) has a lower energy that post-FSR
        # For leptons this is typically from a W decay with an additional hard photon
        # In this case the understanding is that it comes from a W->mu nu gamma matrix element (and weird FSR)
        # -> the photon induces a loss of energy of the W itself as well, even though its mother is associated to the lepton
        # -> in such cases the energy that matches most the reco lepton is the post-FSR one (or at least the one right after satisfying E(pre-FSR)>E(post-FSR))
        # For quarks a similar issue is found when a hard gluon is radiated, which is interpreted as a t->W b g matrix element
        # To circumvent it, two choices :
        #  - just select the cases where E(pre-FSR)>E(post-FSR) (use 99% in case of truncation errors somewhere)
        #  - for leptons use t.GenDressedLepton (that should contain the additional photons)
        # here we choose the former, as it is easier, and would have to be done for quarks anyway (no Dressed counterparts)
        # BUT we only do that for the reco match, as we need the quark gen matches kept for fatjet prong identification
        # See discussion here : https://cms-talk.web.cern.ch/t/nanoaod-genparticles-gaining-energy-with-fsr/33136/15
        # Additionally, this is attached to very big DR between pre-FSR and post-FSR that do not match the reco
        # Note : this requirement is moved to later on in a lambda, to keep the problematic cases in the tree, but not in the plots

        #----- Gen Matching -----#
        gen_match_e = op.combine(
            rng = (t.GenPart,gen_e),
            # We do not apply the eta cut on electrons yet,
            pred = lambda p,e: op.AND(
                # p = all other gen particles (to find electron pre FSR)
                # e = last copy electron (Matching after possible FSR)
                # Mathing last copy e with the initial (before FSR) e
                p.pdgId == e.pdgId,                 # Must be the same particle
                p.pt >= 2,                          # Reco electrons have Pt >= 2 GeV in majotiry of events
                op.deltaR(e.p4,p.p4) < 0.2,         # Make sure no big kick during FSR
                p.statusFlags & (0x1 << 12),        # isFirstCopy
                op.OR(
                    p.statusFlags & (0x1 << 0),     # isPrompt
                    p.statusFlags & (0x1 << 5),     # isDirectPromptTauDecayProduct
                ),
                op.OR(                              # Must be right after tau/Z/W decay
                    op.abs(p.parent.pdgId) == 15,
                    op.abs(p.parent.pdgId) == 23,
                    op.abs(p.parent.pdgId) == 24,
                ),
                op.OR(                              # Must be an ancestor or same particle
                    op.rng_any(
                        e.ancestors,
                        lambda a: a.idx == p.idx,
                    ),
                    e.idx == p.idx,
                ),
            ),
            samePred = lambda p,e: p.idx <= e.idx,  # Precursor always has idx before the last copy
        )
        gen_match_m = op.combine(
            rng = (t.GenPart,gen_m),
            pred = lambda p,m: op.AND(
                # p = all other gen particles (to find muon pre FSR)
                # m = last copy muon (after possible FSR)
                # Mathing last copy m with the initial (before FSR) m
                p.pdgId == m.pdgId,                 # Must be the same particle
                p.pt >= 3,                          # Reco muons have Pt >= 3 GeV
                op.deltaR(m.p4,p.p4) < 0.2,         # Make sure no big kick during FSR
                p.statusFlags & (0x1 << 12),        # isFirstCopy
                op.OR(
                    p.statusFlags & ( 0x1 << 0),    # isPrompt
                    p.statusFlags & ( 0x1 << 5),    # isDirectPromptTauDecayProduct
                ),
                op.OR(                              # Must be right after tau/Z/W decay
                    op.abs(p.parent.pdgId) == 15,
                    op.abs(p.parent.pdgId) == 23,
                    op.abs(p.parent.pdgId) == 24,
                ),
                op.OR(                              # Must be an ancestor or same particle
                    op.rng_any(
                        m.ancestors,
                        lambda a: a.idx == p.idx,
                    ),
                    m.idx == p.idx,
                ),
            ),
            samePred = lambda p,m: p.idx <= m.idx,  # Precursor always has idx before the last copy
        )
        gen_match_t = op.combine(
            rng = (t.GenPart,gen_t),
            pred = lambda p,t: op.AND(
                # p = all other gen particles (to find tau pre FSR)
                # t = last copy tau (after possible FSR)
                # Mathing last copy t with the initial (before FSR) t
                p.pdgId == t.pdgId,                 # Must be the same particle
                p.pt >= 10,                         # Reco taus have Pt >= 18 GeV
                op.deltaR(t.p4,p.p4) < 0.2,         # Make sure no big kick during FSR
                p.statusFlags & (0x1 << 0),         # isPrompt
                p.statusFlags & (0x1 << 12),        # isFirstCopy
                op.OR(                              # Must be right after tau/Z/W decay
                    op.abs(p.parent.pdgId) == 23,
                    op.abs(p.parent.pdgId) == 24,
                ),
                op.OR(                              # Must be an ancestor or same particle
                    op.rng_any(t.ancestors,
                        lambda a: a.idx == p.idx,
                    ),
                    t.idx == p.idx,
                ),
            ),
            samePred = lambda p,t: p.idx <= t.idx,  # Precursor always has idx before the last copy
        )


        #----- Reco Matching -----#
        reco_match_e = op.combine(
            rng = (gen_match_e,t.Electron),
            pred = lambda e,r: op.AND(
                # e = gen match
                #   e[0] : pre-FSR
                #   e[1] : post-FSR
                #e[0].p4.E() >= 0.99*e[1].p4.E(),          # Particle cannot lose energy through FSR
                # Minimum energy for electrons
                e[0].p4.E() > 5,
                # Matching the last copy e to the reco e
                op.deltaR(e[1].p4,r.p4) < 0.2, # Reco should match post-FSR
                # Kill outliers for deltaE > e_gen / 2 (< 0.01% of events in DY)
                (r.p4.E()-e[0].p4.E()) < e[0].p4.E() / 2,
                # Sometimes a hard photon is emitted, and the reco electron
                # gen index is set on the photon, need to handle that
                op.switch(
                    t.GenPart[r.genPart.idx].pdgId == 22, # reco electron matched to gen photon
                    # True case : matched to photon
                    # Sometimes the post-FSR electron is hard, and photon is soft, need distinction
                    op.switch(
                        # Check is most of the pre-FSR energy is still in post-FSR
                        e[1].p4.E() > 0.5 * e[0].p4.E(),
                        # True case : only base the matching on the electron
                        # But need to make sure the photon comes from the pre-FSR considered here
                        # And that the electron is still on the pre-FSR electron axis
                        t.GenPart[r.genPart.idx].parent.idx == e[0].idx,
                        # False case : there is a photon, but the electron is still hard
                        # Make sure the electron and photon have same parent
                        t.GenPart[r.genPart.idx].parent.idx == e[1].parent.idx,
                    ),
                    # False case : reco matched to electron
                    r.genPart.idx == e[1].idx, # reco electron matched to gen electron
                ),

            )
        )
        reco_match_m = op.combine(
            rng = (gen_match_m,t.Muon),
            pred = lambda m,r: op.AND(
                # m = gen match
                #   m[0] : pre-FSR
                #   m[1] : post-FSR
                #m[0].p4.E() >= 0.99*m[1].p4.E(),          # Particle cannot lose energy through FSR
                # Minimum energy for muons
                m[0].p4.E() > 5,
                # Matching the last copy e to the reco e
                r.genPart.idx == m[1].idx,
                op.deltaR(m[1].p4,r.p4) < 0.2, # Reco should match post-FSR
            )
        )
        reco_match_t = op.combine(
            rng = (gen_match_t,t.Jet),
            pred = lambda t,r: op.AND(
                # t = gen match
                #   t[0] : pre-FSR
                #   t[1] : post-FSR
                # Matching the last copy hadronic tau to the jet genjet
                # Since jets are ordered in Pt, if several jets are within DR cone
                # then the first one (leading) is taken, which is what we want
                #t[0].p4.E() >= 0.99*t[1].p4.E(),          # Particle cannot lose energy through FSR
                # Minimum energy fortaus
                t[0].p4.E() > 10,
                # If genjet available use it, otherwise default to reco jet
                op.switch(
                    op.AND(
                        r.genJet.idx >= 0, # jet has a genjet
                        op.deltaR(r.p4,r.genJet.p4) < 0.4, # avoid weird cases
                    ),
                    op.deltaR(t[1].p4,r.genJet.p4) < 0.2,
                    op.deltaR(t[1].p4,r.p4) < 0.2,
                )
            )
        )
        # The method for jets goes in three steps as below
        # pre-FSR quark <−history (1)− post-fSR quark −DR-matching (2)−> genJet −index (3)−>recoJet
        # (1) go back through history to find pre-FSR quark
        # (2) use the post-FSR quark and genJets, match them using DR selection
        # (3) from the index of the genJet, we can get the recoJet
        # Notes :
        # - DR matching is done assuming most of the energy is carried out by the post-FSR quark
        #   there are usually a bunch of hadrons and partons produced between pre and post FSR
        #   This is why we do DR matching with post-FSR, because direction may have changed due to these emissions
        # - What sometimes happen is that more than half the pre-FSR energy is radiated out
        #   In this case the post-FSR quark is relatively soft, and much in line with the jet anymore
        #   So we need a logic to fetch the hardest radiation from that pre-FSR quark, and use that for DR-matching
        # - The AK4 jets use R = 0.4 (so diameter of the jet is 0.8)
        #   From what is observed, the hardest component in the jet (most likely the quark) is mostly within DR<0.2
        #   This also makes sense, as the jet is usually aligned with the largest deposits,
        #   So the hardest component is unlikely to be in 0.2<DR<0.4 from the jet axis
        #   -> DR for jet matchin is DR<0.2
        # - For tau the genJet and jet matching is similar, but since we use the hadronic tau, don't expect big radiations
        # - In some very weird cases, the reco and gen jets are way off, add a cut to remove those weird cases

        #----- Gen Matching -----#
        gen_match_q = op.combine(
            rng = (t.GenPart,gen_q),
            pred = lambda p,q: op.AND(
                # p = all other gen particles (to find pre FSR)
                # q = last copy quark (after possible FSR)
                # Mathing last copy with the initial (before FSR)
                p.pdgId == q.pdgId,                 # Must be the same particle
                p.pt >= 10,                         # Reco (gen) jets have Pt >= (10)15 GeV
                p.statusFlags & (0x1 << 0),         # isPrompt
                # we do not require isFirstCopy, because quarks with parent -1 (ISR) can be lastCopy without firstCopy
                op.OR(
                    # Must be right after Z,H,W,top
                    op.abs(p.parent.pdgId) == 6,
                    op.abs(p.parent.pdgId) == 23,
                    op.abs(p.parent.pdgId) == 24,
                    op.abs(p.parent.pdgId) == 25,
                    # Sometimes the beam can generate quarks directly, in which case quark has parent in the beam
                    # Also can be a gluon generated from the beam, that then generates quarks
                    # Or from ISR (not in the beam, but idx = -1)
                    # In all cases the chain is the same quark until the beam or gluon
                    op.AND(
                        op.OR(
                            # Can be a gluon/quark from the beam (in the beam so eta ~>20000), that produces quarks, in which case get quark right after
                            # Must get the quark right after as pre FSR
                            op.AND(
                                p.parent.idx >= 0,              # particle is not ISR
                                p.parent.parent.idx < 0,        # gluon/quark has no parent
                                op.abs(p.parent.eta) > 10000,   # and from the beam (initial state)
                            ),
                            # Can be a quark in ISR (idx = -1) but not initial state, in which case get that quark as pre FSR
                            op.AND(
                                p.parent.idx < 0,           # has no parent (ISR)
                                op.abs(p.eta) < 10,         # not the beam
                            ),
                        ),
                        # Make sure the whole chain of ancestors is the same type of particles (except beam)
                        # To avoid for example 5 (initial state in beam) -> 6 (top) -> 5 (bottom from top decay pre FSR) -> ... -> 5 (last copy)
                        # However should veto cases with q1(idx=-1, not beam) -> q1 + q2 (ie, FSR radiation of quarks) for q2
                        op.NOT(
                            op.rng_any(
                                q.ancestors,
                                lambda a: op.AND(
                                    op.abs(a.eta) < 10, # Not the beam/initial state
                                    a.pdgId != q.pdgId,
                                )
                            )
                        ),
                    ),
                ),
                op.OR(                              # Must be an ancestor or same particle
                    op.rng_any(q.ancestors,
                        lambda a: a.idx == p.idx,
                    ),
                    q.idx == p.idx,
                ),
            ),
            samePred = lambda p,q: p.idx <= q.idx,  # Precursor always has idx before the last copy
        )

        gen_match_g = op.combine(
            rng = (t.GenPart,gen_g),
            pred = lambda p,g: op.AND(
                # p = all other gen particles (to find initial ISR gluon)
                # g = last copy gluon
                # r = reco jet
                # Mathing last copy with the initial (before FSR)
                p.pdgId == g.pdgId,                 # Must be the same particle
                p.pt >= 10,                         # Reco (gen) jets have Pt >= (10)15 GeV
                p.statusFlags & (0x1 << 0),         # isPrompt
                # Must be right after a quark or gluon that is initial state
                # Or an additional gluon that is not the initial state but is ISR
                op.OR(
                    # If the quark/gluon was the initial state (in the beam axis), take the next one
                    op.AND(
                        p.parent.idx >= 0, # is not initial state
                        p.parent.parent.idx < 0, # the parent has no parent
                        op.abs(p.parent.eta) > 10000, # the parent is the beam (initial state)
                    ),
                    # The gluon is ISR, but not the beam (initial state)
                    op.AND(
                        p.parent.idx < 0, # the parent has no parent (is ISR)
                        op.abs(p.eta) < 10, # ISR but not beam initial state
                    ),
                ),
                # Must be an ancestor or same particle
                op.OR(
                    op.rng_any(g.ancestors,
                               lambda a: a.idx == p.idx,
                               ),
                    g.idx == p.idx,
                ),
            ),
            samePred = lambda p,g: p.idx <= g.idx,  # Precursor always has idx before the last copy
        )
        # Split quarks by flavour (for tree) #
        gen_match_b = op.select(gen_match_q,lambda q : op.abs(q[1].pdgId)==5)
        gen_match_l = op.select(gen_match_q,lambda q : op.abs(q[1].pdgId)!=5)

        #----- Reco Matching -----#
        reco_match_q = op.combine(
            rng = (gen_match_q,t.Jet),
            pred = lambda q,r: op.AND(
                # q = gen match
                #   q[0] : pre-FSR
                #   q[1] : post-FSR
                # Matching the last copy q to the jet genjet
                # Since jets are ordered in Pt, if several jets are within DR cone
                # then the first one (leading) is taken, which is what we want
                #q[0].p4.E() >= 0.99 * q[1].p4.E(),        # Particle cannot lose energy through FSR
                # Use the parton flavour of the reco jet
                r.partonFlavour == q[0].pdgId,
                # Need to deal with cases where a lot of energy is radiated off the quark
                op.switch(
                    q[1].p4.E() >= 0.5 * q[0].p4.E(),
                    # True case : most of the energy is still in the quark (>50% of pre-FSR)
                    # Also handles the cases where post-FSR = pre-FSR (because E1 >= 0.5 * E2 if E1=E2)
                    # Logic of false case below would also work on this case,
                    # but lots of computation can be avoided if indeed post-FSR is hardest
                    op.AND(
                        op.switch(
                            # If genjet available use it, otherwise default to reco jet
                            op.AND(
                                r.genJet.idx >= 0, # jet has a genjet
                                op.deltaR(r.p4,r.genJet.p4) < 0.4, # avoid weird cases
                            ),
                            op.deltaR(q[1].p4,r.genJet.p4) < 0.2,
                            op.deltaR(q[1].p4,r.p4) < 0.2,
                        ),
                        op.deltaR(q[1].p4,q[0].p4) < 0.6,     # Still make sure no large kick
                    ),
                    # False case : need to find the hardest descendant of the pre-FSR quark (that has isLastCopy())
                    # Ideally would want to use an op.select of all potential descendants then op.sort to get the leading
                    # ... but would not work here ...
                    # Instead make two loops (with op.rng_any) of all candidates, the second loop to check that we only look at the leading one
                    # Limit ourselves where hardest is not too far (sometimes meson goes far, DR>1, problematic cases)
                    # -> only gluons and other quarks (veto isDecayedLeptonHadron flag, eg 111 or 531 hadrons) and DR<0.4 of post-FSR quark
                    # In case no other harder constituent that is not a hadron : check below
                    op.switch(
                        # Check if there are other daughters of the pre-FSR that are not the post-FSR
                        op.NOT(
                            op.rng_any(
                                t.GenPart,
                                lambda q1: op.AND(
                                    # Cannot be our post-FSR quark
                                    q1.idx != q[1].idx,
                                    # Descendant, so after in index (not needed, but can make the loop faster)
                                    q1.idx > q[0].idx,
                                    # Prompt and not a decayed hadron
                                    q1.statusFlags & (0x1 << 0),
                                    op.NOT(q1.statusFlags & (0x1 << 1)), # not isDecayedLeptonHadron
                                    # Hardest still around pre-FSR quark
                                    op.deltaR(q1.p4,q[0].p4) < 0.6,
                                    # Make sure it also come from the pre-FSR
                                    op.rng_any(
                                        q1.ancestors,
                                        lambda a : a.idx == q[0].idx
                                    ),
                                )
                            ),
                        ),
                        # True case : only post-FSR quark detected (except for decayed hadrons)
                        # - post-FSR is in DR<0.6 cone of pre-FSR, loosen matching to DR<0.4 with reco jets
                        # - else, only look for jet in pre-FSR cone
                        op.switch(
                            op.deltaR(q[1].p4,q[0].p4)<0.6,
                            # true case : post-FSR is within pre-FSR cone : use loosen DR cone
                            op.switch(
                                # If genjet available use it, otherwise default to reco jet
                                op.AND(
                                    r.genJet.idx >= 0, # jet has a genjet
                                    op.deltaR(r.p4,r.genJet.p4) < 0.4, # avoid weird cases
                                ),
                                op.deltaR(q[1].p4,r.genJet.p4) < 0.4,
                                op.deltaR(q[1].p4,r.p4) < 0.4,
                            ),
                            # false case : look at pre-FSR for matching
                            op.switch(
                                # If genjet available use it, otherwise default to reco jet
                                op.AND(
                                    r.genJet.idx >= 0, # jet has a genjet
                                    op.deltaR(r.p4,r.genJet.p4) < 0.4, # avoid weird cases
                                ),
                                op.deltaR(q[0].p4,r.genJet.p4) < 0.2,
                                op.deltaR(q[0].p4,r.p4) < 0.2,
                            ),
                        ),
                        # False case : there is harder component
                        op.rng_any(
                            # first inner loop to get all the candidate descendants of the pre-FSR quark
                            t.GenPart,
                            lambda q1: op.AND(
                                # Same requirements as before
                                q1.idx > q[0].idx,
                                q1.statusFlags & (0x1 << 0),
                                op.NOT(q1.statusFlags & (0x1 << 1)),
                                op.deltaR(q1.p4,q[0].p4) < 0.6,
                                op.rng_any(
                                    q1.ancestors,
                                    lambda a : a.idx == q[0].idx
                                ),
                                # Need to check all other potential candidates to make sure none has higher energy than q1
                                op.NOT(
                                    op.rng_any(
                                        t.GenPart,
                                        lambda q2: op.AND(
                                            # Obviously not same particle
                                            q2.idx != q1.idx,
                                            # Same requirements as q1
                                            q2.idx > q[0].idx,
                                            q2.statusFlags & (0x1 << 0),
                                            op.NOT(q2.statusFlags & (0x1 << 1)),
                                            op.deltaR(q2.p4,q[0].p4) < 0.6,
                                            op.rng_any(
                                                q2.ancestors,
                                                lambda a : a.idx == q[0].idx
                                            ),
                                            # Higher energy than q1
                                            q2.p4.E() > q1.p4.E()
                                        ),
                                    )
                                ),
                                # Finally check the jet matching
                                op.switch(
                                    # If genjet available use it, otherwise default to reco jet
                                    op.AND(
                                        r.genJet.idx >= 0, # jet has a genjet
                                        op.deltaR(r.p4,r.genJet.p4) < 0.4, # avoid weird cases
                                    ),
                                    op.deltaR(q1.p4,r.genJet.p4) < 0.2,
                                    op.deltaR(q1.p4,r.p4) < 0.2,
                                )
                            )
                        ),
                    ),
                ),
            ),
        )

        reco_match_g = op.combine(
            rng = (gen_match_g,t.Jet),
            pred = lambda g,r: op.AND(
                # g = gen match
                #   g[0] : pre-FSR
                #   g[1] : post-FSR
                # Matching the last copy g to the jet genjet
                # Since jets are ordered in Pt, if several jets are within DR cone
                # then the first one (leading) is taken, which is what we want
                #g[0].p4.E() >= 0.99 * g[1].p4.E(),        # Particle cannot lose energy through FSR
                # Use the parton flavour of the reco jet
                r.partonFlavour == g[0].pdgId,
                # No need for the same hoops as quarks, since we only have selected hardest gen_g
                op.switch(
                    # If genjet available use it, otherwise default to reco jet
                    op.AND(
                        r.genJet.idx >= 0, # jet has a genjet
                        op.deltaR(r.p4,r.genJet.p4) < 0.4, # avoid weird cases
                    ),
                    op.deltaR(g[1].p4,r.genJet.p4) < 0.2,
                    op.deltaR(g[1].p4,r.p4) < 0.2,
                )
            )
        )


        # Must resolve the cases where a gen particle is matched to two reco jets (maybe more?)
        # In which case only select the case with smallest DR
        # Can be a problem for electrons when a photon creates a second reco electron
        # Has to be done for hadronic tau as well though
        # We want to this before checking for uniqueness
        def remove_multiple_reco_matching(match):
            return op.select(
                match,
                lambda m1: op.NOT(
                    op.rng_any(
                        match,
                        lambda m2: op.AND(
                            m1[0][1].idx == m2[0][1].idx,     # same postFSR particle
                            m1[1].idx != m2[1].idx,           # different reco jet
                            op.deltaR(m1[0][1].p4,m1[1].p4) > op.deltaR(m2[0][1].p4,m2[1].p4)
                            # We only want to keep the matched triplet
                            # that has smallest DR between postFSR gen and reco jet
                        )
                    )
                )
            )
        # Must resolve the cases where a pre-FSR gen is matched to several post-FSR
        # This happens mostly for gluons, that can radiate several FSR gluons
        # This means we have the initial energy split into several jets
        # Can also happen for quarks, but also pair creation (eg b -> b (b + bbar))
        # If DR < 0.4 : remove the softest
        # -> Should be treated as AK8
        # If DR >= 0.4 : keep them
        # This should be done after check for uniqueness
        # (fragment can coincide with other jets and overestimate its energy)
        def remove_multiple_preFSR_matching(match):
            return op.select(
                match,
                lambda m1: op.NOT(
                    op.rng_any(
                        match,
                        lambda m2: op.AND(
                            m1[0][0].idx == m2[0][0].idx,               # Same pre-FSR
                            m1[0][1].idx != m2[0][1].idx,               # Different post-FSR
                            op.deltaR(m1[0][1].p4,m2[0][1].p4) < 0.4,   # DR<0.4 between post-FSR
                            m2[0][1].p4.E() > m1[0][1].p4.E(),          # m2 has higher post-FSR energy
                        )
                    )
                )
            )

        reco_match_m = remove_multiple_reco_matching(reco_match_m)
        reco_match_e = remove_multiple_reco_matching(reco_match_e)
        reco_match_t = remove_multiple_reco_matching(reco_match_t)
        reco_match_q = remove_multiple_reco_matching(reco_match_q)
        reco_match_g = remove_multiple_reco_matching(reco_match_g)

        # Split quark flavours #
        reco_match_b = op.select(reco_match_q,lambda q : op.abs(q[0][1].pdgId)==5)
        reco_match_l = op.select(reco_match_q,lambda q : op.abs(q[0][1].pdgId)!=5)
            # contains lightjets (u,d,s) and c (because we do not care about c tagging)

        # Remove cases where matched to the same reco object
        # First we count for each reco object how many matches are attached
        # (with a map to avoid segfault and recomputing several time)
        # Reco unique = exactly 1 matches attached
        map_electron_unique_count = op.map(
            t.Electron,
            lambda e: op.rng_count(
                rng = reco_match_e,
                pred = lambda match: match[1].idx == e.idx,
            )
        )
        map_muon_unique_count = op.map(
            t.Muon,
            lambda m: op.rng_count(
                rng = reco_match_m,
                pred = lambda match: match[1].idx == m.idx,
            )
        )
        map_jet_unique_count = op.map(
            t.Jet,
            lambda j:
                op.rng_count(
                    rng = reco_match_t,
                    pred = lambda match: match[1].idx == j.idx,
                ) + \
                op.rng_count(
                    rng = reco_match_q,
                    pred = lambda match: match[1].idx == j.idx,
                ) + \
                op.rng_count(
                    rng = reco_match_g,
                    pred = lambda match: match[1].idx == j.idx,
                )
        )

        lambda_match_unique_e   = lambda match: map_electron_unique_count[match[1].idx] == 1
        lambda_match_unique_m   = lambda match: map_muon_unique_count[match[1].idx] == 1
        lambda_match_unique_jet = lambda match: map_jet_unique_count[match[1].idx] == 1


        reco_match_e_final = op.select(reco_match_e,lambda_match_unique_e)
        reco_match_m_final = op.select(reco_match_m,lambda_match_unique_m)
        reco_match_t_final = op.select(reco_match_t,lambda_match_unique_jet)
        reco_match_b_final = op.select(reco_match_b,lambda_match_unique_jet)
        reco_match_l_final = op.select(reco_match_l,lambda_match_unique_jet)
        reco_match_g_final = op.select(reco_match_g,lambda_match_unique_jet)
        # we can apply lambda_match_unique_parton to match_b and match_l separately
        # because the lambda looks at match_q = match_b+match_l content (+match_g)

        # Check for fragmentation #
        reco_match_b_final = remove_multiple_preFSR_matching(reco_match_b_final)
        reco_match_l_final = remove_multiple_preFSR_matching(reco_match_l_final)
        reco_match_g_final = remove_multiple_preFSR_matching(reco_match_g_final)

        # Apply the cleaning lambda #
        lambda_match_clean = lambda match : map_clean[match[0][1].idx]

        reco_match_e_final = op.select(reco_match_e_final,lambda_match_clean)
        reco_match_m_final = op.select(reco_match_m_final,lambda_match_clean)
        reco_match_t_final = op.select(reco_match_t_final,lambda_match_clean)
        reco_match_b_final = op.select(reco_match_b_final,lambda_match_clean)
        reco_match_l_final = op.select(reco_match_l_final,lambda_match_clean)
        reco_match_g_final = op.select(reco_match_g_final,lambda_match_clean)


        # Also veto non resolved partons
        # Could use the previous lambda_resolved_parton
        # But this will not consider cases where a jet is not perfectly aligned with the quark,
        # and may miss the DR resolved cone around other partons
        # instead compare the reco jet, with all pre-FSR quarks and gluons
        # But only in cases where the pre-FSR quark pT >= 10 GeV, otherwise little impact
        map_jet_resolved_count = op.map(
            t.Jet,
            lambda j:
                op.rng_count(
                    rng = gen_match_q,
                    pred = lambda q : op.AND(
                        q[0].pt >= 10,                   # pre-FSR pt >= 10 GeV
                        op.deltaR(j.p4,q[1].p4) < 0.4,   # DR(reco jet, post-FSR)
                    )
                ) + \
                op.rng_count(
                    rng = gen_match_g,
                    pred = lambda g : op.AND(
                        g[0].pt >= 10,                   # pre-FSR pt >= 10 GeV
                        op.deltaR(j.p4,g[1].p4) < 0.4,   # DR(reco jet, post-FSR)
                    )
                )
        )

        lambda_match_resolved = lambda m: map_jet_resolved_count[m[1].idx] == 1
        reco_match_b_final = op.select(reco_match_b_final,lambda_match_resolved)
        reco_match_l_final = op.select(reco_match_l_final,lambda_match_resolved)
        reco_match_g_final = op.select(reco_match_g_final,lambda_match_resolved)

        # Remove cases where pre-FSR produces particles with higher energies
        # This can be for example the post-FSR having higher energy
        # or in some cases a meson or gluon is emitted with more energy
        # So make sure that sum of energies of all direct descendants is less or equal to pre-FSR
        map_sum_FSR_energy = op.map(
            t.GenPart,
            lambda p: op.rng_sum(
                t.GenPart,
                lambda d: op.switch(
                    # if particle is direct descendant of pre-FSR, include energy, otherwise a 0
                    d.parent.idx == p.idx,
                    d.p4.E(),
                    op.c_float(0.)
                ),
            )
        )
        lambda_match_check_FSR = lambda match: map_sum_FSR_energy[match[0][0].idx] <= 1.01 * match[0][0].p4.E() # allowing some truncation errors

        reco_match_e_final = op.select(reco_match_e_final,lambda_match_check_FSR)
        reco_match_m_final = op.select(reco_match_m_final,lambda_match_check_FSR)
        reco_match_t_final = op.select(reco_match_t_final,lambda_match_check_FSR)
        reco_match_b_final = op.select(reco_match_b_final,lambda_match_check_FSR)
        reco_match_l_final = op.select(reco_match_l_final,lambda_match_check_FSR)
        reco_match_g_final = op.select(reco_match_g_final,lambda_match_check_FSR)


        # Select eta regions #
        lambda_central_e = lambda e : op.abs(e[1].eta) <= 2.5
        lambda_central_m = lambda m : op.abs(m[1].eta) <= 2.4
        # We do not select forward leptons
        reco_match_e_central = op.select(reco_match_e_final,lambda_central_e)
        reco_match_m_central = op.select(reco_match_m_final,lambda_central_e)

        lambda_central_jet = lambda q : op.abs(q[1].eta) <= 2.4
        lambda_forward_jet = lambda q : op.in_range(2.4,op.abs(q[1].eta),4.7)

        reco_match_b_central = op.select(reco_match_b_final,lambda_central_jet)
        reco_match_l_central = op.select(reco_match_l_final,lambda_central_jet)
        reco_match_g_central = op.select(reco_match_g_final,lambda_central_jet)

        reco_match_b_forward = op.select(reco_match_b_final,lambda_forward_jet)
        reco_match_l_forward = op.select(reco_match_l_final,lambda_forward_jet)
        reco_match_g_forward = op.select(reco_match_g_final,lambda_forward_jet)

        ##########################################################
        #                    Matching trees                      #
        ##########################################################


        # Define variable dicts #
        varsElMatch = {
            'event': None,
            'n_e_gen_match': op.static_cast("UInt_t",op.rng_len(gen_match_e)),
            'n_e_reco_match': op.static_cast("UInt_t",op.rng_len(reco_match_e)),
            'n_e_reco_final': op.static_cast("UInt_t",op.rng_len(reco_match_e_final)),
            'n_e_reco_central': op.static_cast("UInt_t",op.rng_len(reco_match_e_central)),
        }
        varsMuMatch = {
            'event': None,
            'n_m_gen_match': op.static_cast("UInt_t",op.rng_len(gen_match_m)),
            'n_m_reco_match': op.static_cast("UInt_t",op.rng_len(reco_match_m)),
            'n_m_reco_final': op.static_cast("UInt_t",op.rng_len(reco_match_m_final)),
            'n_m_reco_central': op.static_cast("UInt_t",op.rng_len(reco_match_m_central)),
        }
        varsTauMatch = {
            'event': None,
            'n_t_gen_match': op.static_cast("UInt_t",op.rng_len(gen_match_t)),
            'n_t_reco_match': op.static_cast("UInt_t",op.rng_len(reco_match_t)),
            'n_t_reco_final': op.static_cast("UInt_t",op.rng_len(reco_match_t_final)),
        }
        varsBMatch = {
            'event': None,
            'n_b_gen_match': op.static_cast("UInt_t",op.rng_len(gen_match_b)),
            'n_b_reco_match': op.static_cast("UInt_t",op.rng_len(reco_match_b)),
            'n_b_reco_final': op.static_cast("UInt_t",op.rng_len(reco_match_b_final)),
            'n_b_reco_central': op.static_cast("UInt_t",op.rng_len(reco_match_b_central)),
            'n_b_reco_forward': op.static_cast("UInt_t",op.rng_len(reco_match_b_forward)),
        }
        varsLMatch = {
            'event': None,
            'n_l_gen_match': op.static_cast("UInt_t",op.rng_len(gen_match_l)),
            'n_l_reco_match': op.static_cast("UInt_t",op.rng_len(reco_match_l)),
            'n_l_reco_final': op.static_cast("UInt_t",op.rng_len(reco_match_l_final)),
            'n_l_reco_central': op.static_cast("UInt_t",op.rng_len(reco_match_l_central)),
            'n_l_reco_forward': op.static_cast("UInt_t",op.rng_len(reco_match_l_forward)),
        }
        varsGMatch = {
            'event': None,
            'n_g_gen_match': op.static_cast("UInt_t",op.rng_len(gen_match_g)),
            'n_g_reco_match': op.static_cast("UInt_t",op.rng_len(reco_match_g)),
            'n_g_reco_final': op.static_cast("UInt_t",op.rng_len(reco_match_g_final)),
            'n_g_reco_central': op.static_cast("UInt_t",op.rng_len(reco_match_g_central)),
            'n_g_reco_forward': op.static_cast("UInt_t",op.rng_len(reco_match_g_forward)),
        }

        # Loop to fill trees #
        flavours        = ['e','m','t','b','l','g']
        matches         = [reco_match_e,reco_match_m,reco_match_t,reco_match_b,reco_match_l,reco_match_g]
        varDicts        = [varsElMatch,varsMuMatch,varsTauMatch,varsBMatch,varsLMatch,varsGMatch]
        ranges          = [4,4,4,6,6,8]
        for flav, match, varDict, rng in zip(flavours,matches,varDicts,ranges):
            for i in range(rng):
                # Matching flags #
                varDict[f"{flav}_{i}_match_clean"]    = op.switch(op.rng_len(match)>i,op.static_cast("Int_t",lambda_match_clean(match[i])),op.c_int(-9999))
                varDict[f"{flav}_{i}_match_DR"]       = op.switch(op.rng_len(match)>i,op.deltaR(match[i][0][1].p4,match[i][1].p4),op.c_float(-9999))
                varDict[f"{flav}_{i}_match_relE"]     = op.switch(op.rng_len(match)>i,match[i][0][1].p4.E()/match[i][0][0].p4.E(),op.c_float(-9999))
                varDict[f"{flav}_{i}_match_sum_FSR"]  = op.switch(op.rng_len(match)>i,map_sum_FSR_energy[match[i][0][0].idx],op.c_float(-9999))
                varDict[f"{flav}_{i}_match_check_FSR"]= op.switch(op.rng_len(match)>i,op.static_cast("Int_t",lambda_match_check_FSR(match[i])),op.c_float(-9999))
                if flav == 'e':
                    varDict[f"{flav}_{i}_match_unique"]   = op.switch(op.rng_len(match)>i,op.static_cast("Int_t",lambda_match_unique_e(match[i])),op.c_int(-9999))
                    varDict[f"{flav}_{i}_match_unique_count"]    = op.switch(op.rng_len(match)>i,op.static_cast("Int_t",map_electron_unique_count[match[i][1].idx]),op.c_int(-9999))
                    varDict[f"{flav}_{i}_match_central"]  = op.switch(op.rng_len(match)>i,op.static_cast("Int_t",lambda_central_e(match[i])),op.c_int(-9999))
                if flav == 'm':
                    varDict[f"{flav}_{i}_match_unique"]   = op.switch(op.rng_len(match)>i,op.static_cast("Int_t",lambda_match_unique_m(match[i])),op.c_int(-9999))
                    varDict[f"{flav}_{i}_match_unique_count"]    = op.switch(op.rng_len(match)>i,op.static_cast("Int_t",map_muon_unique_count[match[i][1].idx]),op.c_int(-9999))
                    varDict[f"{flav}_{i}_match_central"]  = op.switch(op.rng_len(match)>i,op.static_cast("Int_t",lambda_central_m(match[i])),op.c_int(-9999))
                if flav in ['t','b','l','g']:
                    varDict[f"{flav}_{i}_match_unique"]   = op.switch(op.rng_len(match)>i,op.static_cast("Int_t",lambda_match_unique_jet(match[i])),op.c_int(-9999))
                    varDict[f"{flav}_{i}_match_unique_count"]    = op.switch(op.rng_len(match)>i,op.static_cast("Int_t",map_jet_unique_count[match[i][1].idx]),op.c_int(-9999))
                    varDict[f"{flav}_{i}_match_resolved"] = op.switch(op.rng_len(match)>i,op.static_cast("Int_t",lambda_match_resolved(match[i])),op.c_int(-9999))
                    varDict[f"{flav}_{i}_match_resolved_count"] = op.switch(op.rng_len(match)>i,op.static_cast("Int_t",map_jet_resolved_count[match[i][1].idx]),op.c_int(-9999))
                    if flav in ['b','l','g']:
                        varDict[f"{flav}_{i}_match_central"]  = op.switch(op.rng_len(match)>i,op.static_cast("Int_t",lambda_central_jet(match[i])),op.c_int(-9999))
                        varDict[f"{flav}_{i}_match_forward"]  = op.switch(op.rng_len(match)>i,op.static_cast("Int_t",lambda_forward_jet(match[i])),op.c_int(-9999))
                # Gen information #
                for j,label in enumerate(['preFSR','postFSR']):
                    # Meta data #
                    varDict[f"{flav}_{i}_{label}_idx"]   = op.switch(op.rng_len(match)>i,op.static_cast("Int_t",match[i][0][j].idx),op.c_int(-9999))
                    varDict[f"{flav}_{i}_{label}_pdgId"] = op.switch(op.rng_len(match)>i,match[i][0][j].pdgId,op.c_int(-9999))
                    # Kinematic variables #
                    varDict[f"{flav}_{i}_{label}_E"]     = op.switch(op.rng_len(match)>i,match[i][0][j].p4.E(),op.c_float(-9999))
                    varDict[f"{flav}_{i}_{label}_pt"]    = op.switch(op.rng_len(match)>i,match[i][0][j].pt,op.c_float(-9999))
                    varDict[f"{flav}_{i}_{label}_eta"]   = op.switch(op.rng_len(match)>i,match[i][0][j].eta,op.c_float(-9999))
                    varDict[f"{flav}_{i}_{label}_phi"]   = op.switch(op.rng_len(match)>i,match[i][0][j].phi,op.c_float(-9999))
                    varDict[f"{flav}_{i}_{label}_M"]     = op.switch(op.rng_len(match)>i,match[i][0][j].mass,op.c_float(-9999))

                # Reco information #
                # Meta data #
                varDict[f"{flav}_{i}_reco_idx"] = op.switch(op.rng_len(match)>i,op.static_cast("Int_t",match[i][1].idx),op.c_int(-9999))
                # Kinematic variables #
                varDict[f"{flav}_{i}_reco_E"]   = op.switch(op.rng_len(match)>i,match[i][1].p4.E(),op.c_float(-9999))
                varDict[f"{flav}_{i}_reco_pt"]  = op.switch(op.rng_len(match)>i,match[i][1].pt,op.c_float(-9999))
                varDict[f"{flav}_{i}_reco_eta"] = op.switch(op.rng_len(match)>i,match[i][1].eta,op.c_float(-9999))
                varDict[f"{flav}_{i}_reco_phi"] = op.switch(op.rng_len(match)>i,match[i][1].phi,op.c_float(-9999))
                varDict[f"{flav}_{i}_reco_M"]   = op.switch(op.rng_len(match)>i,match[i][1].mass,op.c_float(-9999))

                # Gen jet information #
                if flav in ['t','b','l','g']: # for hadronic tau the reco is a jet
                    # Meta data #
                    varDict[f"{flav}_{i}_genjet_idx"] = op.switch(op.rng_len(match)>i,op.static_cast("Int_t",match[i][1].genJet.idx),op.c_int(-9999))
                    # Kinematic variables #
                    varDict[f"{flav}_{i}_genjet_E"]   = op.switch(op.AND(op.rng_len(match)>i,match[i][1].genJet.idx>=0),match[i][1].genJet.p4.E(),op.c_float(-9999))
                    varDict[f"{flav}_{i}_genjet_pt"]  = op.switch(op.AND(op.rng_len(match)>i,match[i][1].genJet.idx>=0),match[i][1].genJet.pt,op.c_float(-9999))
                    varDict[f"{flav}_{i}_genjet_eta"] = op.switch(op.AND(op.rng_len(match)>i,match[i][1].genJet.idx>=0),match[i][1].genJet.eta,op.c_float(-9999))
                    varDict[f"{flav}_{i}_genjet_phi"] = op.switch(op.AND(op.rng_len(match)>i,match[i][1].genJet.idx>=0),match[i][1].genJet.phi,op.c_float(-9999))
                    varDict[f"{flav}_{i}_genjet_M"]   = op.switch(op.AND(op.rng_len(match)>i,match[i][1].genJet.idx>=0),match[i][1].genJet.mass,op.c_float(-9999))

            # Add skim to plots #
            plots.append(Skim(f'match_{flav}',varDict,noSel))

        #----- Plotting -----#
        plots.extend(
            plotMatching(
                contName     = "e",
                cont         = reco_match_e_central,
                sel          = noSel,
                type         = 'lepton',
                Egen_binning = (3000,0,3000),
                dE_binning   = (12000,-1500,1500),
            )
        )
        plots.extend(
            plotMatching(
                contName     = "m",
                cont         = reco_match_m_central,
                sel          = noSel,
                type         = 'lepton',
                Egen_binning = (3000,0,3000),
                dE_binning   = (12000,-1500,1500),
            )
        )

        plots.extend(
            plotMatching(
                contName     = "b_central",
                cont         = reco_match_b_central,
                sel          = noSel,
                type         = 'jet',
                Egen_binning = (5000,0,5000),
                dE_binning   = (6000,-3000,3000),
            )
        )
        plots.extend(
            plotMatching(
                contName     = "l_central",
                cont         = reco_match_l_central,
                sel          = noSel,
                type         = 'jet',
                Egen_binning = (5000,0,5000),
                dE_binning   = (6000,-3000,3000),
            )
        )
        plots.extend(
            plotMatching(
                contName     = "g_central",
                cont         = reco_match_g_central,
                sel          = noSel,
                type         = 'jet',
                Egen_binning = (5000,0,5000),
                dE_binning   = (6000,-3000,3000),
            )
        )


        plots.extend(
            plotMatching(
                contName     = "b_forward",
                cont         = reco_match_b_forward,
                sel          = noSel,
                type         = 'jet',
                Egen_binning = (1000,0,5000),
                dE_binning   = (2000,-5000,5000),
            )
        )
        plots.extend(
            plotMatching(
                contName     = "l_forward",
                cont         = reco_match_l_forward,
                sel          = noSel,
                type         = 'jet',
                Egen_binning = (2000,0,10000),
                dE_binning   = (2000,-5000,5000),
            )
        )
        plots.extend(
            plotMatching(
                contName     = "g_forward",
                cont         = reco_match_g_forward,
                sel          = noSel,
                type         = 'jet',
                Egen_binning = (2000,0,10000),
                dE_binning   = (2000,-5000,5000),
            )
        )

        ##########################################################
        #                        FatJets                         #
        ##########################################################

        lambda_valid_subjets = lambda fatj : op.AND(
                fatj.subJet1.isValid,
                fatj.subJet1.isValid,
        )

        # Match to subjets #
        # Initially wanted to recycle the not-unique and not-resolved quarks
        # But fatjets are larger, and can have two quarks with eg DR=0.8
        # -> show up as two valid resolved AK4 jets
        # -> but also as two subjets of an AK8
        # So here need to start with the gen match pairs
        # Need to do the same procedure as for AK4 jets
        # But for the two subjets, so define a lambda externally
        lambda_subjet_matching = lambda  q,subjet: op.switch(
            # Check energy loss through FSR
            q[1].p4.E() >= 0.5 * q[0].p4.E(),
            # True case
            op.AND(
                op.deltaR(q[1].p4,subjet.p4) < 0.2,
                op.deltaR(q[1].p4,q[0].p4) < 0.6,
            ),
            # False case
            op.switch(
                # Check if other genparticles around post-FSR
                op.NOT(
                    op.rng_any(
                        t.GenPart,
                        lambda q1: op.AND(
                            q1.idx != q[1].idx,
                            q1.idx > q[0].idx,
                            q1.statusFlags & (0x1 << 0),
                            op.NOT(q1.statusFlags & (0x1 << 1)),
                            op.deltaR(q1.p4,q[0].p4) < 0.6,
                            op.rng_any(
                                q1.ancestors,
                                lambda a : a.idx == q[0].idx
                            ),
                        )
                    )
                ),
                # No harder constituent
                op.switch(
                    op.deltaR(q[1].p4,q[0].p4)<0.6,
                    # true case : post-FSR is within pre-FSR cone
                    op.deltaR(q[1].p4,subjet.p4) < 0.4,
                    # false case : look at pre-FSR for matching
                    op.deltaR(q[0].p4,subjet.p4) < 0.2,
                ),
                # Use other constituent
                op.rng_any(
                    t.GenPart,
                    lambda q1: op.AND(
                        q1.idx > q[0].idx,
                        q1.statusFlags & (0x1 << 0),
                        op.NOT(q1.statusFlags & (0x1 << 1)),
                        op.deltaR(q1.p4,q[0].p4) < 0.6,
                        op.rng_any(
                            q1.ancestors,
                            lambda a : a.idx == q[0].idx
                        ),
                        # Need to check all other potential candidates to make sure none has higher energy than q1
                        op.NOT(
                            op.rng_any(
                                t.GenPart,
                                lambda q2: op.AND(
                                    q2.idx != q1.idx,
                                    q2.idx > q[0].idx,
                                    q2.statusFlags & (0x1 << 0),
                                    op.NOT(q2.statusFlags & (0x1 << 1)),
                                    op.deltaR(q2.p4,q[0].p4) < 0.6,
                                    op.rng_any(
                                        q2.ancestors,
                                        lambda a : a.idx == q[0].idx
                                    ),
                                    q2.p4.E() > q1.p4.E()
                                ),
                            )
                        ),
                        # Finally check the jet matching
                        op.deltaR(q1.p4,subjet.p4) < 0.2,
                    )
                )
            )
        )

        fatjet_match = op.combine(
            rng = (gen_match_q,gen_match_q,t.FatJet),
            pred = lambda m1,m2,fatj: op.AND(
                # Make sure there are two subjets in the fatjet
                lambda_valid_subjets(fatj),
                # fatjet minimum energy
                # offline cut at pt 200 GeV (gen level is 170), no need for a cut here
                # Subjet minimum energy in gen level is 13 GeV, offline cut is 20
                fatj.subJet1.p4.E() > 15.,
                fatj.subJet2.p4.E() > 15.,
                # Safety check (should be handled in samePred, but just to be sure)
                m1[1].idx != m2[1].idx,
                m1[0].idx != m2[0].idx,
                # Weird FSR cases where gain energy
                m1[0].p4.E() > 0.99 * m1[1].p4.E(),
                m2[0].p4.E() > 0.99 * m2[1].p4.E(),
                # Make sure fatjet and subjets have genjets
                fatj.genJetAK8.idx >= 0,
                #fatj.subJet1.genJet.idx >= 0, # not in nanov7 apparently
                #fatj.subJet2.genJet.idx >= 0,
                # Match m1 to subJet1, and m2 to subJet2
                lambda_subjet_matching(m1,fatj.subJet1),
                lambda_subjet_matching(m2,fatj.subJet2),
                # Use the gen level number of B/C hadrons as proxy for the parton flavour
                op.multiSwitch(
                    (fatj.subJet1.nBHadrons>0, op.abs(m1[0].pdgId) == 5),   # b quark
                    (fatj.subJet1.nCHadrons>0, op.abs(m1[0].pdgId) == 4),   # c quark
                    op.AND(op.abs(m1[0].pdgId)>=1, op.abs(m1[0].pdgId)<=3) # u,d,s quarks
                ),
                op.multiSwitch(
                    (fatj.subJet2.nBHadrons>0, op.abs(m2[0].pdgId) == 5),   # b quark
                    (fatj.subJet2.nCHadrons>0, op.abs(m2[0].pdgId) == 4),   # c quark
                    op.AND(op.abs(m2[0].pdgId)>=1, op.abs(m2[0].pdgId)<=3) # u,d,s quarks
                ),
                # In case of very collimated subjets,
                # the two reco jets can be matched to both subjets
                # -> Resolve that by only getting the closest combination
                #   m1 <-> subjet1
                #   m2 <-> subjet2
                op.deltaR(m1[1].p4,fatj.subJet1.p4) < op.deltaR(m1[1].p4,fatj.subJet2.p4),
                op.deltaR(m2[1].p4,fatj.subJet2.p4) < op.deltaR(m2[1].p4,fatj.subJet1.p4),
            ),
            # We want to test all combinations, except for same post-FSR quarks
            samePred = lambda m1,m2: m1[1].idx != m2[1].idx,
        )

        # Count number of quarks and gluon ending up in the fatjet
        # Similar to resolved in AK4
        # Count is made on the number of gen match with
        #  - pre-FSR pt>=10
        #  - DR(q,fatjet) < 0.8
        # Note : when the post-FSR lost a lot of energy, maybe large kick
        # if that is the case, do the DR based on the pre-FSR
        map_prong_counts = op.map(
            t.FatJet,
            lambda fat:
                op.rng_count(
                    rng = gen_match_q,
                    pred = lambda q: op.AND(
                        q[0].pt >= 10,                   # pre-FSR pt >= 10 GeV
                        op.switch(
                            q[1].p4.E() > 0.5 * q[0].p4.E(),
                            op.deltaR(fat.p4,q[1].p4) < 0.8,   # DR(fat jet, post-FSR)
                            op.deltaR(fat.p4,q[0].p4) < 0.8,   # DR(fat jet, pre-FSR)
                        )
                    )
                ) + \
                op.rng_count(
                    rng = gen_match_g,
                    pred = lambda g: op.AND(
                        g[0].pt >= 10,                   # pre-FSR pt >= 10 GeV
                        op.switch(
                            g[1].p4.E() > 0.5 * g[0].p4.E(),
                            op.deltaR(fat.p4,g[1].p4) < 0.8,   # DR(fat jet, post-FSR)
                            op.deltaR(fat.p4,g[0].p4) < 0.8,   # DR(fat jet, pre-FSR)
                        )
                    )
                )
        )

        lambda_two_prongs = lambda match: map_prong_counts[match[2].idx] == 2
        fatjet_match_final = op.select(fatjet_match,lambda_two_prongs)

        # Exclude fatjets with leptons in the cone
        # Cannot use the previous map_clean, as we want DR>0.8 cleaning from leptons
        lambda_fatjet_clean = return_lambda_clean_overlap([gen_e,gen_m,gen_t],0.8)
        map_fatjet_clean = op.map(
            t.FatJet,
            lambda fatj : lambda_fatjet_clean(fatj)
        )
        fatjet_match_final = op.select(fatjet_match_final,lambda match: map_fatjet_clean[match[2].idx])

        # Apply the same FSR check
        lambda_fatjet_match_check_FSR = lambda match: op.AND(
            map_sum_FSR_energy[match[0][0].idx] <= 1.01 * match[0][0].p4.E(),
            map_sum_FSR_energy[match[1][0].idx] <= 1.01 * match[1][0].p4.E(),
        )
        fatjet_match_final = op.select(fatjet_match_final,lambda_fatjet_match_check_FSR)


        # select central fatjets only
        lambda_central_fatjet = lambda match : op.abs(match[2].eta) <= 2.4
        fatjet_match_central = op.select(fatjet_match_final,lambda_central_fatjet)

        # Monitoring tree #
        varsFatMatch = {
            'event': None,
            'n_fatjets': op.static_cast("UInt_t",op.rng_len(t.FatJet)),
            'n_fatjets_match': op.static_cast("UInt_t",op.rng_len(fatjet_match)),
            'n_fatjets_match_final': op.static_cast("UInt_t",op.rng_len(fatjet_match_final)),
            'n_fatjets_match_central': op.static_cast("UInt_t",op.rng_len(fatjet_match_central)),
        }

        for i in range(4):
            # Clean flag #
            varsFatMatch[f'fat_{i}_clean'] = op.switch(op.rng_len(fatjet_match)>i,op.static_cast("Int_t",map_fatjet_clean[fatjet_match[i].idx]),op.c_int(-9999))
            varsFatMatch[f'fat_{i}_two_prongs'] = op.switch(op.rng_len(fatjet_match)>i,op.static_cast("Int_t",lambda_two_prongs(fatjet_match[i])),op.c_int(-9999))
            varsFatMatch[f'fat_{i}_prong_count'] = op.switch(op.rng_len(fatjet_match)>i,op.static_cast("Int_t",map_prong_counts[fatjet_match[i][2].idx]),op.c_int(-9999))
            varsFatMatch[f'fat_{i}_check_FSR'] = op.switch(op.rng_len(fatjet_match)>i,op.static_cast("Int_t",lambda_fatjet_match_check_FSR(fatjet_match[i])),op.c_int(-9999))
            varsFatMatch[f'fat_{i}_sub1_sum_FSR'] = op.switch(op.rng_len(fatjet_match)>i,map_sum_FSR_energy[fatjet_match[i][0][0].idx],op.c_int(-9999))
            varsFatMatch[f'fat_{i}_sub2_sum_FSR'] = op.switch(op.rng_len(fatjet_match)>i,map_sum_FSR_energy[fatjet_match[i][1][0].idx],op.c_int(-9999))
            for j in range(2): # Subjets
                for k,label in enumerate(['preFSR','postFSR']):
                    varsFatMatch[f'fat_{i}_sub{j+1}_{label}_idx']  = op.switch(op.rng_len(fatjet_match)>i,op.static_cast("Int_t",fatjet_match[i][j][k].idx),op.c_int(-9999))
                    varsFatMatch[f'fat_{i}_sub{j+1}_{label}_pdgId']= op.switch(op.rng_len(fatjet_match)>i,fatjet_match[i][j][k].pdgId,op.c_int(-9999))
                    varsFatMatch[f'fat_{i}_sub{j+1}_{label}_E']    = op.switch(op.rng_len(fatjet_match)>i,fatjet_match[i][j][k].p4.E(),op.c_int(-9999))
                    varsFatMatch[f'fat_{i}_sub{j+1}_{label}_pt']   = op.switch(op.rng_len(fatjet_match)>i,fatjet_match[i][j][k].pt,op.c_int(-9999))
                    varsFatMatch[f'fat_{i}_sub{j+1}_{label}_eta']  = op.switch(op.rng_len(fatjet_match)>i,fatjet_match[i][j][k].eta,op.c_int(-9999))
                    varsFatMatch[f'fat_{i}_sub{j+1}_{label}_phi']  = op.switch(op.rng_len(fatjet_match)>i,fatjet_match[i][j][k].phi,op.c_int(-9999))
                    varsFatMatch[f'fat_{i}_sub{j+1}_{label}_M']    = op.switch(op.rng_len(fatjet_match)>i,fatjet_match[i][j][k].mass,op.c_int(-9999))

            varsFatMatch[f'fat_{i}_sub1_subjet_idx'] = op.switch(op.rng_len(fatjet_match)>i,op.static_cast("Int_t",fatjet_match[i][2].subJet1.idx),op.c_int(-9999))
            varsFatMatch[f'fat_{i}_sub1_subjet_E']   = op.switch(op.rng_len(fatjet_match)>i,fatjet_match[i][2].subJet1.p4.E(),op.c_int(-9999))
            varsFatMatch[f'fat_{i}_sub1_subjet_pt']  = op.switch(op.rng_len(fatjet_match)>i,fatjet_match[i][2].subJet1.pt,op.c_int(-9999))
            varsFatMatch[f'fat_{i}_sub1_subjet_eta'] = op.switch(op.rng_len(fatjet_match)>i,fatjet_match[i][2].subJet1.eta,op.c_int(-9999))
            varsFatMatch[f'fat_{i}_sub1_subjet_phi'] = op.switch(op.rng_len(fatjet_match)>i,fatjet_match[i][2].subJet1.phi,op.c_int(-9999))
            varsFatMatch[f'fat_{i}_sub1_subjet_M']   = op.switch(op.rng_len(fatjet_match)>i,fatjet_match[i][2].subJet1.mass,op.c_int(-9999))
            varsFatMatch[f'fat_{i}_sub1_subjet_nB']  = op.switch(op.rng_len(fatjet_match)>i,fatjet_match[i][2].subJet1.nBHadrons,op.c_int(-9999))
            varsFatMatch[f'fat_{i}_sub1_subjet_nC']  = op.switch(op.rng_len(fatjet_match)>i,fatjet_match[i][2].subJet1.nCHadrons,op.c_int(-9999))

            varsFatMatch[f'fat_{i}_sub2_subjet_idx'] = op.switch(op.rng_len(fatjet_match)>i,op.static_cast("Int_t",fatjet_match[i][2].subJet2.idx),op.c_int(-9999))
            varsFatMatch[f'fat_{i}_sub2_subjet_E']   = op.switch(op.rng_len(fatjet_match)>i,fatjet_match[i][2].subJet2.p4.E(),op.c_int(-9999))
            varsFatMatch[f'fat_{i}_sub2_subjet_pt']  = op.switch(op.rng_len(fatjet_match)>i,fatjet_match[i][2].subJet2.pt,op.c_int(-9999))
            varsFatMatch[f'fat_{i}_sub2_subjet_eta'] = op.switch(op.rng_len(fatjet_match)>i,fatjet_match[i][2].subJet2.eta,op.c_int(-9999))
            varsFatMatch[f'fat_{i}_sub2_subjet_phi'] = op.switch(op.rng_len(fatjet_match)>i,fatjet_match[i][2].subJet2.phi,op.c_int(-9999))
            varsFatMatch[f'fat_{i}_sub2_subjet_M']   = op.switch(op.rng_len(fatjet_match)>i,fatjet_match[i][2].subJet2.mass,op.c_int(-9999))
            varsFatMatch[f'fat_{i}_sub2_subjet_nB']  = op.switch(op.rng_len(fatjet_match)>i,fatjet_match[i][2].subJet2.nBHadrons,op.c_int(-9999))
            varsFatMatch[f'fat_{i}_sub2_subjet_nC']  = op.switch(op.rng_len(fatjet_match)>i,fatjet_match[i][2].subJet2.nCHadrons,op.c_int(-9999))

            varsFatMatch[f'fat_{i}_idx'] = op.switch(op.rng_len(fatjet_match)>i,op.static_cast("Int_t",fatjet_match[i][2].idx),op.c_int(-9999))
            varsFatMatch[f'fat_{i}_E']   = op.switch(op.rng_len(fatjet_match)>i,fatjet_match[i][2].p4.E(),op.c_int(-9999))
            varsFatMatch[f'fat_{i}_pt']  = op.switch(op.rng_len(fatjet_match)>i,fatjet_match[i][2].pt,op.c_int(-9999))
            varsFatMatch[f'fat_{i}_eta'] = op.switch(op.rng_len(fatjet_match)>i,fatjet_match[i][2].eta,op.c_int(-9999))
            varsFatMatch[f'fat_{i}_phi'] = op.switch(op.rng_len(fatjet_match)>i,fatjet_match[i][2].phi,op.c_int(-9999))
            varsFatMatch[f'fat_{i}_M']   = op.switch(op.rng_len(fatjet_match)>i,fatjet_match[i][2].mass,op.c_int(-9999))
            varsFatMatch[f'fat_{i}_nB']  = op.switch(op.rng_len(fatjet_match)>i,fatjet_match[i][2].nBHadrons,op.c_int(-9999))
            varsFatMatch[f'fat_{i}_nC']  = op.switch(op.rng_len(fatjet_match)>i,fatjet_match[i][2].nCHadrons,op.c_int(-9999))
            varsFatMatch[f'fat_{i}_hadronFlavour']  = op.switch(op.rng_len(fatjet_match)>i,fatjet_match[i][2].hadronFlavour,op.c_int(-9999))


        plots.append(Skim(f'match_fatjet',varsFatMatch,noSel))


        fatjet_match_bb = op.select(
            fatjet_match_central,
            lambda match: op.AND(
                op.abs(match[0][1].pdgId) == 5,

                op.abs(match[1][1].pdgId) == 5,
            )
        )
        fatjet_match_ll = op.select(
            fatjet_match_central,
            lambda match: op.AND(
                op.abs(match[0][1].pdgId) != 5,
                op.abs(match[1][1].pdgId) != 5,
            )
        )

        plots.extend(
            plotFatjetMatching(
                cont            = fatjet_match_bb,
                sel             = noSel,
                sub1            = 'b',
                sub2            = 'b',
                Egen_binning    = (3000,0,3000),
                dE_binning      = (4000,-2000,2000),
            )
        )
        plots.extend(
            plotFatjetMatching(
                cont            = fatjet_match_ll,
                sel             = noSel,
                sub1            = 'l',
                sub2            = 'l',
                Egen_binning    = (3000,0,3000),
                dE_binning      = (4000,-2000,2000),
            )
        )

        return plots


    def postProcess(self, taskList, config=None, workdir=None, resultsdir=None):
        pass # No need to do plotting

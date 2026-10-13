# Author: Shubhodeep Chatterjee
# Report text. Executed by build_report.py, which provides the tables (t_*) and result variables.

_v = vs.loc[order_d]
_adequate = [d for d in order_d if vs.loc[d, "verdict_sim"] == "adequate"]
_inadequate = [d for d in order_d if vs.loc[d, "verdict_sim"] != "adequate"]
_pkx = pk[pk.kind == "exclusive"]; _pkc = pk[pk.kind == "cooccur"]; _pkb = pk[pk.kind == "before"]
_cbn_n = int(co.cbn_orders_compared.sum()); _cbn_ok = int(round((co.cbn_agreement * co.cbn_orders_compared).sum()))
_mhn_agree = ", ".join(f"{d} {co.loc[d, 'direction_agreement']:.2f}" for d in order_d)
_tm = tmv.loc[order_d]
_hw, _hr = int(_tm.later_has_higher_base_wrong.sum()), int(_tm.later_has_higher_base_right.sum())
_lw, _lr = int(_tm.later_has_lower_base_wrong.sum()), int(_tm.later_has_lower_base_right.sum())
from scipy.stats import fisher_exact as _fe
_diag_p = _fe([[_hw, _hr], [_lw, _lr]])[1]
_bt_m_n = int(_tm.mhn_vs_bt_pairs.sum()); _bt_m_ok = int(round((_tm.mhn_vs_bt_agreement * _tm.mhn_vs_bt_pairs).sum()))
_bt_c_n = int(_tm.cbn_vs_bt_orders.sum()); _bt_c_ok = int(round((_tm.cbn_vs_bt_agreement.fillna(0) * _tm.cbn_vs_bt_orders).sum()))
_and = lambda xs: ", ".join(xs[:-1]) + " and " + xs[-1] if len(xs) > 1 else xs[0]

ABSTRACT = (
    f"We built per-disease models of cytogenetic progression for seven haematological malignancies from "
    f"{n_total:,} non-duplicate karyotypes in the Mitelman Database. Clone structure was parsed from ISCN "
    f"notation, events were frozen per disease by a power calculation, and three existing methods were fitted "
    f"to the same data: a co-occurrence baseline with a count-preserving permutation null, a conjunctive "
    f"Bayesian network (CBN) and a Mutual Hazard Network (MHN), with a model-free Bradley-Terry timing from the "
    f"within-patient clone orders as reference. Five diseases had at least 200 cases with a "
    f"within-patient clone order; APL and myeloma were limited to the baseline. The MHN fitted held-out data "
    f"better than independence and than the CBN in every disease and every one of {CFG['mhn']['cv_folds']} folds, its interactions were stable under "
    f"{CFG['mhn']['n_bootstrap']} bootstrap refits, and they held in unselected series and in post-2000 publications. "
    f"Simulation recovery was adequate in {_and(_adequate)} and failed in {_and(_inadequate)}. The MHN's event orders, however, "
    f"agreed with the within-patient clone orders only in CML; elsewhere they were close to chance, and they were wrong mainly "
    f"when the later event has the higher base rate ({_hw} of {_hw + _hr} such pairs), so events often seen alone, such as "
    f"trisomy 8, are placed first. The sparse CBN posets agreed with the clone orders in {_cbn_ok} of {_cbn_n} testable cases "
    f"and with the clone timing in {_bt_c_ok} of {_bt_c_n}. The supported conclusions are therefore the interaction "
    f"structure from the MHN and the orders backed by the clone data and the CBN, not the MHN's orders on their own. "
    f"Code, event tables and figures are released.")


RESULTS_BLOCKS = [
    ("h1", "Results"),
    ("h2", "Co-occurrence baseline"),
    ("p", f"{{ref:tab:cooc}} counts the significant pairs. Mutual exclusivity dominates wherever founding "
          "rearrangements are in the vocabulary: in AML, t(8;21), inv(16), KMT2A rearrangements and "
          "inv(3)/t(3;3) exclude one another and the secondary changes of complex karyotypes, and in ALL the "
          "founding translocations exclude the trisomies of high hyperdiploidy. Requiring Fisher's test and "
          "the permutation null to agree on the direction mattered: some pairs, 31 of them in ALL, co-occur "
          "more often than under independence but less often than busy karyotypes alone would predict, and "
          "they are no longer called."),
    t_cooc,
    ("h2", "Model-free clone orders"),
    ("p", "The clone trees give direct within-patient evidence ({ref:tab:cloneorders}). The orders are "
          "lopsided and match the established biology: 5q loss precedes 17p loss (96:5) and -7 in AML, "
          "t(8;21) precedes -Y and -X, inv(16) precedes +22, t(9;22) precedes the extra Ph and +8 in CML, "
          "and 5q loss precedes 17p loss in MDS. Trisomy 8 is a secondary event in nearly every case where it "
          "shares a karyotype with a founding change (KMT2A 61:4, inv(16) 44:1)."),
    t_cloneorders,
    ("p", "Pooling all pairs per disease into a Bradley-Terry timing gives a full order with uncertainty "
          "({ref:tab:timing}, {ref:fig:timing}). Founding rearrangements come first wherever they are in the "
          "vocabulary: t(8;21), inv(3)/t(3;3), inv(16) and KMT2A in AML (with 5q loss among them), t(9;22) in CML, t(9;22) and KMT2A in ALL, "
          "and 5q loss in MDS. Losses of 7, 7q and 17p follow, and the trisomies (+8, +11, +13, +19, +21, +22) "
          "come last. In CLL +12 ranks first, consistent with its role as a founding lesion. Sex-chromosome loss "
          "ranks early in CML and MDS, which fits its age-related origin outside the leukaemic process."),
    t_timing,
    ("figure", "fig:timing", str(FIGURES / "clone_timing.png"),
     "Event timing from within-patient clone orders: Bradley-Terry earliness with 90% case-bootstrap intervals "
     "(orange; higher = earlier), and the MHN's mean probability that each event comes before the others (blue).", True),
    ("figure", "fig:aml", str(FIGURES / "AML.png"),
     "AML. (a) MHN log θ (row: affected event, column: acting event); dots mark stable edges. (b) Strongest "
     "stable promotions; edge labels give bootstrap sign stability. (c) Co-occurrence baseline, log odds ratio; "
     "* significant. (d) For the pairs with the most clone evidence, the MHN's P(a before b) with its bootstrap "
     "10-90% interval (blue) against the share of clone pairs ordered a→b (orange).", True),
    ("figure", "fig:cml", str(FIGURES / "CML.png"),
     "CML, panels as in {ref:fig:aml}. Here the MHN's orders and the clone orders agree.", True),
    ("figure", "fig:overview", str(FIGURES / "overview.png"),
     "Held-out log-likelihood gain over independence (MHN on all model events, CBN on ten) and share of "
     "stable MHN edges per disease; the dashed line is the 50% stability rule.", False),
    ("h2", "Mutual Hazard Networks"),
    ("p", f"{{ref:tab:mhn}} summarises the MHN fits. Held-out log-likelihood beat the independence model in all "
          f"five diseases, by {vs.loc['AML', 'heldout_mhn'] - vs.loc['AML', 'heldout_independence']:.2f} nats per case in AML "
          f"and {vs.loc['ALL', 'heldout_mhn'] - vs.loc['ALL', 'heldout_independence']:.2f} in ALL, and in every fold (smallest "
          f"fold gain {_v.heldout_mhn_gain_fold_min.min():.3f}, in {_v.heldout_mhn_gain_fold_min.idxmin()}). The held-out "
          "likelihood is flat across the penalty grid near its optimum, so the choice of λ matters little. Between "
          f"{_v.mhn_frac_stable.min():.0%} and {_v.mhn_frac_stable.max():.0%} of edges kept their sign in at "
          f"least {TH['stable_sign_agreement']:.0%} of bootstrap refits, so every disease passes the stability rule. The refits on "
          f"unselected series (r = {_v.unselected_corr.min():.2f} to {_v.unselected_corr.max():.2f}) and on "
          f"cases published from 2000 (r = {_v.post2000_corr.min():.2f} to {_v.post2000_corr.max():.2f}) "
          "reproduce the full-data networks, so publication bias and era do not drive them."),
    ("p", "Most interactions are negative. This is expected: log θ is measured against the positive "
          "association that a shared observation time alone produces, so mild raw co-occurrence still carries "
          "a negative log θ. In CLL, for example, every interaction except 11q- on 13q loss is negative, yet "
          "the fitted model reproduces every observed pairwise odds ratio. The strongest stable promotions "
          "are biologically familiar: t(8;21) on -X and -Y, 5q loss on 17p loss, +22 on +19 and +21 in AML; "
          "the extra Ph on +8 and +21 on +19 in CML; +21 on the other trisomies of high hyperdiploidy in ALL "
          "[@moorman2012]; and 18 loss on 17p loss in MDS. Two strong edges are definitional: between 17p loss "
          "and i(17q) in CML, where the events are nested, and t(1;19) on 1q gain in ALL, where the unbalanced "
          "der(19)t(1;19) carries the gain."),
    ("p", f"Simulation recovery was adequate in {_and(_adequate)} (mean precision "
          f"{min(vs.loc[d, 'mhn_sim_precision'] for d in _adequate):.2f} or more, recall "
          f"{min(vs.loc[d, 'mhn_sim_recall'] for d in _adequate):.2f} or more over {CFG['mhn']['n_simulations']} data sets) and "
          f"failed in {_and(_inadequate)}: refitted networks recovered only {vs.loc['MDS', 'mhn_sim_precision']:.2f} "
          f"(MDS) and {vs.loc['CML', 'mhn_sim_precision']:.2f} (CML) of their own edges. By the pre-specified "
          "rule this failure is the finding for those two diseases: at their sample sizes the MHN's individual "
          "edges are not identifiable, even though the networks are stable under resampling."),
    t_mhn,
    ("h2", "Conjunctive Bayesian networks"),
    ("p", "The CBNs are sparse ({ref:tab:cbn}). In AML the poset has 5q loss before -7, 7q- and 17p loss, "
          "+8 before +21 and t(8;21) before -Y; in CML t(9;22) precedes the extra Ph, 17p loss, +19 and 3q26 "
          "rearrangements; in ALL +21 precedes the other hyperdiploid trisomies; in CLL 13q loss precedes 11q- "
          "and 17p loss. Simulation recovery was perfect in every disease and the noise rate ε stayed between "
          f"{_v.cbn_eps.min():.3f} and {_v.cbn_eps.max():.3f}. On the same {CFG['cbn']['n_events']} events the MHN fitted held-out "
          "data better than the CBN in every disease and every fold, so the CBN's hard precedence constraints cost likelihood."),
    t_cbn,
    ("h2", "Validation against prior knowledge and clone orders"),
    ("p", f"Against the written test set ({{ref:tab:pk}}), the baseline confirmed all {len(_pkx)} expected "
          f"exclusivities and {int(_pkc.baseline_pass.sum())} of {len(_pkc)} expected co-occurrences (the "
          "exception is 5q loss with -7 in MDS, which is not significant). The MHN confirmed every exclusivity "
          f"it could test, but {int((_pkc.mhn_pass == True).sum())} of {int(_pkc.mhn_pass.notna().sum())} co-occurrences, "
          "for the shared-time reason above. Every expected order "
          "with at least five clone pairs was confirmed by the clone data. The MHN confirmed "
          f"{int((_pkb.mhn_pass == True).sum())} of {int(_pkb.mhn_pass.notna().sum())} expected orders and the CBN "
          f"implied {int((_pkb.cbn_pass == True).sum())} of {int(_pkb.cbn_pass.notna().sum())}; the CBN's misses "
          "are orders its poset leaves unconstrained, while the MHN's misses mostly reverse the expected order."),
    t_pk,
    ("p", "The decisive check is the comparison with clone orders, which the cross-sectional models never "
          f"see ({{ref:tab:clone}}). The MHN's order probabilities agreed in direction with the clone pairs in "
          f"{_mhn_agree}: clearly above chance only in CML (CLL rests on six pairs). The CBN's implied orders agreed in {_cbn_ok} of "
          f"{_cbn_n} testable cases. The MHN's errors follow one pattern ({{ref:fig:aml}}d): it places "
          "first the events with high base rates that are often seen alone. Trisomy 8 is the only model event "
          "in 43% of AML cases that carry it and -7 in 34%, so the MHN puts +8 and -7 before 5q loss, KMT2A "
          "and inv(16), the reverse of what the clone trees show in the same patients. In CML, where t(9;22) "
          "is both the most common event and the founder, the two signals coincide."),
    t_clone,
    ("p", f"Against the Bradley-Terry timing the picture is the same: on event pairs with a bootstrap-stable "
          f"timing order, the MHN agreed in {_bt_m_ok} of {_bt_m_n} and the CBN in {_bt_c_ok} of {_bt_c_n}. "
          "The disagreement has a measurable cause. Across the five diseases, when the event that comes later in "
          f"the clone trees has the higher MHN base rate, the MHN orders the pair wrongly in {_hw} of {_hw + _hr} "
          f"cases; when it has the lower base rate, in {_lw} of {_lw + _lr} (Fisher p = {_diag_p:.1e}). The MHN's "
          "order is set by how often an event appears, not by when it appears."),
    ("h2", "Verdict per disease"),
    ("list", [
        "**AML and ALL:** stable, recoverable and bias-robust interaction networks; orders supported where "
        "the clone data and the CBN agree (for example 5q loss before 17p loss and -7, t(8;21) before -Y, "
        "inv(16) before +22, +21 before the other trisomies). The MHN's own orders are not supported.",
        "**CML:** the MHN's orders agree with the clone data (t(9;22) first, then the extra Ph, +19 and "
        "i(17q)), but its simulation recovery failed, so only the stable, clone-backed subset is reported.",
        "**MDS:** simulation recovery failed and the MHN's orders disagree with the clone data; only the "
        "clone-backed orders (5q loss before 17p loss, 18 loss and 13q loss) and the baseline stand.",
        "**CLL:** recoverable network of mostly negative interactions; +12 is exclusive of 13q loss, and "
        "13q loss precedes 11q- [@dohner2000].",
        "**APL and MM:** co-occurrence only under the stop rule (64 and 55 cases with a clone order).",
    ]),
]

CONCLUSION = (
    "Fitted per disease on Mitelman karyotypes, the MHN describes which cytogenetic events go together and "
    "which avoid each other better than independence or a CBN, and that description is stable and robust to "
    "publication bias. Its event orders are a different matter: checked against the clone trees inside the same "
    "patients, they hold only in CML. Cross-sectional order inference is pulled toward events that are common "
    "as sole abnormalities, and the within-patient clone structure in Mitelman is the stronger ordering signal. "
    "Models that use the clone trees directly, and external validation on sequencing-era cohorts, are the "
    "natural next steps.")


META = dict(
    title="Per-Disease Cytogenetic Progression Models from the Mitelman Database: "
          "Co-occurrence, Conjunctive Bayesian Networks and Mutual Hazard Networks",
    author="Shubhodeep Chatterjee",
    affiliation=["Independent research", "github.com/shubhodeepcx/MitelmanProgression"],
    keywords="Cancer cytogenetics, karyotype, ISCN, Mitelman Database, tumour progression, "
             "Mutual Hazard Networks, conjunctive Bayesian networks, haematological malignancy",
    abstract=ABSTRACT,
)

BLOCKS = [
    ("h1", "Introduction"),
    ("p", "Chromosomal aberrations in haematological malignancies are not acquired at once. A founding "
          "rearrangement such as t(8;21) or t(9;22) is followed by secondary gains and losses, and the "
          "secondary changes that accompany a given founder are characteristic rather than random "
          "[@johansson2002, @grimwade2010]. Knowing the usual order of events, and which events promote or "
          "suppress one another, is a description of clonal evolution that complements outcome-based "
          "classifications."),
    ("p", "The Mitelman Database of Chromosome Aberrations and Gene Fusions in Cancer [@mitelman] collects "
          "every published cancer karyotype in ISCN notation [@iscn2020]. Two kinds of ordering signal are "
          "available in it. Within a patient, a karyotype that lists a stemline and sidelines states which "
          "aberrations are shared by all abnormal cells and which were added later. Across patients, the "
          "joint distribution of events carries the information that cancer progression models such as "
          "oncogenetic trees [@desper1999], conjunctive Bayesian networks (CBN) [@beerenwinkel2007, @gerstung2009] "
          "and Mutual Hazard Networks (MHN) [@schill2020] turn into an order. Earlier multivariate analyses of "
          "Mitelman data described karyotypic pathways in solid tumours [@hoglund2001]; here the aim is a "
          "per-disease model for the haematological malignancies with explicit validation."),
    ("p", "**Claim and scope.** For each disease with enough multi-clone cases, the study asks which "
          "pairwise orderings and interactions between recurrent cytogenetic events are supported by the "
          "data, with an uncertainty estimate for each, and which are not. It fits existing methods only, "
          "and it makes no clinical claims: the orders describe the published case mix."),
    ("p", "The contributions are: (i) an open pipeline with one entry point, from the Mitelman bulk files to "
          "per-disease event tables, models, validation tables and figures; (ii) a released tidy event table "
          "with inferred clone trees; (iii) a co-occurrence baseline, a CBN and an MHN fitted per disease on "
          "the same data; and (iv) six validation checks with pre-specified pass thresholds and stop rules."),

    ("h1", "Data and Biases"),
    ("h2", "Source and cohorts"),
    ("p", f"The Mitelman bulk download (md5 {frozen['mitelman_zip_md5']}) holds 81,148 cases from 12,342 "
          "references. Diseases were assigned from morphology codes: acute myeloid leukaemia (AML, excluding "
          "acute promyelocytic leukaemia), acute promyelocytic leukaemia (APL, morphology code or t(15;17) "
          "in a myeloid case), myelodysplastic neoplasms (MDS), chronic myeloid leukaemia (CML), acute "
          "lymphoblastic leukaemia (ALL), chronic lymphocytic leukaemia (CLL) and multiple myeloma (MM). "
          f"After removing likely duplicates, {n_total:,} cases remain ({{ref:tab:data}})."),
    t_data,
    ("h2", "Biases"),
    ("list", [
        "**Publication and selection bias.** Mitelman marks each series as Unselected (consecutive cases) or "
        "Selected (published for an unusual feature). Selected series over-represent rare rearrangements: "
        "in AML, for example, t(8;21) and KMT2A rearrangements are about 1.3 to 1.4 times as frequent in "
        "Selected as in Unselected series. Every MHN was refitted on Unselected cases only.",
        "**Era.** Banding resolution and reporting practice changed over five decades; MHNs were "
        "refitted on cases published from 2000.",
        f"**Duplicate patients.** There is no usable patient identifier. Cases from different references "
        f"with the same disease, sex, age and an identical karyotype string seen at most three times in "
        f"the database were treated as one patient; {n_dup} cases were flagged and excluded.",
        "**Cryptic lesions.** Karyotyping misses cryptic rearrangements (for example t(12;21) in ALL) and "
        "small copy-number changes, so frequencies understate them.",
        "**Cross-sectional sampling.** Most cases contribute one investigation; ordering is inferred "
        "from patterns across cases plus the clone structure within cases.",
    ]),

    ("h1", "Method"),
    ("h2", "Parsing clones and building events"),
    ("p", "Each karyotype is split into clones, and each clone into aberration tokens. The idem shorthand is expanded to "
          "the stemline's aberrations (and multiplied for idemx2 polyploid sidelines), opposite numerical "
          "changes cancel, sex-chromosome complements and uncounted markers are ignored. A clone's parent "
          "is the largest abnormal clone in the same investigation whose aberrations are a proper subset of "
          "its own (ties go to the earlier-listed clone); every event in a parent precedes every event the "
          "child adds. Copy-number events are the union of what the ISCN text states and arm-level calls "
          "from the CytoConverter [@cytoconverter] segments shipped with Mitelman, computed relative to each "
          "clone's own baseline ploidy (≥ 50% of an arm, ≥ 90% of a chromosome). The complex-karyotype flag "
          "follows ELN 2022 [@dohner2022]."),
    ("p", "The parser was tested on 50 karyotypes drawn by stratified sampling (single clone, idem "
          "sidelines, three or more clones, idemx2, unrelated clones, constitutional and uncertain "
          "aberrations, KMT2A variants, CML and MDS multi-clone cases) and annotated by hand in an "
          "independent review pass. The parser reproduces every annotated feature set, parent and ordered "
          f"pair; the full suite has {N_TESTS} tests, run with a linter on every push."),
    ("h2", "Event vocabulary and power"),
    ("p", "Candidate events per disease were listed from the cytogenetics literature before any ordering "
          "result was looked at. An event enters the vocabulary if at least N non-duplicate cases carry it, "
          "where N is the smallest count giving 80% power for Fisher's exact test at a Bonferroni level of "
          "0.05/190 (20 events) to detect both an odds ratio of 3 and one of 0.2 against a partner of 10% "
          "frequency. This gives N ≈ 200 in every disease. At most 20 events are kept per disease and "
          "defining events are always kept. The vocabulary was then frozen. Count-derived events (complex, "
          "high hyperdiploid, hypodiploid) are functions of the others and enter the baseline only."),
    ("h2", "Co-occurrence baseline"),
    ("p", "For every event pair, Fisher's exact test is compared with a permutation null that keeps each "
          "case's number of events fixed and the event frequencies approximately matched (weighted sampling "
          f"without replacement through Gumbel top-k draws, weights calibrated by iterative scaling, {CFG['cooccurrence']['n_permutations']:,} "
          "permutations). This removes the busy-karyotype effect, where complex cases make every pair look "
          "co-occurring. Both p-value sets are corrected by Benjamini-Hochberg [@bh1995]; a pair is "
          "significant if both q < 0.05 and both tests agree on the direction."),
    ("h2", "Conjunctive Bayesian network"),
    ("p", "A CBN is a partial order on events: an event can occur only after all its parents [@beerenwinkel2007]. "
          "Hidden genotypes compatible with the poset are generated with per-event probabilities θ_i, and "
          "each observed bit is flipped with probability ε, as in H-CBN [@gerstung2009]. θ and ε are fitted by "
          "expectation maximisation over the enumerated genotypes, and the poset by greedy BIC search with "
          "add, remove and reverse moves. The reference C implementations could not be obtained in the "
          "analysis environment, so the published discrete-time model was implemented and checked on "
          f"simulated posets. Because genotypes are enumerated, the CBN uses the {CFG['cbn']['n_events']} most frequent model "
          "events per disease."),
    ("h2", "Mutual Hazard Network"),
    ("p", "An MHN [@schill2020] models progression as a continuous-time Markov chain on the 2ⁿ event "
          "states. With event set x already present, event i occurs at rate"),
    ("eq", r"\lambda_i(x) = \Theta_{ii} \prod_{j \in x} \Theta_{ij},",
     "λ_i(x) = Θ_ii · ∏_{j ∈ x} Θ_ij"),
    ("p", "so log Θ_ij > 0 means event j, once present, promotes event i, and log Θ_ij < 0 means it "
          "inhibits it. Tumours are observed at an Exp(1) time. Because every event shares that observation "
          "time, events with no interaction are still positively associated in the observed data, so log θ "
          "is relative to this shared-time association, not to independence. The likelihood and its "
          "gradient were computed with state-space restriction (only states below an observed state are "
          "visited). The reference mhn package returned non-finite scores for states with five or six events "
          "in this environment, so the likelihood was reimplemented and checked against the package on "
          "smaller states, against a closed form for diagonal Θ and against finite differences. Θ is fitted "
          "by L-BFGS-B on the mean log-likelihood minus λ times a smooth L1 penalty on the off-diagonal "
          f"entries, with λ = c/N and c ∈ {{{', '.join(f'{c:g}' for c in CFG['mhn']['lambda_c_grid'])}}} chosen by "
          f"{CFG['mhn']['cv_folds']}-fold held-out likelihood. Order "
          f"probabilities P(a before b | both present) come from {CFG['mhn']['order_trajectories']:,} simulated trajectories."),
    ("h2", "Model-free timing from clone orders"),
    ("p", "The clone trees also support an ordering that needs no progression model. Pooling the parent-to-child "
          "pairs of one disease, the Bradley-Terry paired-comparison model [@bradleyterry1952] gives every event "
          "an earliness score s, with P(a before b) = 1/(1 + exp(-(s_a - s_b))), fitted by penalised maximum "
          "likelihood (ridge 0.01, centred scores). Only events in at least "
          f"{CFG['clone_timing']['min_pairs_per_event']} ordered pairs enter. Intervals come from "
          f"{CFG['clone_timing']['n_bootstrap']} bootstrap resamples of cases, so a case that contributes many "
          "pairs is resampled as a unit. These scores are the reference against which the MHN and CBN orders "
          "are compared."),
    ("h2", "Validation and stop rules"),
    ("p", f"Six checks were specified before the fits: (1) bootstrap stability ({CFG['mhn']['n_bootstrap']} MHN and "
          f"{CFG['cbn']['n_bootstrap']} CBN refits); "
          f"(2) {CFG['mhn']['cv_folds']}-fold held-out likelihood against an independence model, with the spread "
          f"across folds; (3) simulation recovery ({CFG['mhn']['n_simulations']} data sets "
          "drawn from each fitted model at the real sample size and refitted); (4) concordance with a "
          "written test set of established relationships; (5) agreement with the within-patient clone "
          "orders and the clone timing, which the cross-sectional models never see; and (6) publication-bias "
          f"and era sensitivity refits. An MHN interaction is an edge if |log θ| ≥ {TH['edge_abs_log_theta']} and "
          f"stable if it keeps its sign in ≥ {TH['stable_sign_agreement']:.0%} "
          f"of bootstrap refits. Simulation recovery is adequate if mean precision ≥ {TH['sim_recovery_precision']} "
          f"and recall ≥ {TH['sim_recovery_recall']}. All thresholds, counts and seeds are fixed in one "
          "configuration file. "
          f"Stop rules: a disease with fewer than {TH['min_cases_with_clone_order']} cases with a clone order gets the co-occurrence "
          "baseline only; failed simulation recovery is reported as the finding; when fewer than half the "
          "edges are stable, only the stable subset is presented."),
] + RESULTS_BLOCKS + [
    ("h1", "What the Models Cannot Say"),
    ("list", [
        "**Not patient trajectories.** The MHN describes the published case mix; a high P(a before b) "
        "does not mean that every patient acquires a before b.",
        "**Interactions are relative.** Negative log θ often means less co-occurrence than shared "
        "observation time would produce, not raw mutual exclusivity (CLL shows this most clearly).",
        "**Vocabulary limits.** Events below N, cryptic lesions and sub-arm changes are invisible. "
        "Nested events (17p loss and i(17q) in CML) give definitional edges.",
        "**Residual bias.** The Unselected and post-2000 refits bound but do not remove publication bias.",
        "**Pooled timing.** The Bradley-Terry timing pools clone orders across patients and assumes one "
        "shared order per disease; subgroups with different routes (for example core-binding-factor and "
        "complex-karyotype AML) are averaged.",
        "**No cross-dataset check.** TCGA and BeatAML could not be reached from the analysis environment "
        "(no outbound access to their portals), so external validation on sequencing-era cohorts remains "
        "to be done.",
        "**No clinical use.** None of these results should inform diagnosis or treatment.",
    ]),
    ("h1", "Conclusion"),
    ("p", CONCLUSION),
    ("h1", "Data and Code Availability"),
    ("p", "Code, the frozen vocabulary, the released event tables, model outputs and figures are at "
          "github.com/shubhodeepcx/MitelmanProgression under the MIT licence. Running run_all.sh reproduces "
          "every number in this report from the Mitelman bulk download. A command-line tool runs any step, "
          "every threshold, resampling count and seed sits in one configuration file, a Docker image pins the "
          "exact package versions, and continuous integration runs the linter and tests on every push."),
    ("h1", "Appendix: Supplementary Figures"),
    ("p", "Figures for the remaining diseases, with the panels of {ref:fig:aml} for the diseases with enough "
          "clone orders and the co-occurrence baseline for APL and MM."),
] + supp_figs

REFS = [
    ("mitelman", "F. Mitelman, B. Johansson, and F. Mertens (Eds.), \"Mitelman Database of Chromosome Aberrations and Gene Fusions in Cancer,\" 2026. [Online]. Available: https://mitelmandatabase.isb-cgc.org"),
    ("iscn2020", "J. McGowan-Jordan, R. J. Hastings, and S. Moore (Eds.), ISCN 2020: An International System for Human Cytogenomic Nomenclature. Basel: Karger, 2020."),
    ("johansson2002", "B. Johansson, T. Fioretos, and F. Mitelman, \"Cytogenetic and molecular genetic evolution of chronic myeloid leukemia,\" Acta Haematologica, vol. 107, no. 2, pp. 76-94, 2002."),
    ("grimwade2010", "D. Grimwade et al., \"Refinement of cytogenetic classification in acute myeloid leukemia,\" Blood, vol. 116, no. 3, pp. 354-365, 2010."),
    ("desper1999", "R. Desper et al., \"Inferring tree models for oncogenesis from comparative genome hybridization data,\" Journal of Computational Biology, vol. 6, no. 1, pp. 37-51, 1999."),
    ("beerenwinkel2007", "N. Beerenwinkel, N. Eriksson, and B. Sturmfels, \"Conjunctive Bayesian networks,\" Bernoulli, vol. 13, no. 4, pp. 893-909, 2007."),
    ("gerstung2009", "M. Gerstung, M. Baudis, H. Moch, and N. Beerenwinkel, \"Quantifying cancer progression with conjunctive Bayesian networks,\" Bioinformatics, vol. 25, no. 21, pp. 2809-2815, 2009."),
    ("schill2020", "R. Schill, S. Solbrig, T. Wettig, and R. Spang, \"Modelling cancer progression using Mutual Hazard Networks,\" Bioinformatics, vol. 36, no. 1, pp. 241-249, 2020."),
    ("hoglund2001", "M. Höglund et al., \"Multivariate analyses of genomic imbalances in solid tumors reveal distinct and converging pathways of karyotypic evolution,\" Genes, Chromosomes and Cancer, vol. 31, no. 2, pp. 156-171, 2001."),
    ("cytoconverter", "J. Wang and T. LaFramboise, \"CytoConverter: a web-based tool to convert karyotypes to genomic coordinates,\" BMC Bioinformatics, vol. 20, art. 467, 2019."),
    ("dohner2022", "H. Döhner et al., \"Diagnosis and management of AML in adults: 2022 recommendations from an international expert panel on behalf of the ELN,\" Blood, vol. 140, no. 12, pp. 1345-1377, 2022."),
    ("bradleyterry1952", "R. A. Bradley and M. E. Terry, \"Rank analysis of incomplete block designs: I. The method of paired comparisons,\" Biometrika, vol. 39, no. 3/4, pp. 324-345, 1952."),
    ("bh1995", "Y. Benjamini and Y. Hochberg, \"Controlling the false discovery rate: a practical and powerful approach to multiple testing,\" Journal of the Royal Statistical Society: Series B, vol. 57, no. 1, pp. 289-300, 1995."),
    ("moorman2012", "A. V. Moorman, \"The clinical relevance of chromosomal and genomic abnormalities in B-cell precursor acute lymphoblastic leukaemia,\" Blood Reviews, vol. 26, no. 3, pp. 123-135, 2012."),
    ("dohner2000", "H. Döhner et al., \"Genomic aberrations and survival in chronic lymphocytic leukemia,\" New England Journal of Medicine, vol. 343, no. 26, pp. 1910-1916, 2000."),
    ("schanz2012", "J. Schanz et al., \"New comprehensive cytogenetic scoring system for primary myelodysplastic syndromes (MDS) and oligoblastic acute myeloid leukemia after MDS derived from an international database merge,\" Journal of Clinical Oncology, vol. 30, no. 8, pp. 820-829, 2012."),
]

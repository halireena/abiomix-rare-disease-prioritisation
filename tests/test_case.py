"""Phenotype-aware fold (acmg.case.classify_case): case-level inheritance + gated agent facts on top of the
phenotype-free points. Hermetic."""
import pandas as pd
from acmg import case, agent


def _free(rows):
    # rows: (variant_key, total_points, criteria)
    return pd.DataFrame(rows, columns=["variant_key", "total_points", "criteria"]).assign(gene="G", acmg_class="VUS")


def test_de_novo_lifts_vus_to_lp():
    # trio: phenotype-free VUS at 5 pts; de novo (PS2 +4) -> 9 -> Likely Pathogenic
    df = case.classify_case(_free([("1-100-A-G", 5, "PM2_Supporting, PP3_Strong")]), de_novo=["1-100-A-G"], case_structure="trio")
    r = df.iloc[0]
    assert r.case_total_points == 9 and r.case_class == "Likely Pathogenic" and "PS2" in r.case_criteria


def test_comp_het_and_segregation_points():
    df = case.classify_case(_free([("2-1-C-T", 4, "PM2_Supporting"), ("2-2-G-A", 2, "PM2_Supporting")]),
                            comp_het=["2-1-C-T"], segregating=["2-2-G-A"], case_structure="trio")
    d = df.set_index("variant_key")
    assert d.loc["2-1-C-T", "case_points"] == 2 and "PM3" in d.loc["2-1-C-T", "case_criteria"]
    assert d.loc["2-2-G-A", "case_points"] == 1 and "PP1" in d.loc["2-2-G-A", "case_criteria"]


def test_singleton_does_not_award_family_criteria():
    # de novo passed, but a SINGLETON can't confirm it -> no PS2 (stays a candidate, not points)
    df = case.classify_case(_free([("1-100-A-G", 5, "PM2_Supporting")]), de_novo=["1-100-A-G"], case_structure="singleton")
    assert df.iloc[0].case_points == 0 and df.iloc[0].case_class == "VUS"


def test_genotype_fit_singleton_recessive():
    cv = _free([("6-1-A-T", 5, "PM2"), ("6-2-C-G", 4, "PM2"), ("7-1-A-T", 5, "PM2"), ("8-1-A-T", 6, "PVS1")])
    cv["gene"] = ["ARG", "ARG", "ARG", "ADG"]     # ARG=recessive gene (2 hets -> comp-het candidate), ADG=dominant
    zyg = {"6-1-A-T": "het", "6-2-C-G": "het", "7-1-A-T": "hom", "8-1-A-T": "het"}
    moi = {"ARG": "Autosomal recessive", "ADG": "Autosomal dominant"}
    d = case.classify_case(cv, zygosity=zyg, gene_moi=moi, case_structure="singleton").set_index("variant_key")
    assert d.loc["6-1-A-T", "genotype_fit"] == "biallelic_candidate"     # two hets in a recessive gene (unphased)
    assert d.loc["7-1-A-T", "genotype_fit"] == "biallelic_homozygous"    # hom in a recessive gene -> reportable
    assert d.loc["8-1-A-T", "genotype_fit"] == "monoallelic_dominant"    # het in a dominant gene -> reportable
    # a lone het in a recessive gene would be a carrier
    d2 = case.classify_case(_free([("9-1-A-T", 5, "PM2")]).assign(gene="ARG"),
                            zygosity={"9-1-A-T": "het"}, gene_moi={"ARG": "recessive"}).set_index("variant_key")
    assert d2.loc["9-1-A-T", "genotype_fit"] == "monoallelic_carrier"


def test_agent_fact_direct_points():
    df = case.classify_case(_free([("3-1-A-T", 6, "PS1")]),
                            agent_facts=[{"variant_key": "3-1-A-T", "criterion": "PS3", "points": 4}])
    assert df.iloc[0].case_total_points == 10 and df.iloc[0].case_class == "Pathogenic"


def test_agent_fact_gated_only_when_approved(tmp_path):
    ledger = agent.AgentLedger(str(tmp_path / "runs.jsonl"))
    prop = agent.propose({"variant_key": "4-1-A-T", "gene": "G"}, llm=lambda p: "PS3 supported",
                         kind="literature-criterion", model="test", ledger=ledger)
    fact = {"variant_key": "4-1-A-T", "criterion": "PS3", "points": 4, "proposal": prop, "ledger": ledger}
    # pending -> contributes 0
    pend = case.classify_case(_free([("4-1-A-T", 4, "PS1")]), agent_facts=[fact])
    assert pend.iloc[0].case_points == 0
    # approved -> contributes 4
    agent.approve(prop, "curator", ledger=ledger)
    appr = case.classify_case(_free([("4-1-A-T", 4, "PS1")]), agent_facts=[fact])
    assert appr.iloc[0].case_points == 4


def test_ba1_stand_alone_not_overridden():
    df = case.classify_case(_free([("5-1-A-T", -8, "BA1")]), de_novo=["5-1-A-T"])  # de novo can't rescue BA1
    assert df.iloc[0].case_class == "Benign"


def test_agent_fact_family_criterion_gated_on_singleton():
    """A gated agent fact asserting a family criterion (PS2/PM3/PP1) must NOT count on a singleton — only the
    pedigree-supported structures. Closes the propose->approve gate-bypass."""
    v = pd.DataFrame([{"variant_key": "1-1-A-G", "gene": "G", "total_points": 5, "criteria": "PM2", "acmg_class": "VUS"}])
    fact = [{"variant_key": "1-1-A-G", "criterion": "PS2", "points": 4}]
    singleton = case.classify_case(v, case_structure="singleton", agent_facts=fact)
    assert singleton.iloc[0]["case_points"] == 0          # PS2 dropped: no parents on a singleton
    trio = case.classify_case(v, case_structure="trio", agent_facts=fact)
    assert trio.iloc[0]["case_points"] == 4               # PS2 counts where the pedigree supports it
    # a non-family agent fact (PP4 phenotype) still counts on a singleton
    pp4 = case.classify_case(v, case_structure="singleton", agent_facts=[{"variant_key": "1-1-A-G", "criterion": "PP4", "points": 1}])
    assert pp4.iloc[0]["case_points"] == 1


def test_dual_moi_lone_het_stays_reportable_dominant():
    """A gene curated BOTH AD and AR (moi 'AD,AR'): a lone het must be reportable (monoallelic_dominant), not
    demoted to monoallelic_carrier — the dominant path applies. Pure-AR lone het stays a carrier."""
    v = pd.DataFrame([{"variant_key": "1-1-A-G", "gene": "DUAL", "total_points": 0, "criteria": "", "acmg_class": "VUS"},
                      {"variant_key": "2-1-C-T", "gene": "PUREAR", "total_points": 0, "criteria": "", "acmg_class": "VUS"}])
    fit = case.genotype_fit(v, {"1-1-A-G": "het", "2-1-C-T": "het"}, {"DUAL": "AD,AR", "PUREAR": "AR"})
    g = dict(zip(fit["variant_key"], fit["genotype_fit"]))
    assert g["1-1-A-G"] == "monoallelic_dominant"   # AD+AR lone het -> reportable, not penalized
    assert g["2-1-C-T"] == "monoallelic_carrier"    # pure-AR lone het -> carrier

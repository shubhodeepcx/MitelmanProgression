"""Build clone-level and event-level tables from the Mitelman download."""
import pandas as pd
from . import io
from .iscn import parse_investigation, relate, ordered_pairs
from .events import events_of, is_complex


def iter_investigations(clones):
    for (ref, case, inv), g in clones.groupby(["RefNo", "CaseNo", "InvNo"], sort=False):
        yield ref, case, inv, list(g.clone)


def build_tables(data_dir=io.DEFAULT_DIR):
    cases = io.load_cases(data_dir)
    clones = io.load_clones(data_dir)
    inv_rows, event_rows, pair_rows = [], [], []
    for ref, case, inv, texts in iter_investigations(clones):
        cl = parse_investigation(texts)
        parent, depth = relate(cl)
        ab = [c for c in cl if not c.is_normal]
        n_edges = sum(p is not None for p in parent.values())
        ab_pairs = ordered_pairs(cl)
        ev_pairs = ordered_pairs(cl, key=events_of)
        ev_pairs_impl = ordered_pairs(cl, key=events_of, implied=True)
        n_unrelated_founders = sum(p is None for p in parent.values())
        inv_rows.append(dict(
            RefNo=ref, CaseNo=case, InvNo=inv, n_clones=len(cl), n_abnormal=len(ab),
            n_edges=n_edges, n_founders=n_unrelated_founders,
            uses_idem=any(c.uses_idem for c in cl),
            n_ab_pairs=len(ab_pairs), n_ev_pairs=len(ev_pairs),
            n_ev_pairs_implied=len(ev_pairs_impl),
            complex=any(is_complex(c.abn) for c in ab),
            n_aberrations_max=max((len(c.abn) for c in ab), default=0),
        ))
        for c in ab:
            for e in sorted({e for a in c.abn for e in events_of(a)}):
                event_rows.append(dict(RefNo=ref, CaseNo=case, InvNo=inv, clone_id=c.index,
                                       event=e, inferred_depth=depth[c.index],
                                       parent_clone=parent[c.index]))
        for e, l in ev_pairs:
            pair_rows.append(dict(RefNo=ref, CaseNo=case, InvNo=inv, earlier=e, later=l))
    inv = pd.DataFrame(inv_rows).merge(cases, on=["RefNo", "CaseNo"], how="left")
    ev = pd.DataFrame(event_rows).merge(cases[["RefNo", "CaseNo", "disease", "Series", "Year"]],
                                        on=["RefNo", "CaseNo"], how="left")
    pairs = pd.DataFrame(pair_rows).merge(cases[["RefNo", "CaseNo", "disease"]],
                                          on=["RefNo", "CaseNo"], how="left")
    for df in (inv, ev, pairs):
        df.insert(0, "case_id", df.RefNo.astype(str) + ":" + df.CaseNo.astype(str))
    return cases, inv, ev, pairs

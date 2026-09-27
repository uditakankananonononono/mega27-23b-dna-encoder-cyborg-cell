#!/usr/bin/env python3
"""Verdict foldback: #7 evolutionary stability (mutation accumulation, event-driven),
#1/#12 realistic error profiles, #11 injection-recovery proxy.
Event-driven: each mutation event = one subst/indel at a random position;
events-per-generation ~ Poisson(L*(mu+indel)); we recover after each event, so
generations-to-loss = events-to-loss / (L*(mu+indel)).
"""
import json, pathlib, sys, os
sys.path.insert(0, os.path.dirname(__file__))
import numpy as np
from rs_inner import rsns_encode, rsns_decode

BASES = "ACGT"

def one_event(a, rng, indel_frac=0.1):
    i = rng.integers(0, len(a))
    r = rng.random()
    if r < indel_frac / 2:
        return a[:i] + a[i+1:]                       # deletion
    if r < indel_frac:
        return a[:i] + rng.choice(list(BASES)) + a[i:]     # insertion
    return a[:i] + rng.choice([b for b in BASES if b != a[i]]) + a[i+1:]

def recovers(msg, dna):
    try:
        return rsns_decode(dna) == msg
    except Exception:
        return False

def mutate_profile(seq, sub_p, indel_p, rng, homo_mult=1.0):
    out = []
    for i, ch in enumerate(seq):
        homo = i > 0 and seq[i-1] == ch
        m = homo_mult if homo else 1.0
        r = rng.random()
        if r < sub_p * m:
            out.append(rng.choice([b for b in BASES if b != ch]))
        elif r < (sub_p + indel_p) * m:
            if rng.random() < 0.5: pass
            else: out.extend((ch, rng.choice(list(BASES))))
        else:
            out.append(ch)
    return "".join(out)

def main():
    rng = np.random.default_rng(11)
    # #7 event-driven mutation accumulation, coded payloads, 12 trials per mu
    acc = []
    for mu in (1e-6, 1e-5, 1e-4):
        ev = []
        for t in range(8):
            msg = rng.bytes(256)
            dna = rsns_encode(msg); L = len(dna)
            n = 0
            while n <= 200:
                dna = one_event(dna, rng); n += 1
                if not recovers(msg, dna): break
            lam = L * (mu + mu / 10)
            ev.append({"events_to_loss": n if n <= 200 else ">200",
                       "expected_generations_to_loss": (n / lam) if n <= 200 else None})
        gens = [e["expected_generations_to_loss"] for e in ev if e["expected_generations_to_loss"]]
        acc.append({"mu_per_base_gen": mu, "trials": len(ev),
                    "median_events_to_loss": float(np.median([e["events_to_loss"] if isinstance(e["events_to_loss"], int) else 200 for e in ev])),
                    "median_expected_generations_to_loss": float(np.median(gens)) if gens else ">200/lambda",
                    "censored_trials": sum(1 for e in ev if not isinstance(e["events_to_loss"], int))})
        print(acc[-1], flush=True)
    # uncoded intact-probability reference at mu=1e-5, L=2048
    curve = {str(g): round((1 - 1e-5) ** (2048 * g), 6) for g in (10, 100, 500, 1000)}
    # #1/#11/#12 profile injection-recovery
    profs = {
        "illumina_like": dict(sub_p=1e-3, indel_p=1e-5, homo_mult=2.0),
        "nanopore_like": dict(sub_p=2e-2, indel_p=3e-2, homo_mult=3.0),
        "synthesis_like": dict(sub_p=5e-3, indel_p=1e-3, homo_mult=2.5),
    }
    inj = []
    for name, pr in profs.items():
        ok, dens = 0, []
        for t in range(16):
            msg = rng.bytes(256)
            dna = rsns_encode(msg)
            dens.append(8.0 * 256 / len(dna))
            if recovers(msg, mutate_profile(dna, pr["sub_p"], pr["indel_p"], rng, pr["homo_mult"])):
                ok += 1
        inj.append({"profile": name, **pr, "recovery_256B": round(ok / 16, 4),
                    "density_bits_per_base": round(float(np.mean(dens)), 4)})
        print(inj[-1], flush=True)
    out = {"mutation_accumulation": acc,
           "uncoded_intact_probability_mu1e-5_L2048": curve,
           "injection_recovery": inj,
           "note": "orders-of-magnitude error profiles from published Illumina/Nanopore/synthesis literature; computational proxy, not wet-lab validation"}
    pathlib.Path("results/evolution_error_models.json").write_text(json.dumps(out, indent=1))
    print("saved results/evolution_error_models.json")

if __name__ == "__main__":
    main()

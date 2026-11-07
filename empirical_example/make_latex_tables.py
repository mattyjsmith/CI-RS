#!/usr/bin/env python3
"""Builds the LaTeX table bodies for the empirical example from the CSV outputs of the R scripts.
Run from this folder:  python3 make_latex_tables.py   (rerun after any new run of 03_ or 05_ scripts).
Writes paper/tab_*_body.tex (tabular environments only; captions and labels live in the manuscript)."""
import csv, os

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "paper"); os.makedirs(OUT, exist_ok=True)
def read(f): return list(csv.DictReader(open(os.path.join(HERE, f), encoding="utf-8")))
def write(name, text): open(os.path.join(OUT, name), "w", encoding="utf-8").write(text)

def num(x, d=1, scale=1.0):
    v = round(float(x) * scale, d)
    if v == 0: v = 0.0
    s = f"{abs(v):.{d}f}"
    return f"$-${s}" if v < 0 else s
def est_ci(e, lo, hi, d=1, scale=1.0): return f"{num(e, d, scale)} ({num(lo, d, scale)}, {num(hi, d, scale)})"
DOSES = ["0.2 mg estrogen", "1.0 mg estrogen", "5.0 mg estrogen"]
DOSE_LAB = {d: d.replace(" estrogen", "") + " vs placebo" for d in DOSES}
METHODS = ["Aalen-Johansen (non-parametric)", "Cause-specific g-formula", "Relative survival g-formula"]
METH_LAB = {METHODS[0]: "Aalen--Johansen", METHODS[1]: "Cause-specific", METHODS[2]: "Relative survival"}

# ---- 1. weights against the observed causes of death (by comparison and arm) ----
w = read("weights_diagnostic.csv")
rows = []
for d in DOSES:
    for r in [x for x in w if x["dose"] == d]:
        rows.append(" & ".join([DOSE_LAB[d] if r["arm"] == "Placebo" else "", r["arm"], r["deaths"], r["observed_prostate"],
                                num(r["sum_wE"]), r["observed_other"], num(r["sum_wP"]), num(r["expected_other_life_table"]),
                                num(r["observed_over_expected_other"])]) + r"\\")
    rows[-1] += r" \hline"
write("tab_weights_body.tex", "\\begin{tabular}{llrrrrrrr}\n\\hline\n"
      "Comparison & Arm & Deaths & \\multicolumn{2}{c}{Prostate cancer deaths} & \\multicolumn{2}{c}{Other-cause deaths} & Expected other-cause & Observed / \\\\\n"
      " & & & Observed & Sum of weights & Observed & Sum of weights & (life table) & expected \\\\ \\hline\n" + "\n".join(rows) + "\n\\end{tabular}\n")

# ---- 2. observed vs expected other-cause deaths by subgroup (appendix) ----
ov = {r["Trt"]: r for r in read("observed_vs_expected_overall.csv")}
sub = read("observed_vs_expected_by_subgroup.csv")
def cell(r): return f"{r['n']} & {r['observed_other']} & {num(r['expected_other'])} & {num(r['observed_over_expected'])}"
order = [("CVD history", ["no", "yes"]), ("ECG", ["normal", "abnormal"]), ("Performance status", ["normal activity", "limited"]),
         ("Stage", ["stage 3", "stage 4"]), ("Age", ["<60", "60-74", "75+"])]
rows = [f"Overall & & {ov['placebo']['n']} & {ov['placebo']['observed_other']} & {num(ov['placebo']['expected_other'])} & {num(ov['placebo']['observed_over_expected'])}"
        f" & {ov['All four arms']['n']} & {ov['All four arms']['observed_other']} & {num(ov['All four arms']['expected_other'])} & {num(ov['All four arms']['observed_over_expected'])}\\\\ \\hline"]
for g, levels in order:
    for i, lv in enumerate(levels):
        p = [x for x in sub if x["Trt"] == "placebo" and x["subgroup"] == g and x["level"] == lv][0]
        a = [x for x in sub if x["Trt"] == "All four arms" and x["subgroup"] == g and x["level"] == lv][0]
        rows.append(f"{g if i == 0 else ''} & {lv.replace('-', '--').replace('<', '$<$')} & {cell(p)} & {cell(a)}\\\\")
    rows[-1] += r" \hline"
write("tab_subgroups_body.tex", "\\begin{tabular}{llrrrrrrrr}\n\\hline\n"
      " & & \\multicolumn{4}{c}{Placebo arm} & \\multicolumn{4}{c}{All four arms} \\\\\n"
      "Subgroup & Level & $n$ & Observed & Expected & Ratio & $n$ & Observed & Expected & Ratio \\\\ \\hline\n" + "\n".join(rows) + "\n\\end{tabular}\n")

# ---- 3. 50-month risks by arm (appendix) ----
rk = read("risks_50m.csv")
rows = []
for d in DOSES:
    for i, m in enumerate(METHODS):
        p = [x for x in rk if x["dose"] == d and x["method"] == m and x["arm"] == "Placebo"][0]
        e = [x for x in rk if x["dose"] == d and x["method"] == m and x["arm"] == "Estrogen"][0]
        rows.append(f"{DOSE_LAB[d] if i == 0 else ''} & {METH_LAB[m]} & {num(p['prostate_cancer_death'], 1, 100)} & {num(e['prostate_cancer_death'], 1, 100)}"
                    f" & {num(p['other_cause_death'], 1, 100)} & {num(e['other_cause_death'], 1, 100)}\\\\")
    rows[-1] += r" \hline"
write("tab_risks_body.tex", "\\begin{tabular}{llrrrr}\n\\hline\n & & \\multicolumn{2}{c}{Prostate cancer death (\\%)} & \\multicolumn{2}{c}{Other-cause death (\\%)} \\\\\n"
      "Comparison & Method & Placebo & Estrogen & Placebo & Estrogen \\\\ \\hline\n" + "\n".join(rows) + "\n\\end{tabular}\n")

# ---- 4. risk differences and risk ratios (main text) ----
res = read("results_empirical.csv")
def get(d, m, e):
    r = [x for x in res if x["dose"] == d and x["method"] == m and x["estimand"] == e][0]
    return r["estimate"], r["lower"], r["upper"]
rows = []
for d in DOSES:
    for i, m in enumerate(METHODS):
        c = [est_ci(*get(d, m, "Prostate cancer death: risk difference"), 1, 100), est_ci(*get(d, m, "Prostate cancer death: risk ratio"), 2),
             est_ci(*get(d, m, "Other-cause death: risk difference"), 1, 100), est_ci(*get(d, m, "Other-cause death: risk ratio"), 2)]
        rows.append(f"{DOSE_LAB[d] if i == 0 else ''} & {METH_LAB[m]} & " + " & ".join(c) + r"\\")
    rows[-1] += r" \hline"
write("tab_results_body.tex", "\\begin{tabular}{llcccc}\n\\hline\n & & \\multicolumn{2}{c}{Prostate cancer death} & \\multicolumn{2}{c}{Other-cause death} \\\\\n"
      "Comparison & Method & Risk difference (\\%) & Risk ratio & Risk difference (\\%) & Risk ratio \\\\ \\hline\n" + "\n".join(rows) + "\n\\end{tabular}\n")

# ---- 5. sensitivity to the level of the expected rates (appendix) ----
se = read("sensitivity_lifetable_estimates.csv"); sw = read("sensitivity_lifetable_weights.csv")
rows = []
for d in DOSES:
    factors = sorted({float(x["mx_multiplier"]) for x in se if x["dose"] == d})
    for i, f in enumerate(factors):
        def pick(e): return [x for x in se if x["dose"] == d and abs(float(x["mx_multiplier"]) - f) < 1e-9 and x["estimand"] == e][0]
        wr = [x for x in sw if x["dose"] == d and abs(float(x["mx_multiplier"]) - f) < 1e-9][0]
        cal = pick("RD")["calibrated"].strip().upper() == "TRUE"
        lab = f"{f:g}" + (r"$^{\dagger}$" if cal else "")
        a, b = pick("RD"), pick("RD_comp")
        rows.append(f"{DOSE_LAB[d] if i == 0 else ''} & {lab} & {est_ci(a['estimate'], a['lower'], a['upper'], 1, 100)} & {est_ci(b['estimate'], b['lower'], b['upper'], 1, 100)}"
                    f" & {num(wr['sum_wE'])} ({wr['observed_prostate']}) & {num(wr['sum_wP'])} ({wr['observed_other']})\\\\")
    rows[-1] += r" \hline"
write("tab_sensitivity_body.tex", "\\begin{tabular}{llcccc}\n\\hline\n & & \\multicolumn{2}{c}{Risk difference (\\%)} & \\multicolumn{2}{c}{Sum of weights (observed deaths)} \\\\\n"
      "Comparison & Factor & Prostate cancer death & Other-cause death & Prostate cancer & Other-cause \\\\ \\hline\n" + "\n".join(rows) + "\n\\end{tabular}\n")
print("written to", OUT, ":", sorted(os.listdir(OUT)))

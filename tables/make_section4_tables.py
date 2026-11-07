#!/usr/bin/env python3
"""Builds the LaTeX table bodies for Section 4 (simulation results) from the performance CSV files.
Run:  python3 make_section4_tables.py   (rerun after any new simulation run). Writes section4/tab_s4_*_body.tex."""
import csv, os
U = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "simulations") + os.sep
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "section4"); os.makedirs(OUT, exist_ok=True)
def rd(f): return list(csv.DictReader(open(U + f, encoding="utf-8")))
def get(rows, **kw):
    r = [x for x in rows if all(x[k] == v for k, v in kw.items())]; assert len(r) == 1, (kw, len(r)); return r[0]
def num(x, d=1, scale=1.0, plus=False):
    v = round(float(x) * scale, d)
    if v == 0: v = 0.0
    s = f"{abs(v):.{d}f}"
    if v < 0: return f"$-${s}"
    return ("+" + s) if plus else s
def write(name, text): open(os.path.join(OUT, name), "w", encoding="utf-8").write(text)
ESTS = ["Total RD", "Total logRR", "Competing RD", "Competing logRR"]

# ---- A. main simulation: relative bias (MCSE) ----
M = rd("main/performance_bias.csv"); rows = []
for est in ESTS:
    for i, n in enumerate(["1000", "3000", "5000"]):
        cells = []
        for m in "12345":
            r = get(M, stat="rbias", method=m, estimand=est, SampleSize=n)
            cells.append(f"{num(r['est'], 1, 100)} ({num(r['mcse'], 1, 100).replace('$-$', '')})")
        rows.append(f"{est if i == 0 else ''} & {n} & " + " & ".join(cells) + r"\\")
    rows[-1] += r" \hline"
write("tab_s4_main_body.tex", "\\begin{tabular}{llccccc}\n\\hline\nEstimand & $n$ & Method 1 & Method 2 & Method 3 & Method 4 & Method 5 \\\\ \\hline\n" + "\n".join(rows) + "\n\\end{tabular}\n")

# ---- B. treatment effect on other-cause death ----
O = rd("other_cause/performance_othercause_spline.csv")
SC = [("base_HR1", "Usual age, HR 1.0 (main run)"), ("base_HR1.3", "Usual age, HR 1.3"), ("base_HR1.65", "Usual age, HR 1.65"),
      ("older_HR1", "Older, HR 1.0"), ("older_HR1.3", "Older, HR 1.3"), ("older_HR1.65", "Older, HR 1.65")]
rows = []
for sc, lab in SC:
    for j, meth in enumerate(["Cause-specific", "Relative survival"]):
        c = []
        for est in ESTS:
            b = get(O, stat="bias", method=meth, estimand=est, Scenario=sc)
            t = f"{num(b['est'], 4, 1, True)} ({float(b['mcse']):.4f})"
            if est == "Total RD": t += f" [{num(get(O, stat='rbias', method=meth, estimand=est, Scenario=sc)['est'], 1, 100)}\\%]"
            c.append(t)
        rows.append(f"{lab if j == 0 else ''} & {meth} & " + " & ".join(c) + r"\\")
    rows[-1] += r" \hline"
write("tab_s4_other_body.tex", "\\begin{tabular}{llcccc}\n\\hline\nScenario & Method & Total RD & Total logRR & Competing RD & Competing logRR \\\\ \\hline\n" + "\n".join(rows) + "\n\\end{tabular}\n")

# ---- C. cancer aggressiveness ----
C = rd("cancer_type/performance_cancertype.csv"); T = rd("cancer_type/true_values_cancertype.csv"); rows = []
for sc in ["aggressive", "medium", "indolent"]:
    share = float(get(T, Scenario=sc)["cancer_share_of_deaths_untreated"])
    for j, meth in enumerate(["Cause-specific", "Relative survival"]):
        rb = get(C, stat="rbias", method=meth, estimand="Total RD", Scenario=sc); es = get(C, stat="empse", method=meth, estimand="Total RD", Scenario=sc)
        rl = get(C, stat="rbias", method=meth, estimand="Total logRR", Scenario=sc)
        cb = get(C, stat="bias", method=meth, estimand="Competing RD", Scenario=sc); ce = get(C, stat="empse", method=meth, estimand="Competing RD", Scenario=sc)
        rows.append(f"{sc.capitalize() + ' (' + num(share, 0, 100) + chr(92) + '%)' if j == 0 else ''} & {meth} & "
                    f"{num(rb['est'], 1, 100)} ({num(rb['mcse'], 1, 100).replace('$-$', '')}) & {float(es['est']):.4f} & "
                    f"{num(rl['est'], 1, 100)} ({num(rl['mcse'], 1, 100).replace('$-$', '')}) & "
                    f"{num(cb['est'], 4, 1, True)} ({float(cb['mcse']):.4f}) & {float(ce['est']):.4f}" + r"\\")
    rows[-1] += r" \hline"
write("tab_s4_cancer_body.tex", "\\begin{tabular}{llccccc}\n\\hline\n & & \\multicolumn{2}{c}{Total RD} & Total logRR & \\multicolumn{2}{c}{Competing RD} \\\\\n"
      "Scenario (cancer share of deaths) & Method & Relative bias (\\%) & Empirical SE & Relative bias (\\%) & Bias & Empirical SE \\\\ \\hline\n" + "\n".join(rows) + "\n\\end{tabular}\n")

# ---- D. time model (n = 5000) ----
Tm = rd("time_model/performance_timemodel.csv"); rows = []
for meth in ["Cause-specific", "Relative survival"]:
    for j, (tm, tl) in enumerate([("linear", "Linear"), ("factor", "Interval-specific intercepts"), ("spline", "Spline (4 df)")]):
        c = []
        for est in ESTS:
            r = get(Tm, stat="rbias", method=meth, estimand=est, TimeModel=tm)
            c.append(f"{num(r['est'], 2, 100)} ({num(r['mcse'], 2, 100).replace('$-$', '')})")
        rows.append(f"{meth if j == 0 else ''} & {tl} & " + " & ".join(c) + r"\\")
    rows[-1] += r" \hline"
write("tab_s4_time_body.tex", "\\begin{tabular}{llcccc}\n\\hline\nMethod & Time model & Total RD & Total logRR & Competing RD & Competing logRR \\\\ \\hline\n" + "\n".join(rows) + "\n\\end{tabular}\n")
print("written:", sorted(os.listdir(OUT)))

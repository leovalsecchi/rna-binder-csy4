import io
ra="run_all.py"
r=io.open(ra,encoding="utf-8").read()
old='    "plot_19_refolder_iptm_pairs", "plot_20_refolder_dockq_pairs",\n]'
new='    "plot_19_refolder_iptm_pairs", "plot_20_refolder_dockq_pairs",\n    "plot_21_violin_iptm_refolder",\n]'
if "plot_21_violin_iptm_refolder" not in r:
    assert r.count(old)==1, r.count(old)
    io.open(ra,"w",encoding="utf-8").write(r.replace(old,new))
    print("registered plot_21 in run_all.PLOTS")
else:
    print("plot_21 already registered")

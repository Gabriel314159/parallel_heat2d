#!/usr/bin/env python3
"""Run actual time-dependent Model B, keeping raw results, logs and median timings."""
import argparse
import json
import statistics
from pathlib import Path
import sys
from runner_common import run_case, write_csv


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--build-dir", default="build")
    p.add_argument("--output", default="results/model_b")
    p.add_argument("--launcher", default="mpirun --bind-to none")
    p.add_argument("--suite", choices=("smoke","benchmark"), default="smoke")
    p.add_argument("--layouts", default="serial,1x1,4x1,4x2", help="serial or MPI-ranks x threads")
    for name,value in (("nx",32),("ny",32),("nz",16),("steps",100),("spinup",0),("progress",0),("max_iter",10000)):
        p.add_argument("--"+name, type=int, default=value)
    for name,value in (("lx",1.),("ly",1.),("lz",.5),("nu",1e-3),("kappa",1e-4),("force",.8),("tol",1e-10),("dt",0.)):
        p.add_argument("--"+name, type=float, default=value)
    p.add_argument("--profile", choices=("surface","uniform"), default="surface")
    p.add_argument("--decomp", choices=("slabs","pencils","blocks"), default="blocks")
    p.add_argument("--repeats", type=int, default=5)
    p.add_argument("--warmups", type=int, default=1)
    p.add_argument("--timeout", type=int, default=1800)
    args = p.parse_args()
    if args.repeats<1 or args.warmups<0:
        p.error("invalid repetition count")
    out = Path(args.output).resolve()
    out.mkdir(parents=True,exist_ok=True)
    if (out/"raw.csv").exists():
        p.error("output already contains raw.csv; choose a new output directory")
    (out/"status.json").write_text(json.dumps(dict(status="RUNNING"))+"\n")
    base = {key:getattr(args,key) for key in ("nx","ny","nz","lx","ly","lz","steps","spinup","nu","kappa","force","tol","dt","profile","decomp","progress","max_iter")}
    reps = 1 if args.suite=="smoke" else args.repeats
    warmups = 0 if args.suite=="smoke" else args.warmups
    metadata = dict(vars(args), effective_repeats=reps, effective_warmups=warmups,
                    timing="loop_s excludes setup, spinup and field/CSV writes; includes per-step stability/divergence checks",
                    trace_enabled=args.suite=="smoke", device="CPU", pressure_boundary="homogeneous Neumann, zero mean gauge")
    (out/"setup.json").write_text(json.dumps(metadata,indent=2)+"\n")
    rows=[]
    for layout in args.layouts.split(","):
        layout=layout.strip()
        if layout=="serial":
            ranks,threads=0,1
        else:
            ranks,threads=map(int,layout.split("x"))
            if ranks<1 or threads<1:
                raise ValueError("invalid layout "+layout)
        for phase,count in (("warmup",warmups),("measured",reps)):
            for repeat in range(1,count+1):
                name="{}_{}_{:02d}".format(layout,phase,repeat)
                values=dict(base)
                if args.suite=="smoke":
                    values.update(trace=out/(name+"_trace.csv"), output=out/(name+"_field.csv"))
                print("Running "+name, flush=True)
                r=run_case(args.build_dir,args.launcher,out,name,values,ranks,threads,args.timeout)
                r.update(layout=layout,phase=phase,repeat=repeat)
                rows.append(r)
                write_csv(out/"raw.csv",rows)
                if r["cg_solves"]!=args.steps:
                    raise RuntimeError("expected one pressure solve per timed step")
                if r["mass_abs_error"]>1e-10*max(abs(r["mass_initial"]),1e-12):
                    raise RuntimeError("pigment mass conservation failed")
                if args.profile=="surface" and args.force!=0 and r["w_max"]<=1e-12:
                    raise RuntimeError("expected nonzero vertical velocity")
                if r["c_min"] < -1e-9 or r["c_max"] > 1+1e-9:
                    raise RuntimeError("pigment bounds failed")
    summaries=[]
    for layout in dict.fromkeys(r["layout"] for r in rows):
        measured=[r for r in rows if r["layout"]==layout and r["phase"]=="measured"]
        row=dict(layout=layout,repetitions=len(measured),ranks=measured[0]["ranks"],threads=measured[0]["threads"],nodes=measured[0]["nodes"])
        for key in ("loop_s","poisson_s","spinup_s","end_to_end_s","cg_iterations","max_divergence","mass_abs_error","w_max"):
            row[key+"_median"]=statistics.median(r[key] for r in measured)
        row["loop_s_min"]=min(r["loop_s"] for r in measured)
        row["loop_s_max"]=max(r["loop_s"] for r in measured)
        summaries.append(row)
    reference=next((s["loop_s_median"] for s in summaries if s["layout"]=="serial"),None)
    if reference is not None and args.suite=="benchmark":
        for s in summaries:
            s["speedup_vs_serial"]=reference/s["loop_s_median"]
            s["efficiency_vs_serial"]=s["speedup_vs_serial"]/(s["ranks"]*s["threads"])
    write_csv(out/"summary.csv",summaries)
    (out/"status.json").write_text(json.dumps(dict(status="PASS", runs=len(rows)),indent=2)+"\n")
    print("PASS: "+str(out/"summary.csv"),flush=True)


if __name__=="__main__":
    try:
        main()
    except Exception as exc:
        print("ERROR: "+str(exc),file=sys.stderr)
        sys.exit(1)

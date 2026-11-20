#!/usr/bin/env bash
# Job status for Snakemake's cluster-generic executor: prints running | success | failed.
# sacct is not available on this cluster, so use squeue, then scontrol for finished jobs.
jobid="$1"

state=$(squeue -h -j "$jobid" -o %T 2>/dev/null)
case "$state" in
    PENDING|CONFIGURING|RUNNING|COMPLETING|SUSPENDED|REQUEUED|RESIZING|REQUEUE_HOLD|REQUEUE_FED)
        echo running; exit 0 ;;
esac

info=$(scontrol show job "$jobid" 2>/dev/null)
st=$(sed -n 's/.*JobState=\([A-Z_]*\).*/\1/p' <<< "$info" | head -n 1)
ec=$(sed -n 's/.* ExitCode=\([0-9]*\):.*/\1/p' <<< "$info" | head -n 1)

case "$st" in
    COMPLETED) [ "${ec:-0}" = "0" ] && echo success || echo failed ;;
    PENDING|CONFIGURING|RUNNING|COMPLETING|SUSPENDED|REQUEUED) echo running ;;
    "")
        # Job already purged from the controller: report success and let Snakemake
        # verify the outputs (missing outputs still fail the job).
        echo success ;;
    *) echo failed ;;   # FAILED, CANCELLED, TIMEOUT, OUT_OF_MEMORY, NODE_FAIL, ...
esac

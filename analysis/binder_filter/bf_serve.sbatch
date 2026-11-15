#!/bin/bash
# bf_serve.sbatch — run the binder_filter web app on a compute node.
#
# Not a login node: the app holds ~1 GB resident and lives for hours, which is exactly what
# login nodes kill. A compute node also means several people can each run their own instance
# without fighting over a port.
#
# The app binds 127.0.0.1 on the compute node and is reached by an ssh tunnel, so it is never
# exposed to the cluster network.
#SBATCH --job-name=bf_serve
#SBATCH --partition=normal.4h
#SBATCH --time=03:59:00
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=2G
#SBATCH --output=%u_bf_serve_%j.log
set -uo pipefail

APP=/cluster/work/wicky/lvalsecchi/rna_binder/tools/binder_filter
STORE=/cluster/work/wicky/lvalsecchi/rna_binder/analyses/cross_campaign/filter_store_v2v8
INFO="$HOME/.binder_filter_serve"

# a free port on THIS node, so two instances on the same node cannot collide
PORT=$(python3 - <<'EOF'
import socket
s = socket.socket(); s.bind(("127.0.0.1", 0)); print(s.getsockname()[1]); s.close()
EOF
)

# values are quoted: bf-serve sources this file, and an unquoted timestamp (which contains a
# space) would be run as a command
cat > "$INFO" <<EOF
node="$(hostname)"
port="$PORT"
job="$SLURM_JOB_ID"
user="$(whoami)"
started="$(date '+%F %T')"
EOF

echo "=========================================================================="
echo " binder_filter serving on $(hostname):$PORT   (job $SLURM_JOB_ID)"
echo " store: $STORE"
echo ""
echo " From your OWN machine, open a tunnel:"
echo "   ssh -o StrictHostKeyChecking=accept-new -J euler \\"
echo "       -L 8010:127.0.0.1:$PORT $(whoami)@$(hostname)"
echo " then browse to  http://127.0.0.1:8010"
echo "=========================================================================="

cd "$APP" || exit 1
export BF_STORE="$STORE"
export BF_PORT="$PORT"
export BF_HOST=127.0.0.1
export BF_REMOTE=1          # ChimeraX is handled browser-side; this node cannot reach it
exec "$APP/venv/bin/python" "$APP/app.py"

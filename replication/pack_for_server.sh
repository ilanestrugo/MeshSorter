#!/bin/sh
# Build a self-contained bundle for the Linux server.
#
#   ./pack_for_server.sh
#
# Writes meshsorter-server.tar.gz containing the sources, the table scripts,
# geometry.py and the order-derived CSV, and nothing that has to be rebuilt:
# no compiled binaries, no results, no caches.  On the server:
#
#   tar xzf meshsorter-server.tar.gz && cd replication && ./run_server.sh
set -e
cd "$(dirname "$0")"

STAGE=$(mktemp -d)
DEST="$STAGE/replication"
mkdir -p "$DEST"

# C++ sources
cp meshsorter.cpp meshsorter_core.hpp certify_rep.cpp sweep_rep.cpp "$DEST"/

# every Python module: the expensive scripts import common.py and geometry.py,
# and approx_eval.py imports approx_model.py.  The cheap table scripts are
# included too so the bundle can reproduce anything, at no cost in size.
cp ./*.py "$DEST"/

# runners and documentation
cp run_server.sh run_all.sh README.md "$DEST"/ 2>/dev/null || true
chmod +x "$DEST"/run_server.sh

# the order-derived destination sequence.  order_derived.py looks for it
# beside itself first, so putting it inside the bundle makes it portable.
if [ -f olist_orders_ForRun.csv ]; then
    cp olist_orders_ForRun.csv "$DEST"/
elif [ -f ../olist_orders_ForRun.csv ]; then
    cp ../olist_orders_ForRun.csv "$DEST"/
else
    echo "WARNING: olist_orders_ForRun.csv not found; Table S12 will not run" >&2
fi

rm -rf "$DEST/__pycache__" "$DEST/_geometry_backup" "$DEST/results"
rm -f "$DEST/_apply_geometry.py" "$DEST/pack_for_server.sh"

tar czf meshsorter-server.tar.gz -C "$STAGE" replication
rm -rf "$STAGE"

echo "wrote meshsorter-server.tar.gz"
echo
echo "contents:"
tar tzf meshsorter-server.tar.gz | sed 's/^/  /'
echo
echo "on the server:"
echo "  tar xzf meshsorter-server.tar.gz"
echo "  cd replication"
echo "  ./run_server.sh --quick      # smoke test, a few minutes"
echo "  ./run_server.sh              # the real run"

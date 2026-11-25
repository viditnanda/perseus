#!/bin/sh
# Finite groups: integral homology, mod-p cohomology rings, and persistence along a series of groups.
# The inputs are written from GAP by tools/export_pcgroup.g, tools/export_sylow.g and tools/export_pcseries.g:
#   gap> Read("tools/export_pcgroup.g");  WritePerseusPcGroup(SmallGroup(16, 9), "q16_pcgroup.txt");
#   gap> Read("tools/export_sylow.g");    WritePerseusSylow(AlternatingGroup(5), "a5_sylow.txt");
#   gap> Read("tools/export_pcseries.g"); G := DihedralGroup(64);;
#   gap> WritePerseusQuotientSeries(G, LowerCentralSeries(G), "d64_lower_central.txt");
# Run:  sh examples/cli/09_finite_groups.sh
set -e
P=${PERSEUS:-perseus}
cd "$(dirname "$0")/../data"

echo "== Q_16 from its polycyclic presentation: H_n(G; Z) and the mod-2 cohomology ring through degree 6"
"$P" pcgroup q16_pcgroup.txt --maxdim 6 --ring 2

echo
echo "== A_5 through its Sylow subgroups"
"$P" sylow a5_sylow.txt --maxdim 6 --ring 2

echo
echo "== D_64 -> D_32 -> D_16 -> D_8 -> C_2 x C_2: the induced maps on H_n(-; F_2) and their persistence"
"$P" pcseries d64_lower_central.txt --maxdim 2

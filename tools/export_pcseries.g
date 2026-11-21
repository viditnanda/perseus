# export_pcseries.g -- write a chain of homomorphisms of finite solvable groups for `perseus pcseries` /
# perseus.pc_series: the maps induced on mod-p homology and their persistence (Ellis-King).
#
#   gap> Read("tools/export_pcseries.g");
#   gap> WritePerseusQuotientSeries(DihedralGroup(64), LowerCentralSeries(DihedralGroup(64)), "d64_L.txt");
#   gap> WritePerseusSubgroupSeries(SmallGroup(64, 138), LowerCentralSeries(SmallGroup(64, 138)), "s.txt");
#   gap> WritePerseusHoms([hom1, hom2, ...], "chain.txt");   # Source(hom_{i+1}) = Range(hom_i)
#
# A normal series S gives the chain of quotients G/N, N in S, from the largest quotient to the smallest, joined by
# the natural surjections (the trivial quotient G/G is left out, as in HAP's PersistentHomologyOfQuotientGroupSeries);
# a subgroup series gives the chain of subgroups from the smallest to G, joined by inclusions (the trivial subgroup
# is left out). Every group is written with a pcgs, so all must be solvable.
#
# File format (whitespace separated; lines starting with # are comments):
#   groups m
#   m groups, each as in a pcgroup file: n k / f_1 ... f_k / r_1 ... r_k / e / n rows of n entries, each followed
#   optionally by "alternative k' / pcgs / relative orders" (another pcgs; Perseus keeps the cheaper)
#   maps
#   m - 1 lines: the images in group i + 1 of the n_i elements of group i
# (elements of each group are numbered by position in AsSortedList)

if not IsBoundGlobal("PERSEUS_DIR_OF_") then
    BindGlobal("PERSEUS_DIR_OF_", function(f) local i; i := Length(f); while i > 0 and f[i] <> '/' do i := i - 1; od; if i = 0 then return "./"; fi; return f{[1..i]}; end);
fi;
if not IsBoundGlobal("PerseusShortPcgs") then      # tools/perseus_pcgs.g, next to this file
    Read(Concatenation(PERSEUS_DIR_OF_(INPUT_FILENAME()), "perseus_pcgs.g"));
fi;

PerseusWriteGroupBlock_ := function(out, G, el)
    local pc, n, row, ids, alt;
    pc := PerseusPcgs(G); n := Length(el);
    ids := l -> JoinStringsWithSeparator(List(l, x -> String(Position(el, x) - 1)), " ");
    AppendTo(out, n, " ", Length(pc), "\n", ids(pc), "\n", JoinStringsWithSeparator(List(PerseusRelativeOrders(G, pc), String), " "), "\n", Position(el, One(G)) - 1, "\n");
    for row in el do AppendTo(out, JoinStringsWithSeparator(List(el, y -> String(Position(el, row * y) - 1)), " "), "\n"); od;
    alt := PerseusAlternativePcgs(G, pc);
    if alt <> fail then AppendTo(out, "alternative ", Length(alt), "\n", ids(alt), "\n", JoinStringsWithSeparator(List(PerseusRelativeOrders(G, alt), String), " "), "\n"); fi;
end;

WritePerseusHoms := function(homs, file)
    local groups, els, out, i, G, h;
    groups := Concatenation([Source(homs[1])], List(homs, Range));
    for i in [1..Length(homs) - 1] do if Range(homs[i]) <> Source(homs[i + 1]) then Error("the homomorphisms do not compose"); fi; od;
    els := List(groups, AsSortedList);
    out := OutputTextFile(file, false); SetPrintFormattingStatus(out, false);
    PrintTo(out, "# perseus pcseries: ", JoinStringsWithSeparator(List(groups, StructureDescription), " -> "), "\n");
    AppendTo(out, "groups ", Length(groups), "\n");
    for i in [1..Length(groups)] do PerseusWriteGroupBlock_(out, groups[i], els[i]); od;
    AppendTo(out, "maps\n");
    for i in [1..Length(homs)] do
        h := homs[i];
        AppendTo(out, JoinStringsWithSeparator(List(els[i], x -> String(Position(els[i + 1], Image(h, x)) - 1)), " "), "\n");
    od;
    CloseStream(out);
end;

WritePerseusQuotientSeries := function(G, S, file)
    local N, quots, homs, i, nat, hs, Q;
    N := Filtered(SortedList(ShallowCopy(S), function(a, b) return Size(a) < Size(b); end), x -> Size(x) < Size(G));
    N := Set(N);  N := SortedList(N, function(a, b) return Size(a) < Size(b); end);
    nat := List(N, x -> NaturalHomomorphismByNormalSubgroup(G, x));
    quots := List(nat, h -> Image(h));
    homs := [];
    for i in [1..Length(N) - 1] do
        # G/N_i -> G/N_{i+1}, N_i < N_{i+1}
        hs := GeneratorsOfGroup(quots[i]);
        Add(homs, GroupHomomorphismByImagesNC(quots[i], quots[i + 1], hs, List(hs, q -> Image(nat[i + 1], PreImagesRepresentative(nat[i], q)))));
    od;
    if Length(homs) = 0 then Error("the series gives fewer than two nontrivial quotients"); fi;
    WritePerseusHoms(homs, file);
end;

WritePerseusSubgroupSeries := function(G, S, file)
    local N, homs, i;
    N := Filtered(Set(S), x -> Size(x) > 1);
    N := SortedList(N, function(a, b) return Size(a) < Size(b); end);
    homs := List([1..Length(N) - 1], i -> GroupHomomorphismByImagesNC(N[i], N[i + 1], GeneratorsOfGroup(N[i]), GeneratorsOfGroup(N[i])));
    if Length(homs) = 0 then Error("the series gives fewer than two nontrivial subgroups"); fi;
    WritePerseusHoms(homs, file);
end;

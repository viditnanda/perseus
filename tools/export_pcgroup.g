# export_pcgroup.g -- write a finite solvable group for `perseus pcgroup` / perseus.pc_group.
#
#   gap> Read("tools/export_pcgroup.g");
#   gap> WritePerseusPcGroup(SmallGroup(64, 138), "g64_138.txt");
#   gap> WritePerseusPcGroup(G, "g.txt", PcgsPCentralSeriesPGroup(G));   # any pcgs of G
#
# File format (whitespace separated; lines starting with # are comments):
#   n k                      order and pcgs length
#   f_1 ... f_k              pcgs elements (0-based ids)
#   r_1 ... r_k              relative orders
#   e                        identity
#   n rows of n entries      multiplication table: row a, column b holds a * b
# Elements are numbered by their position in AsSortedList(G).

WritePerseusPcGroup := function(arg)
    local G, file, pcgs, el, n, out, row;
    G := arg[1]; file := arg[2];
    if Length(arg) >= 3 then pcgs := arg[3]; else pcgs := Pcgs(G); fi;
    if pcgs = fail then Error("the group is not solvable: it has no pcgs"); fi;
    el := AsSortedList(G); n := Length(el);
    out := OutputTextFile(file, false); SetPrintFormattingStatus(out, false);
    PrintTo(out, "# perseus pcgroup: ", StructureDescription(G), "\n");
    AppendTo(out, n, " ", Length(pcgs), "\n");
    AppendTo(out, JoinStringsWithSeparator(List(pcgs, x -> String(Position(el, x) - 1)), " "), "\n");
    AppendTo(out, JoinStringsWithSeparator(List(RelativeOrders(pcgs), String), " "), "\n");
    AppendTo(out, Position(el, One(G)) - 1, "\n");
    for row in el do
        AppendTo(out, JoinStringsWithSeparator(List(el, y -> String(Position(el, row * y) - 1)), " "), "\n");
    od;
    CloseStream(out);
end;

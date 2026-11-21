# export_pcgroup.g -- write a finite solvable group for `perseus pcgroup` / perseus.pc_group.
#
#   gap> Read("tools/export_pcgroup.g");
#   gap> WritePerseusPcGroup(SmallGroup(64, 138), "g64_138.txt");
#   gap> WritePerseusPcGroup(G, "g.txt", Pcgs(G));   # any pcgs of G, as a pcgs or a list of elements
#
# By default the pcgs is a short one with cyclic steps (tools/perseus_pcgs.g), which can make Perseus much faster;
# PerseusUseShortPcgs := false; gives GAP's Pcgs(G) instead.
#
# File format (whitespace separated; lines starting with # are comments):
#   n k                      order and pcgs length
#   f_1 ... f_k              pcgs elements (0-based ids)
#   r_1 ... r_k              relative orders
#   e                        identity
#   n rows of n entries      multiplication table: row a, column b holds a * b
#   (optional, repeated) k' / f'_1 ... f'_k' / r'_1 ... r'_k'   an alternative pcgs; Perseus keeps the cheaper
# Elements are numbered by their position in AsSortedList(G).

if not IsBoundGlobal("PERSEUS_DIR_OF_") then
    BindGlobal("PERSEUS_DIR_OF_", function(f) local i; i := Length(f); while i > 0 and f[i] <> '/' do i := i - 1; od; if i = 0 then return "./"; fi; return f{[1..i]}; end);
fi;
if not IsBoundGlobal("PerseusShortPcgs") then      # tools/perseus_pcgs.g, next to this file
    Read(Concatenation(PERSEUS_DIR_OF_(INPUT_FILENAME()), "perseus_pcgs.g"));
fi;

WritePerseusPcGroup := function(arg)
    local G, file, pcgs, el, n, out, row, alt;
    G := arg[1]; file := arg[2];
    if Length(arg) >= 3 then pcgs := arg[3]; else pcgs := PerseusPcgs(G); fi;
    if pcgs = fail then Error("the group is not solvable: it has no pcgs"); fi;
    el := AsSortedList(G); n := Length(el);
    out := OutputTextFile(file, false); SetPrintFormattingStatus(out, false);
    PrintTo(out, "# perseus pcgroup: ", StructureDescription(G), "\n");
    AppendTo(out, n, " ", Length(pcgs), "\n");
    AppendTo(out, JoinStringsWithSeparator(List(pcgs, x -> String(Position(el, x) - 1)), " "), "\n");
    AppendTo(out, JoinStringsWithSeparator(List(PerseusRelativeOrders(G, AsList(pcgs)), String), " "), "\n");
    AppendTo(out, Position(el, One(G)) - 1, "\n");
    for row in el do
        AppendTo(out, JoinStringsWithSeparator(List(el, y -> String(Position(el, row * y) - 1)), " "), "\n");
    od;
    alt := PerseusAlternativePcgs(G, AsList(pcgs));
    if alt <> fail and Length(arg) < 3 then                  # an alternative pcgs: k / elements / relative orders
        AppendTo(out, "# alternative pcgs (GAP's Pcgs): Perseus keeps the cheaper\n", Length(alt), "\n");
        AppendTo(out, JoinStringsWithSeparator(List(alt, x -> String(Position(el, x) - 1)), " "), "\n");
        AppendTo(out, JoinStringsWithSeparator(List(PerseusRelativeOrders(G, alt), String), " "), "\n");
    fi;
    CloseStream(out);
end;

# export_sylow.g -- write any finite group for `perseus sylow` / perseus.finite_group, through its Sylow subgroups.
#
#   gap> Read("tools/export_sylow.g");
#   gap> WritePerseusSylow(AlternatingGroup(5), "a5.txt");
#   gap> WritePerseusSylow(MathieuGroup(11), "m11.txt", [2, 3]);   # only these primes
#
# Perseus never sees G itself. For each prime p it gets a Sylow p-subgroup P (a pc group: table and pcgs) and, for
# each double coset P g P with g outside P, the subgroup Q = P meet g P g^-1 of P with the conjugation map
# Q -> P, y -> g^-1 y g. By Cartan-Eilenberg, H^n(G; M) for p-local M is the set of x in H^n(P; M) whose restrictions
# to Q along the inclusion and along the conjugation agree, for every such g; dually, the p-part of H_n(G; Z) is a
# cokernel of the maps H_n(Q) -> H_n(P). G need not be solvable; only its Sylow subgroups are.
#
# File format (whitespace separated; lines starting with # are comments):
#   order N
#   then, for each prime p:
#     prime p
#     n k / f_1 ... f_k / r_1 ... r_k / e / n rows of n entries    P, exactly as in a pcgroup file
#     cosets m
#     (optional, repeated) alternative k' / pcgs / relative orders     another pcgs of P; Perseus keeps the cheaper
#     m times:  q kq / the q elements of Q / a pcgs of Q (kq elements) / the q images g^-1 y g, in the same order
#               (optional, repeated) alternative k' / pcgs                another pcgs of Q
#   (elements of Q and their images are ids of P's elements, numbered by position in AsSortedList(P))

if not IsBoundGlobal("PERSEUS_DIR_OF_") then
    BindGlobal("PERSEUS_DIR_OF_", function(f) local i; i := Length(f); while i > 0 and f[i] <> '/' do i := i - 1; od; if i = 0 then return "./"; fi; return f{[1..i]}; end);
fi;
if not IsBoundGlobal("PerseusShortPcgs") then      # tools/perseus_pcgs.g, next to this file
    Read(Concatenation(PERSEUS_DIR_OF_(INPUT_FILENAME()), "perseus_pcgs.g"));
fi;

WritePerseusSylow := function(arg)
    local G, file, primes, out, p, P, pc, el, n, row, reps, g, Q, qel, qpc, ids, alt;
    G := arg[1]; file := arg[2];
    if Length(arg) >= 3 then primes := arg[3]; else primes := Set(Factors(Size(G))); fi;
    out := OutputTextFile(file, false); SetPrintFormattingStatus(out, false);
    PrintTo(out, "# perseus sylow: ", StructureDescription(G), "\n");
    AppendTo(out, "order ", Size(G), "\n");
    for p in primes do
        if Size(G) mod p <> 0 then Error("the prime does not divide |G|"); fi;
        P := SylowSubgroup(G, p);
        pc := PerseusPcgs(P);
        el := AsSortedList(P); n := Length(el);
        ids := l -> JoinStringsWithSeparator(List(l, x -> String(Position(el, x) - 1)), " ");
        AppendTo(out, "prime ", p, "\n", n, " ", Length(pc), "\n", ids(pc), "\n");
        AppendTo(out, JoinStringsWithSeparator(List(PerseusRelativeOrders(P, pc), String), " "), "\n");
        AppendTo(out, Position(el, One(G)) - 1, "\n");
        for row in el do AppendTo(out, JoinStringsWithSeparator(List(el, y -> String(Position(el, row * y) - 1)), " "), "\n"); od;
        alt := PerseusAlternativePcgs(P, pc);
        if alt <> fail then AppendTo(out, "alternative ", Length(alt), "\n", ids(alt), "\n", JoinStringsWithSeparator(List(PerseusRelativeOrders(P, alt), String), " "), "\n"); fi;
        reps := Filtered(List(DoubleCosetRepsAndSizes(G, P, P), x -> x[1]), x -> not x in P);
        AppendTo(out, "cosets ", Length(reps), "\n");
        for g in reps do
            Q := Intersection(P, P^(g^-1));              # y in P with g^-1 y g in P
            qel := AsSortedList(Q);
            qpc := PerseusPcgs(Q);
            AppendTo(out, Length(qel), " ", Length(qpc), "\n", ids(qel), "\n", ids(qpc), "\n", ids(List(qel, y -> y^g)), "\n");
            alt := PerseusAlternativePcgs(Q, qpc);
            if alt <> fail then AppendTo(out, "alternative ", Length(alt), "\n", ids(alt), "\n"); fi;
        od;
    od;
    CloseStream(out);
end;

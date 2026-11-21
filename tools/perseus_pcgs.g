# perseus_pcgs.g -- a short pcgs with cyclic steps, for the Perseus exporters (read by them; you need not read it).
#
# Perseus resolves a solvable group through a pcgs f_1, ..., f_k: a sequence with G_i = <f_i, ..., f_k>, G_{i+1}
# normal in G_i and G_i / G_{i+1} cyclic of order r_i (any r_i, not only primes). Its Anick resolution has
# C(n + k - 1, k - 1) generators in degree n, so a short pcgs is much cheaper: the dihedral group of order 128 has
# a GAP pcgs with 7 steps of order 2, but f_1 = a reflection, f_2 = a rotation of order 64 is a pcgs with 2 steps.
#
# PerseusShortPcgs(G) builds one greedily, from the top: at each step, with H the current subgroup, it takes the
# largest cyclic quotient of H, of order e = exponent of H / H', and among the normal subgroups N of H with H / N
# cyclic of order e (preimages of complements of cyclic subgroups of order e of H / H') it keeps the one whose
# abelianisation has fewest cyclic factors, then the one of largest exponent. It is a heuristic, not a proven
# minimum; Perseus checks the pcgs it is given. Set PerseusUseShortPcgs := false; to have the exporters use GAP's
# Pcgs (prime steps) instead.


BindGlobal("PerseusShortPcgs", function(G)
    local seq, H, D, nat, A, e, Zs, seen, a, Z, C, N, sc, best, bestsc, bestx, comps, cnt;
    seq := []; H := G;
    while Size(H) > 1 do
        D := DerivedSubgroup(H);
        if Size(D) = Size(H) then Error("PerseusShortPcgs: the group is not solvable"); fi;
        nat := NaturalHomomorphismByNormalSubgroup(H, D);
        A := Image(nat); e := Exponent(A);
        Zs := []; seen := [];
        for a in A do
            if Order(a) = e then
                Z := Subgroup(A, [a]);
                if not ForAny(seen, z -> z = Z) then Add(seen, Z); Add(Zs, a); fi;
                if Length(Zs) >= 64 then break; fi;
            fi;
        od;
        best := fail; bestsc := fail; bestx := fail;
        for a in Zs do
            comps := ComplementClassesRepresentatives(A, Subgroup(A, [a])); cnt := 0;
            for C in comps do
                cnt := cnt + 1; if cnt > 64 then break; fi;
                N := PreImage(nat, C);
                sc := [Length(AbelianInvariants(N)), -Exponent(N)];
                if bestsc = fail or sc < bestsc then bestsc := sc; best := N; bestx := PreImagesRepresentative(nat, a); fi;
            od;
        od;
        Add(seq, bestx); H := best;
    od;
    return seq;
end);

# the relative orders |G_i : G_{i+1}| of a sequence of elements of G
BindGlobal("PerseusRelativeOrders", function(G, seq)
    local k, sizes, i;
    k := Length(seq);
    sizes := List([1..k + 1], i -> Size(Subgroup(G, seq{[i..k]})));
    if sizes[1] <> Size(G) then Error("PerseusRelativeOrders: the sequence does not generate G"); fi;
    return List([1..k], i -> sizes[i] / sizes[i + 1]);
end);

# the alternative the exporters also write: GAP's own Pcgs(G), when the short pcgs is in use and differs from it.
# Perseus tries both (and their prime refinements and coarsenings) on a short trial resolution and keeps the cheaper:
# neither is better for every group.
BindGlobal("PerseusAlternativePcgs", function(G, used)
    local pc;
    if Size(G) = 1 then return fail; fi;
    pc := AsList(Pcgs(G));
    if pc = used then return fail; fi;
    return pc;
end);

# the pcgs the exporters use: short by default, GAP's Pcgs when PerseusUseShortPcgs is false
BindGlobal("PerseusPcgs", function(G)
    if Size(G) = 1 then return []; fi;
    if not IsBoundGlobal("PerseusUseShortPcgs") or ValueGlobal("PerseusUseShortPcgs") <> false then return PerseusShortPcgs(G); fi;
    return AsList(Pcgs(G));
end);

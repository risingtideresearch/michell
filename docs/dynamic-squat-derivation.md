# Thin-ship dynamic sinkage/trim: independent derivation and adjudication

Raw output of a 7-agent research workflow (“Independently derive thin-ship sinkage force and trim moment, find validation anchors, map code integration, then adjudicate”), run to independently re-derive the thin-ship dynamic vertical force and pitch moment used by `crates/michell/src/squat.rs`, cross-check them against published literature/experimental anchors, map the derivation onto the crate's existing quadrature machinery, and adjudicate any disagreements between routes. Content is preserved as returned (only reformatted from JSON to Markdown), so the reasoning and validation anchors are exact, not paraphrased.

Workflow: 7 agents, 206 tool calls, 1030328 total tokens across the run.

## Adjudicated verdict

**Confidence:** High (~0.95) that the force formula F_up = F_DB + F_wPV + F_res as stated is the correct second-order thin-ship vertical force in the crate's conventions, including sign, prefactors and the radiation prescription: three independent derivations reduce to it term by term, it reproduces the crate's R_w exactly from the same kernel, I verified the Rankine cap constant independently (2-D deep-strut limit), and two independent numerical implementations (the scratchpad oracle and my judge_check.py) agree with each other and with Havelock/Kajitani (sinkage proportional to Fn^2, C = 0.0206 vs 0.019-0.025). High (~0.9) for the x-arm moment M_v and its bow-up sign (the residue trim reproduces the measured Wigley bow-up trim and crossover; the PV/Rankine moment terms follow from the same verified block but have no numerical anchor yet beyond parity). Medium for two modelling choices rather than derivations: whether to include the z-arm term M_h (recommend opt-in, default off) and how to weight the transom appendage (recommend staged: composite first, hull-only second). The orchestrator's derivation is refuted (missing waterline cap term); the previous "least sure" point (rise at low Fn) is resolved in favour of sinkage with no conflict against experiment.

### Force consensus

Conventions: bow at high x, stream -U x-hat, z down, nu = g/U^2, hull closed at both x-ends, S = {x_s<=x<=x_b, 0<=z<=T}. Transforms in the e^{-i kx x} convention (the crate's InnerIntegral kernel is e^{+i kx x} = the conjugate; every Re[...] term below is unaffected, every Im[A* F] term must be formed as Im[A_c F_c*] = -Im[A_c* F_c] from crate-convention values or its sign flips):
  F(kx,kappa) = iint_S f e^{-i kx x} e^{-kappa z} dx dz            (Q := i kx F is the transform of f_x; crate I+iJ = Q(-nu lam, nu lam^2) = conj Q(nu lam, nu lam^2))
  W(kx)       = int f(x,0) e^{-i kx x} dx                          (waterline half-beam, 1-D osc moments at z=0)
  B(kx)       = int f(x,T) e^{-i kx x} dx                          (keel line; 0 for a hull closed at the keel)
  Z(kx,kappa) := W - kappa F + e^{-kappa T} B                      (transform of -f_z including the flat-keel delta; = W - kappa F for a closed keel)
Polar wavenumber coordinates kx = k cos th, ky = k sin th, k = sqrt(kx^2+ky^2), th in (0, pi/2), all integrands folded x4 (even in kx and ky). F = F(k cos th, k), W = W(k cos th) unless stated.

F_up (positive upward, hydrostatic buoyancy excluded, second order in f) = F_DB + F_wPV + F_res with

  F_DB  = -(2 rho U^2/pi^2) int_0^{pi/2} dth cos^2 th int_0^inf k dk [ |W|^2 - |W - k F|^2 ]
          (Rankine + Neumann image = rigid-lid double body; speed-independent apart from rho U^2, cache per attitude; strictly < 0 for sections that narrow with depth; for an open keel replace |W|^2 - |W-kF|^2 by 2k Re(Z* F) - k^2|F|^2 ... i.e. use the general form -(2 rho U^2/pi^2) int dth cos^2 th int k^2 dk { Re[Z* F] + Re[W* F] } )

  F_wPV = (4 rho U^2/pi^2) int_0^{pi/2} dth cos^4 th  PV int_0^inf k^2 dk  Re[Z* F](k cos th, k) / (k cos^2 th - nu)
          (free-surface near-field; pole at k0(th) = nu sec^2 th; vanishes like 1/nu as Fn -> 0)

  F_res = (4 rho U^2 nu^3/pi) int_0^{pi/2} dth sec^4 th  Im[Z* F](nu sec th, nu sec^2 th)
        = (4 rho U^2 nu^3/pi) int_1^inf lam^3/sqrt(lam^2-1) Im[W* F](nu lam, nu lam^2) dlam   (+ e^{-kT} Im[B* F] if the keel is open)
          (wave/radiating part; identically zero for fore-aft symmetric or separable hulls)

Equivalent single expression (fourier route = scratchpad oracle wigley_oracle.py, verified identical to the havelock and pressure routes term by term):
  F_up = -(rho U^2/(2 pi^2)) iint d^2k [ A |Q|^2 - ((A-1)/k) kx^2 Re(W* F) ],  A = (kx^2 + nu k)/(kx^2 - nu k + i0 kx),  Q = i kx F,
whose polar PV form is (2 rho U^2/pi^2) int dth cos^2 th PV int k^2 dk [2 nu Re(W*F) - k(nu + k cos^2 th)|F|^2]/(k cos^2 th - nu) plus F_res.

General block used for moments/transoms: for any centreplane weight a(x,z) with x-transform a~(kx,z) and full transform A(kx,kappa), L[a] := 2 rho U iint a phi_x dx dz =
  L_R[a]   = -(rho U^2/2pi^2) iint d^2k (kx^2/k) Re int_0^T int_0^T e^{-k|z-zeta|} a~*(kx,z) f~_s(kx,zeta) dz dzeta   (Rankine, non-separable; = cap term W*F_s when a = -f_z and f_s = f)
  L_I[a]   = -(rho U^2/2pi^2) iint d^2k (kx^2/k) Re[A*(kx,k) F_s(kx,k)]
  L_PV[a]  = (rho U^2/pi^2) PV iint d^2k kx^4 Re[A* F_s] /(k (kx^2 - nu k))
  L_res[a] = (4 rho U^2 nu^3/pi) int_0^{pi/2} sec^4 th Im[A* F_s](nu sec th, nu sec^2 th) dth
F_up = L[-f_z] + keel line; L[f_x] reproduces the crate's R_w exactly and positive (this pins the -i pi sgn(kx) radiation prescription). Transom hulls: F_s = composite hull + virtual appendage in the source slot; phase 1 uses the composite for the field weights too (then L_R reduces to the cap term exactly); phase 2 restricts W, Z, X, Xw, B to the real hull, which makes the Rankine volume cross term nonzero and requires the ordered z-pair moments (exp_pair_moments_ordered, already in moments.rs).

### Moment consensus

M (positive bow-up, about the station x_ref and depth z_ref; right-hand rule with x forward and z down gives y starboard, and positive rotation about +y moves the bottom forward and the bow up, so M = int (x - x_ref) dF_up + int (z - z_ref) dF_x with dF_up = -2 p1 f_z dx dz, dF_x = +2 p1 f_x dx dz). Additional transforms:
  X(kx,kappa)  = iint (x - x_ref) f e^{-i kx x} e^{-kappa z}     (= (i d/dkx - x_ref) F; build as (x - x_center) weight + (x_center + place.x - x_ref) F)
  Xw(kx)       = int (x - x_ref) f(x,0) e^{-i kx x} dx
  Bx(kx)       = int (x - x_ref) f(x,T) e^{-i kx x} dx          (keel line, 0 if closed)
  Z_x(kx,kappa) := Xw - kappa X + e^{-kappa T} Bx                (transform of -(x - x_ref) f_z incl. keel delta)
  Zf(kx,kappa) = iint (z - z_ref) f e^{..} = (-d/dkappa - z_ref) F
  S_x(kx,k)    = int_0^T int_0^T sgn(z - zeta) e^{-k|z-zeta|} conj(xf~(kx,z)) f~(kx,zeta) dz dzeta,  xf~ = x-transform of (x - x_ref) f
                 (ordered pair kernel: same z-span T_{bb'} - T_{b'b} from exp_pair_moments_ordered; distinct spans separable N_b N+_b' with the ordering sign)

M = M_v + [M_h, optional, default OFF]

M_v = L[-(x - x_ref) f_z] = M_v,DB + M_v,wPV + M_v,res:
  M_v,DB  = -(2 rho U^2/pi^2) int_0^{pi/2} dth cos^2 th int_0^inf k^2 dk { Re[Z_x* F] + Re[Xw* F] - k Re S_x }(k cos th, k)
            (image + Rankine cap + Rankine volume; all three vanish identically for a hull fore-aft symmetric about x_ref)
  M_v,wPV = (4 rho U^2/pi^2) int_0^{pi/2} dth cos^4 th PV int_0^inf k^2 dk Re[Z_x* F](k cos th, k)/(k cos^2 th - nu)   (zero for symmetric hulls)
  M_v,res = (4 rho U^2 nu^3/pi) int_0^{pi/2} dth sec^4 th Im[Z_x* F](nu sec th, nu sec^2 th)
          = (4 rho U^2 nu^3/pi) int_1^inf lam^3/sqrt(lam^2-1) Im[F Xw* - k F X*](nu lam, nu lam^2) dlam   (closed keel)
            = the ENTIRE trim moment of a fore-aft symmetric hull (numerically: Wigley +0.73 deg bow-up at Fn 0.40, zero crossing Fn 0.339).
  M_v(x_ref) = M_v(0) - x_ref F_up exactly.

M_h = L[(z - z_ref) f_x] (couple of the longitudinal pressure force about depth z_ref), A = i kx Zf:
  residue part = -(4 rho U^2 nu^4/pi) int_0^{pi/2} sec^5 th Re[F Zf*](nu sec th, nu sec^2 th) dth = -iint (z - z_ref) dR_w (the z_ref part is exactly +z_ref R_w);
  DB/PV parts: L_I, L_PV with A = i kx Zf plus the Rankine pair term -(rho U^2/2pi^2) iint (kx^2/k) Re int int e^{-k|z-zeta|} conj((z-z_ref) f~_x) f~_s.
  It is second order like M_v, but classical thin-ship trim (Yeung 1972, Baar 1986, Noblesse) omits it because without a thrust/tow-force model the couple depends on the unknown thrust line; keep it as an opt-in with z_ref = tow-point/thrust-line depth, default off, and document.

Keel line for open-keel lofts: the Z, Z_x definitions above already include 2 int p1(x,T) f(x,T) (x - x_ref) dx via the e^{-kappa T} B, Bx terms; equivalently use the volume form 2 rho U [iint (x-x_ref) f phi_xz + int (x-x_ref) f(x,0) phi_x(x,0,0) dx] which needs only F, X, W, Xw and is automatically keel-robust.

### Rankine term verdict

FORCE: YES, the Rankine -1/(4 pi r) term contributes to F_z. Decisive argument: the volume piece iint f phi^R_xz vanishes by the (x,z)<->(xi,zeta) antisymmetry of d_x^2 d_z (1/r) against the symmetric weight f(P)f(Q) (Fourier: kernel kx^2 kz/|K|^2 odd in kz), but the wetted hull is open at the waterplane, so the divergence theorem gives F(sides+bottom) = int_V p_z dV + int_WP p1 dA and the cap term survives: F^R = 2 rho U int f(x,0) phi^R_x(x,0,0) dx = -(rho U^2/2pi^2) iint (kx^2/k) Re[W* F] != 0, negative for ordinary hulls (Bernoulli suction on the missing deck; d'Alembert holds for hull+deck, so the hull alone carries minus the deck force; equals 37% of the double-body force for the pressure route's test hull; verified by me in the 2-D deep-strut limit, constant -(rho U^2/2pi) int |kx||W|^2 dkx, by two independent arguments). The orchestrator's "Rankine gives zero force" came from dropping this cap term (its check p = rho g z passes only because rho g z = 0 on the cap). With hull-only field weights on a transom hull (f_h != f_s) the volume cross term also becomes nonzero and needs the pair kernel.
MOMENT: YES in general. Cap term Re[Xw* F] plus the non-separable volume term -k Re S_x (sgn(z-zeta) e^{-k|z-zeta|} kernel) — the Munk-type moment of the lidded half body, nonzero for hulls asymmetric fore-aft about x_ref (raked keel, deep forefoot, LCB off x_ref). Identically zero for fore-aft symmetric hulls with x_ref at midship (every integrand odd in kx). The orchestrator was right about the moment.

### Wave term verdict

FORCE: the residue (radiating) part is F_res = (4 rho U^2 nu^3/pi) int_0^{pi/2} sec^4 th Im[W* F](nu sec th, nu sec^2 th) dth: only the waterline x volume cross term survives because the |F|^2 self term is real. Nonzero in general, but identically ZERO for fore-aft symmetric hulls (F, W real with midship origin) and for separable hulls f = b(x) g(z) (same phase), so it is a small "depth-warp" correction (e.g. raked bow with vertical stern). For the Wigley it is exactly zero; sinkage there is entirely double-body + PV near-field. The orchestrator's "residue gives zero force" is right only in these special cases.
MOMENT: YES and dominant. M_v,res = (4 rho U^2 nu^3/pi) int sec^4 th Im[F Xw* - k F X*] dth is nonzero for symmetric hulls (X, Xw purely imaginary) and is the ENTIRE trim of a fore-aft symmetric hull (PV/image/Rankine moments vanish there by parity); flow reversal flips it, as trim must. Confirmed numerically: Wigley trim +0.73 deg bow-up at Fn 0.40 (oracle) / +0.71 deg (my independent code) vs experiment +0.010 rad ~ 0.57 deg and NM linear theory +0.0095 rad; principal zero crossing Fn 0.339 vs ~0.36 measured / ~0.35 NM. M_h's residue part is -iint (z - z_ref) dR_w. Both PV and residue therefore enter F_z and M; only R_w is residue-only.

### Low froude sign verdict

LINEAR THEORY: as nu -> infinity at fixed hull the wave kernel -> 0 like 1/nu and the pole k0 = nu sec^2 th runs to wavenumbers where the transforms vanish, so F_up -> F_DB = -(rho U^2/2pi^2) iint d^2k (kx^2/k^2) [|W|^2 - |W - kF|^2] (Rankine + Neumann image = rigid-lid double body). For sections that narrow with depth (f = b(x) g(z), 0 <= g <= 1, g(0) = 1) the bracket is |b^|^2 k g^ (2 - k g^) > 0, so F_DB < 0: the hull SINKS, force proportional to rho U^2, sinkage/L = C Fn^2 with C speed-independent. Wigley: C = 0.020564 (oracle, converged 1e-8; my check 0.02056). Trim -> 0 for fore-aft symmetric hulls; for asymmetric hulls the double-body moment gives an Fn^2 trim of hull-dependent sign.
EXPERIMENT/LITERATURE: Havelock 1939 (double-body, downward, proportional to U^2; half-ellipsoid L/D=16, B/D=1.6 interpolates to C ~ 0.019), Tuck & Taylor 1970 / Gourlay & Tuck 2001 (F_inf < 0, M_inf = 0 for symmetric hulls), Kajitani 1983 Wigley sigma = 2 s/Fn^2 ~ 0.04-0.05 flat for 0.1 < Fn < 0.3 (C ~ 0.020-0.025), IWWWFB28 linear Neumann-Michell within ~10% of experiment, Ma/Noblesse empirical C ~ 0.022. THEY DO NOT CONFLICT: linear theory and every experiment say sinkage proportional to Fn^2 at low Fn in deep water. The only "rise" claim (orchestrator) is an artefact of the dropped waterline cap term; with it the sign is unambiguous. Quantitatively the thin-ship Wigley sinkage runs ~10-25% above the measured values (s/L 0.00147 vs 0.0011-0.0016 at Fn 0.25; 0.00225 vs 0.0020-0.0023 at Fn 0.30-0.316; 0.00515 vs 0.0047 at Fn 0.40), consistent with thin-ship overprediction at B/L = 0.1.
WHAT TO PRINT: sinkage positive DOWN (draft increase), trim positive BOW-UP, both as the increment over the hydrostatic attitude; a one-line note that the low-speed sinkage is the speed-independent double-body suction (proportional to Fn^2, present even when R_w is negligible) and that trim at low Fn is zero for fore-aft symmetric hulls; a validity note "linear thin-ship on the at-rest hull; sinkage reliable to Fn ~ 0.45, trim to Fn ~ 0.40; trim values below ~0.05 deg at Fn 0.25-0.35 are within linear-theory noise (the residue moment oscillates in sign there, 35 crossings for the Wigley)"; and a diagnostic when |F_z| exceeds ~0.5 x weight (planing regime, outside thin-ship theory). No sign-conflict disclaimer is needed.

### Errors found

1. ORCHESTRATOR (decisive): F_up = int_V dp/dz dV is incomplete for the wetted hull, which is open at the waterplane; the divergence theorem gives F(sides+bottom) = int_V p_z dV + int_WP p1 dA with int_WP p1 dA = 2 rho U int f(x,0) phi_x(x,0,0) dx (the wave-hollow loss of buoyancy). Fix: F_up = -(rho U^2/2pi^2) iint [A|Q|^2 - ((A-1)/k) kx^2 Re(W*F)] (the orchestrator's term is exactly the first piece). Consequences corrected: Rankine contributes (cap term), the residue contributes for asymmetric hulls, and the low-Fn limit is sinkage, not rise.
2. CODE MAP numerics: the claim that the crate's lambda quiet-window truncation applies unchanged to the force integral is wrong for the non-residue (double-body) part, whose theta-integrand is finite/log-singular at theta -> pi/2; grade panels in u = pi/2 - theta instead (oracle recipe). Also "cache the Rankine part" should be "cache the whole double-body part F_DB (Rankine + Neumann image)", which is the entire speed-independent piece. The code map's residue-only reuse of eval(lambda) is fine.
3. IMPLEMENTATION SIGN TRAP (not an error in any route, but load-bearing): the crate's InnerIntegral kernel is e^{+i kx x} while all three routes define transforms with e^{-i kx x}; the Re[...] terms are invariant, but the residue terms Im[W*F], Im[F Xw* - k F X*] flip sign under conjugation. Form them as Im[A_c F_c*] = -Im[A_c* F_c] from crate-convention values, then verify with test 2 (L[f_x] == R_w > 0) and the Wigley bow-up trim at Fn 0.40.
4. HAVELOCK route: (1.1) assumes f(x,T) = 0; for lofted hulls with an open keel the keel-line term 2 int p1(x,T) f(x,T) dx (transform e^{-kappa T} B) must be added, or the keel-robust volume form used. Not an error under its stated assumption.
5. Minor: the pressure route's F^R constant was verified only through R_w and a two-representation agreement; I added an independent analytic constant check (deep wall-sided strut, 2-D Hilbert-transform limit) which confirms -(rho U^2/2pi) int |kx||W|^2, so no factor-2 slip exists. No pi, factor-2, image-sign, or radiation-sign errors were found in fourier, havelock or pressure; all three reduce to the same F and M term by term including prefactors and polar forms (the pressure route's s = k/k0 substitution and the havelock route's F1+F2+F3 grouping were checked explicitly).
6. Literature report 1's transcription F_up = -2 rho U iint phi_x f_z (surface form) is correct and equals the routes' formulas; its warning about the factor 2 and f_z sign is resolved (both sides counted, n_z dS = -f_z dx dz per side).

### Disagreements between routes (resolved)

**Does the Rankine -1/r term contribute to F_z?**

- Positions: fourier/havelock/pressure: yes, through the waterline (cap) term only, volume term zero by antisymmetry. Orchestrator: zero, by the same antisymmetry applied to the whole force.
- Resolution: Routes are right. The divergence theorem over the open wetted hull leaves int_WP p1 dA; the orchestrator's F_up = int_V p_z dV omits it (its hydrostatic check passes only because rho g z vanishes on z=0). Independently confirmed: 2-D deep-strut limit gives -(rho U^2/2pi) int |kx||W|^2 from both the Fourier formula and a direct Hilbert-transform calculation.
- Test to settle: Compute F^R two ways on one hull: -(rho U^2/2pi^2) iint (kx^2/k) Re[W*F] versus the (kx,kz) plane form -(rho U^2/2pi^2) iint (kx^2/|K|) Re[fhat W*] (pressure route: -0.2856 vs -0.2858 for f = b(1-x^2)(1-z/T), T=0.125, rho=U=1); and check that the pair-kernel volume term returns 0 when f_h = f_s.

**Sign of F_z as Fn -> 0**

- Positions: Routes: sinks (F_DB < 0). Orchestrator: rises (+(rho U^2/2pi^2) iint k|Q|^2 with A = -1).
- Resolution: Sinks. The orchestrator's expression is exactly the volume term alone; the dropped cap term -(rho U^2/pi^2) iint (kx^2/k) Re(W*F) is larger and negative. Wigley C = s/L/Fn^2 = 0.02056 agrees with Havelock's ellipsoid (0.019) and Kajitani's sigma/2 (0.020-0.025).
- Test to settle: Evaluate F_up at Fn = 0.02, 0.05, 0.10 for the crate's Wigley and compare s/L/Fn^2 with 0.020564 (oracle) and with the rigid-lid formula A = -1.

**Does the residue term contribute to F_z?**

- Positions: Orchestrator: never (odd in kx against |Q|^2). Routes: yes via Im[W*F], zero only for symmetric or separable hulls.
- Resolution: Routes are right; the orchestrator lacked the cross term. Practically small.
- Test to settle: Hull with raked bow and vertical stern and non-affine sections: F_res != 0; reflect x -> -x: F_res flips sign while F_DB + F_PV unchanged; separable hull: F_res = 0 to roundoff.

**Include the z-arm (resistance depth) term M_h in the pitching moment?**

- Positions: fourier and pressure: include with free z_ref (second order like the rest). havelock: same. Literature (Yeung, Baar, Noblesse): omitted in classical thin-ship trim. Code map: pivot at lcg, no thrust model.
- Resolution: Implement M_v as the default; expose M_h as an opt-in with z_ref = tow-point/thrust-line depth. Without a thrust or tow-force model the drag couple is not closed and the term would bias trim by an arbitrary reference choice; the Wigley towed-model data implicitly contain the tow-line couple.
- Test to settle: Compare Wigley trim vs Kajitani/IWWWFB with M_v alone (expected: crossing ~0.34, +0.73 deg at Fn 0.40) and with M_h at z_ref = 0 and z_ref = T/2; the M_v-alone result should be the closest.

**Transom appendage on the pressure-weighting side**

- Positions: All three routes assumed one common closed S (composite in both slots). havelock unsure (4) and the code map argue the appendage is a virtual source only and the field weights should cover the real hull.
- Resolution: Physically the pressure acts on the real hull only, so hull-only field weights are more correct, but this reintroduces the non-separable Rankine volume cross term (f_h != f_s). Stage it: phase 1 composite in both slots (all Rankine force terms reduce to the cap, exactly as derived), phase 2 hull-only weights using exp_pair_moments_ordered; the difference is transom-local.
- Test to settle: For a transom hull evaluate F_up, M with composite vs hull-only weights at several L_v; the difference must be continuous in L_v, vanish for a hull that closes aft, and be localized (insensitive to bow shape changes).

**theta truncation of the outer integral**

- Positions: Code map and fourier numerics: reuse the crate's Michell quiet-window truncation past lambda ~ 2 for the force integrals.
- Resolution: Valid only for the residue and the pole-subtraction value. The double-body (Rankine + image) theta-integrand does not vanish at theta -> pi/2 (kx -> 0 at fixed k keeps |W(k cos th)|^2 O(1) out to k ~ 1/(L cos th)); the oracle found log(1/cos th) growth and needed geometric panels in u = pi/2 - theta down to 1e-10. Split F = F_DB (no truncation, cached per attitude) + F_wPV + F_res (truncatable).
- Test to settle: Compute F_DB with theta cutoffs at lambda_max = 2, 10, 100 and with the u-graded full range; the cutoff versions must converge to the graded value only slowly (log), confirming truncation is unsafe.

### Validation plan

- 1. Unit: eval_at(nu*lam, nu*lam^2) == InnerIntegral::eval(lam) to 1e-12; transforms F, W, X, Xw (and B for open keels) vs brute-force 2-D quadrature on a spline hull; W(0) = A_w/2 and F(0,0) = V/2 to 1e-10.
- 2. Constant and radiation-sign pin: run the generic block with weight a = f_x; its residue part must equal the crate's R_w (positive) to 1e-10 and its DB/PV parts must vanish (d'Alembert, kx-odd).
- 3. Antisymmetry: the Rankine pair-kernel volume force term (built from exp_pair_moments_ordered) must return 0 to roundoff when field and source hulls coincide; for a fore-aft symmetric hull about x_ref, M_DB = M_PV = 0 and F_res = 0 to roundoff.
- 4. Analytic limits: nu -> infinity gives F_up -> F_DB with O(1/nu) difference; s/L/Fn^2 flat as Fn -> 0 (Fn = 0.02, 0.03, 0.05); wall-sided strut closed form F_DB = -(rho U^2/2pi^2) iint (kx^2/k^2)|W|^2 (1 - e^{-2kT}); deep-strut 2-D limit -(rho U^2/2pi) int |kx||W|^2 dkx.
- 5. Symmetry: reflecting the hull x -> -x leaves F_up unchanged and flips M; M(x_ref) = M(0) - x_ref F_up to roundoff.
- 6. Wigley oracle (scratchpad wigley_oracle.json, 13 Froude numbers, converged to 1e-8): F_up (e.g. -986.5 N at Fn 0.25, -3450.7 N at Fn 0.40 for L=10, B=1, T=0.625, rho=1025) and M_res (4267 N m at Fn 0.40, trim +0.7297 deg) to ~1e-3 relative (B-spline loft error); low-Fn C = s/L/Fn^2 = 0.020564; principal trim zero crossing Fn 0.339.
- 7. Literature anchors, sinkage: Havelock 1939 ellipsoid C ~ 0.019 (expect thin-ship +8%); Kajitani 1983 sigma = 2s/Fn^2 ~ 0.04-0.05 for 0.1 < Fn < 0.3, rising to ~0.06 by 0.40; Wigley s/L 0.0011-0.0012 (Fn 0.25), 0.0020-0.0023 (0.30-0.316), 0.0047 (0.405); accept theory within +10-25%.
- 8. Literature anchors, trim: Wigley ~0 below Fn 0.28, slight bow-down -0.0005 to -0.0015 at 0.30-0.35 (theory only marginally reproduces this: oscillating +/-0.03 deg), bow-up crossing 0.34-0.36 (theory 0.339, NM 0.35), +0.010 rad at 0.405 (theory +0.0127 rad); NM linear theory (IWWWFB28) +0.0095 rad at 0.40.
- 9. Asymmetric anchor (Series 60 Cb 0.6, qualitative): sinkage 0.0018 (Fn 0.25), 0.0032 (0.32); trim slightly bow-down at low Fn turning bow-up near 0.37 — checks the sign of the DB/PV moment terms that vanish for the Wigley.
- 10. Transom hulls: F, M continuous as L_v -> 0+ (Fixed{length} sweep), finite for every positive hollow, inert for a hull that closes aft (bit-identical to the closed-hull path); phase-2 hull-only weights differ from composite weights by a transom-local amount insensitive to bow changes.
- 11. Equilibrium closure: for a symmetric waterplane s = -F_up/(rho g A_w), theta = M/(rho g I_L) (Gourlay 2003; Yang 2000 eqs 8-9); general 2x2 with H0, H1, H2 (Tarafder 2006 5.3-5.4); Picard/Newton loop converges for Fn <= 0.45 and reproduces the report-only values at convergence.

### Numerics plan

Coordinates: polar (theta, k) in the (kx, ky) plane at y = 0, kx = k cos th, kappa = k, theta in (0, pi/2), x4 folding. Split every load into DB (speed-independent, no pole), wPV (pole) and res (Kelvin curve).
DOUBLE-BODY PART (F_DB, M_v,DB): no pole. Its theta-integrand does NOT vanish at theta -> pi/2 (kx -> 0 at fixed k keeps |W(k cos th)|^2 = O(1) out to k ~ 1/(L cos th); the k-integral grows like log(1/cos th) and cos^2 th only tames it): integrate the full range with geometric panels in u = pi/2 - theta on [1e-10, 0.1] plus uniform panels on [0.1, pi/2] (oracle recipe, converged 1e-8), GL 16 per panel. k-axis: [0, ~0.05/L] + geometric panels (~40/decade) to k_max = max(20/(L cos th) x few, 2e4/(a cos th)); integrand decays like k^-2 (|W|^2 ~ kx^-4 for a C^0 waterline) so the tail ~ 1/k_max; optionally subtract the large-k asymptote F ~ W/k - W_z/k^2 analytically. Compute once per attitude and cache across the speed loop (scales exactly as rho U^2). Rankine pair term (asymmetric-hull moments, hull-only transom weights): same nodes, per (kx,k) contract the x-reduced vectors with the ordered pair moments (T_{bb'} - T_{b'b} same span; N_b N+_b' with sign for distinct spans), O(n_z^2 (q+1)^2) per node.
WAVE PV PART: per theta, k0 = nu sec^2 th; PV int_0^inf g(k)/(k cos^2 th - nu) dk = sec^2 th { int_0^{2k0} [g(k) - g(k0)]/(k - k0) dk + int_{2k0}^inf g(k)/(k - k0) dk } (the log term vanishes on the symmetric window; the oracle's [k0/2, 3k0/2] window is equivalent). g is entire in k for a bounded hull so the subtracted integrand is smooth; g(k0) comes from eval(lambda = sec th) / eval_at(nu lam, nu lam^2) — the same node the residue uses; even GL rules never hit k = k0. The kernel 2 k cos^2 th/(k cos^2 th - nu) x Re[Z*F] decays like k^-2 x |W|^2 beyond k0 (geometric octave panels, ~8 nodes each, ~12 octaves) and the whole PV integrand vanishes like cos^4 th as theta -> pi/2, so the crate's marching theta panels sized by rate(theta) with the quiet-window truncation ARE applicable here (regular part ~24-32 GL nodes on [0, 2k0]).
RESIDUE PART: 1-D theta (or lambda) integral at Michell's nodes (kx, kappa) = (nu sec th, nu sec^2 th), weight sec^4 th (force, moment) — dtheta sec^4 th = lam^3 dlam/sqrt(lam^2 - 1); integrand ~ cos^6 th |W|^2/nu, faster than R_w's own; reuse integrate_outer verbatim with amp_sq replaced by Im[Z*F], Im[Z_x*F].
TRANSFORMS per node: one fill_zm(kappa = k), then nets f (F), (x - x_center) f (X), 1-D waterline nets W, Xw (+ keel line B, Bx if the loft leaves f(x,T) != 0); transom appendage only in the source slot (F_s) via transom_term at independent (kx, kappa). Sign traps: crate kernel e^{+i kx x} conjugates every transform (Im terms flip; Re terms invariant); x-arm = (x - x_center) weight + (x_center + place.x - x_ref) F; multihull loads use cross products F_s,i A_j* with placement phases, not |sum|^2.
DECAY/TRUNCATION: TransomClosure::None leaves a step (W ~ 1/kx) making the DB integrand ~ 1/k and the Rankine pressure log-singular at the step — require a positive hollow for loads. Large kappa underflow: skip pairs whose e^{-kappa (z0_t - z1_t')} underflows (mirror fill_zm's decay == 0 branch).
EXPECTED COST (release, exact transforms ~1.7 us per fill+accumulate at full resolution): F_DB ~ 3e2 theta x 6e2 k = 2e5 nodes ~ 0.3-0.5 s once per attitude; per speed F_wPV + F_res ~ (Michell theta nodes ~ 6e3 at frac 1, or ~3e2 with coarser panels) x ~150 k nodes ~ 5e4-1e6 nodes ~ 0.1-2 s; moments add extra nets on the same nodes (~x1.5-2.5). One refinement pass (frac 0.5) for the reported row only; rel_tol 1e-3 for the equilibrium loop. Coarse fleet (7x4 spans) is ~3x cheaper per node. Well under the code map's 3-5 s per full evaluation.

### Suggested implementation order

- 1. InnerIntegral::eval_at(kx, kappa) + eval_nets_at with PolyNet (f, x-weighted f, 1-D waterline/keel line nets); unit tests: eval_at(nu lam, nu lam^2) == eval(lam), transforms vs quadrature, W(0) = A_w/2, F(0,0) = V/2.
- 2. Generic block L[a] with the three pieces (DB on a u-graded theta grid, wPV with pole subtraction, residue on Michell nodes); first weight a = f_x: must reproduce R_w exactly and positive, with DB/PV parts zero. This pins constants and the residue sign (and exposes the e^{+ikx x} conjugation trap).
- 3. F_up for closed hulls (composite in both slots): F_DB (cached per attitude) + F_wPV + F_res; tests: nu -> inf limit, wall-sided strut closed form, Wigley oracle table (C = 0.020564 at low Fn; -3450.7 N at Fn 0.40), symmetric hull F_res = 0.
- 4. M_v residue (x-arm) on the same Michell nodes; test: Wigley trim +0.7297 deg at Fn 0.40, crossing 0.339, M(x_ref) = M(0) - x_ref F_up, reflection flips M.
- 5. M_v DB + wPV separable terms (image, Rankine cap, PV) for asymmetric hulls; test: exactly zero for symmetric hulls, Series 60 qualitative trim sign.
- 6. Rankine pair-kernel volume moment term via exp_pair_moments_ordered; test: force analogue returns 0 for f_h = f_s, moment zero for symmetric hulls, convergence vs a brute-force (kx,kz) plane-transform evaluation on one hull.
- 7. Report-only CLI/manifest mode (dynamic = report): print F_z, M, sinkage_dyn (positive down), trim_dyn (positive bow-up) at the hydrostatic attitude with the validity/low-Fn notes.
- 8. Equilibrium hook (dynamic = solve): Picard with lazy re-evaluation, continuation in speed, freeze transom detection during the dynamic solve, cap |lift| at 0.5 x weight with a diagnostic; tests: Wigley free-to-sink-and-trim vs Kajitani/IWWWFB within stated tolerances.
- 9. Transom phase 2: hull-only field weights with the pair-kernel Rankine cross term; tests: continuity in L_v, inert for closed hulls, transom-local difference.
- 10. Optional: M_h (z-arm) opt-in with z_ref; multihull cross terms (superpose_force); heel approximation documented.

## Independent derivations (three routes)

Three agents independently derived the same force/moment starting from different formulations (Fourier Green's function, Havelock source distribution, direct pressure integration), then were cross-checked against each other term by term.

### Route: fourier

**Approach**

Start from steady Bernoulli in the ship frame (Phi = -Ux + phi, z down): p = rho g z + rho U phi_x. Thin-ship body condition gives a centreplane source sheet sigma = -2U f_x. Build the Green function by 2-D Fourier transform in (x,y): FT of -1/(4 pi r) is -e^{-k|z-zeta|}/(2k); add A(kx,k) e^{-k(z+zeta)} and fix A from phi_z = phi_xx/nu on z=0 (the z-down form of U^2 phi_xx + g phi_Z = 0). Then write F_z and M as hull integrals of phi_x, phi_xz (volume/cap form), insert the Fourier Green function, and reduce every term to products of exact hull transforms F, W, X, Xw, Zf (and Phi, X3, Z3 in (kx,kz) for the Rankine moment) integrated over wavenumber space, with the Kelvin pole split into PV + residue. Full working (5 stages) at /private/tmp/claude-501/-Users-avi-michell/2fdb50af-0291-45b8-bdde-e7c422fd969b/scratchpad/derive-A.md.

**Pressure and force form**

p = rho g z + p1 + O(f^2), p1 = rho U phi_x (from |grad Phi|^2 = U^2 - 2U phi_x + ...; z down, datum p=0 at Z=0). Body condition phi_y(x,0+,z) = -U f_x, so sigma = -2U f_x (positive = source at the bow, where f_x < 0). Linear free-surface elevation at the hull: eta = U phi_x(x,0,0)/g = p1/(rho g) (positive up).

Vertical force (positive up), excluding rho g V0:
(a) Volume form. Close the wetted hull with a cap on the actual free surface (p=0 there): F_up = -∫_V dp/dZ dV = +∫_V dp/dz dV over the hull volume below the actual free surface. Splitting off the sliver -eta<z<0 of width 2f(x,0): F_z = rho g ∫ 2 f(x,0) eta dx + ∬_S 2 f dp1/dz dx dz = ∫ 2 f(x,0) p1(x,0,0) dx + ∬_S 2 f dp1/dz dx dz.
(b) Surface form. Sides (n_Z = +f_z per dx dz, both sides) plus flat bottom (n_Z = -1): F_z = -2 ∬_S p1 f_z dx dz + 2 ∫ f(x,T) p1(x,T) dx; the sliver contributes O(f^3).
Equality: integrate (a) by parts in z: ∬ 2f dp1/dz = 2∫f(x,T)p1(x,T) - 2∫f(x,0)p1(x,0) - ∬ 2 f_z p1; the -2∫f(x,0)p1(x,0) term cancels the wave-elevation term exactly. So the "waterline/wave-elevation term" is precisely the pressure on the cap z=0 of the closed lower body; both forms agree station by station. I use (a) because it needs only transforms of f and of the waterline f(x,0):
  F_z = 2 rho U [ ∬_S f phi_xz(x,0,z) dx dz + ∫ f(x,0) phi_x(x,0,0) dx ]                         (F1)
Moment, positive bow-up (bow at high x), about the axis through (x_ref, z_ref): a +x force below the axis lifts the bow, so
  M = ∫ (x-x_ref) dF_up + ∫ (z-z_ref) dF_x,  dF_x = +2 p1 f_x dx dz (hydrostatic part drops since ∫ f_x dx = 0):
  M = 2 rho U [ ∬_S (x-x_ref) f phi_xz + ∫ (x-x_ref) f(x,0) phi_x(x,0,0) dx + ∬_S (z-z_ref) f_x phi_x ]   (M1)
The third (z-arm) term is O(f^2) like the others; Tuck's classical trim formula omits it. I keep it with z_ref free (z_ref = KG for the equilibrium problem; set the term to zero to recover the classical form). Same pressure gives R = -F_x = -2 rho U ∬ f_x phi_x.

**Green function**

Convention: phi = (2pi)^-2 ∬ phi_hat e^{i(kx x + ky y)} dkx dky, k = sqrt(kx^2+ky^2), source at (xi,0,zeta), zeta>0, grad^2 G = delta.
Fourier: G_hat = -e^{-k|z-zeta|}/(2k) + A e^{-k(z+zeta)}.  Free-surface condition: kinematic phi_Z = -U eta_x, dynamic eta = U phi_x/g, so U^2 phi_xx + g phi_Z = 0 on Z=0, i.e. with z down phi_z = phi_xx/nu on z=0, in Fourier dG_hat/dz = -(kx^2/nu) G_hat at z=0. Solving: -1/2 - kA = -(kx^2/nu)(-1/(2k) + A)  =>
  A(kx,k) = (nu k + kx^2) / (2k (kx^2 - nu k)) = -1/(2k) + kx^2/(k (kx^2 - nu k)).
Checks: nu->inf gives A = -1/(2k) (positive image, rigid lid); nu->0 gives +1/(2k) (phi = 0 on z=0).
  G(x,y,z; xi,zeta) = -1/(4 pi r) - 1/(4 pi r') + (2pi)^-2 ∬ dkx dky kx^2 e^{-k(z+zeta)} e^{i(kx(x-xi) + ky y)} / ( k (kx^2 - nu k + i0 kx) ),
  r^2 = (x-xi)^2 + y^2 + (z-zeta)^2, r'^2 = (x-xi)^2 + y^2 + (z+zeta)^2.  Call the three pieces G_R, G_I, G_W.
Pole: kx^2 = nu k; with kx = k cos th, k0(th) = nu sec^2 th (Kelvin dispersion).
Radiation condition: Rayleigh damping following the fluid (material derivative -U d/dx here): (-U d/dx + mu)^2 phi + g phi_Z = 0 => kx^2 -> kx^2 + i eps kx, eps>0, both roots kx ≈ ±sqrt(nu k) - i eps/2 in the LOWER half kx-plane. For x - xi > 0 (upstream) close upward: no waves; for x - xi < 0 the poles contribute: waves trail toward -x. In polar form kx^2/(k(kx^2 - nu k + i0 kx)) = cos^2 th /(D + i0 sgn cos th), D = k cos^2 th - nu, and 1/(D + i0 sgn cos th) = PV(1/D) - i pi sgn(cos th) delta(D). This sign is verified independently: it reproduces the crate's Michell R_w exactly and positive (the opposite choice gives R_w < 0).

**Force formula**

Cartesian form (all integrals over the full (kx,ky) plane; A carries the +i0 kx prescription):
  F_z = (rho U^2/pi^2) ∬ dkx dky kx^2 { A(kx,k) [ F(kx,k) W*(kx) - k |F(kx,k)|^2 ] - (1/(2k)) F(kx,k) W*(kx) }
(first brace = image + wave acting on volume + cap terms; last term = Rankine, cap only).
Polar, PV + residue, explicit (kx = k cos th, everything evaluated at F = F(k cos th, k), W = W(k cos th)):
  F_z = (2 rho U^2/pi^2) ∫_0^{pi/2} dth cos^2 th  PV∫_0^inf dk  k^2 [ 2 nu Re(F W*) - k (nu + k cos^2 th) |F|^2 ] / (k cos^2 th - nu)
      + (4 rho U^2 nu^3/pi) ∫_0^{pi/2} dth sec^4 th  Im[ F(nu sec th, nu sec^2 th) W*(nu sec th) ].
In the crate's lambda = sec th: residue = (4 rho U^2 nu^3/pi) ∫_1^inf lambda^3/sqrt(lambda^2-1) Im[F(nu lambda, nu lambda^2) W*(nu lambda)] dlambda, where the crate's I + iJ = -i nu lambda F*(nu lambda, nu lambda^2).
Limits: nu->inf: F_z -> F_z^N = -(rho U^2/(2pi^2)) ∬ dkx dky (kx^2/k^2) [ |W|^2 - |W - kF|^2 ] (double-body, < 0 for ordinary hulls). nu->0: F_z -> -(rho U^2/(2pi^2)) ∬ kx^2 |F|^2.
Cross-check from the same kernel: R = -2 rho U ∬ f_x phi_x = (2 rho U^2 nu^4/pi) ∫_{-pi/2}^{pi/2} sec^5 th |F(nu sec th, nu sec^2 th)|^2 dth = (4 rho g^2/(pi U^2)) ∫_1^inf |I+iJ|^2 lambda^2/sqrt(lambda^2-1) dlambda — Michell exactly, PV/Rankine/image giving zero.

**Moment formula**

Positive bow-up about (x_ref, z_ref). Define Gm(kx,k) = F Xw* - k F X* - i kx F Zf*  (Gm(-kx) = Gm*), with all of F, X, Zf at (kx,k) and Xw at kx.
Cartesian:
  M = (rho U^2/pi^2) ∬ dkx dky kx^2 { A(kx,k) Gm(kx,k) - (1/(2k)) F Xw* }
      - (rho U^2/(2pi^2)) ∬ dkx dkz (kx^2/kap) Im[ kx Phi Z3* - kz Phi X3* ],   kap = sqrt(kx^2 + kz^2)
(first line: image + wave on all three terms of (M1), plus Rankine cap; second line: Rankine volume (x-arm) and z-arm terms, which need the (kx,kz) transforms because the |z-zeta| kernel is not separable off z=0).
Polar PV + residue for the first line:
  M_FS = (2 rho U^2/pi^2) ∫_0^{pi/2} dth cos^2 th PV∫_0^inf dk k^2 [ (nu + k cos^2 th) Re Gm - (k cos^2 th - nu) Re(F Xw*) ] / (k cos^2 th - nu)
       + (4 rho U^2 nu^3/pi) ∫_0^{pi/2} dth sec^4 th Im Gm(nu sec th, nu sec^2 th).
Drop the -i kx F Zf* piece of Gm and the kx Phi Z3* piece of the Rankine line to recover the classical "moment of the vertical force distribution only" (z-arm omitted). For a fore-aft symmetric hull with x_ref at midship the PV and Rankine parts vanish identically and M is entirely the residue term (Im Gm = -F Im Xw + k F Im X - kx F Zf).

**Hull transforms**

S = centreplane region x_s <= x <= x_b, 0 <= z <= T (composite closed body, f = 0 at both x-ends). All are exact per knot span via the crate's osc/exp moments (complex kappa for the (kx,kz) ones).
  F(kx,kappa)  = ∬_S f(x,z) e^{-i kx x} e^{-kappa z} dx dz            (F(-kx,kappa) = F*)
  Q(kx,kappa)  = ∬_S f_x e^{-i kx x} e^{-kappa z} dx dz = i kx F      (source transform; the crate's I+iJ = Q(-nu lambda, nu lambda^2))
  W(kx)        = ∫ f(x,0) e^{-i kx x} dx                             (waterline; 1-D, z = 0)
  X(kx,kappa)  = ∬_S (x - x_ref) f e^{-i kx x} e^{-kappa z} dx dz = (i d/dkx - x_ref) F
  Xw(kx)       = ∫ (x - x_ref) f(x,0) e^{-i kx x} dx
  Zf(kx,kappa) = ∬_S (z - z_ref) f e^{-i kx x} e^{-kappa z} dx dz = (-d/dkappa - z_ref) F
  Phi(kx,kz)   = ∬_S f e^{-i(kx x + kz z)} dx dz = F(kx, i kz)
  X3(kx,kz)    = ∬_S (x - x_ref) f e^{-i(kx x + kz z)} dx dz = (i d/dkx - x_ref) Phi
  Z3(kx,kz)    = ∬_S (z - z_ref) f e^{-i(kx x + kz z)} dx dz = (i d/dkz - z_ref) Phi
  (B(kx) = ∫ f(x,T) e^{-i kx x} dx and Fz = ∬ f_z e^{..} appear only in the identity W - kF = B e^{-kT} - Fz, i.e. the surface form.)
In the FS integrals kappa = k = sqrt(kx^2+ky^2) = k; on the residue locus kx = nu sec th, k = nu sec^2 th. Gf := F W* - k|F|^2; Gm := F Xw* - k F X* - i kx F Zf*.

**Rankine direct term**

- **Contributes to force:** True
- **Contributes to moment:** True
- **Argument:** Volume term ∬ f phi^R_xz: the kernel d_x d_z(-1/4 pi r) = -3(x-xi)(z-zeta)/(4 pi r^5) on y=0 is symmetric under (x,z)<->(xi,zeta); integrating f f_xi K by parts in xi, relabelling, and by parts in x gives V^R = -V^R = 0 (Fourier: -(iU/4pi^2) ∬ kx^2 kz |Phi|^2/kap, odd in kz). Same antisymmetry gives zero Rankine resistance (d'Alembert). BUT the thin hull is cut open at the waterplane, so d'Alembert does not apply to it: the force on sides+bottom = -∫_V grad p dV (zero) PLUS the cap pressure on z=0, and for observation points on z=0 the Rankine kernel is separable (|z-zeta| = zeta), giving F_z^R = -(rho U^2/(2pi^2)) ∬ dkx dky (kx^2/k) Re[F(kx,k) W*(kx)] — identical in form to the image cap term, nonzero, and negative (suction: phi_x<0 alongside) for ordinary hulls. So Rankine contributes to F_z through the waterline/cap term only. Moment: the cap term with W -> Xw, plus the volume and z-arm terms -(rho U^2/2pi^2) ∬ dkx dkz (kx^2/kap) Im[kx Phi Z3* - kz Phi X3*], whose kernels lose the symmetry because of the (x - x_ref) and (z - z_ref) weights. For a fore-aft symmetric hull (x_ref midship) every Rankine moment integrand is odd in kx and M^R = 0 (reversibility); for an asymmetric hull (raked keel, deep forefoot) M^R != 0 — the Munk moment of the doubled body's lower half. So Rankine contributes to M in general.

**Wave residue term**

- **Contributes to force:** True
- **Contributes to moment:** True
- **Argument:** The residue of 1/(k cos^2 th - nu + i0 cos th) at k0 = nu sec^2 th picks out Im of the hull bilinears (the -i pi sgn(cos th) delta term times the kx-odd imaginary part survives; the kx-even real part cancels between ±kx). For F_z the residue is (4 rho U^2 nu^3/pi) ∫_0^{pi/2} sec^4 th Im[F W*] dth: only the cross term between the immersed hull transform F and the waterline transform W survives, because the volume piece k|F|^2 is real. Nonzero for a generic hull; identically zero for a fore-aft symmetric hull (F and W real with midship origin), whose sinkage is then a pure near-field (PV) effect. For M the residue is (4 rho U^2 nu^3/pi) ∫ sec^4 th Im Gm dth, Gm = F Xw* - k F X* - i kx F Zf*; for a symmetric hull with x_ref midship Re Gm = 0 (so PV and Rankine parts vanish) while Im Gm = -F Im Xw + k F Im X - kx F Zf != 0: the trim moment of a symmetric hull is entirely the wave residue. Validation: the same residue machinery applied to R = -2 rho U ∬ f_x phi_x gives exactly the crate's Michell integral, positive, confirming both the normalisation and the pole-displacement sign.

**Low froude limit**

nu -> inf at fixed hull: the wave kernel kx^2/(k(kx^2 - nu k)) -> -kx^2/(nu k^2) -> 0 and the pole k0 = nu sec^2 th moves to wavenumbers where |F| ~ |W|/k has decayed, so wave PV and residue parts vanish like 1/nu. What remains is Rankine + positive image = Neumann double-body flow:
  F_z -> F_z^N = (rho U^2/pi^2) ∬ dkx dky kx^2 [ |F|^2/2 - Re(F W*)/k ] = -(rho U^2/(2pi^2)) ∬ dkx dky (kx^2/k^2) [ |W|^2 - |W - kF|^2 ].
With sectional transform F~(kx,z) = ∫ f e^{-i kx x} dx: W - kF = F~(0) e^{-kT} + ∫_0^T [F~(0) - F~(z)] k e^{-kz} dz, a convex combination of sectional transforms. For sections affine to the waterline, f(x,z) = f(x,0) g(z) with 0 <= g <= 1 (wall-sided or narrowing with depth), |W - kF| = |W| [e^{-kT} + ∫(1-g) k e^{-kz} dz] < |W|, hence F_z^N < 0 strictly: the hull SINKS. Wall-sided strut: F_z^N = -(rho U^2/2pi^2) ∬ (kx^2/k^2)|W|^2 (1 - e^{-2kT}). Physics: the double-body flow speeds up alongside the hull, phi_x < 0, p1 < 0, the linear free surface eta = U phi_x/g is depressed along the waterline and buoyancy is lost; this waterline (cap) term -Re(F W*)/k outweighs the upward volume term |F|^2/2 (suction weakening with depth). This is the classical Fn^2 low-speed squat, with zero trim for a symmetric hull (all M-integrands odd in kx). Caveat: a hull much fuller below the waterline than at it (bulb) can make |W|^2 - |W-kF|^2 negative at some k; the sign theorem is for g <= 1. Note also nu -> 0 gives A -> +1/(2k), the cap term cancels exactly, and F_z -> -(rho U^2/2pi^2) ∬ kx^2 |F|^2 < 0 too; any linear-theory rise lives at intermediate Froude numbers via the PV part near the pole.

**Numerics**

PV: for each th, k0 = nu sec^2 th; PV∫_0^inf g(k)/(k cos^2 th - nu) dk = sec^2 th [ ∫_0^{2k0} (g(k) - g(k0))/(k - k0) dk + ∫_{2k0}^inf g(k)/(k - k0) dk ] — the log term vanishes on the symmetric interval, g is entire in k for a bounded hull so the subtracted integrand is smooth; Gauss-Legendre panels sized to the oscillation rate (as in the crate's outer integral) suffice. Alternative: shift k -> k - i delta (transforms already accept complex kappa) and take Re. Because the ±kx halves are conjugate and the integrand is even in th, integrate th over (0, pi/2) only with the factor 4 already included. As th -> pi/2, k0 -> inf and e^{-k0 z} kills g(k0): truncate th (or lambda) exactly as the crate does for Michell. PV and residue share the same F(k cos th, k) nodes.
Decay: for large k at fixed th, F ≈ W(k cos th)/k (e^{-kz} collapses onto the waterline), so the PV integrand ~ -k |W(k cos th)|^2; with u = k cos th, ∫ dk k |W|^2 = cos^-2 th ∫ u |W(u)|^2 du, finite because |W|^2 ~ u^-4 (C^0 waterline with end kinks) or faster for smoother splines; the cos^2 th prefactor cancels cos^-2 th, so the th integral is finite. Residue integrands ~ sec^4 th |F||W| ~ cos^2 th |W(nu sec th)|^2/nu ~ cos^6 th: converge faster than Michell's own sec^5 |F|^2. Rankine (kx,kz) moment integrals: Phi ~ W/(i kz), X3 ~ Xw/(i kz) for large kz, integrand ~ kx^2 W Xw/(kz kap): convergent; the |Phi|^2 force analogue is exactly zero (odd in kz) and need not be computed. No k -> 0 trouble: kx^2/k and kx^2/k^2 are integrable in 2-D. Transom hulls: use the composite (virtually closed) f so that Q = i kx F holds and the endpoint terms vanish; a bare step would make W decay only like 1/kx and slow the k-integrals (same issue the crate documents for Michell).

**Unsure about**

(1) The z-arm (resistance x depth) term of M: it is second order like the rest and I include it with a free z_ref, but classical thin-ship trim formulas omit it; I am confident of its derivation but not of whether the crate wants it. (2) The thin-ship replacement of p1 inside the hull volume by its centreplane value in the volume form; this is the standard thin-ship approximation and consistent with Michell's order, but phi_x has the usual logarithmic behaviour at the hull ends, which the surface and volume forms share. (3) The sign theorem F_z^N < 0 is proven only for sections affine to the waterline (g <= 1); for bulbous hulls I only claim "typically negative". (4) I did not numerically evaluate any formula; the analytic checks are: Neumann/Dirichlet limits of A, d'Alembert (zero Rankine and PV resistance), exact recovery of the crate's Michell R_w including its sign (which fixes the radiation prescription), and the reversibility symmetries for fore-aft symmetric hulls. (5) The residue-term overall sign for F_z and M (Im[F W*] etc.) depends on the transform sign convention e^{-i kx x} and on the +i0 kx prescription; both are stated explicitly and the R_w check pins the latter, but a numerical test against a known trim result (e.g. Tuck's symmetric-strut trim) would be prudent before trusting the bow-up/bow-down direction.

### Route: havelock

**Approach**

Route B, first principles. Linearised steady Bernoulli in the ship frame (stream -U xhat, z down) gives p_d = rho U phi_x; the thin-ship body condition gives a centreplane source sheet sigma = -2U f_x, and phi = INT sigma G with the Havelock Green function derived by Fourier transform in (x,y) with the free-surface condition phi_xx - nu phi_z = 0 and a Rayleigh (e^{eps t}) radiation condition. F_z and M are written as centreplane integrals of f_z phi_x and (z f_x - (x-x_ref) f_z) phi_x; each of the three parts of G (Rankine, image, free-surface) is then reduced by Fourier-space separability plus integration by parts in x and z to products of exact hull transforms integrated against a kernel over (kx,ky), with the dispersion pole split into a principal value and a residue on the Kelvin curve. Normalisation of the residue term was checked by re-deriving the crate's R_w exactly. Full working (5 stages) is in /private/tmp/claude-501/-Users-avi-michell/2fdb50af-0291-45b8-bdde-e7c422fd969b/scratchpad/derive-B.md.

**Pressure and force form**

Bernoulli with z down: p = p_atm + rho g z + rho U phi_x - (rho/2)|grad phi|^2; linear dynamic pressure p_d = rho U phi_x (phi_x>0 = flow slowed = higher pressure). Body condition on y=f: phi_y(x,0+,z) = -U f_x, source density sigma = -2U f_x (source at bow where f_x<0). Free-surface condition on z=0 (derived): phi_xx - nu phi_z = 0; elevation eta = U phi_x/g (upward).

Pressure form: outward normal on y=+f is (-f_x, 1, -f_z)/N, so upward force on both sides F_up = -2 INT_S p f_z dx dz. Errors from evaluating p at y=0 instead of y=f, from |grad phi|^2, and from the wetted sliver -eta<z<0 (where p = rho g(z+eta) = O(f) over height O(f)) are all O(f^3). Hence to Michell order
  F_up = -2 rho U INT_S f_z(x,z) phi_x(x,0,z) dx dz.        (1.1)

Divergence form: INT_V p_z dV over the volume below z=0 equals INT_wetted p n_z dS - INT_WP p dA, so F_up = INT_V p_z dV + INT_WP p_d dA. The bare '-INT_V dp/dZ' is INCOMPLETE by the waterplane term INT_WP p_d dA = rho g INT_WP eta dA (it is complete only if V is closed at the actual free surface where p=0). Thin-ship reduction INT_V p_z dV = 2 rho U INT_S f phi_xz, integrated by parts in z (f=0 at the keel) gives -2 rho U INT f(x,0) phi_x(x,0,0) dx - 2 rho U INT_S f_z phi_x; the first term cancels the waterplane term exactly, reproducing (1.1). I use (1.1) (fewest terms, no waterplane bookkeeping).

Moment: with x forward, z down, y is starboard; positive rotation about +y lifts the +x (bow) end, so bow-up M = M_y = INT[(x-x_ref) dF_up + z dF_x] about (x_ref,0,0). dF_x = 2 p f_x dx dz (both sides), dF_up = -2 p f_z dx dz:
  M = 2 rho U INT_S [ z f_x - (x - x_ref) f_z ] phi_x(x,0,z) dx dz.   (1.2)
(Check: drag F_x<0 at depth z>0 gives M<0, bow-down, like a box dragged by friction at its base. For a reference at depth z_ref add +z_ref R_w.) Hydrostatic rho g z part excluded (gives rho g Vol exactly).

**Green function**

Unit-outflow source at (x',0,z'), z'>0; X = x-x', Y = y, r^2 = X^2+Y^2+(z-z')^2, r'^2 = X^2+Y^2+(z+z')^2, k = sqrt(kx^2+ky^2), nu = g/U^2:
  G = -1/(4 pi r) + 1/(4 pi r') + G_F,
  G_F = (nu/4pi^2) INT INT dkx dky e^{i(kx X + ky Y)} e^{-k(z+z')} / ((kx + i0)^2 - nu k),
equivalently
  G = -(1/8pi^2) INT INT (dkx dky/k) e^{i(kx X+ky Y)} [ e^{-k|z-z'|} - e^{-k(z+z')} - 2 nu k e^{-k(z+z')}/((kx+i0)^2 - nu k) ]
  = -(1/8pi^2) INT INT (dkx dky/k) e^{ik.X} e^{-k|z-z'|} + (1/8pi^2) INT INT (dkx dky/k) e^{ik.X} e^{-k(z+z')} R(kx,k),
  R = ((kx+i0)^2 + nu k)/((kx+i0)^2 - nu k)  (R -> +1 as nu->0: pressure release; R -> -1 as nu->inf: rigid lid; both limits correct).
Derived by imposing phi_xx - nu phi_z = 0 at z=0 on the (x,y)-Fourier form of -1/(4 pi r).

Radiation condition: Rayleigh/Lighthill growth e^{eps t} with D/Dt = d_t - U d_x turns kx^2 into (kx + i eps/U)^2; both poles kx = +/- sqrt(nu k) move into the lower half kx-plane. For X<0 (field point aft of the source = downstream, since the fluid runs toward -x) the kx contour closes below and picks up the residues; for X>0 none. Waves trail toward -x. In polar form kx = k cos th, ky = k sin th:
  1/(k cos^2 th - nu + i0 sgn cos th) = PV - i pi sgn(cos th) sec^2 th delta(k - nu sec^2 th),
and the residue (wave) part of G_F is the classical Havelock term
  G_W = (nu/2pi) INT_{-pi/2}^{pi/2} sec^2 th e^{-nu sec^2 th (z+z')} sin(nu sec th (X + Y tan th)) dth,
whose PV partner makes the pattern vanish upstream and double downstream. Normalisation check: this G reproduces the crate's R_w = (4 rho g^2/(pi U^2)) INT_0^{pi/2} sec^3 th (I^2+J^2) dth exactly from F_x = 2 rho U INT f_x phi_x (residue part only; Rankine and PV parts of F_x vanish by kx-oddness).

**Force formula**

F_up (positive upward, hydrostatics excluded) = F1 + F2 + F3, with k = sqrt(kx^2+ky^2):

  F1 = -(rho U^2 / (2 pi^2)) INT_{-inf}^{inf} INT_{-inf}^{inf} dkx dky  kx^2 |Fo(kx,k)|^2
       [Rankine -1/4pi r plus pressure-release image +1/4pi r'; their waterline cross terms cancel identically; F1 < 0 always]

  F2 = (rho U^2 nu / pi^2) PV INT INT dkx dky  kx^2 [ Re( conj(f_wl(kx)) Fo(kx,k) ) - k |Fo(kx,k)|^2 ] / (kx^2 - nu k)
       [free-surface term, principal value across the Kelvin curve kx^2 = nu k]

  F3 = -(4 rho g nu / pi) INT_1^inf (I0 I + J0 J) lam^2 dlam / sqrt(lam^2 - 1)
     = (4 rho g nu^2/pi) INT_0^{pi/2} sec^4 th Im[ conj(f_wl(nu sec th)) Fo(nu sec th, nu sec^2 th) ] dth
       [wave/residue term; I+iJ = the crate's Michell amplitude, I0+iJ0 = INT f(x,0) e^{i nu lam x} dx its waterline analogue]

Both 2-D integrands are even in kx and depend on ky only through k, so
  INT INT dkx dky g(kx,k) = 4 INT_0^inf dkx INT_{kx}^inf g(kx,k) k dk / sqrt(k^2 - kx^2).
Equivalent grouping: Rankine alone F_up^R = -(rho U^2/2pi^2) INT INT (dkx dky/k) kx^2 Re[conj(f_wl) Fo(kx,k)] = INT_WP p_d^R dA (the missing-deck pressure); image (R=1) part = -F_up^R - F1... i.e. F1 = F_up^R + F_up^{image}. Low-Froude (rigid-lid) limit: F_up -> -(rho U^2/2pi^2) INT INT dkx dky kx^2 [ (2/k) Re(conj(f_wl) Fo) - |Fo|^2 ].

**Moment formula**

M (positive bow-up, about (x_ref, 0, 0), hydrostatics excluded) = M1 + M2 + M3 + M4, with the moment weight transform
  w^(kx,kap) = i kx Zo(kx,kap) - kap Xo'(kx,kap) + x_wl'(kx)   [transform of z f_x - (x-x_ref) f_z, valid for any complex kap]:

  M1 = -(rho U^2 / (2 pi^3)) INT_{R^3} d^3k (kx^2 / (kx^2+ky^2+kz^2)) Re[ conj(w^(kx, i kz)) Fo(kx, i kz) ]
       [Rankine; kap = i kz imaginary. Equivalent 2-D form: -(rho U^2/2pi^2) INT INT (dkx dky/k) kx^2 Re INT_0^T INT_0^T e^{-k|z-z'|} conj(w~(kx,z)) a(kx,z') dz dz', with w~, a the x-only transforms (polynomials in z per span): separable exponential moments for distinct z-spans, a nested polynomial-exponential moment on diagonal span pairs]

  M2 = (rho U^2 / (2 pi^2)) INT INT (dkx dky / k) kx^2 Re[ conj(w^(kx,k)) Fo(kx,k) ]            [pressure-release image]

  M3 = (rho U^2 nu / pi^2) PV INT INT dkx dky kx^2 Re[ conj(w^(kx,k)) Fo(kx,k) ] / (kx^2 - nu k)   [free-surface PV]

  M4 = (4 rho g nu^2 / pi) INT_1^inf lam^3 Im[ conj(w^(nu lam, nu lam^2)) Fo(nu lam, nu lam^2) ] dlam / sqrt(lam^2 - 1)
     = (4 rho g nu^2/pi) INT_0^{pi/2} sec^4 th Im[conj(w^) Fo](nu sec th, nu sec^2 th) dth              [wave/residue]

with, for real kap = k,
  Re[conj(w^) Fo] = kx Im[conj(Zo) Fo] - k Re[conj(Xo') Fo] + Re[conj(x_wl') Fo],
  Im[conj(w^) Fo] = -kx Re[conj(Zo) Fo] - k Im[conj(Xo') Fo] + Im[conj(x_wl') Fo].
Dependence on x_ref: M(x_ref) = M(0) - x_ref F_up (exact, term by term). Reference at depth z_ref: add + z_ref R_w. For a hull fore-aft symmetric about x_ref, M1 = M2 = M3 = 0 and M = M4 (trim is carried entirely by the trailing wave pattern).

**Hull transforms**

All over S = {x0<=x<=x1, 0<=z<=T}, f = 0 at both x-ends and at the keel (closed hull; a transom must be closed by the crate's appendage or end terms appear):
  Fo(kx,kap)  = INT_S f(x,z) e^{-i kx x} e^{-kap z} dx dz
  F(kx,kap)   = INT_S f_x e^{-i kx x - kap z} dx dz = i kx Fo            (crate: I + iJ = F(-nu lam, nu lam^2) = conj F(nu lam, nu lam^2))
  Fz(kx,kap)  = INT_S f_z e^{..} = -f_wl(kx) + kap Fo(kx,kap)
  Zo(kx,kap)  = INT_S z f e^{..} = -d/dkap Fo
  Xo'(kx,kap) = INT_S (x - x_ref) f e^{..} = i d/dkx Fo - x_ref Fo
  f_wl(kx)    = INT f(x,0) e^{-i kx x} dx   (waterline half-beam; I0 + iJ0 = conj f_wl(nu lam))
  x_wl'(kx)   = INT (x - x_ref) f(x,0) e^{-i kx x} dx = i d/dkx f_wl - x_ref f_wl
  w^(kx,kap)  = i kx Zo - kap Xo' + x_wl'   (transform of the moment weight z f_x - (x-x_ref) f_z)
  a_z(kx,0)   = INT f_z(x,0) e^{-i kx x} dx (only used for asymptotics)
Evaluation points: kap = k = sqrt(kx^2+ky^2) >= 0 for the image/PV/residue terms (on the residue curve kx = nu lam, kap = nu lam^2, exactly Michell's); kap = i kz (imaginary) for the 3-D Rankine moment form, or the per-span-pair e^{-k|z-z'|} kernel instead. All are the crate's moments INT t^a e^{i k t} dt, INT t^b e^{-kap t} dt with polynomial weights (x f, z f are polynomials on each span).

**Rankine direct term**

- **Contributes to force:** True
- **Contributes to moment:** True
- **Argument:** Force: F_up^R = -2 rho U INT f_z phi^R_x = -(rho U^2/2pi^2) INT INT (dkx dky/k) kx^2 Re[conj(f_wl(kx)) Fo(kx,k)]. The volume part INT_S f phi^R_xz vanishes identically (by parts in x, then P<->Q antisymmetry of d_x^2 d_z (1/r) against the symmetric weight f(P)f(Q)), so the entire Rankine vertical force equals the waterplane term 2 rho U INT f(x,0) phi^R_x(x,0,0) dx = INT_WP p_d^R dA: it is exactly the pressure the missing 'deck' of the half-body would have carried. d'Alembert applies to the closed body hull+deck in unbounded flow, not to the hull surface alone, so there is no cancellation; fore-aft symmetry does not help either (phi^R is odd in x, phi^R_x even, f_z even). Counterexample/sign: for f = b(x) g(z) with g(0)=1, F_up^R = -(rho U^2/2pi^2) INT INT (kx^2/k) |b^|^2 g^(k) dkx dky < 0 (downward; midbody Bernoulli suction on the sloping bottom). The x-force from the Rankine term is zero (conj(F)F real, odd kx kernel), which is the d'Alembert statement that does survive. Moment: M^R = -(rho U^2/2pi^3) INT d^3k (kx^2/K^2) Re[conj(w^(kx,ikz)) Fo(kx,ikz)], nonzero in general (it is the half-body analogue of the Munk moment: the double body has zero total moment, but the lower half's share of z dF_x and (x-x_ref) dF_up does not vanish). It vanishes for a hull fore-aft symmetric about x_ref (phi^R_x even, both weights odd). The Rankine kernel e^{-k|z-z'|} is not separable in (z,z'), so the moment does not collapse to waterline transforms as the force does.

**Wave residue term**

- **Contributes to force:** True
- **Contributes to moment:** True
- **Argument:** The residue at the dispersion pole is the odd-in-kx (imaginary) half of the i0 prescription, so it picks out Im of the kx-integrand P = (i kx) conj(w^) F. Force: with w = f_z, P = kx^2 conj(f_wl) Fo - kx^2 k |Fo|^2; the |Fo|^2 self-term is real and gives no residue, only the waterline cross term survives: F3 = -(4 rho g nu/pi) INT_1^inf (I0 I + J0 J) lam^2 dlam/sqrt(lam^2-1), the cross-correlation of Michell's amplitude with the waterline's. It is identically zero for fore-aft symmetric hulls (all amplitudes real) and for separable hulls f = b(x)g(z) (amplitudes differ by the factor -i nu lam g^), so it is a small 'depth-warp' correction in general; this matches the flow-reversal argument: F_z of a symmetric hull cannot depend on which way the waves trail, and the residue term is the only piece odd under reversal. Moment: M4 = (4 rho g nu^2/pi) INT lam^3 Im[conj(w^) Fo] dlam/sqrt(lam^2-1) with Im[conj(w^)Fo] = -kx Re[conj(Zo)Fo] - k Im[conj(Xo')Fo] + Im[conj(x_wl')Fo]; the first term is nonzero even for symmetric hulls, so the wave term carries the whole trim of a symmetric hull (flow reversal flips bow-up to bow-down, so the PV/image/Rankine parts of M vanish there). Normalisation of the residue is confirmed by reproducing the crate's R_w exactly from the same machinery.

**Low froude limit**

nu -> inf at fixed hull: the pole k0 = nu sec^2 th runs off to k = inf where the transforms are negligible, so F3, M4 -> 0 and the PV kernel nu/(kx^2 - nu k) -> -1/k pointwise, giving the rigid-lid (double-body, G = -1/4pi r - 1/4pi r') result
  F_up -> F_DB = -(rho U^2/2pi^2) INT INT dkx dky kx^2 [ (2/k) Re(conj(f_wl(kx)) Fo(kx,k)) - |Fo(kx,k)|^2 ] = 2 F_up^R + (rho U^2/2pi^2) INT INT kx^2 |Fo|^2.
Sign: Fo = INT_0^T e^{-kz} a(kx,z) dz with a(kx,0) = f_wl, so |Fo| <= max_z|a|/k. For sections that shrink with depth without strong warp (|a(kx,z)| <= |f_wl(kx)| and phases within 60 deg) the bracket is >= |Fo|^2 >= 0 and F_DB < 0: the hull SINKS (positive sinkage, force ~ rho U^2 L^2, i.e. squat proportional to Fn^2). Exact for f = b(x)g(z), g(0)=1, 0<=g<=1: bracket = |b^|^2 g^(k)(2/k - g^(k)) > 0. Physics: at low speed the free surface acts as a rigid lid; the double-body flow accelerates round the hull and lowers the pressure on the sloping bottom (Bernoulli suction); the mirror 'deck' that would cancel this in unbounded flow is absent, so the net linear force is downward. Equivalent divergence-form view: F_DB = INT_WP p_d dA + INT_V p_z dV with p_d = rho U phi_x^DB < 0 over most of the waterplane. Moment in the same limit: for a fore-aft symmetric hull M -> 0 (the double-body flow is symmetric); low-speed trim comes only from fore-aft asymmetry of the sections and the wave-induced trim vanishes with R_w. I commit to: linear theory predicts sinking as nu -> inf for ordinary hulls; a pathological hull with |a(kx,z)| >> |f_wl| (e.g. a bulb far wider than the waterline) could flip the sign of F_DB and the formula shows exactly how.

**Numerics**

Principal value: use polar (k, th) at fixed th (kx = k cos th, ky = k sin th, dkx dky = k dk dth); the pole in k at k0 = nu sec^2 th is simple with a smooth numerator N(k) built from exact transforms, so PV INT N(k)/(k - k0) dk = INT [N(k) - N(k0)]/(k - k0) dk + N(k0) ln|(k_max - k0)/(k0 - k_min)| (or symmetric Gauss panels about k0). The residue half is F3/M4 and must not be added again. The residue integrals have Michell's weight lam^2/sqrt(lam^2-1) (force) or lam^3/sqrt(lam^2-1) (moment): th = arcsec(lam) removes the integrable lam=1 endpoint singularity and the crate's adaptive Gauss-Legendre outer quadrature applies unchanged; on the curve the transforms are evaluated at exactly Michell's (kx, kap) = (nu lam, nu lam^2), and f_wl, x_wl' are the same oscillatory moments on the z=0 boundary. Rankine parts: no singularity (kx^2/K^2 bounded, 1/k integrable at k=0 in 2-D).

Decay: large-k behaviour is set by the waterline strip, Fo(kx,k) = f_wl(kx)/k - a_z(kx,0)/k^2 + O(k^-3). F1: kx^2|Fo|^2 ~ kx^2|f_wl|^2/k^2, ky-integral gives pi|kx||f_wl(kx)|^2, absolutely convergent iff f(x,0) is continuous (f_wl decays faster than 1/kx): a bare transom step diverges logarithmically (the Rankine pressure at a blunt step is log-singular), so the transom closure must be applied to f before forming the transforms. F2: the 1/k terms in the bracket cancel exactly, leaving Re(conj(f_wl) a_z(kx,0))/k^2 + O(k^-3) against a kernel -> nu, i.e. ~|f_wl||a_z0|/|kx| after ky, ~|kx|^-5 for a hull with an entrance-angle discontinuity. F3/M4 along the Kelvin curve: Im[conj(f_wl) Fo] ~ lam^-6, so the integrand decays like lam^-4 or better (faster than R_w's, since |F|^2 ~ lam^-2 . lam^{-4}... comparable or better). M1-M3: the x-weights do not change the kx-decay (same endpoint smoothness), the z-weights gain one power of k. Cost: F1/F2/M2/M3 are 2-D quadratures of exact transforms (evaluate on a (kx, k) grid with the k-inner integral 4 INT_kx^inf ... k dk/sqrt(k^2-kx^2), or in polar); M1 needs either an imaginary-kap 3-D quadrature (kz-decay 1/kz^2 . 1/K^2, fine) or, preferably, the per-z-span-pair e^{-k|z-z'|} evaluation, which is exact and O(n_z^2) per (kx,k) node.

**Unsure about**

(1) Orientation of the moment sign: I derived y = starboard and 'positive rotation about +y lifts the bow' twice, and checked with the drag-at-depth heuristic, but a right-hand-rule slip here would flip the sign of every M term uniformly; verify numerically with a hull having a single obvious asymmetry. (2) The Rankine vertical force F_up^R reduces to a waterline transform by an antisymmetry argument (INT f phi^R_xz = 0); I checked it two independent ways (Fourier by-parts and the P<->Q symmetry of d_x^2 d_z (1/r)) but it is the least 'standard' step, worth a brute-force numerical check on one hull. (3) The wave residue part of F_z vanishing for separable and symmetric hulls is a strong structural claim that follows from the algebra; if a numerical experiment shows a nonzero residue contribution for f = b(x)g(z), the Fz = -f_wl + k Fo by-parts step (which needs f = 0 at the keel) is the place to look. (4) For a transom hull, the same integration domain S is used for the source sheet (hull + virtual appendage) and for the pressure integration; physically the pressure on the appendage acts on water, not on the hull, so the weights f_z, w in (1.1)/(1.2) should arguably be restricted to the real hull while F/Fo in the phi factor include the appendage. All formulas above assume one common closed S. (5) The low-Froude sign argument is rigorous only under the stated 'no strong depth-warp' condition on a(kx,z); it is exact for separable hulls. (6) I did not attempt to reduce the Rankine moment M1 beyond the 3-D/per-span-pair forms; it may admit a further identity I did not find.

### Route: pressure

**Approach**

Start from the linearised Bernoulli pressure p = rho g z + rho U phi_x on the actual wetted surface y = +/-f(x,z) (z downward, bow at high x, stream -U x-hat), integrate p n_z over both sides including the strip up to the wave elevation zeta = (U/g) phi_x|_{z=0}, and show that at second order the hydrostatic part gives exactly rho g V_0 while the dynamic part gives F_z = 2 rho U int int phi_x (-f_z) dx dz; integration by parts in z turns this into the volume form plus the explicit waterline (Bernoulli-hollow buoyancy) term 2 rho U int f(x,0) phi_x(x,0,0) dx. Represent phi by the centreplane source sheet sigma = -2U f_x with the Havelock Green function written as a 2-D Fourier integral in (kx,ky) (Rankine + image + wave), which makes every x,z derivative algebraic and separates the image/wave kernels into products of hull transforms; the non-separable Rankine part is handled by a 2-D Fourier transform on the centreplane. The resistance drops out of the same machinery and reproduces the crate's Michell integral exactly, which fixes normalisation and the radiation sign.

**Pressure and force form**

Conventions: z downward (gravity +z), bow at high x, u = -U x-hat + grad phi, y = +/-f, f>=0, f=0 at both x-ends and at z=T (flat keel = distributional f_z at z=T). Bernoulli with constant fixed far upstream: p = rho g z + rho U phi_x + O(f^2) (stagnation at the bow: phi_x>0, p_d>0). Elevation (upward) zeta = (U/g) phi_x|_{z=0}; FS condition phi_xx - nu phi_z = 0 on z=0 (minus sign because z is down; check: e^{-kz}cos kx needs k = nu). Body condition phi_y(x,0+,z) = -U f_x; source sheet sigma = -2U f_x.
Force: dF = -p n dS; F_up = closed-int p n_z dS; on y=+f, n_z dS = -f_z dx dz, so F_up = 2 int int_{wetted} p (-f_z) dx dz. Split p = rho g z + p_d over the ACTUAL wetted surface (to z=-zeta): (a) hydrostatic: closing with the cap on z=-zeta (n=-z-hat) gives int_wetted z n_z dS = V - int zeta dA = V_0, i.e. exactly the still-water buoyancy, no strip correction at 2nd order; (b) dynamic: the strip -zeta<z<0 contributes O(zeta f_z phi_x) = O(f^3). Hence, buoyancy removed,
  F_z (up) = 2 rho U int_{x_s}^{x_b} dx int_0^T phi_x(x,0,z) (-f_z) dz   [surface form: Bernoulli pressure on the sloped bottom; NO separate waterline term].
Integrating by parts in z (f(x,T)=0): int_0^T phi_x(-f_z) dz = f(x,0) phi_x(x,0,0) + int_0^T f phi_xz dz, so
  F_z = rho U int_{V_0} phi_xz dV + 2 rho U int f(x,0) phi_x(x,0,0) dx = int_{V_0}(-dp_d/dZ_up) dV + 2 rho g int f(x,0) zeta(x) dx   [volume form + explicit loss-of-buoyancy of the strip 2 f(x,0) zeta].
The two forms agree identically; the elevation term is the boundary term of the z integration by parts. In wavenumber space the same split is: transform of (-f_z) = W(kx) - kappa F(kx,kappa) (waterline minus kappa x volume). I use the surface form as master and the volume form for interpretation and the low-Froude sign.
Moment (bow-up = rotation about +y-hat, taking +x toward -z) about (x_ref, 0, z_ref): M = closed-int p[(x-x_ref) n_z - (z-z_ref) n_x] dS = 2 int int p[(x-x_ref)(-f_z) + (z-z_ref) f_x] dx dz; hydrostatic parts are the buoyancy moment (excluded) and zero (int f_x dx = 0). So M = M_v + M_h with M_v = 2 rho U int int (x-x_ref) phi_x (-f_z) dx dz (moment of the vertical-force distribution) and M_h = 2 rho U int int (z-z_ref) phi_x f_x dx dz = -int int (z-z_ref) dR (resistance acting below the reference point gives bow-down). M_h is the same order as M_v; drop it only if the caller wants the vertical distribution alone. Resistance from the same pressure: R = -2 rho U int int phi_x f_x dx dz.

**Green function**

Laplacian_P G = delta(P-Q), (d_xx - nu d_z)G = 0 on z=0, G -> 0 as z -> inf, waves only for X = x-x' < 0 (downstream).
G = G_R + G_I + G_w,
 G_R = -1/(4 pi r),  r = sqrt(X^2+Y^2+(z-z')^2);  G_I = -1/(4 pi r'), r' = sqrt(X^2+Y^2+(z+z')^2);
 G_w = -(1/(4 pi^2)) int int dkx dky  kx^2 e^{i(kx X + ky Y)} e^{-k(z+z')} / [ k (k nu - kx^2 - 2 i mu kx) ],  k = sqrt(kx^2+ky^2), mu -> 0+.
Derivation: 2-D Fourier of -1/(4 pi r) = -(1/(8 pi^2)) int int (1/k) e^{ik.X} e^{-k|z-z'|} dk; add -(1/(8 pi^2)) int int (A/k) e^{ik.X} e^{-k(z+z')} dk; the FS condition at z=0 gives (-kx^2 - nu k) + A(-kx^2 + nu k) = 0, A = (k nu + kx^2)/(k nu - kx^2) = 1 + 2kx^2/(k nu - kx^2). Limits: nu -> inf: A -> 1 (rigid lid, double body -1/r - 1/r'); nu -> 0: A -> -1 (G = 0 on z=0). Dispersion pole k nu = kx^2, i.e. k = nu sec^2(theta) = nu lambda^2 with kx = k cos theta — the crate's exp(-nu lambda^2 z).
Radiation condition: Rayleigh damping (flow built up as e^{eps t}) turns the FS condition into (d_x - mu)^2 phi - nu phi_z = 0, mu = eps/U > 0, replacing kx^2 by (kx + i mu)^2 ~ kx^2 + 2 i mu kx. With k in (-inf,inf), theta in (-pi/2,pi/2) both k-poles have Im k < 0, so residues appear only where X cos theta + Y sin theta < 0: waves trail toward -x. Plemelj: 1/(k nu - kx^2 - 2 i mu kx) = PV 1/(k nu - kx^2) + i pi sgn(kx) delta(k nu - kx^2). Havelock polar form: G_w = -(1/(4 pi^2)) int_{-pi}^{pi} dtheta cos^2 theta int_0^inf k dk e^{ik(X cos th + Y sin th)} e^{-k(z+z')} / (nu - k cos^2 th - i0 sgn(cos th)).
On y=0: G = int int dk e^{i kx (x-x')} g^(k;z,z') with g^_R = -(1/(8 pi^2 k)) e^{-k|z-z'|} (non-separable), g^_I = -(1/(8 pi^2 k)) e^{-k(z+z')}, g^_w = -(kx^2/(4 pi^2 k)) e^{-k(z+z')} [PV 1/(k nu - kx^2) + i pi sgn(kx) delta(k nu - kx^2)]. Derivatives: d_x -> i kx, d_z -> -k (image/wave), -k sgn(z-z') (Rankine). Verification of sign and constant: R = 4 rho U^2 int int dk (i kx) h^(k) |S|^2 gets zero from image and PV (odd in kx) and from the residue (i kx)(i pi sgn kx) = -pi|kx| gives R = (rho U^2/pi) int int (|kx|^3/k) delta(k nu - kx^2)|S|^2 = (4 rho g^2/(pi U^2)) int_1^inf (I^2+J^2) lambda^2/sqrt(lambda^2-1) dlambda > 0 — the crate's Michell integral exactly (the opposite sgn would give negative resistance).

**Force formula**

With nu = g/U^2, lambda = sec theta, kx = k cos theta, and transforms S, F, W, Z = W - kappa F defined below:
F_z (positive up, buoyancy excluded) = F^R + F^I + F^PV + F^res, where
 F^R  = -(2 rho U^2/pi^2) int_0^{pi/2} cos th dth int_0^inf k dk Im[ conj(W(k cos th)) S(k cos th, k) ]
 F^I  = -(2 rho U^2/pi^2) int_0^{pi/2} cos th dth int_0^inf k dk Im[ conj(Z(k cos th, k)) S(k cos th, k) ]
 (F^R + F^I = -(2 rho U^2/pi^2) int_0^{pi/2} cos th dth int_0^inf k dk Im[ (2 conj W - k conj F) S ](k cos th, k) = the double-body / low-Froude force)
 F^PV = -(4 rho g^2/(pi^2 U^2)) int_0^{pi/2} lambda^3 dth  PV int_0^inf s^2 Im[ conj(Z) S ](nu lambda s, nu lambda^2 s) ds/(1 - s)
 F^res = -(4 rho g^2/(pi U^2)) int_0^{pi/2} lambda^3 Re[ conj(Z) S ](nu lambda, nu lambda^2) dth
       = -(4 rho g^2/(pi U^2)) int_1^inf Re[conj(Z) S] lambda^2/sqrt(lambda^2-1) dlambda   (Michell's integral with |S|^2 -> Re[conj(Z) S]; in crate notation conj S(nu lambda, nu lambda^2) = I + iJ).
Unfolded (kx,ky)-plane form (theta in (-pi,pi), dk = dkx dky): F^I = -(rho U^2/(2 pi^2)) int int (kx/k) Im[conj(Z)S]; F^PV = -(rho U^2/pi^2) int int (kx^3/k) PV[1/(k nu - kx^2)] Im[conj(Z)S]; F^res = -(rho U^2/pi) int int (|kx|^3/k) delta(k nu - kx^2) Re[conj(Z)S]; F^R same as F^I with Z -> W. Equivalent plane-transform form of F^R: -(rho U^2/(2 pi^2)) int int dkx dkz (kx^2/|K|) Re[fhat(K) conj W(kx)] (checked numerically equal, -0.2856 vs -0.2858 for the test hull). In every line, Z = W - k F splits into the waterline (loss-of-buoyancy) part W and the sloped-bottom volume part -kF.

**Moment formula**

M (positive bow-up) about (x_ref, 0, z_ref), z_ref measured downward, = M_v + M_h.
Define the operator, for any P-side transform A(kx,kappa):
 Phi[A] = -(2 rho U^2/pi^2) int_0^{pi/2} cos th dth int_0^inf k dk Im[conj(A) S](k cos th, k)
          -(4 rho g^2/(pi^2 U^2)) int_0^{pi/2} lambda^3 dth PV int_0^inf s^2 Im[conj(A) S](nu lambda s, nu lambda^2 s) ds/(1-s)
          -(4 rho g^2/(pi U^2)) int_0^{pi/2} lambda^3 Re[conj(A) S](nu lambda, nu lambda^2) dth      (image + wave PV + wave residue).
Then F_z = F^R + Phi[Z], and
 M_v = M_v^R + Phi[Z_x],  Z_x = W_x - kappa F_x  (transform of (x - x_ref)(-f_z)),
   M_v^R = -(rho U^2/(2 pi^2)) int int dkx dkz (kx^2/|K|) Re{ fhat(K) conj[ W_x(kx) - i kz fhat_x(K) ] },  |K| = sqrt(kx^2+kz^2);
 M_h = M_h^R + Phi[S_z],  S_z = -dS/dkappa - z_ref S  (transform of (z - z_ref) f_x),
   M_h^R = +(rho U^2/(4 pi^2)) int int dkx dkz |fhat(K)|^2 kx^3 kz/|K|^3   (z_ref drops out: total Rankine drag is zero).
Physical content: M_v is the moment of the vertical (sloped-bottom + waterline) pressure distribution, M_h = -int int (z - z_ref) dR is the bow-down moment of the resistance distribution acting below the reference point (for a fore-aft symmetric hull it reduces to the depth moment of R_w, from the residue only). M_v^R and M_h^R vanish identically for fore-aft symmetric hulls about midships; for such hulls M_v comes ONLY from the residue (trim is pure wave radiation) while sinkage comes from double body + PV.

**Hull transforms**

All over the closed hull region x_s <= x <= x_b, 0 <= z <= T (f = 0 at both x-ends and at z = T; a flat keel is included via the distributional f_z at z = T):
 S(kx,kappa)  = int int f_x(x,z) e^{-i kx x} e^{-kappa z} dx dz   (the crate's amplitude: I + iJ = S(-nu lambda, nu lambda^2) = conj S(nu lambda, nu lambda^2))
 F(kx,kappa)  = int int f e^{-i kx x} e^{-kappa z} dx dz
 W(kx)        = int f(x,0) e^{-i kx x} dx                          (waterline transform; 1-D oscillatory moments at z = 0)
 Z(kx,kappa)  = int int (-f_z) e^{-i kx x} e^{-kappa z} dx dz = W(kx) - kappa F(kx,kappa)
 F_x, W_x     = F, W with weight (x - x_ref);  Z_x = W_x - kappa F_x = transform of (x - x_ref)(-f_z)
 S_z(kx,kappa)= int int (z - z_ref) f_x e^{-i kx x} e^{-kappa z} dx dz = -dS/dkappa - z_ref S  (exp moments of one higher power in z)
 fhat(K)      = int int f e^{-i(kx x + kz z)} dx dz,  K = (kx,kz)  (plane transform; oscillatory in z too: exp_moments_complex with kappa = i kz, Re kappa = 0 >= 0)
 fhat_x(K)    = int int (x - x_ref) f e^{-i(kx x + kz z)} dx dz
 hhat_z(K)    = int int (z - z_ref) f e^{-i K.P} = i d_kz fhat - z_ref fhat (only needed to derive M_h^R; the final M_h^R uses |fhat|^2).
All are evaluated at (kx,kappa) = (k cos th, k) for the image/Rankine lines, (nu lambda s, nu lambda^2 s) for the PV line, (nu lambda, nu lambda^2) for the residue line. Useful identities: S = i kx F (f closes in x); S(0,kappa) = 0; S = O(kx) as kx -> 0.

**Rankine direct term**

- **Contributes to force:** True
- **Contributes to moment:** True
- **Argument:** Generic block: 2 rho U int int a phi^R_x = -(rho U^2/(2 pi^2)) int int d^2K (kx^2/|K|) Re[fhat(K) conj A^(K)] with A^ the plane transform of the weight a. For F_z, a = -f_z, A^ = W(kx) - i kz fhat(K). The |fhat|^2 piece has kernel kx^2 kz/|K|, ODD under K -> -K, so it integrates to zero: rho U int_V phi^R_xz dV = 0, the linearised d'Alembert theorem (direct space: kernel f(P) f(Q) d_X^2 d_Z (1/r) is antisymmetric under P <-> Q). What survives is exactly the waterline/lid term F^R = 2 rho U int f(x,0) phi^R_x(x,0,0) dx = -(rho U^2/(2 pi^2)) int int d^2K (kx^2/|K|) Re[fhat conj W] = -(2 rho U^2/pi^2) int_0^{pi/2} cos th dth int k dk Im[conj(W) S](k cos th, k), which is generically NON-zero and negative (phi_x < 0 along the midbody where f(x,0) is largest). Interpretation: the source sheet on 0<z<T in unbounded fluid is a closed body with a flat lid on z=0; d'Alembert says lid + sides carry zero net force, so the sides carry minus the lid's Bernoulli-suction force. Counterexample to 'Rankine gives nothing': the separable hull f = b(1-x^2)(1-z/T), T=0.125, gives F^R = -0.286 rho U^2 b^2 (37% of the total double-body force -0.455). Moment: the lidded half body is a top/bottom-asymmetric (cambered) body, so unlike the classical Munk moment (proportional to sin 2 alpha, zero at zero incidence for a symmetric body) it has a zero-incidence moment (m_13-type coupling). Explicitly M_v^R = -(rho U^2/(2 pi^2)) int int d^2K (kx^2/|K|) Re{fhat conj[W_x - i kz fhat_x]} (the x-weight breaks the P<->Q antisymmetry; symmetrised kernel proportional to kx kz (kx^2 + 2 kz^2)/|K|^3 on |fhat|^2) and M_h^R = +(rho U^2/(4 pi^2)) int int d^2K |fhat|^2 kx^3 kz/|K|^3 (depth distribution of a drag whose total is zero). Both vanish identically for fore-aft symmetric hulls about midships (kernels odd in kx, |fhat|^2 even in kx), as a symmetric body at zero incidence must have a symmetric pressure field.

**Wave residue term**

- **Contributes to force:** True
- **Contributes to moment:** True
- **Argument:** The residue i pi sgn(kx) delta(k nu - kx^2) inserted in the generic block gives F^res = -(rho U^2/pi) int int dk (|kx|^3/k) delta(k nu - kx^2) Re[conj(Z) S] = -(4 rho g^2/(pi U^2)) int_1^inf Re[conj(Z) S] lambda^2/sqrt(lambda^2-1) dlambda: Michell's integrand |S|^2 replaced by the interference Re[conj(Z) S] between the (-f_z) spectrum and the f_x spectrum on the free-wave curve. It is real, not sign-definite, and generically non-zero. It vanishes IDENTICALLY for fore-aft symmetric hulls: W, F real and S = i kx F make conj(Z) S purely imaginary; in direct space d_X G_res is even in X while sigma = -2U f_x is odd and -f_z is even. For those hulls the whole free-surface effect on sinkage is the PV (local, non-radiating) part F^PV = -(rho U^2/pi^2) int int (kx^3/k) PV[1/(k nu - kx^2)] Im[conj(Z) S], which is non-zero for all hulls (test hull: F^PV = -0.25 rho U^2 b^2 at Fn 0.41, +0.25 at Fn 0.71). For the moment the residue contributes through Phi_res[Z_x] and Phi_res[S_z]; for a symmetric hull about midships it is the ONLY contribution to M_v (the (x - x_ref) weight flips the parity: PV, image and Rankine give zero M_v), so trim of a symmetric hull is pure wave radiation, and Phi_res[S_z] = -int (z - z_ref) dR_w is the depth moment of the wave resistance. Both PV and residue therefore enter F_z and M; only R_w is residue-only.

**Low froude limit**

nu -> inf at fixed hull: h^_w = -(kx^2/(4 pi^2 k))/(k nu - kx^2) -> 0 like 1/nu at every fixed k and the residue moves to k = nu lambda^2 -> inf where the transforms decay, so F^PV, F^res -> 0 (numerically O(1/nu): -0.030, -0.011 at nu = 20, 50 for the test hull). G -> -1/(4 pi r) - 1/(4 pi r'), the rigid-lid double-body flow, and F_z -> F^DB = -(2 rho U^2/pi^2) int_0^{pi/2} cos th dth int k dk Im[(2 conj W - k conj F) S] = 4 rho U int f(x,0) phi^R_x(x,0,0) dx + rho U int_V phi^I_xz dV = 2 rho g int 2 f(x,0) zeta^DB dx + (sloped-bottom Bernoulli suction, volume form), zeta^DB = (U/g) phi^DB_x|_{z=0}. SIGN: NEGATIVE — linear theory predicts the hull SINKS (squats) at low Froude number, with force of order rho U^2 x (waterline area-like) that does not vanish as nu -> inf (it is the U^2-scaling double-body force; zeta^DB itself -> 0 like U^2/g). Physics: the double-body flow accelerates along the midbody (phi_x < 0, Bernoulli hollow, depressed waterline) where the waterline half-beam f(x,0) is largest, and decelerates only briefly at the stagnation ends; the loss of buoyancy of the hollow strip dominates, the sloped-bottom term (phi_xz > 0 where phi_x < 0 since the disturbance decays with depth) offsets at most part of it. Formula check: for any hull with separable sections f = f0(x) g(z), g(0)=1, g(T)=0, g' <= 0, the integrand is kx |f0^|^2 ghat (2 - k ghat) with ghat > 0 and 2 - k ghat = 1 + int(-g') e^{-kz} dz > 0, so F^DB < 0 strictly; the waterline part carries weight 2 and the volume part -k ghat in (-1,0), i.e. the elevation term is at least twice the magnitude of the sloped-bottom term and of opposite sign. Numerical test (f = b(1-x^2)(1-z/T), T = L/16, rho = U = 1): F^R = -0.286, F^I = -0.170, F^DB = -0.455 b^2; total F_z = -0.466 (Fn 0.10), -0.524 (0.22), -0.707 (0.41), -0.542 (0.50), -0.205 (0.71): sinkage at all speeds, maximum near Fn 0.4 and rising toward high Fn, matching the classical thin-ship/Wigley pattern (Tuck 1966). Low-Froude trim M^DB: generically non-zero, sign-indefinite, identically zero for fore-aft symmetric hulls; no sign commitment.

**Numerics**

Radiation/PV: the Rayleigh prescription fixes the split 1/(k nu - kx^2 - i0 sgn kx) = PV + i pi sgn(kx) delta. Residue line: regular theta (or lambda) integral, identical quadrature to the crate's Michell R_w with |S|^2 -> Re[conj(A) S], same lambda cutoff. PV line: per theta, substitute k = nu lambda^2 s so the pole sits at s = 1, and subtract it: PV int_0^inf g(s)/(1-s) ds = int_0^2 [g(s) - g(1)]/(1-s) ds + int_2^inf g(s)/(1-s) ds with g(s) = s^2 Im[conj(A) S](nu lambda s, nu lambda^2 s) (PV int_0^2 ds/(1-s) = 0); g is entire in s because the transforms are entire in (kx,kappa), so the subtracted integrand is smooth (implemented and verified in the test script). As theta -> pi/2 the pole goes to k -> inf where the transforms are negligible: truncate theta as for R_w. Image and Rankine lines: no singularity; kx/k and kx^2/|K| are bounded at the origin. Near theta -> pi/2, S = O(kx) because the hull closes (S(0,kappa) = 0) and the explicit cos theta factors keep the integrands bounded. Large-wavenumber decay: e^{-kz} localises the z integrals at the waterline, S ~ (i kx/k) W(kx), Z ~ (1/k) int(-f_z(x,0)) e^{-i kx x} dx, F ~ W/k, so conj(A) S = O(k^-2) times products of 1-D x-transforms, which for spline hulls decay like kx^-2 (slope discontinuity where f -> 0 at the ends or C^1 knots) up to kx^-4; the image integrand k dk (kx/k) conj(A) S decays like k^-4 at fixed theta, the PV integrand like k^-1 x (transform decay), the residue like R_w's. Plane-transform Rankine moment pieces: fhat(kx,kz) ~ W(kx)/(i kz) from the waterline jump and the kernel kx^2/|K| grows like |K|, so the polar integrand is O(|K|^-3): absolutely convergent but slow; use polar coordinates with a tail estimate/Richardson on the cutoff or subtract the waterline asymptote analytically. Prefer the (kx,ky) route for F^R (same S, W as everything else); the plane route was used only as a cross-check (agreement -0.2856 vs -0.2858) and is required only for M_v^R, M_h^R. All 2-D integrals converge absolutely; no regularisation beyond the single PV is needed. Transform evaluation uses the crate's closed-form moments: osc_moments in x (with weight x for the x_ref terms), exp_moments in z (one extra power for S_z), exp_moments_complex with kappa = i kz for fhat, fhat_x.

**Unsure about**

(1) The magnitude of the Rankine/image constants was verified only through (a) exact recovery of the crate's Michell R_w from the same kernel and (b) numerical agreement of two independent representations of F^R; I did not independently verify F^I or the PV line against a direct-space computation, so a factor-of-2 slip in F^I or F^PV, while unlikely (they share the identical algebra with R_w), is the most likely residual error. (2) The vanishing of the extra wetted strip's dynamic pressure at second order assumes the linearised zeta and a hull with finite f_z at the waterline; for a hull with a horizontal flat at z=0 (f_z -> infinity there) the O(zeta^2) strip term is not negligible and the classical thin-ship expansion is non-uniform. (3) The interpretation of the Rankine moment as an m_13-type camber moment is by analogy; the formulas themselves follow directly. (4) M_h requires a z_ref; if the caller's 'moment about station x_ref' means the vertical-force distribution only, use M_v alone. (5) The claimed decay rates assume the crate's usual spline continuity; a bare transom step (Fixed{length: 0}) degrades the kx decay to kx^-1 and makes the image/PV tails converge only like k^-2 to k^-3, requiring much larger cutoffs, exactly as the crate already notes for R_w.

## Literature and experimental anchors

### Literature set 1

**Sources**

- **Note on the sinkage of a ship at low speeds**
  - *Authors:* T. H. Havelock
  - *Year:* 1939
  - *Where:* Z. Angew. Math. Mech. (ZAMM) 19, 202-205; reprinted in 'Collected Papers of Sir Thomas Havelock on Hydrodynamics' (ONR/ACR-103, 1965) pp. 458-461 — full text read from https://archive.org/download/collectedpaperso00have/collectedpaperso00have.pdf (pdf pages 472-475)
  - *What it gives:* The canonical DEEP-WATER low-speed result. Assumes the flow is the double-body (rigid free surface) potential flow, so the vertical force is the 'defect of vertical pressure' Q = -rho ∬(U φ_x + ½|∇φ|²) n dS over the immersed half-ellipsoid; Q is downward and 'should be proportional to the square of the speed'. Equivalent sinkage h defined by Q = rho g (pi a b) h (waterplane area × h). Exact closed forms for a half-ellipsoid (eqs 12,14,16) and Table I of gh/U²: L/D=10: B/D=1 -> 0.0253, 2 -> 0.0453, 3 -> 0.0612, 4 -> 0.0735; L/D=16: B/D=1 -> 0.0138, 2 -> 0.0231, 3 -> 0.0318, 4 -> 0.0397. Thus h/L = (gh/U²)·Fn². Compared with Horn's empirical model formula (h = 0.0283 U²/g vs 0.0231 for L/B=8, B/D=2) and with Amtsberg's body of revolution (Q ≈ 0.0284 rho U² × section area).
- **Shallow-water flows past slender bodies**
  - *Authors:* E. O. Tuck
  - *Year:* 1966
  - *Where:* J. Fluid Mech. 26, 81-95 (not accessed directly; formulas transcribed from Gourlay 2011 review, eqs 4-10)
  - *What it gives:* SHALLOW-water slender-body sinkage/trim (squat). Ship = line of sources with strength ∝ U S'(x)/h in the (x,y) plane; hull b.c. φ_y = ±(U/2h) dS/dx on y=0±; subcritical potential φ = (U/(4π h sqrt(1-Fh²))) ∫ S'(ξ) ln[(x-ξ)² + (1-Fh²) y²] dξ. Pressure from Bernoulli, sinkage & trim by hydrostatics. For a fore-aft symmetric hull: sinkage ≠ 0, trim = 0 at subcritical speed; supercritical version gives trim ≠ 0, sinkage = 0. Sinkage ∝ Fh²/sqrt(1-Fh²) — this is the depth-Froude squat law and is NOT the deep-water law.
- **A brief history of mathematical ship-squat prediction, focussing on the contributions of E.O. Tuck**
  - *Authors:* T. P. Gourlay
  - *Year:* 2011
  - *Where:* J. Eng. Math. 70, 5-16; author preprint http://www.perthhydro.com/pdf/Gourlay2011HistorySquatTuck.pdf (read in full)
  - *What it gives:* Transcribes Tuck's shallow-water formulas with conventions (x positive AFT from midships, z positive UP, θ positive BOW-DOWN): s_LCF = c_s (∇/L²) Fh²/sqrt(1-Fh²) (eq 6); c_s = (L²/(2π ∇ A_WP)) ∫∫ (dS/dξ) B(x)/(x-ξ) dξ dx (eq 7), c_s ≈ 1.3-1.5 for all hulls, 1.5 recommended; low-Fh approx s_LCF = 1.5 (∇/L²) U²/(g h) (eq 8); θ = c_θ (∇/L³) Fh²/sqrt(1-Fh²) (eq 9); c_θ = -(L³/(2π ∇ I_LCF)) ∫∫ (dS/dξ)(x-LCF) B(x)/(x-ξ) dξ dx (eq 10), c_θ = 0 for fore-aft symmetric hulls. Confirms Tuck 1964 (JSR 8:15-23) as the deep-water slender-ship asymptotic paper and Havelock 1939 as the deep-water low-speed sinkage reference; 1960s model tests showed both bow and stern sink at moderate speed.
- **The maximum sinkage of a ship**
  - *Authors:* T. P. Gourlay & E. O. Tuck
  - *Year:* 2001
  - *Where:* J. Ship Research 45(1), 50-58; https://cmst.curtin.edu.au/wp-content/uploads/sites/4/2016/05/gourlay-2001-the_maximum_sinkage_of_a_ship.pdf (read in full, eq page rendered)
  - *What it gives:* Finite-depth slender-body theory (Tuck & Taylor 1970 corrected): F = F_∞ + F_d, M = M_∞ + M_d, where F_∞, M_∞ are the force/moment on the lower half of the equivalent DOUBLE BODY in unbounded fluid (the deep-water, zero-Froude limit) and F_d, M_d are finite-depth corrections. F_∞ approximated by Havelock's (1939) spheroid: F_∞ = rho U² A_W ε² ( ln(ε/2) + 3/2 − ε ) (eq 23), ε = sqrt(12 C_V/π) = beam/length of equivalent spheroid, C_V = V/L³; F_∞ < 0 for small ε, i.e. DOWNWARD, ∝ U². 'For fore-aft symmetric ships M_∞ is identically zero.' Depth corrections: F_d = -(rho U²/4π²) ∫ k² S̄(k) B̄*(k) A(k) dk, M_d = (rho U²/4π²) ∫ k² S̄(k) (xB)̄*(k) A(k) dk (eq 25), A(k) = -2 ∫_{|k|}^∞ [1 + q/(Fh² k² h − q tanh(qh))] dq/sqrt(q²−k²) (eq 26). Also notes 'the calculation of F_∞ and M_∞ is a well known but difficult problem'. Finite-width channel version eq 20 with coth(λw) kernel. Sign: F positive up (s = -F/(rho g A_W) in Gourlay 2003), M positive bow-up (M = rho g I_W θ, θ bow-up, Gourlay 2003).
- **Sinkage and trim in first-order thin-ship theory**
  - *Authors:* R. W. Yeung
  - *Year:* 1972
  - *Where:* J. Ship Research 16(1), 47-59 (abstract only, via TRID https://trid.trb.org/View/6630 and OnePetro listing; full text paywalled)
  - *What it gives:* The classical Michell-theory sinkage/trim paper. 'Sinkage and trim of a ship moving with constant speed into still water can be obtained from a pair of linear equations associated with force and moment. The right-hand sides are triple integrals of the Havelock source function over the undisturbed underwater profile', evaluated with piecewise-linear hull approximations. Computed for two mathematical hulls and five Series 60 models; 'agreement between theoretical and experimental values considered satisfactory, especially for sinkage and for sufficiently small beam/length ratio.' Formulas NOT transcribed (paywalled).
- **Thin-ship theory and influence of rake and flare**
  - *Authors:* F. Noblesse, G. Delhommeau, H. Y. Kim, C. Yang
  - *Year:* 2009
  - *Where:* J. Eng. Math. 64, 49-80, doi 10.1007/s10665-008-9247-x (abstract only; Semantic Scholar reports openAccessPdf CLOSED)
  - *What it gives:* Gives 'a straightforward method for evaluating the pressure and the wave profile at a ship hull (the wave drag, hydrodynamic lift and pitch moment, and sinkage and trim are also considered) in accordance with Michell's thin-ship theory', using a simple analytical approximation to the local-flow (near-field) part of the Michell-Kelvin Green function; 'the hydrodynamic lift and pitch moment cause the ship hull to sink and trim, with nondimensional sinkage and trim angle given by formulas involving hydrodynamic coefficients'; rake/flare effects 'can be significant, especially at low Froude numbers'. Explicit formulas NOT obtainable.
- **The Summary of the Cooperative Experiment on Wigley Parabolic Model in Japan (17th ITTC Resistance Committee, Varna 1983; Proc. 2nd DTNSRDC Workshop on Ship Wave-Resistance Computations, Nov 1983)**
  - *Authors:* H. Kajitani, H. Miyata, M. Ikehata, H. Tanaka, H. Adachi, M. Namimatsu, S. Ogiwara
  - *Year:* 1983
  - *Where:* DTIC ADP003037; OCR text https://archive.org/stream/DTIC_ADP003037/DTIC_ADP003037_djvu.txt; PDF https://archive.org/download/DTIC_ADP003037/DTIC_ADP003037.pdf (Fig. 4 on p.13 rendered and read)
  - *What it gives:* THE primary Wigley experimental sinkage/trim dataset. Hull y = (B/2)(1-(2x/L)²)(1-(z/D)²), L/B = 10, B/D = 1.6 (IHI 6.0 m: B 0.6, D 0.375; SRI 4.0 m; UT 2.5 m; YNU 2.0 m). Definitions: s = (Δd_F + Δd_A)/(2L) POSITIVE DOWN; t = (d_A − d_F)/L POSITIVE BOW-UP; plotted as σ = 2 k0 L s = 2 s/Fn² and trim in %. Fig. 4: σ for IHI (6 m) and SRI (4 m) free-to-sink-and-trim models is roughly CONSTANT ≈ 0.04-0.05 for 0.1 < Fn < 0.3, rising to ≈ 0.055-0.065 by Fn 0.4 (direct experimental confirmation of s/L ∝ Fn² at low Fn); trim ≈ 0 for Fn < 0.28, slightly bow-down (≈ -0.1%) around Fn 0.32-0.34, then bow-up rising to ≈ +0.75-1.0% at Fn 0.40. Data only as plots (no table).
- **Calculation of ship sinkage and trim using unstructured grids**
  - *Authors:* C. Yang, R. Löhner, F. Noblesse, T. T. Huang
  - *Year:* 2000
  - *Where:* ECCOMAS 2000, Barcelona; http://congress2.cimne.com/eccomas/proceedings/eccomas2000/pdf/367.pdf (read in full; Fig. 4a rendered)
  - *What it gives:* Nonlinear (Euler, free-surface FE) computations of Wigley sinkage & trim free to sink and trim at Fr = 0.177, 0.25, 0.316, 0.374, 0.408, compared with University of Tokyo experiments. Same sign conventions as Kajitani (s positive down, t = (d_A−d_F)/L positive bow-up). Hydrostatic update eqs (8)-(9): ΔH = L/(rho g A_w0), Δα = M/(rho g A_w2) with A_w0 waterplane area, A_w2 its moment of inertia about y. Numbers read from Fig. 4a are listed in wigley_data.
- **Illustrative applications of the Neumann-Michell theory of ship waves**
  - *Authors:* F. Huang, X. Li, F. Noblesse, C. Yang, W. Duan
  - *Year:* 2013
  - *Where:* 28th IWWWFB, http://www.iwwwfb.org/Abstracts/iwwwfb28/iwwwfb28_22.pdf (read in full; Fig. 1 rendered)
  - *What it gives:* LINEAR potential-flow (Neumann-Michell, i.e. Hogner slender-ship + wave correction, no waterline integral) predictions of Wigley sinkage and trim vs F, from two independent codes (GMU, HEU), overlaid on IHHI and SRI experiments (from the 1983 Japanese cooperative experiments). Linear theory tracks the measured sinkage to within ~10% over 0.15<F<0.40 and reproduces the trim sign change (≈0 below F 0.3, slight dip, steep bow-up rise above F≈0.36). Numbers read from the figure are in wigley_data. Also states NM predictions in these plots are for hull in fixed position.
- **Calculation of ship sinkage and trim in deep water using a potential based panel method**
  - *Authors:* M. S. Tarafder & G. M. Khalil
  - *Year:* 2006
  - *Where:* Int. J. Applied Mechanics and Engineering 11(2), 401-414; http://www.ijame.uz.zgora.pl/ijame_files/archives/v11PDF/n2/401-414_Article_14.pdf (read)
  - *What it gives:* Deep-water Rankine-panel (Dawson-type) sinkage/trim with explicit hydrostatic equilibrium equations (5.1)-(5.4): F3 = -rho g s ∫ f_w dx + rho g t ∫ x f_w dx, F5 = -rho g s ∫ x f_w dx ... solved with H0 = rho g ∫ f_w dx, H1 = rho g ∫ x f_w dx, H2 = rho g ∫ x² f_w dx (f_w = waterline width). Sinkage defined positive DOWN at x=0, trim positive BOW-UP about y=0 (their nomenclature list contradicts this, saying 'sinkage positive upward, trim by stern positive' — inconsistent within the paper). Sinkage/trim results are for SERIES 60 only (Figs 5-6), NOT Wigley.
- **Ship squat in water of varying depth**
  - *Authors:* T. P. Gourlay
  - *Year:* 2003
  - *Where:* Int. J. Maritime Eng.; https://www.perthhydro.com/pdf/Gourlay2003VaryingDepth.pdf (read)
  - *What it gives:* States the pioneering constant-depth theories: 'Havelock [1939] for a slender ship in open water of infinite depth', Constantine (narrow channel), Tuck 1966 (shallow open water), Tuck 1967 (channel). Gives the hydrostatic conversions used throughout the Tuck/Gourlay literature: for a fore-aft symmetric ship Z = -rho g A_W s (Z vertical force, s midship sinkage) and M = rho g I_W θ with θ the BOW-UP trim angle.
- **Practical estimation of sinkage and trim for common generic monohull ships**
  - *Authors:* C. Ma, C. Zhang, F. Huang, C. Yang, X. Chen, F. Noblesse
  - *Year:* 2016
  - *Where:* Ocean Engineering 126 (Sept 2016); also ISOPE 2016 paper ISOPE-I-16-021 (abstracts only; ScienceDirect/OnePetro 403)
  - *What it gives:* Two practical approaches: an 'experimental approach' from 22 ship models giving 'particularly simple relations' for sinkage and trim vs Froude number, and a 'numerical approach' using Neumann-Michell linear flow computations for the hull at rest; both reasonable for F ≤ 0.45. The explicit empirical relations could NOT be retrieved.

**Deep water low fn sinkage**

SIGN: every deep-water source found (Havelock 1939; Tuck & Taylor 1970 / Gourlay & Tuck 2001; Yeung 1972 & Noblesse et al. 2009 by description; all Wigley experiments) says the hull SINKS (net hydrodynamic force is DOWNWARD) at low Froude number in deep water. Physical mechanism (Havelock 1939, §1 and §7): at low speed the flow is the double-body streamline flow with an undisturbed (rigid) free surface; the speed-up along the sides lowers the pressure (Bernoulli suction) and the 'defect of vertical pressure' pulls the hull down. Havelock: 'the sinkage is then due to the defect of vertical pressure caused by the fluid motion and should be proportional to the square of the speed.' In the finite-depth decomposition F = F_∞ + F_d (Tuck & Taylor 1970; Gourlay & Tuck 2001 eq 22), F_∞ is exactly this deep-water double-body force; for a slender spheroid F_∞ = rho U² A_W ε² (ln(ε/2) + 3/2 − ε) < 0 (downward) with ε ≈ B/L (eq 23), and 'F_∞ increases in proportion to U²'. Experiments (Graff et al 1964, cited by Gourlay 2011) show 'both the bow and stern normally sink further downwards as the speed increases' at displacement speeds; a RISE of the hull appears only near/above critical depth-Froude number in shallow water (supercritical planing-type lift) — never at low Fn in deep water.

Fn SCALING: F_z / (rho U² L²) tends to a NONZERO NEGATIVE CONSTANT as Fn → 0 (double-body limit), so sinkage/L = −F_z/(rho g A_W L) = C · Fn² with C = O(B/L)²·log(L/B) for slender forms. Havelock's exact half-ellipsoid values give C = gh/U² directly: for L/D = 16 (Wigley-like L/T) C = 0.0138 (B/D=1), 0.0231 (B/D=2), 0.0318 (3), 0.0397 (4); interpolating to the Wigley B/T = 1.6 gives C ≈ 0.019, i.e. s/L ≈ 0.019 Fn² (0.0012 at Fn 0.25, 0.0019 at Fn 0.316), which is within ~10-15% of the measured Wigley sinkage (0.0011-0.0012 at Fn 0.25; 0.0020-0.0023 at Fn 0.30-0.316). The Kajitani 1983 plot is drawn as σ = 2 s/Fn² precisely because it is nearly flat (≈0.04-0.05, i.e. C ≈ 0.020-0.025) over 0.1 < Fn < 0.3; it rises to ≈0.06 by Fn 0.4 as wave effects add to the sinkage. Slender-body asymptotics (Gourlay & Tuck 2001 eq 23) reproduce Havelock's exact spheroid value to ~5% (ε = 1/8: 0.0218 vs exact 0.0231).

TRIM at low Fn, deep water: the double-body moment M_∞ vanishes identically for fore-aft symmetric hulls (Gourlay & Tuck 2001), so for the Wigley the trim is purely wave-induced and is o(Fn²): measured trim is ≈ 0 for Fn < 0.28, slightly BOW-DOWN (t ≈ −0.0005 to −0.0015) for 0.30 < Fn < 0.35, then strongly BOW-UP (t ≈ +0.003 at Fn 0.375, +0.010 at Fn 0.405). Linear NM theory reproduces this qualitative shape including the sign change near Fn ≈ 0.35 (IWWWFB28 Fig 1). For non-symmetric hulls M_∞ ≠ 0 and trim ∝ Fn² at low speed with sign set by hull shape (Tuck's c_θ: 'positive [bow-down in Gourlay's convention] for modern bulk carriers with centre of buoyancy well forward, either sign for containerships').

SHALLOW vs DEEP: Tuck 1966 shallow-water squat scales as Fh²/sqrt(1−Fh²) with Fh = U/sqrt(gh) and is O(∇/L²) with an O(1) coefficient c_s ≈ 1.5 — it is a depth effect (line sources in the horizontal plane with strength ∝ S'(x)/h) and diverges at Fh → 1; the deep-water sinkage is the separate F_∞ term ∝ U² with a much smaller coefficient ∝ (B/L)² ln(L/B). Do not use Tuck's c_s, c_θ formulas for deep water.

**Standard formulas**

1. THIN-SHIP (Michell) — structure per Yeung 1972 / Noblesse et al. 2009 (full formulas paywalled; structure confirmed from abstracts): the linear pressure p = −rho U_∞·∇φ = rho U φ_x (for the crate's free stream −U x̂; p is O(f)) is integrated over the undisturbed centreplane against the hull-surface slope, giving a force and moment O(f²), then sinkage and trim from a 2×2 hydrostatic system. Derived transcription in the crate's conventions (z DOWN, bow at high x, outward normal on y = +f is ∝ (−f_x, 1, −f_z)), FLAGGED AS MY DERIVATION, not copied from a source:
   F_up = −2 rho U ∬_{0<z<T} φ_x(x, 0, z) f_z(x,z) dx dz,   M_bow-up(x_ref) = −2 rho U ∬ (x − x_ref) φ_x f_z dx dz,
   with φ the Michell (Havelock-source) centreplane potential driven by f_x. All other second-order terms (−½ rho |∇φ|² × f_z, and the wetted strip between z = 0 and z = −ζ carrying pressure ~rho g ζ) are O(f³), so to O(f²) F_z is exactly the linear pressure on the mean wetted centreplane. Sanity: at midbody φ_x < 0 (speed-up) and f_z < 0 for a Wigley (beam shrinks with depth), so F_up < 0: the hull sinks. In Fourier–Kochin form (Noblesse's near-field approach: 'space integration over hull panels first, Fourier integration subsequently'; 'the panel integration merely involves an exponential-trigonometric function') this is a wavenumber-space integral of the product of the transform of f_x·e^{−κz} and the transform of f_z·e^{−κz}, the latter obtainable from the crate's Q(kx, κ) of f by parts (boundary terms at z=0 need the waterline half-beam f(x,0), and f(x,T)=0). Yeung's phrase: 'triple integrals of the Havelock source function over the undisturbed underwater profile'.
   Hydrostatic closure (Yang et al. 2000 eqs 8-9; Tarafder & Khalil 2006 eqs 5.1-5.4; Gourlay 2003): for a symmetric hull ΔH = F_up/(rho g A_W) (sinkage = −F_up/(rho g A_W) downward), Δα = M/(rho g I_yy,W); general case solve [H0, −H1; −H1, H2]·(s, t) = (−F3, −F5) with H0 = rho g ∫ B(x) dx, H1 = rho g ∫ x B(x) dx, H2 = rho g ∫ x² B(x) dx (Tarafder eq 5.3-5.4).

2. DEEP-WATER LOW-SPEED (double-body) — Havelock 1939: Q = −rho ∬ (U φ_x + ½|∇φ|²) n_z dS over the immersed half-ellipsoid (his eq 7; body moving +x in fluid at rest, so his U φ_x has the opposite sign convention to the ship-fixed −U free stream), equivalent sinkage h from Q = rho g (π a b) h (eq 11). Exact results eqs 12, 14, 16; Table I of gh/U² listed under sources. Slender-spheroid asymptote (Gourlay & Tuck 2001 eq 23, from Havelock): F_∞ = rho U² A_W ε² (ln(ε/2) + 3/2 − ε), F positive UP, ε = sqrt(12 C_V/π); ⇒ sinkage s = (U²/g) ε² (ln(2/ε) − 3/2 + ε). M_∞ ≡ 0 for fore-aft symmetric hulls.

3. SHALLOW-WATER slender-body (Tuck 1966, as transcribed by Gourlay 2011; x positive AFT from midships, z positive UP, θ positive BOW-DOWN): hull condition φ_y = ±(U/2h) dS/dx on y = 0±; φ = (U/(4π h sqrt(1−Fh²))) ∫_{−L/2}^{L/2} S'(ξ) ln[(x−ξ)² + (1−Fh²) y²] dξ (subcritical); s_LCF = c_s (∇/L²) Fh²/sqrt(1−Fh²), c_s = (L²/(2π ∇ A_WP)) ∫∫ S'(ξ) B(x)/(x−ξ) dξ dx ≈ 1.3-1.5; θ = c_θ (∇/L³) Fh²/sqrt(1−Fh²), c_θ = −(L³/(2π ∇ I_LCF)) ∫∫ S'(ξ) (x−LCF) B(x)/(x−ξ) dξ dx, zero for symmetric hulls. Low-Fh: s_LCF ≈ 1.5 (∇/L²) U²/(gh).

4. FINITE-DEPTH slender-body (Tuck & Taylor 1970; Gourlay & Tuck 2001 eqs 22, 25, 26): F = F_∞ + F_d, M = M_∞ + M_d; F_d = −(rho U²/4π²) ∫_{−∞}^{∞} k² S̄(k) B̄*(k) A(k) dk; M_d = (rho U²/4π²) ∫ k² S̄(k) (xB)̄*(k) A(k) dk; A(k) = −2 ∫_{|k|}^{∞} [1 + q/(Fh² k² h − q tanh(qh))] dq/sqrt(q² − k²); S̄, B̄, (xB)̄ are Fourier transforms of section area, waterline beam and x×beam; F positive up, M positive bow-up (Gourlay 2003: Z = −rho g A_W s, M = rho g I_W θ_bow-up). This is the closest published analogue to the crate's 'transforms integrated against a wavenumber kernel' form, but for finite depth and slender-body (not thin-ship) theory.

5. Sign conventions in the Wigley data: s = (Δd_F + Δd_A)/(2L) positive DOWN (draft increase), t = (d_A − d_F)/L positive BOW-UP (Kajitani 1983; Yang et al. 2000). Kajitani additionally plots σ = 2 k0 L s = 2 s/Fn² and trim in percent.

**Wigley data points**

| Fn | kind | sinkage/L | trim | source |
|---|---|---|---|---|
| 0.15 | experiment | ≈0.0004 (IHHI 6 m / SRI 4 m, free to sink & trim; read from IWWWFB28 Fig.1 replot of 1983 Japanese data) | ≈0.000 (t=(dA−dF)/L) | Kajitani et al. 1983 data as replotted in Huang, Li, Noblesse, Yang, Duan, IWWWFB28 (2013) Fig.1 |
| 0.2 | experiment | ≈0.0007 | ≈ −0.0002 | Kajitani et al. 1983 (IHHI/SRI) via IWWWFB28 Fig.1 |
| 0.25 | experiment | ≈0.0011–0.0012 | ≈ −0.0003 | Kajitani et al. 1983 (IHHI/SRI) via IWWWFB28 Fig.1; UT curve in Yang et al. 2000 Fig.4a gives ≈0.0012, trim ≈0 |
| 0.3 | experiment | ≈0.0020–0.0021 | ≈ −0.0005 (slightly bow-down) | Kajitani et al. 1983 (IHHI/SRI) via IWWWFB28 Fig.1 |
| 0.316 | experiment | ≈0.0023 (UT 2.5 m) | ≈ −0.0003 | UT experiment curve, Yang et al. 2000 Fig.4a |
| 0.33 | experiment | ≈0.0027 | ≈ −0.0015 (most bow-down point) | Kajitani et al. 1983 (IHHI) via IWWWFB28 Fig.1 |
| 0.35 | experiment | ≈0.0033 | ≈ −0.0010 to −0.0005 | Kajitani et al. 1983 (IHHI/SRI) via IWWWFB28 Fig.1 |
| 0.375 | experiment | ≈0.0034–0.0038 | ≈ +0.002 to +0.003 (bow-up) | Kajitani et al. 1983 (IHHI/SRI) via IWWWFB28 Fig.1; UT curve in Yang 2000 |
| 0.405 | experiment | ≈0.0047–0.0048 (SRI) | ≈ +0.010 (bow-up); Kajitani Fig.4 shows ≈ +0.75–1.0 % | Kajitani et al. 1983 (SRI) via IWWWFB28 Fig.1 and Yang 2000 Fig.4a |
| 0.2 | experiment | σ = 2 s/Fn² ≈ 0.04–0.05 ⇒ s/L ≈ 0.0008–0.0010 (approximate, read from scanned plot) | ≈ 0 % | Kajitani et al. 1983 Fig.4, IHI 6 m and SRI 4 m models free/free (σ roughly constant 0.1<Fn<0.3) |
| 0.15 | linear-theory | ≈0.0004 | ≈0.000 | Neumann–Michell linear theory (GMU & HEU codes), IWWWFB28 Fig.1; hull in fixed position |
| 0.2 | linear-theory | ≈0.0008 | ≈0.000 | Neumann–Michell linear theory, IWWWFB28 Fig.1 |
| 0.25 | linear-theory | ≈0.0012 | ≈0.000 to +0.0002 | Neumann–Michell linear theory, IWWWFB28 Fig.1 |
| 0.3 | linear-theory | ≈0.0020 | ≈ +0.0002 (small bump), dipping to ≈ −0.0005 at Fn≈0.32 | Neumann–Michell linear theory, IWWWFB28 Fig.1 |
| 0.35 | linear-theory | ≈0.0030 | ≈0.000 | Neumann–Michell linear theory, IWWWFB28 Fig.1 |
| 0.375 | linear-theory | ≈0.0035 | ≈ +0.003 | Neumann–Michell linear theory, IWWWFB28 Fig.1 |
| 0.4 | linear-theory | ≈0.0046 | ≈ +0.0095 | Neumann–Michell linear theory, IWWWFB28 Fig.1 |
| 0.45 | linear-theory | ≈0.0056 | ≈ +0.015 | Neumann–Michell linear theory, IWWWFB28 Fig.1 |
| 0.177 | nonlinear-numerics | 0.0010 (exp UT ≈0.0007) | +0.0002 (exp ≈0) | Yang, Löhner, Noblesse, Huang, ECCOMAS 2000, Fig.4a (Euler + nonlinear free surface, free to sink and trim) |
| 0.25 | nonlinear-numerics | 0.0017 (exp UT ≈0.0012) | +0.0004 (exp ≈0) | Yang et al. 2000 Fig.4a |
| 0.316 | nonlinear-numerics | 0.0026 (exp UT ≈0.0023) | +0.0008 (exp ≈ −0.0003) | Yang et al. 2000 Fig.4a |
| 0.374 | nonlinear-numerics | 0.0040 (exp UT ≈0.0034) | +0.0060 (exp ≈ +0.002) | Yang et al. 2000 Fig.4a |
| 0.408 | nonlinear-numerics | 0.0047 (exp UT ≈0.0047) | +0.0123 (exp ≈ +0.011) | Yang et al. 2000 Fig.4a |
| 0.25 | linear-theory | Havelock 1939 double-body half-ellipsoid, L/D=16, B/D=1.6 (interpolated gh/U² ≈ 0.019): s/L ≈ 0.019 Fn² = 0.0012; analytic zero-Fn limit, not a Wigley computation | 0 (symmetric double body) | Havelock 1939 Table I (interpolated between B/D=1: 0.0138 and B/D=2: 0.0231) |

**Other validation anchors**

(a) ZERO-FROUDE (double-body) LIMIT: as Fn → 0 the Michell Green function tends to the rigid-wall (positive image) Rankine source, so the crate's F_z/(rho U²) must tend to a finite negative constant (hull sinks) and M → 0 for any fore-aft symmetric hull (Gourlay & Tuck 2001: M_∞ ≡ 0). Compare the constant with Havelock 1939 Table I for a half-ellipsoid of the same L, B, T: gh/U² = 0.0253/0.0453/0.0612/0.0735 (L/D = 10, B/D = 1,2,3,4) and 0.0138/0.0231/0.0318/0.0397 (L/D = 16); equivalently F_up = −rho g A_W h with A_W = π L B/4. Expect agreement to O(B/L) (thin-ship vs exact ellipsoid) — the Havelock number is exact potential flow for that body, not thin-ship, so a 10-20% gap at B/L = 0.1 is acceptable; the asymptotic form ε²(ln(2/ε) − 3/2 + ε) (Gourlay & Tuck eq 23) is itself 5% below Havelock's exact value at ε = 1/8. (b) Fn² SCALING CHECK: s/L·Fn⁻² should be flat as Fn → 0 (Kajitani's σ = 2 s/Fn² is experimentally ≈0.04-0.05 for 0.1 < Fn < 0.3); any Fn⁴ or log-Fn behaviour at low Fn indicates a missing double-body term. (c) SYMMETRY: for f(x) even about midships, M(x_ref = midships) must vanish in the Fn → 0 limit and be purely wave-generated; reflecting the hull (x → −x) must flip the sign of M and leave F_z unchanged. (d) TRIM SIGN CHANGE for the Wigley: trim ≈ 0 (slightly bow-down, t ≈ −0.0005 to −0.0015) for 0.30 < Fn < 0.35 and bow-up for Fn > 0.36, reaching t ≈ +0.010 at Fn 0.405 — linear NM theory reproduces this (IWWWFB28 Fig 1), so a thin-ship code should too. (e) SLENDER-BODY CROSS-CHECK for non-symmetric hulls (deep water, low Fn): trim direction determined by double-body M_∞; for a hull with LCB forward Tuck's c_θ > 0 (bow-down in Gourlay's convention) — shallow-water only, use qualitatively. (f) HYDROSTATIC CLOSURE: s = −F_up/(rho g A_W), θ = M/(rho g I_L) for a symmetric waterplane; general 2×2 system with H0, H1, H2 (Tarafder & Khalil 2006 eqs 5.3-5.4; Yang et al. 2000 eqs 8-9). (g) Series 60 (Cb = 0.6) has parallel experimental sinkage/trim data (IHHI, SRS, UT) plotted in Yang et al. 2000 Fig 4b: s/L ≈ 0.0008 (Fr 0.18), 0.0018 (0.25), 0.0032 (0.32), 0.0048 (0.368), 0.0055 (0.388); trim ≈ −0.0005 (0.18-0.32), +0.0065 (0.368), +0.012 (0.388) — nonlinear computations; useful only as a second, non-symmetric anchor.

**Notes**

NOT FOUND / NOT ACCESSIBLE: (1) Full text of Yeung 1972 (JSR 16:47-59) — only the abstract; the explicit Michell force/moment integrals could not be transcribed. (2) Noblesse, Delhommeau, Kim, Yang 2009 'Thin-ship theory and influence of rake and flare' (J Eng Math 64:49-80) — paywalled (Springer/ResearchGate 403; Semantic Scholar: no OA PDF); this is the paper most likely to contain the explicit Michell-theory hydrodynamic lift and pitch-moment formulas with Fourier-Kochin near-field evaluation, and Wigley thin-ship sinkage/trim curves. (3) Ma et al. 2016 (Ocean Eng.) empirical sinkage/trim relations and NM Wigley numbers — 403. (4) Tuck 1964 JSR deep-water slender-ship paper — not accessible; Gourlay & Tuck 2001 confirm the deep-water double-body force is 'a well known but difficult problem' handled via Havelock's spheroid. (5) Doctors & Day: no paper specifically on thin-ship sinkage/trim located; their transom virtual-appendage work (Doctors & Day 2000) is resistance-only. (6) Wehausen & Laitone 1960 — no sinkage/trim section located online. (7) Lazauskas 2009 Adelaide thesis (reportedly compares thin-ship 'squat' predictions with DTMB 5415 data) — repository not reachable. (8) Journée 1992 Delft report 909 (four Wigley hullforms) — repository 403. (9) Kajitani 1983 gives sinkage/trim only as Fig. 4 plots (no table); the numbers above are read from a scanned figure (σ scale ±0.01) and from the clean vector replots in IWWWFB28 Fig 1 and Yang 2000 Fig 4a (±0.0001 in s/L, ±0.0005 in t). (10) Tarafder & Khalil 2006 turned out to give Series 60, not Wigley, sinkage/trim. (11) Chen & Noblesse 1983 JSR (Wigley wave-resistance comparison 'corrected for sinkage and trim') — abstract only.

CAVEATS: The F_up / M formulas in item 1 of standard_formulas are my derivation in the crate's conventions from the thin-ship pressure p = rho U φ_x, consistent with Yeung's description; verify the factor 2 (two sides) and the sign of f_z (z down ⇒ f_z < 0 for a Wigley) before implementing. Tarafder & Khalil's nomenclature list contradicts their §5 text on the sign of s and t; I used §5 (s down, t bow-up), which also matches Kajitani/Yang. Gourlay 2011 uses θ positive BOW-DOWN and x positive AFT for Tuck's shallow-water formulas, opposite to the crate's bow-up/bow-at-high-x convention. Havelock's ellipsoid has vertical sides at the waterline and elliptic sections, so its gh/U² is only an order-of-magnitude anchor for the parabolic-section Wigley (though numerically it lands within ~15% of experiment at Fn 0.25).

LOCAL FILES (scratchpad, absolute paths): /private/tmp/claude-501/-Users-avi-michell/2fdb50af-0291-45b8-bdde-e7c422fd969b/scratchpad/havelock.pdf (collected papers; 1939 note at pdf pages 472-475), gourlay2001.pdf, gourlay2001_eq23.pdf.png, kajitani.pdf, kaj_p9.pdf.png (Fig. 4 sinkage/trim), iwwwfb28_22.pdf, iwwwfb_p3.pdf.png (Fig. 1 Wigley NM vs exp), eccomas_p11.pdf.png (Yang 2000 Fig. 4a), webfetch-1789574174409-1rriek.txt (Gourlay 2011 text), tarafder.pdf.

### Literature set 2

**Sources**

- **Sinkage and Trim in First-Order Thin-Ship Theory**
  - *Authors:* R. W. Yeung
  - *Year:* 1972
  - *Where:* J. Ship Research 16(1):47-59, doi 10.5957/jsr.1972.16.1.47 (abstract only reachable via TRID; full text paywalled)
  - *What it gives:* THE canonical thin-ship sinkage/trim paper. Abstract: sinkage and trim follow from a pair of linear equations (vertical force + moment equilibrium) whose right-hand sides are 'triple integrals of the Havelock source function over the undisturbed underwater profile' (centreplane). Hull approximated piecewise-linearly so the integrals reduce to two closed-form types. Computed for two mathematical hulls and five Series 60 models; agreement with experiment 'satisfactory', better for sinkage than trim and for modest B/L. Formulas themselves NOT retrievable here.
- **The Maximum Sinkage of a Ship**
  - *Authors:* T. P. Gourlay & E. O. Tuck
  - *Year:* 2001
  - *Where:* J. Ship Research 45(1):50-58; PDF cmst.curtin.edu.au/wp-content/uploads/sites/4/2016/05/gourlay-2001-the_maximum_sinkage_of_a_ship.pdf (text extracted to scratchpad gourlay_tuck2001_clean.txt)
  - *What it gives:* (i) Tuck & Taylor (1970) decomposition F = F_inf + F_d, M = M_inf + M_d, where F_inf, M_inf are the force/moment on the lower half of the equivalent DOUBLE BODY in unbounded fluid (the deep-water, zero-Fn limit) and F_d, M_d are finite-depth corrections. (ii) Havelock (1939) closed form used for F_inf on an equivalent spheroid: F_inf = rho U^2 A_W eps^2 [ln(eps/2) + 3/2 - eps] (decoded from garbled PDF text; bracket sign robust: negative, i.e. downward, for slender eps), eps = sqrt(12 C_vol/pi), C_vol = Vol/L^3, A_W waterplane area; M_inf = 0 for fore-aft symmetric hulls. (iii) Shallow-water Tuck 1966 form F = (rho U^2 /(2 pi h sqrt(1-Fh^2))) double-int B'(x) S'(xi) ln|x-xi| dx dxi, 'usually negative, i.e. downward'; sinkage coefficient C_S almost universal. (iv) Fourier-space forms: F ∝ ∫ k S(k) B*(k) ... dk, M ∝ ∫ ... S(k) xB(k) dk, i.e. force = wavenumber integral of (transform of section-area) x conj(transform of beam or x*beam) x kernel.
- **Slender-Body Methods for Predicting Ship Squat**
  - *Authors:* T. P. Gourlay
  - *Year:* 2008
  - *Where:* Ocean Engineering 35(2):191-200; PDF cmst.curtin.edu.au (scratchpad gourlay2008slender.txt)
  - *What it gives:* General Fourier-transform statement (eq. 6/29): upward vertical force on a ship held at static draft/trim  Z = (rho U^2 /(2 pi h sqrt(1-Fh^2))) * i ∫ k S'(k) B*(k) K(k) dk  (shallow water; S'(k) = FT of dS/dx, B(k) = FT of waterline beam B(x), * = conjugate, K depends on transverse geometry); 'the bow-down trim moment is found by replacing B(k) by xB(k)'; 'steady sinkage and trim then follow hydrostatically'. Written for transom sterns using dS/dx rather than S (transom modelled as an infinitely long cylinder downstream). Pressure = excess above hydrostatic from linearised Bernoulli. Non-singular integrand is highlighted as the computational advantage - the same structure the crate wants (products of exact transforms integrated against a kernel).
- **Sinkage and Trim of a Fast Displacement Catamaran in Shallow Water**
  - *Authors:* T. P. Gourlay
  - *Year:* 2008
  - *Where:* J. Ship Research 52(3):175-183; PDF perthhydro.com/pdf/Gourlay2008Catamaran.pdf (scratchpad gourlay2008cat.txt)
  - *What it gives:* Catamaran sinkage/trim by LINEAR SUPERPOSITION of the two demihull slender-body potentials (no Kutta condition at sterns). NPL round-bilge demihull stretched to L/B = 14.0, B/T = 1.5 (cB 0.397, cP 0.693, cM 0.573, LCB 6.4% L aft, LCF 8.4% L aft), catamaran s/L = 0.2-0.5, h/L = 0.10. Max (transcritical) sinkage coefficients Table 3: NPL demihull C_max_mid 0.31, theta_max 2.29, C_max_stern 1.30; Wigley 0.56 / 2.98 / 1.95; Taylor A3 0.56/3.18/2.03. Notes deep-water cross-flow under demihulls of 5-7% U (Miyazawa 1979) and Insel 1990 breaking waves between hulls in deep water. Results are plots only (Fig 2-6); shallow water, not deep.
- **A brief history of mathematical ship-squat prediction, focussing on the contributions of E.O. Tuck**
  - *Authors:* T. P. Gourlay
  - *Year:* 2011
  - *Where:* J. Eng. Math. 70:5-16; PDF perthhydro.com/pdf/Gourlay2011HistorySquatTuck.pdf (scratchpad gourlay2011.txt)
  - *What it gives:* Tuck 1963 thesis / 1964 JSR deep-water slender-ship matched asymptotics as the parent of the 1966 shallow-water squat theory. Shallow-water results: s = c_s (Vol/L^2) Fh^2/sqrt(1-Fh^2), c_s = (L^2/(2 pi A_WP Vol)) ∬ S'(xi) B(x)/(x-xi)... ≈ 1.3-1.5 (nearly universal, 1.5 recommended); trim theta = c_theta (Vol/L^3) Fh^2/sqrt(1-Fh^2), c_theta zero for fore-aft symmetric hulls, positive (bow-down) for bulk carriers with LCB forward, either sign for containerships. 'For a fore-aft symmetric hull, sinkage is non-zero and trim is zero at subcritical speeds.' Cites Havelock 1939 ZAMM 19:458-461 'Note on the sinkage of a ship at low speeds' (ellipsoid, potential flow, decimetre-order sinkage) as the first deep-water sinkage calculation.
- **Calculation of Ship Sinkage and Trim Using Unstructured Grids**
  - *Authors:* C. Yang, R. Lohner, F. Noblesse, T. T. Huang
  - *Year:* 2000
  - *Where:* ECCOMAS 2000, Barcelona; PDF congress2.cimne.com/eccomas/proceedings/eccomas2000/pdf/367.pdf (scratchpad eccomas2000.txt)
  - *What it gives:* Sign conventions matching the crate's need: s = (dDF + dDA)/Lpp positive DOWNWARD, t = (dA - dF)/Lpp positive BOW-UP; sinkage/trim corrections dH = L/(rho g A_w0), dalpha = M/(rho g A_w2) with A_w2 waterplane inertia. Wigley (B/L 0.1, D/L 0.0625) at Fr 0.177, 0.25, 0.316, 0.374, 0.408 vs Univ. Tokyo experiments; Series 60 CB 0.6 at Fr 0.18, 0.25, 0.32, 0.3682, 0.388 vs IHHI, SRS(Korea), UT experiments. Plots only: Wigley s in 0-0.006, t in -0.002-0.014; S60 s in 0-0.007, t in -0.002-0.014 (Euler/FE, nonlinear-numerics).
- **Sinkage, trim, drag of a common freely floating monohull ship (and: Practical estimation of sinkage and trim for common generic monohull ships, Ocean Eng. 126:203-216, 2016)**
  - *Authors:* C. Ma, Y. Zhu, H. Wu, W. Li, H. Fu, F. Noblesse
  - *Year:* 2017 (2016)
  - *Where:* MARINE 2017 (UPC upcommons PDF; scratchpad noblesse_sinkage.txt); Ocean Engineering 126 (2016)
  - *What it gives:* Linear (Neumann-Michell) lift/moment: (Cz, Czx) = ∫_Sigma (n_z, n_x z - n_z x) p da, p from linear Bernoulli; equilibrium Hm/L ≈ F^2 (Cz + eps2 Czx)/(a0(1-eps0 eps2)), 2Htau/L ≈ F^2 (Czx + eps0 Cz)/(a2(1-eps0 eps2)), a_k = ∫_W x^k dxdy / L^(k+2), eps0 = a1/a0, eps2 = a1/a2. Positive Hm = downward, positive tau = bow-up. Key finding: sinkage/trim computed from pressure on the AT-REST hull are adequate for F <= 0.45. EMPIRICAL fit to 22 monohull models: Hm ≈ 0.9 sqrt(B D)(Cb - 0.13) F^2 (20-30% accuracy), Hs (stern sinkage) ≈ 0.025 sqrt(BD) F*^2 sqrt(1+F*^8), F* = F/0.33, Hb = 2Hm - Hs, Htau = Hs - Hm = L tau_deg pi/360. 'Midship sinkage increases approximately like F^2 as F <= 0.45'. Applied to Wigley (L 2.5 m), S60 (4 m), DTMB 5415 (5.72 m); plots only.
- **Nonlinear corrections of linear potential-flow theory of ship waves**
  - *Authors:* C. Ma, Y. Zhu, J. He, C. Zhang, D. Wan, C. Yang, F. Noblesse
  - *Year:* 2018
  - *Where:* Eur. J. Mech. B/Fluids 67:1-14 (scratchpad ejmb2018.txt)
  - *What it gives:* Same lift/moment/equilibrium relations (4a-c) with Cz = Fz/(rho L^2 V^2), Czx = Mzx/(rho L^3 V^2); trim tau_deg ≈ (9/pi) sqrt(BD)/L [F*^2 sqrt(1+F*^8) - 36 (Cb-0.13) F^2]. Nonlinear effects on sinkage and trim 'relatively small'; the -|grad phi|^2/2 Bernoulli term gives a small INCREASE of sinkage (bottom n_z ≈ -1) but little change in trim. Wigley, S60, DTMB 5415, KCS results as plots.
- **A Three-Dimensional Linear Analysis of Steady Ship Motion in Deep Water (PhD thesis)**
  - *Authors:* J. J. M. Baar
  - *Year:* 1986
  - *Where:* Brunel University; PDF saved in tool-results (scratchpad thesis3d.txt, 184 pp.)
  - *What it gives:* Neumann-Kelvin (hull-surface source) formulation. Eq. 2.21a-c: d_w, l_w, m_w = -∬_h (phi_x - |grad phi|^2/2 ...)(n_x, n_z, z n_x - x n_z) da (nondimensional by rho V^2 L^2, rho V^2 L^3). Static equilibrium 2.22b-c (Wehausen 1969, Yeung 1972, Gadd 1973, Noblesse & Dagan 1976): i0 s_w - i1 theta_w + Fn^2 l_w = 0, i1 s_w - i2 theta_w + Fn^2 m_w = 0 with i_k = ∫ x^k b(x) dx / L^(k+2); s positive for increasing draft, theta positive bow-up; hull surface taken at the at-rest position (small sinkage assumption); trimming moment from the resistance neglected. Friesland-class destroyer (transom): measured sectional lift is DOWNWARD and roughly uniform along the hull for Fn <= 0.35, peaks aft for Fn > 0.35; total lift peaks at Fn 0.45; trim small up to Fn 0.35; linear theory good for Fn < 0.35.
- **Calculation of ship sinkage and trim in deep water using a potential based panel method**
  - *Authors:* M. S. Tarafder & G. M. Khalil
  - *Year:* 2006
  - *Where:* Int. J. Applied Mechanics and Engineering 11(2):401-414 (scratchpad tarafder.txt)
  - *What it gives:* Series 60 in deep water (Morino panel + Dawson operator). Defines s = downward displacement at x=0, t = bow-up rotation; equilibrium (5.3-5.4): -(F3^s + H0) s + (F3^t + H1) t = -F3^0, (F5^s + H1) s - (F5^t - H2) t = -F5^0 with H_k = rho g ∫ x^k f_w(x) dx (f_w = waterplane width). Compares with IHI free-to-trim experiments (Takeshi et al. 1987): Fig 5 sinkage (s/L)*100 on a 0-1.0 axis, Fig 6 trim 100*t on a 0-3 axis, Fn 0-0.5; 'agreement quite satisfactory'. Cites Doctors & Day (2000) IWWWFB and Suzuki 1979, Yasukawa 1993, Bessho & Sakuma 1992 as earlier sinkage/trim work.
- **The Influence of a Bottom Mud Layer on the Steady-State Hydrodynamics of Marine Vehicles**
  - *Authors:* L. J. Doctors, G. Zilman, T. Miloh
  - *Year:* 1996
  - *Where:* 21st Symposium on Naval Hydrodynamics, nationalacademies.org/read/5870/chapter/50
  - *What it gives:* Doctors' thin-ship generalized-force convention: x forward (bow at +x, as in the crate), y to port, z UP; pressure p = rho U phi_x (linear Bernoulli, hydrostatic part dropped); wave resistance R_W, 'sinkage force' S_W, 'bow-up moment' M_W obtained by integrating p against a 'generalized hull slope' (eq. 27/29, images only). Coefficients C_S = S_W/(1/2 rho U^2 S), C_M = M_W/(1/2 rho U^2 S L) with S wetted area (inferred: for the standard Wigley S/L^2 = 0.14879, A_W = 2BL/3, I_W = BL^3/30 these give exactly the quoted s/L = 1.116 C_S F^2 and t = 22.32 C_M F^2). Finite depth + mud; deep-water numbers not given.
- **The squat of a vessel with a transom stern (IWWWFB-15, 2000); Nonlinear free-surface effects on the resistance and squat of high-speed vessels with a transom stern (24th Symp. Naval Hydrodynamics, 2003); Resistance prediction for transom-stern vessels (FAST'97)**
  - *Authors:* L. J. Doctors & A. H. Day
  - *Year:* 1997-2003
  - *Where:* iwwwfb.org (site unreachable during this search); nationalacademies.org/read/10834/chapter/35 (page images only)
  - *What it gives:* Thin-ship (Michell-type) computation of resistance, sinkage and trim for transom-stern hulls with the virtual-appendage closure the crate already uses; Tarafder describes it as 'an inviscid linearized near-field solution within classical thin-ship theory'. Formulas and numbers NOT retrievable in this session.
- **An investigation into the resistance components of high speed displacement catamarans**
  - *Authors:* M. Insel & A. F. Molland
  - *Year:* 1992
  - *Where:* Trans. RINA 134 (full text not reachable; pdfcoffee/undip mirrors did not serve the paper)
  - *What it gives:* Experiments on a Wigley hull and NPL-derived round-bilge C3/C4/C5 demihulls, Fn 0.2-1.0, s/L 0.2,0.3,0.4,0.5,inf: total resistance, running trim, sinkage, wave-cut analysis; thin-ship theory used for wave resistance/interference (not for sinkage/trim). Sinkage and trim 'change significantly for Fn >= 0.35'. Numbers not retrieved.
- **Resistance experiments on a systematic series of high speed displacement catamaran forms (Ship Science Report 71); Practical evaluation of high-speed round bilge catamaran resistance**
  - *Authors:* A. F. Molland, J. F. Wellicome, P. R. Couser
  - *Year:* 1994
  - *Where:* eprints.soton.ac.uk/46409 (HTTP 403 here)
  - *What it gives:* Tabulated running trim and sinkage vs Fn 0.1-1.0 for NPL-series demihulls (L/B 7-15.1, B/T 1.5-2.5) alone and as catamarans - the best public deep-water catamaran sinkage/trim dataset, but not retrievable in this session.
- **Optimum hull spacing of a family of multihulls**
  - *Authors:* E. O. Tuck & L. Lazauskas
  - *Year:* 1998
  - *Where:* Ship Technology Research 45
  - *What it gives:* Thin-ship (Michell) wave-resistance optimisation of Wigley-demihull multihulls; NO sinkage/trim. Lazauskas states explicitly on boatdesign.net (thread 47953) that 'Michlet cannot calculate sinkage and trim, or any near-field effects'; it can only take a prescribed sunk/trimmed hull.
- **Trim and Sinkage Effects on Wave Resistance with Series 60, CB = 0.60**
  - *Authors:* H. C. Kim & D. Jenkins
  - *Year:* 1981
  - *Where:* DTNSRDC report, DTIC ADA105972 (HTTP 403 on both DTIC paths)
  - *What it gives:* Measured sinkage and trim vs Fn for Series 60 CB 0.60 fixed vs free; the primary source for the S60 sinkage/trim curves used by Yeung-type comparisons. Not retrievable here.
- **Thin-ship theory and influence of rake and flare**
  - *Authors:* F. Noblesse, G. Delhommeau, H. Y. Kim, C. Yang
  - *Year:* 2009
  - *Where:* J. Eng. Math. 64:49-80 (abstract only)
  - *What it gives:* Abstract states a 'straightforward method for evaluating the pressure and the wave profile at a ship hull (the wave drag, hydrodynamic lift and pitch moment, and sinkage and trim are also considered) in accordance with Michell's thin-ship theory', with a practical Green-function approximation; rake/flare effects 'significant especially at low Froude numbers'. Formulas paywalled.

**Deep water low fn sinkage**

Consistent picture from every source reached: (1) SIGN: in deep water at low Fn the linear hydrodynamic vertical force is DOWNWARD (suction), so the hull SINKS (positive sinkage = increase of draft) - Havelock 1939 (spheroid, potential flow), Tuck & Taylor 1970 / Gourlay & Tuck 2001 ('F_inf ... on the lower half of an equivalent double body' with F_inf = rho U^2 A_W eps^2 [ln(eps/2)+3/2-eps] < 0 for slender eps), Baar 1986 measured sectional lift on a destroyer 'downward' and nearly uniform along the hull for Fn <= 0.35, Ma/Noblesse 2016 empirical Hm ≈ 0.9 sqrt(BD)(Cb-0.13)F^2 > 0 (downward) for all 22 models. (2) Fn-SCALING: the force is proportional to rho U^2 x (hull-shape integral) in the zero-Fn (double-body / rigid free surface) limit, hence sinkage/L ∝ Fn^2 with an almost speed-independent coefficient at low Fn; Ma et al. find 'midship sinkage increases approximately like F^2 as F <= 0.45'; the Doctors form s/L = const x C_S Fn^2 and Baar's i0 s - i1 theta = -Fn^2 l_w make the Fn^2 explicit. Wave (Michell-type, e^{-nu ...}) contributions modulate this above Fn ~0.25-0.3 and make sinkage grow faster than Fn^2 toward the hump (Fn 0.4-0.5), where linear theory degrades (Baar: poor above 0.35 for a transom destroyer; Ma: adequate to 0.45 for common monohulls). (3) TRIM at low Fn: the double-body moment vanishes for fore-aft symmetric hulls (Gourlay & Tuck 2001: M_inf = 0; Gourlay 2011: 'for a fore-aft symmetric hull sinkage is non-zero and trim is zero'), so low-Fn trim is set by fore-aft asymmetry and is small; typical hulls (Series 60, Wigley-like with LCB near midships) show slight BOW-DOWN trim at low Fn turning bow-up around Fn 0.30-0.37 (Ma et al. formula: Htau = sqrt(BD)F^2[0.2296 sqrt(1+F*^8) - 0.9(Cb-0.13)], negative (bow-down) at low F when Cb > 0.385, crossing zero near F ≈ 0.37 for Cb = 0.6). Order of magnitude check for the crate's Wigley (B/L 0.1, D/L 0.0625, Cb 0.444): Havelock/Tuck-Taylor equivalent-spheroid double-body estimate gives eps = sqrt(12 C_vol/pi) = 0.103, s/L ≈ 0.0166 Fn^2; Ma et al. empirical fit gives s/L ≈ 0.0224 Fn^2 - same sign and order, the ~30% gap being the free-surface (wave + waterline) part and non-spheroidal shape.

**Standard formulas**

All linear theories reached use the same structure, differing only in the pressure model: (A) Pressure: linearised Bernoulli p = rho U phi_x (Doctors; x forward, free stream -U along x in the ship frame -> in the crate's frame p = -rho * (perturbation velocity along the free stream) * U; Gourlay: 'excess above hydrostatic'; Baar/Noblesse keep the -rho|grad phi|^2/2 term as a second-order correction that slightly increases sinkage). (B) Force/moment by integrating p over the hull with the outward normal: F_z(up) = -∬ p n_z dS, M = -∬ p (x n_z - z n_x) dS ... [Noblesse (4c): (Cz, Czx) = ∫ (n_z, n_x z - n_z x) p da, z up; Baar (2.21b,c): l_w = -∬(phi_x - |grad phi|^2/2) n_z da, m_w = -∬(...)(z n_x - x n_z) da; all normalise by rho U^2 L^2 and rho U^2 L^3]. For a THIN ship with y = ±f(x,z), z downward, n_z dS -> ∓f_z dx dz per side, so to second order in f (same order as Michell's R_w = 2∬ p1 f_x dxdz): F_up = -2 ∬_centreplane p1 (∂f/∂z) dx dz [+ 2∫ p1(x,T) f(x,T) dx if the hull has a flat bottom where f(x,T) != 0], M_bow-up = -2 ∬ p1 (x - x_ref)(∂f/∂z) dx dz, with the z*f_x (resistance) lever-arm term negligible for slender hulls (Baar, Wehausen 1973, Noblesse & Dagan 1976) and the waterline-strip (wave-elevation) term O(f^3). Yeung 1972 is exactly this with p1 from the Havelock (Kelvin) source distribution sigma ∝ U f_x on the centreplane ('triple integrals of the Havelock source function over the undisturbed profile'). (C) Equilibrium (Baar 2.22b-c, Yang et al. 2000 eq 8-9, Tarafder 5.1-5.4, Noblesse 3a/4a): rho g (A_W s - I_1 theta) = -F_up ... i.e. i0 s - i1 theta + Fn^2 l = 0, i1 s - i2 theta + Fn^2 m = 0, i_k = ∫ x^k b(x) dx / L^(k+2) about the chosen reference; s positive downward, theta positive bow-up; decoupled if x_ref = LCF. Doctors' coefficient form for the standard Wigley: s/L = 1.116 C_S Fn^2, trim = 22.32 C_M Fn^2 with C_S = S_W/(1/2 rho U^2 S), C_M = M_W/(1/2 rho U^2 S L), S = wetted area (0.14879 L^2). (D) Wavenumber-space form (Tuck 1966 as restated by Gourlay 2008; Gourlay & Tuck 2001 eq 17-20): Z(up) = (rho U^2/(2 pi h sqrt(1-Fh^2))) i ∫ k S'(k) B*(k) K(k) dk, moment with B(k) -> xB(k): the force is an integral over wavenumber of [transform of the source strength (dS/dx or, in thin-ship, f_x)] x conj[transform of the pressure-weighting hull function (beam B or, in thin-ship, f_z)] x a real kernel - i.e. exactly the product-of-transforms-against-kernel structure the crate prefers. In deep water the same construction with the Havelock Green function gives F_up as (a) a Michell-like wave part ∝ ∫ dtheta sec^k(theta) Re[ (I+iJ)_{f_x}(theta) conj((I+iJ)_{f_z}(theta)) ] over the Kelvin fan (both amplitudes use e^{-nu lambda^2 z} e^{i nu lambda x}, available now) plus (b) a local/double-body part given by a principal-value integral over k of Q_{f_x}(k cos theta, k) conj(Q_{f_z}(k cos theta, k)) times k/(k - nu sec^2 theta); the f_z transform follows from Q_f by parts: ∬ f_z e^{-kappa z} e^{-i kx x} = kappa Q_f(kx,kappa) - ∫ f(x,0) e^{-i kx x} dx (+ e^{-kappa T}∫ f(x,T) e^{-i kx x} dx if the bottom is open). Nobody reached states this deep-water spectral form explicitly for f_z; Yeung 1972 (paywalled) is the reference to check it against.

**Wigley data points**

| Fn | kind | sinkage/L | trim | source |
|---|---|---|---|---|
| 0.177 | nonlinear-numerics | ~0.001 (lowest of five plotted points; axis 0-0.006, values not printed) | ~0 (axis -0.002 to 0.014, positive bow-up) | Yang, Lohner, Noblesse, Huang 2000 (ECCOMAS) Fig 4a, Euler/FE free to sink & trim vs Univ. Tokyo experiment; Wigley B/L 0.1, D/L 0.0625 |
| 0.408 | nonlinear-numerics | largest of the five points, plotted axis tops at 0.006 (so <= 0.006; ~0.004-0.005) | largest, plotted axis tops at 0.014 (bow-up) | Yang et al. 2000 Fig 4a (plot only; Fr set 0.177, 0.25, 0.316, 0.374, 0.408); experiments UT |
| 0.25 | unknown | 0.00140 (midship, downward) from the Ma/Noblesse 2016 EMPIRICAL 22-model fit Hm = 0.9 sqrt(BD)(Cb-0.13)F^2 with B/L 0.1, D/L 0.0625, Cb 0.444 | Htau/L = (Hs-Hm)/L ≈ -0.0002 (slightly bow-down; tau ≈ -0.02 deg) from the same fit | Ma, Zhang, Chen, Jiang, Noblesse 2016 Ocean Eng. eq (6) - empirical, not a measurement; verify against colleague's 1983 Wigley data |
| 0.3 | unknown | 0.00201 (empirical fit as above); Havelock/Tuck-Taylor double-body spheroid estimate gives 0.0166 F^2 = 0.0015 | ≈ -0.0001 L trim-sinkage (near zero; crossing to bow-up ~Fn 0.32-0.35 per the fit) | Ma et al. 2016 empirical fit; Gourlay & Tuck 2001 eq 23-24 for the double-body estimate |
| 0.408 | unknown | 0.00372 (empirical fit) | Htau/L ≈ +0.0032 (bow-up, tau ≈ +0.37 deg) from the fit | Ma et al. 2016 empirical fit (F* = F/0.33 term dominates) |

**Other validation anchors**

SERIES 60 CB 0.60: measured free-to-trim sinkage/trim exist in Kim & Jenkins 1981 (DTNSRDC, DTIC ADA105972), IHHI/Takeshi et al. 1987, SRS Korea, Univ. Tokyo; plotted (not tabulated) in Yang et al. 2000 Fig 4b (Fr 0.18-0.388: s from ~0.001 to ~0.006 L, t from ~0 slightly negative to ~0.012 bow-up), Tarafder & Khalil 2006 Fig 5-6 (s/L up to ~0.5-1% by Fn 0.5, trim 100*t up to ~3), Ma/Noblesse 2017-2018 (S60 L = 4 m). Ma et al. empirical fit with B/L 0.1333, D/L 0.0533, Cb 0.6: Hm/L ≈ 0.0357 F^2 (0.00223 at Fn 0.25, 0.0032 at Fn 0.30), bow-down trim at low Fn crossing to bow-up at Fn ≈ 0.37. Yeung 1972 computed five S60 models with linear thin-ship theory and called sinkage agreement satisfactory, trim less so. FRIESLAND destroyer (Baar 1986 / Andrew 1985): measured sectional vertical force downward and nearly uniform for Fn <= 0.35, aft-peaked above; total lift peaks Fn 0.45; linear NK theory quantitatively good for Fn <= 0.35. NPL CATAMARANS: Insel & Molland 1992 and Molland, Wellicome & Couser 1994 (SSR 71) measured running trim and sinkage Fn 0.1-1.0 for demihulls alone and at s/L 0.2-0.5 ('sinkage and trim change significantly for Fn >= 0.35'); Gourlay 2008 catamaran theory = linear superposition of demihull potentials, so a thin-ship F_z for a catamaran would be the sum of the self term plus a cross term evaluated at the hull offset (the crate's existing multihull phase machinery). Deep-water zero-Fn limit checks: Havelock 1939 spheroid F_inf = rho U^2 A_W eps^2[ln(eps/2)+3/2-eps] (downward), eps = sqrt(12 Vol/(pi L^3)); M_inf = 0 for fore-aft symmetric hulls. Shallow-water limit (if ever wanted): Tuck 1966 c_s ≈ 1.3-1.5 nearly universal, c_theta = 0 for symmetric hulls. Standard-hull numbers I recall but could NOT verify online in this session (all data pages 403/404): KCS Fr 0.26 sinkage ≈ -1.39e-3 L (down), trim ≈ -0.17 deg (bow-down); DTMB 5415 Fr 0.28 sinkage ≈ -1.82e-3 L, trim ≈ -0.11 deg (Longo & Stern 2005) - treat as unverified.

**Notes**

NOT FOUND / NOT REACHABLE: Yeung 1972 full text (OnePetro 403) - the single most relevant formula source; Kim & Jenkins 1981 DTIC report (403); Doctors & Day 2000 IWWWFB-15 and 2003 24th-ONR papers (iwwwfb.org socket errors; NAP page serves images only); Insel & Molland 1992 and Molland/Wellicome/Couser SSR 71 full text (403 / wrong mirror); Havelock 1939 ZAMM original; Tuck 1963 thesis / 1964 JSR deep-water slender-body sinkage formula; Physics of Fluids 37:015122 (2025) 'Scale effects on ship vertical force and trim moment' (AIP 403); any paper stating a deep-water thin-ship F_z in the (I+iJ)-product spectral form. Michlet/Tuck-Lazauskas: confirmed Michlet does NOT compute sinkage/trim (Lazauskas), so no Michlet anchors exist. Eq (23) of Gourlay & Tuck 2001 was decoded from garbled PDF glyphs; the bracket [ln(eps/2) + 3/2 - eps] is my best reading (the '- eps' term is least certain) - the negative sign and eps^2 ln eps scaling are robust. Crate conventions confirmed against crates/michell/src/michell.rs module docs (nu = g/U^2, I+iJ = ∬ f_x e^{-nu lambda^2 z} e^{i nu lambda x}, R_w prefactor 4 rho g^2/(pi U^2), z downward, TransomClosure smoothstep appendage following Doctors & Day / Couser & Molland). Text dumps of every PDF read are in /private/tmp/claude-501/-Users-avi-michell/2fdb50af-0291-45b8-bdde-e7c422fd969b/scratchpad/ (gourlay_tuck2001_clean.txt, gourlay2011.txt, gourlay2008slender.txt, gourlay2008cat.txt, eccomas2000.txt, noblesse_sinkage.txt, ejmb2018.txt, thesis3d.txt [Baar 1986], tarafder.txt, iwwwfb28_22.txt). Sign-convention warning for the crate: Doctors and Noblesse use z UP with p = rho U phi_x for x forward; the crate has z DOWN and the bow at +x with stream -U x-hat, so n_z dS = -f_z dx dz per side and the upward force is F_up = -2∬ p1 f_z dx dz (plus bottom term if f(x,T) != 0); positive sinkage should be reported downward and positive trim bow-up to match every source above.

## Code integration map

**Reuse from the existing crate:**

- {'symbol': 'InnerIntegral::fill_zm / accumulate / transom_term', 'file': '/Users/avi/michell/crates/michell/src/michell.rs (lines ~886-960)', 'how': 'Reuse verbatim. They already take kx and kappa independently: eval_pair (line ~872) is the only place tying kx = nu*lambda to kappa = nu*lambda^2. fill_zm(kappa) fills self.zm (per z-span exp moments times e^{-kappa z0}); accumulate(kx, coeff) contracts any coefficient net with layout ((sx*nsz+sz)*p + a)*(q+1) + b; transom_term(kx) reads self.zm (whatever kappa it was filled at) and self.nu (physical nu, correct for the hollow length). Add one thin method eval_at(kx, kappa) that calls fill_zm then accumulate+transom_term; zero changes inside the three primitives. Note the existing kernel is e^{+i kx x}; the wanted Q(kx,kappa) with e^{-i kx x} is the complex conjugate (all coefficients real).'}
- {'symbol': 'osc_moments, exp_moments, C64, SERIES_THRESHOLD', 'file': '/Users/avi/michell/crates/michell/src/moments.rs', 'how': 'Verbatim. osc_moments(k,h,a_max) is called with a_max = p-1 today; the x-weighted transforms need a_max = p (one more degree), which the routine already supports. exp_moments to degree q+1 if a z-weighted (z*f_x) transform is wanted. The series/recurrence switch at SERIES_THRESHOLD=4 is the pattern to copy for the new diagonal moment.'}
- {'symbol': 'hollow_shape_moment', 'file': '/Users/avi/michell/crates/michell/src/michell.rs (line ~150)', 'how': "Reuse for the appendage's df/dx term. Add sibling hollow_value_moment (phi = 1 - 3s^2 + 2s^3 -> M0 - 3M2 + 2M3, osc_moments to degree 3) for f-type appendage transforms and hollow_arm_moment (s*phi' = 6s^3 - 6s^2) for the x-weighted appendage term."}
- {'symbol': 'compute_fx_coeff, BSplineSurface::corner_partials', 'file': '/Users/avi/michell/crates/michell/src/hull.rs (line ~410), /Users/avi/michell/crates/michell/src/bspline.rs (line 151)', 'how': 'Copy the pattern into a generalised compute_poly_coeff(surface, xs, zs, dx, dz, na, nb): coefficient of X^a Z^b in d^dx_x d^dz_z f is d[a+dx][b+dz]/(a! b!) for a+dx<=p, b+dz<=q, zero-padded to a common (na, nb) box. corner_partials already returns the full (p+1)x(q+1) Taylor table, so f, f_z, x-shifts and the keel line all come from the same call.'}
- {'symbol': 'detect_transom / Transom::coeff / Transom::depth', 'file': '/Users/avi/michell/crates/michell/src/hull.rs (lines ~500-560)', 'how': "tr.coeff (f_T(z) per z-span, layout sz*(q+1)+b) is exactly the z-factor the appendage transforms need; f_T'(z) coefficients are b*c[b] shifted down one degree. tr.depth feeds TransomClosure::hollow_length unchanged."}
- {'symbol': 'InnerIntegral::eval(lambda)', 'file': '/Users/avi/michell/crates/michell/src/michell.rs (line ~860)', 'how': 'Verbatim for the residue (wave) part of the force: the PV pole at k0(theta) = nu sec^2 theta is exactly the point (kx = nu lambda, kappa = nu lambda^2) that eval(lambda) evaluates, so the pole-subtraction value g(k0) reuses eval() for f_x and eval_at(kx0,kappa0) for the new nets.'}
- {'symbol': 'integrate_outer (marching theta panels, rate(), cap, STOP_WINDOW_PHASE truncation), run_outer (frac-halving refinement + est_rel), OuterParams, fleet_outer_params, fleet_phase_refs', 'file': '/Users/avi/michell/crates/michell/src/michell.rs (lines ~560-700, 760-800)', 'how': "Copy the theta-marching structure (GL_N=16 panels sized by rate(theta), quiet-window truncation past lambda=2, LAMBDA_HARD_CAP) into integrate_force_2d; the panel sizing is designed for the phase of products of hull transforms at the Michell dispersion point and is the right resolution for the theta direction of the force integral too. run_outer's refinement loop is the template for error control (optionally one refinement pass only)."}
- {'symbol': 'superpose / MemberWave, SourceMember', 'file': '/Users/avi/michell/crates/michell/src/michell.rs (lines ~700-830)', 'how': 'Pattern for multihull force: per-member transforms with placement phase exp(i(kx dx +/- ky dy)), ky = k sin theta in polar form; average the +/-theta systems as superpose does. Force on member j from all sources needs the cross products Q_fz,j * conj(sum_i Q_fx,i * phase_i), not |sum|^2, so write a sibling superpose_force rather than reusing superpose directly.'}
- {'symbol': 'gauss_legendre', 'file': '/Users/avi/michell/crates/michell/src/quadrature.rs', 'how': 'Verbatim for both the theta panels and the inner k panels (regular part on [0,2k0] after pole subtraction; geometric octave panels on the tail).'}
- {'symbol': 'Body::situate, FrameMap', 'file': '/Users/avi/michell/crates/michell/src/body.rs (line 129)', 'how': "Unchanged. It re-lofts the wetted hull in the water frame (x along the flow, z down from the actual free surface), so the situated Hull's f_x, f_z are exactly the water-frame slopes the force integrals want; placement.x = 0 for bodies."}
- {'symbol': 'solve_equilibrium_with, totals, Equilibrium, FleetState', 'file': '/Users/avi/michell/crates/michell/src/float.rs (lines 73-280)', 'how': 'Core Newton loop to extend (see equilibrium_hook). Keep the adaptive relax logic (lines 197-206), the step clamps, the best-near-miss acceptance (244-258) and the 20-degree abort (235-242) as-is; only the residuals r1 (line 165) and r2 (line 177) gain hydrodynamic terms.'}
- {'symbol': 'solve_equilibrium_bodies / solve_equilibrium_heeled coarse_opts recipe', 'file': '/Users/avi/michell/crates/michell/src/float.rs (lines 384-388, 515-519)', 'how': 'Same coarse_opts derivation (stations/2 >= 31, waterlines/2 >= 11, n_ctrl 10x7) for the force evaluated during iterations; the coarse hull (7x4 spans, cubic) is ~3x cheaper per transform than full resolution.'}
- {'symbol': 'multihull_resistance_with / validate_fleet', 'file': '/Users/avi/michell/crates/michell/src/lib.rs (line 163), michell.rs (line ~740)', 'how': 'validate_fleet verbatim at the top of the new dynamic_loads entry point; multihull_resistance_with is the model for a public fn multihull_dynamic_loads(members, cond, opts) -> DynamicLoads.'}
- {'symbol': 'manifest::run row assembly, header, options parsing, parse_transom', 'file': '/Users/avi/michell/crates/michell-cli/src/manifest.rs (lines 68-102, 364-393, 439-459, 504-573); /Users/avi/michell-cli main.rs line 653', 'how': "options block gains a 'dynamic' key parsed like 'transom'; header/nums chains gain the new columns; the equilibrium call moves inside the speed loop when dynamic is on (see cli_manifest)."}
- {'symbol': 'FreeWaveSpectrum::amp_pair / Member', 'file': '/Users/avi/michell/crates/michell/src/spectrum.rs (line 430)', 'how': "Reference for how a second consumer already holds InnerIntegral per member and evaluates at arbitrary (sec, ky); the force module should follow the same ownership pattern (Vec<Member<'h>> with InnerIntegral + dx + y)."}

### Generalise inner integral

Goal: evaluate Q_net(kx, kappa) = sum over span pairs of coefficient nets contracted with x-oscillatory and z-exponential moments, for arbitrary kx (real) and kappa >= 0, for several nets sharing one fill_zm.

1. New method on InnerIntegral (michell.rs, next to eval_pair ~line 872):
   pub(crate) fn eval_at(&mut self, kx: f64, kappa: f64) -> C64 {
       if !self.fill_zm(kappa) { return C64::ZERO; }
       self.accumulate(kx, self.hull.fx_coeff()) + self.transom_term(kx)
   }
   This is eval_pair with the lambda coupling removed; fill_zm/accumulate/transom_term are untouched. Q(kx,kappa) in the requested e^{-i kx x} convention = conj(eval_at(kx,kappa)); phases are relative to hull.x_center() (constant phase, cancels in Q_a * conj(Q_b) products of the same hull; for cross-hull products the existing dx = x_center + place.x - cx_ref offset phase applies).

2. Multi-net evaluation without re-filling zm: add
   pub(crate) fn eval_nets_at(&mut self, kx: f64, kappa: f64, nets: &[&PolyNet], out: &mut [C64])
   where PolyNet { coeff: Vec<f64>, na: usize, nb: usize, xweight: bool } and a generalised accumulate_net(kx, net) that (a) calls osc_moments(kx, sx.len, na - 1 + xweight as usize, ..) and (b) forms the per-a weight as M_a when xweight=false or (sx.start - x_center) * M_a + M_{a+1} when xweight=true (this yields the transform of (x - x_center) * poly exactly; the caller adds (x_center - x_ref) * Q for the arm about x_ref), and (c) indexes coeff as ((s*nsz+t)*na + a)*nb + b with zrow = zm[t*(q+1) .. t*(q+1)+nb] (nb <= q+1). The existing accumulate is the na=p, nb=q+1, xweight=false special case; leave it in place (bit-for-bit safety for all existing paths) and let accumulate_net be a sibling.

3. x-reduced form for the non-separable |z - zeta| kernel: add
   fn x_reduce(&mut self, kx: f64, net: &PolyNet, out: &mut [C64])  // out[t*nb + b] = sum_s phase_s * sum_a coeff[s,t,a,b] * M_a(kx)
   i.e. the same loops as accumulate but stopping before the contraction with zm. Then the separable transform is Q = sum_{t,b} out[t,b] * zm[t,b] (identical to accumulate up to rounding), and the |z-zeta| kernel term is sum_{t,t',b,b'} A[t,b] * conj(B[t',b']) * K_{t t' b b'}(kappa) with K from the new moments (see new_moments_needed). Cost per node: nsz^2 * nb^2 complex madds (~36^2 = 1300 for 9 z-spans, cubic), comparable to one accumulate.

4. Transom appendage at independent (kx, kappa): transom_term(kx) is already correct as long as fill_zm(kappa) ran first (z_factor = tr.coeff . zm, shape = hollow_shape_moment(-kx*L_v), phase = cis(kx*(tr.x - x_center)), L_v from self.transom.hollow_length(tr.depth, self.nu)). Generalise to fn transom_term_net(&mut self, kx, kind: AppendageKind) with kind in { DfDx (existing: phi', tr.coeff), F (phi, tr.coeff), DfDz (phi, tr.coeff shifted: c'[b] = (b+1) c[b+1]), XDfDx ((tr.x - x_center) * phi'-moment  - L_v * (s phi')-moment, tr.coeff) }, all using osc_moments(-kx*L_v, 1.0, 3). IMPORTANT physics asymmetry to encode: the appendage is a virtual source distribution, not hull surface. It must be included in the source-side net (Q_fx, which generates the pressure) but excluded from the force-side nets (Q_fz, Q_f, Q_xfz that weight the pressure over the real hull); expose a flag include_transom: bool on eval_nets_at, default true for f_x nets and false for the pressure-weighting nets. With TransomClosure::None the source side keeps the bare step (as today), and the 2-D tail then decays only algebraically (same caveat as TransomClosure::Fixed{length: 0} docs).

5. Keel line and waterline line transforms (for the closure/boundary terms the physics derivation may produce): waterline f(x,0) line net = coefficients b=0 of the first z-span (zm not needed: Q_line(kx) = sum_s phase_s sum_a c[s,0,a,0] M_a); keel line f(x,T) = for the last z-span, sum_b d[a][b]/(a! b!) h^b (from corner_partials) then e^{-kappa T} * x-transform. Both are one-dimensional osc_moments sums and reuse x_reduce with nb=1.

### New moments needed

Polynomial nets beyond fx_coeff (all from one corner_partials(sx, sz) call per span pair; fact[] as in compute_fx_coeff, hull.rs ~line 410-435), stored on Hull lazily (OnceCell) or built by the force module from hull.surface():

(a) f itself: coeff[a][b] = d[a][b]/(a! b!), a = 0..=p, b = 0..=q  -> PolyNet { na: p+1, nb: q+1 }. Needed for any integration-by-parts form and for the appendage/waterline closure terms.
(b) f_z = df/dz: coeff[a][b] = d[a][b+1]/(a! b!), a = 0..=p, b = 0..q-1 -> PolyNet { na: p+1, nb: q }. This is the pressure-weighting net of the vertical force (do NOT derive it from f by parts: a lofted hull need not close at the keel, so f(x,T) != 0 in general; see (e)).
(c) x * f_x and x * f_z (pitch moment arms): no new net; set xweight = true, which uses osc_moments to degree na (one higher) and the per-span weight (x0_s - x_center) M_a + M_{a+1}; add (x_center + place.x - x_ref) * Q for the arm about the fleet station x_ref (= load.lcg / pivot_x). Exact.
(d) optional z * f_x (arm of the longitudinal pressure force in the pitch moment, if the derivation keeps it): exp_moments to degree q+1 and per-span weight (z0_t) zm[b] + zm[b+1]; implement as a zweight flag on the contraction, not a new net.
(e) keel line f(x, T) and waterline line f(x, 0): 1-D nets (nb = 1) described in generalise_inner_integral item 5; the keel one matters because thin-ship f_z cannot represent a flat bottom, so a flat-keeled loft loses the bottom pressure unless the line term -2 int p(x,T) f(x,T) dx is added explicitly.
(f) Transom appendage nets: tr.coeff (f_T), shifted tr.coeff for f_T' (c'[b] = (b+1) c[b+1], degree q-1), together with the hollow moments int_0^1 {phi', phi, s phi'} e^{iKs} ds (osc_moments(K, 1.0, 3)).

Ordered z-pair / diagonal moments for the non-separable e^{-kappa |z - zeta|} kernel (moments.rs, new):
(g) Same-span diagonal D_{b b'}(kappa, h) = int_0^h int_0^h Z^b W^{b'} e^{-kappa |Z - W|} dZ dW, b,b' = 0..=q (symmetric, (q+1)(q+2)/2 distinct values per z-span). Small kappa*h (<= SERIES_THRESHOLD): series e^{-kappa|Z-W|} = sum_m (-kappa)^m |Z-W|^m / m! with the closed form int int Z^b W^{b'} |Z-W|^m = h^{b+b'+m+2} m! [ b'!/(b'+m+1)! + b!/(b+m+1)! ] / (b+b'+m+2) (Beta-function identity; entire in kappa, no cancellation). Large kappa*h: split W<Z and W>Z; int_0^Z W^{b'} e^{kappa W} dW = e^{kappa Z} P_{b'}(Z) - P_{b'}(0) with P_{b'}(Z) = sum_j (-1)^j b'!/(b'-j)! Z^{b'-j} kappa^{-(j+1)}, so int_0^h Z^b e^{-kappa Z}[e^{kappa Z} P(Z) - P(0)] dZ = int_0^h Z^b P(Z) dZ (pure polynomial) - P(0) N_b(kappa,h) (exp_moments); add the mirrored term with b <-> b'. Signature: pub fn exp_moments_diag(kappa: f64, h: f64, b_max: usize, out: &mut Vec<f64>) filling out[b*(b_max+1)+b'].
(h) Distinct spans t (deeper, z in span t) and t' (shallower, zeta in span t'): z - zeta = (z0_t - z1_{t'}) + Z + (h_{t'} - W) >= 0 always, so e^{-kappa(z-zeta)} = e^{-kappa (z0_t - z1_{t'})} e^{-kappa Z} e^{-kappa (h_{t'} - W)} is separable: N_b(kappa, h_t) (existing exp_moments) times the top-anchored moment N^+_{b'}(kappa, h_{t'}) = int_0^h W^{b'} e^{-kappa (h - W)} dW = sum_j C(b',j) h^{b'-j} (-1)^j N_j(kappa, h) (binomial in exp_moments; b' <= 3 so the alternating sum is benign; or compute directly with the same series/recurrence on (h-u)). The prefactor e^{-kappa (z0_t - z1_{t'})} <= 1 (no overflow). Provide pub fn exp_moments_top(kappa, h, b_max, out) and let the force module fill three per-span tables per kappa: N (already zm_raw), N^+, D. Assembly: K_{t t' b b'} = D_{bb'}(t) if t == t'; e^{-kappa(z0_t - z1_{t'})} N_b(t) N^+_{b'}(t') if t > t' (z deeper); symmetric counterpart if t < t'. Contract with the x-reduced vectors A[t,b] (f_z net) and conj(B[t',b']) (f_x net).

Note: exp_moments' underflow branch (decay == 0.0 in fill_zm) must be mirrored: skip pairs whose prefactor underflows.

### Outer quadrature

Coordinates: polar (theta, k) in the (kx, ky) wavenumber plane at y = 0 (both source and field points on the centreplane), kx = k cos theta, ky = k sin theta, kappa = k. The free-surface denominator (k nu - kx^2) = nu (k - nu sec^2 theta) has, at fixed theta, a single simple pole at k0(theta) = nu sec^2 theta for every theta in [0, pi/2) — exactly Michell's lambda = sec theta point, so the pole value reuses InnerIntegral::eval(lambda). (Cartesian (kx, ky) is worse: the pole at kx = nu, ky = 0 becomes a double zero in ky.) Use theta in [0, pi/2) with the +/- theta systems averaged as superpose does (for a single hull the integrand is even in theta).

Substitute k = k0(theta) u, u in (0, inf). Then kx = nu sec theta u, kappa = nu sec^2 theta u: as theta -> pi/2 both grow so the transforms decay and the theta integrand dies, which is what makes the theta truncation criteria of integrate_outer (quiet-window past lambda = 2, LAMBDA_HARD_CAP) applicable unchanged.

Theta direction: copy integrate_outer's marching panels verbatim (GL_N = 16, dt from rate(theta) with OuterParams from fleet_outer_params, cap near theta = 0). The integrand per theta node is the inner k-integral, so 'amp_sq(sec)' becomes 'force_density(theta)'.

Inner k (u) direction, per theta node, with g(u) = the product of transforms times the smooth prefactor:
  (i) regular part on [0, 2]: PV int_0^2 g(u)/(u-1) du = int_0^2 [g(u) - g(1)]/(u-1) du exactly (PV int_0^2 du/(u-1) = 0); integrate with one GL rule of ~24-32 nodes (the subtracted integrand is smooth; g(1) comes from eval(lambda) / eval_at(nu lambda, nu lambda^2)). Guard u -> 1 by evaluating the difference quotient with a Taylor fallback when |u-1| < 1e-6 (or simply place GL nodes so none hits u = 1: even GL rules have no node at the centre).
  (ii) tail [2, u_max]: g(u)/(u-1) with g decaying algebraically (~1/u^2: at large kappa each transform tends to f-line(x,0)/kappa, so the integrand ~ 1/u^3 to 1/u^2 depending on the kx factor); use geometric octave panels [2,4],[4,8],... each with GL 8, stopping when a panel contributes < 1e-6 relative or u > 2^12; typical 80-100 nodes. Optional acceleration: subtract the analytic large-kappa asymptote (waterline line transforms / kappa^2) and integrate it in closed form.
  Total ~110-130 (kx,kappa) nodes per theta node, each evaluating fill_zm once and accumulate_net for 3-4 nets (f_x with transom, f_z, x-weighted f_z, optionally f) plus the |z-zeta| contraction.

Residue (wave) part: the i*pi residue at u = 1 is a 1-D theta integral of Q_fz(k0) * conj(Q_fx(k0)) sec^3-type weight — evaluate it in the same theta loop at the same panel nodes (free; reuses eval(lambda)).

Rankine (double-body, speed-independent up to rho U^2) part: no pole; same polar grid without the subtraction, or (cheaper) integrate once per attitude on a fixed (theta, k) grid independent of nu and cache it while the speed loop runs (its kappa values need not track nu; pick k panels geometric from 0.1/L to 1e3/T).

Error control: replicate run_outer — pass 1 at frac = 1, one refinement at frac = 0.5, est_rel = |refined - first|/|refined|; for the equilibrium loop stop at rel_tol_force = 1e-3 (a force fed into a Newton solve with 1e-3 volume tolerance needs no more), and only do the refinement pass on the final full-resolution report. Expose as ForceOptions { rel_tol: 1e-3, max_refinements: 1, k_regular_nodes: 32, k_tail_nodes_per_octave: 8, k_octaves: 12 } beside WaveOptions.

Sanity checks to build in: (1) Q evaluated at theta nodes at u = 1 must equal InnerIntegral::eval(sec theta) to 1e-12 (same code path); (2) the residue-only part of F_z for a Wigley hull must match the known Michell-consistent wave-lift sign/scale; (3) the total F_z must converge as k_octaves grows (tail 1/u^2 decay visible in the panel sums).

### Equilibrium hook

Signature change (float.rs line 107): keep solve_equilibrium_with(situate, load, density) as a thin wrapper and add

pub struct DynamicLoads { pub lift_volume: f64, pub moment_volume: f64, pub fz: f64, pub my: f64 }  // lift_volume = Fz_up/(rho g) [m^3], moment_volume = My_bowup_about_pivot/(rho g) [m^4]; fz [N], my [N m] carried for reporting

pub fn solve_equilibrium_dynamic_with(
    mut situate: impl FnMut(f64, f64, bool) -> Result<FleetState>,
    mut hydro: impl FnMut(&FleetState, f64, f64, bool) -> Result<Option<DynamicLoads>>,  // (fleet, s, tau, coarse); None = hydrostatic only
    load: &LoadCase, density: f64, start: Option<(f64, f64)>,
) -> Result<Equilibrium>

The core divides by rho g nowhere (v_target = mass/density), so the closure returns loads already in volume units; it owns cond (rho, g), WaveOptions/ForceOptions and the moment reference x_ref = load.lcg.unwrap_or(0.0) (= pivot_x, line 128), so no gravity parameter is needed in the core.

Line-level hooks inside the existing loop:
- line 130-131: initialise s, tau from `start` (warm start from the hydrostatic solution or the previous speed) instead of 0.0.
- line 151: after `let fleet = situate(s, tau, coarse)?;` and `let t = totals(&fleet);` (157) add `let h = hydro(&fleet, s, tau, coarse)?.unwrap_or_default();`.
- line 165: `let r1 = t.volume + h.lift_volume - v_target;` (an upward hydrodynamic force reduces the buoyancy required; a suction, lift_volume < 0, sinks the hull).
- line 177: `let r2 = t.moment_x - lcg * t.volume + h.moment_volume;` with h.moment_volume the bow-up moment about x_ref = lcg. (If a different x_ref is used: add (M_ref + Fz (x_ref - lcg))/(rho g).) Sign check: positive r2 means net bow-up moment, which the existing Newton step turns into an increase of tau (positive tau raises the +x = bow end), consistent with lines 191-193 and 230-231.
- lines 179-182 (Jacobian): default leave the waterplane Jacobian untouched — Picard (frozen-load) iteration: the hydrodynamic loads act as a residual offset only. Convergence condition is spectral radius of K_hs^{-1} K_hd < 1, with K_hd = -d(Fz,My)/d(s,tau) ~ rho U^2 B (order) versus K_hs = rho g A_w ~ rho g L B, ratio ~ Fn^2: Picard converges comfortably for Fn <~ 0.6, i.e. the displacement regime thin-ship theory covers. Optional quasi-Newton for robustness: DynamicLoads gets an Option<[[f64;2];2]> jac (d lift_volume/ds, d lift_volume/dtau, d moment_volume/ds, d moment_volume/dtau) that the hook adds to dv_ds, dv_dt, dm_ds, dm_dt; obtain it once per phase by two extra force evaluations at (s + 0.05 z_scale, tau) and (s, tau + 0.005) on the coarse fleet, or by rank-1 Broyden updates from successive (delta residual, delta state) pairs. Do not finite-difference every iteration (loft roughness ~1e-3 of volume, noted at lines 143-147, would make the FD Jacobian noisy).
- lines 198-204 (relaxation): key the sign-flip test on both r1 and r2 (currently r1 only) since the dynamic moment can drive a trim limit cycle at the sinkage hump.
- line 262: final `situate(s, tau, false)` then one full-resolution `hydro(&fleet, s, tau, false)` for the reported fz/my; add to Equilibrium: `dynamic: Option<DynamicLoads>`, `static_sinkage: f64`, `static_trim: f64` (from the hydrostatic pre-solve, so the dynamic increments are reportable).

Lazy re-evaluation (the key cost control, implemented inside the closure in a new float::solve_equilibrium_dynamic_bodies(bodies, water_offset, poses, load, density, cond, wave_opts, force_opts, opts) mirroring solve_equilibrium_bodies lines 368-418): keep (s_last, tau_last, loads_last); recompute only when |s - s_last| > 0.02 * z_scale or |tau - tau_last| > 0.002 rad, or on phase change (coarse -> fine). Typical result: hydrostatic Newton converges in 2-7 cheap iterations between 2-4 force evaluations. Recommended driver order per (point, speed): (1) hydrostatic solve (existing, warm start none); (2) dynamic solve warm-started from (1) or from the previous speed's dynamic solution (continuation in U, which keeps the Picard steps small); (3) full-resolution report.

Heeled path: solve_equilibrium_heeled (line 467) wraps the same core; the upright force kernel on the heel_poses reposition is an approximation (no complex-kappa force kernel yet); either forbid dynamic + heel != 0 initially or document it as the reposition-only approximation.

Convergence risks: (a) Fn 0.45-0.55 (stern-down trim hump): loads change fastest with attitude and the lazy threshold may need tightening; (b) the appendage hollow length L_v depends on tr.depth, which changes with attitude (continuous, but the transom detection threshold TRANSOM_AREA_REL = 1e-3 can switch the transom on/off between iterations -> discontinuous residual; mitigate by freezing transom presence at the hydrostatic attitude for the dynamic solve); (c) a member going dry mid-iteration drops its force contribution (already handled for buoyancy at lines 152-156); (d) the 20-degree abort (235) is still appropriate; (e) hydrostatic Jacobian underestimates when lift is a large fraction of weight — cap |lift_volume| at, say, 0.5 v_target with a diagnostic rather than let Newton wander (planing is out of scope).

### Situate cost

Measured (release, LTO, this machine) on a Wigley 10 m x 1 m x 0.625 m body, default BodyOptions (121x33 samples, 20x12 cubic net) and the coarse Newton options (61x17, 10x7 net):
- Body::situate: 72.6 ms full, 12.7 ms coarse (sample + loft). Negligible against the force.
- Existing wave resistance: 33 ms full (~19k inner evaluations => ~1.7 us per fill_zm+accumulate+transom_term on 17x9 spans), 3.8 ms coarse (~7k, ~0.55 us on 7x4 spans). Theta nodes for the frac=1 pass: ~6400 full / ~2700 coarse.
- Hydrostatic solve_equilibrium_bodies: 157 ms (2 iterations, on-design) to 370 ms (7 iterations, -15% mass, lcg -0.15 m).

Force evaluation estimate (3-4 nets sharing fill_zm ~2.5x one eval, ~120 k-nodes per theta node, frac = 1 only): coarse 2700 x 120 x 0.55 us x 2.5 ~= 0.45 s; full 6400 x 120 x 1.7 us x 2.5 ~= 3.3 s; the |z-zeta| contraction (~1300 complex madds per node) adds ~30-50%. With one refinement pass (frac 0.5) triple those. Per (point, speed): hydrostatic solve 0.2-0.4 s + 3 lazy coarse force evaluations ~1.5-2 s + one full-resolution report ~3-5 s (10 s with refinement) => about 5-8 s single hull, 15-30 s for a trimaran (3 members: transforms per member once per node, cross-member products cheap). A 12-speed x 20-point sweep is therefore 20-90 min single-threaded; points and speeds are independent, so std::thread::scope over (point, speed) gives near-linear speed-up (already README roadmap item 1). Cheapest sane default: solve the attitude entirely on the coarse fleet and evaluate the full-resolution force only for the reported row; optionally cache the speed-independent Rankine part per attitude across the speed loop.

### Cli manifest

Manifest (crates/michell-cli/src/manifest.rs):
- options block (lines 68-102): add `"dynamic": "off" | "report" | "solve"` (default off) parsed beside "transom" via a parse_dynamic(Option<&str>) modelled on main.rs parse_transom (line 653), plus optional `"dynamic_tol"` (ForceOptions.rel_tol, default 1e-3). "report" evaluates (Fz, My) once per row at the hydrostatic/fixed attitude (cheap diagnostic, works without a weight axis); "solve" requires float_mode (weight axis) — reject otherwise next to the existing axis-constraint checks at lines 289-309; "solve" with a heel axis: reject (or document as reposition-only approximation).
- Restructure (lines 439-459 and 504-573): with "solve", keep the hydrostatic solve_equilibrium_heeled per point as the warm start and for the static columns, then move a per-speed `solve_equilibrium_dynamic_bodies(&bodies, 0.0, &poses, &LoadCase{mass,lcg}, density, &cond, &wave_opts, &force_opts, &bopts, start)` inside `for &u in &speeds` and build `members` from that equilibrium's fleet before calling multihull_resistance_with — so rw/rv/rt are computed on the dynamic attitude. Continuation: pass the previous speed's (sinkage, trim) as `start`.
- Columns (header lines 364-386, nums 536-554): keep `sinkage`, `trim_deg` as the attitude the resistance was computed at (dynamic when solving); append after `lcb`: `sinkage_static`, `trim_static_deg` (hydrostatic values), `sinkage_dyn` (= sinkage - sinkage_static), `trim_dyn_deg`, `fz` (N, positive up), `my` (N m, bow-up about lcg), `lift_frac` (= fz/(mass g); NaN/0 in non-float mode), `dyn_iterations`. Emit them when dynamic != off so existing CSVs are unchanged by default.
- CLI flags (main.rs): `--dynamic MODE` in PHYSICS OPTIONS help (line ~124) and parsed where `--transom` is (lines 735, 1106); wire into the `--float` IGES sweep (line 1199, solve_equilibrium -> dynamic variant; the SourceFleet path shares the generic core) and into the single-hull `resistance` command as a `fz`/`my`/`lift_frac` report.
- README: Sweeps paragraph (lines 337-346) drops "All poses are hydrostatic (no speed-dependent squat)" and documents `options.dynamic`; lib.rs assumptions (lines 48-49) likewise; add `multihull_dynamic_loads` and `float::solve_equilibrium_dynamic_bodies` to the API list (~line 505-517).
- Public API (lib.rs re-exports): `pub use michell::{multihull_dynamic_loads, dynamic_loads, DynamicLoads, ForceOptions}`; `float::{solve_equilibrium_dynamic_with, solve_equilibrium_dynamic_bodies}`.

### Risks

1. Physics/sign conventions: the crate's frame (bow at +x, free stream -U x-hat, z down) must be matched in p = rho U phi_x sign, in n_z = -f_z, and in the Havelock Green function's radiation condition (waves trail toward -x). Wrong signs give a bow-down instead of bow-up trim and are easy to miss; validate against Wigley squat data / Tuck's deep-water sinkage curve before trusting.
2. Thin-ship f_z cannot represent a flat bottom: lofted hulls with residual f(x,T) > 0 lose the keel pressure unless the explicit keel-line term is added (see new_moments (e)); similarly the waterline term if the derivation keeps one.
3. Transom: the appendage belongs to the source side only; including it on the pressure-weighting side double counts; TransomClosure::None leaves a step whose 2-D tail decays algebraically (slow k convergence); transom on/off detection (TRANSOM_AREA_REL) can toggle between Newton iterations — freeze it during the dynamic solve.
4. Numerical: the PV subtraction must never evaluate at u = 1 exactly except via the eval(lambda) path; the k-tail decays only algebraically (waterline contributions ~1/kappa), so truncation error must be monitored (panel sums); at theta near pi/2 the k-substitution k = k0 u is what keeps the integrand decaying — do not use a fixed absolute k grid there. The non-separable |z-zeta| diagonal moment needs a series branch (cancellation in kappa^{-(j+1)} terms at small kappa h).
5. Cost: 5-30 s per (point, speed) makes existing sweeps 10-100x slower when dynamic is on; keep it opt-in, coarse-resolution during iterations, lazy force re-evaluation, continuation in speed, and parallelise over rows.
6. Convergence: Picard on frozen loads is safe below Fn ~0.6 but the sinkage/trim hump (Fn 0.45-0.55) and light hulls with large lift fractions can limit-cycle; extend the relax logic to r2 and offer the FD/Broyden Jacobian; cap |lift| relative to weight with a diagnostic (planing is outside thin-ship theory).
7. Loft roughness: the force is evaluated on a re-lofted hull each iteration; ~1e-3 relative roughness in volume (float.rs comment) implies similar noise in Fz; do not build Jacobians from per-iteration finite differences.
8. Heel + dynamic: the upright force kernel on heel_poses is an approximation; a complex-kappa force kernel (mirroring HeelInner) is future work. Asymmetric hulls' dipole contribution to the vertical force is not covered by the source-only plan.
9. Bit-for-bit safety: add new methods (eval_at, accumulate_net, x_reduce) beside accumulate rather than refactoring it, so every existing resistance path stays bit-identical; the multihull force needs a superpose_force (cross products, not |sum|^2) — do not overload superpose.
10. Manifest semantics change: with dynamic solve the equilibrium becomes per-speed; the current per-point `state`/`members` reuse across speeds (manifest.rs 501-504) must move inside the speed loop, and the `volume`/`lcb` columns then also become per-row dynamic values — document the column meanings.

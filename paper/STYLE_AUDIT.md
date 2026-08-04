# Paper A style audit

Frozen source: `paper-a-jsr-v3` (`9cd6c79dd870c20e9b4226a9e93e7636f9394a185719428a9deee128ff60c0ed`).
The mechanical comparison excludes only the marked Nomenclature apparatus,
which repeats symbols and units already defined in the frozen text.

## Mechanical freeze checks

| Check | Result | Before | After |
|---|---|---:|---:|
| Numeral token sequence | PASS | 450 | 450 |
| Unit token sequence | PASS | 16 | 16 |
| Equation environments, byte for byte | PASS | 20 | 20 |
| Citation-key sequence | PASS | 39 | 39 |

## Sentence-length histogram

| Words | Before | After |
|---|---:|---:|
| 0--7 | 19 | 19 |
| 8--14 | 61 | 61 |
| 15--21 | 68 | 68 |
| 22--28 | 36 | 36 |
| 29--40 | 19 | 19 |
| 41+ | 7 | 7 |

## Register-term grep

| Term | Before | After |
|---|---:|---:|
| `Note that` | 0 | 0 |
| `It is worth noting` | 0 | 0 |
| `It should be noted` | 0 | 0 |
| `Importantly` | 0 | 0 |
| `Crucially` | 0 | 0 |
| `Specifically,` | 0 | 0 |
| `In particular,` | 1 | 1 |
| `Moreover` | 0 | 0 |
| `Furthermore` | 0 | 0 |
| `comprehensive` | 0 | 0 |
| `robust` | 0 | 0 |
| `novel` | 0 | 0 |
| `carefully` | 0 | 0 |
| `deliberately` | 0 | 0 |
| `explicitly` | 2 | 2 |
| `we emphasize` | 0 | 0 |
| `we stress` | 0 | 0 |
| `leverages` | 0 | 0 |
| `utilizes` | 0 | 0 |
| `prior to` | 0 | 0 |
| `in order to` | 0 | 0 |
| `is able to` | 0 | 0 |
| `a number of` | 0 | 0 |

## Frozen numeral and unit manifest

### Numerals

- N001 `1960` — line 5: `Michelsen's 1960 dissertation and 1972 sequel developed analytic reductions`
- N002 `1972` — line 5: `Michelsen's 1960 dissertation and 1972 sequel developed analytic reductions`
- N003 `16` — line 8: `B-spline hulls of degree at most 16 in each direction.  Repeated integration`
- N004 `1` — line 13: `K_s(\omega)=\int_1^\infty`
- N005 `2` — line 14: `\e^{\iu\omega\lambda}\lambda^{-s}(\lambda^2-1)^{-1/2}\,d\lambda,`
- N006 `1` — line 14: `\e^{\iu\omega\lambda}\lambda^{-s}(\lambda^2-1)^{-1/2}\,d\lambda,`
- N007 `1` — line 14: `\e^{\iu\omega\lambda}\lambda^{-s}(\lambda^2-1)^{-1/2}\,d\lambda,`
- N008 `2` — line 14: `\e^{\iu\omega\lambda}\lambda^{-s}(\lambda^2-1)^{-1/2}\,d\lambda,`
- N009 `0.02` — line 22: `For the standard Wigley hull at \(\Fn=0.02\), the method uses 1,152 kernel nodes`
- N010 `1,152` — line 22: `For the standard Wigley hull at \(\Fn=0.02\), the method uses 1,152 kernel nodes`
- N011 `511,504` — line 23: `instead of 511,504 inner-transform evaluations, reduces median runtime from`
- N012 `21.541` — line 24: `\SI{21.541}{ms} to \SI{0.030}{ms} (about 700-fold), and differs from an`
- N013 `0.030` — line 24: `\SI{21.541}{ms} to \SI{0.030}{ms} (about 700-fold), and differs from an`
- N014 `700` — line 24: `\SI{21.541}{ms} to \SI{0.030}{ms} (about 700-fold), and differs from an`
- N015 `1.83` — line 25: `independent real-axis reference by \(1.83\times10^{-13}\) relatively.  The`
- N016 `13` — line 25: `independent real-axis reference by \(1.83\times10^{-13}\) relatively.  The`
- N017 `1898` — line 40: `Michell's 1898 thin-ship theory expresses wave resistance as a one-dimensional`
- N018 `2` — line 66: `\(\nu=g/U^2\), where \(U\) is speed and \(g\) gravitational acceleration.`
- N019 `5.123` — line 75: `criterion reported \(5.123\times10^{-9}\) relative error at \(\Fn=0.05\), but`
- N020 `9` — line 75: `criterion reported \(5.123\times10^{-9}\) relative error at \(\Fn=0.05\), but`
- N021 `0.05` — line 75: `criterion reported \(5.123\times10^{-9}\) relative error at \(\Fn=0.05\), but`
- N022 `5.781` — line 76: `the measured error was \(5.781\times10^{-7}\), about 113 times larger.  The`
- N023 `7` — line 76: `the measured error was \(5.781\times10^{-7}\), about 113 times larger.  The`
- N024 `113` — line 76: `the measured error was \(5.781\times10^{-7}\), about 113 times larger.  The`
- N025 `10` — line 80: `claims that the requested \(10^{-8}\) tolerance was achieved.`
- N026 `8` — line 80: `claims that the requested \(10^{-8}\) tolerance was achieved.`
- N027 `83` — line 88: `\citep[pp.~83--85]{Gotman2002}.  Huybrechs and Vandewalle develop numerical steepest descent`
- N028 `85` — line 88: `\citep[pp.~83--85]{Gotman2002}.  Huybrechs and Vandewalle develop numerical steepest descent`
- N029 `9.80665` — line 118: `All calculations below use \(g=\SI{9.80665}{m.s^{-2}}\), freshwater density`
- N030 `2` — line 118: `All calculations below use \(g=\SI{9.80665}{m.s^{-2}}\), freshwater density`
- N031 `999.1` — line 119: `\(\rho=\SI{999.1}{kg.m^{-3}}\), and the length Froude number`
- N032 `3` — line 119: `\(\rho=\SI{999.1}{kg.m^{-3}}\), and the length Froude number`
- N033 `0` — line 126: `rectangle \(S=[x_0,x_1]\times[z_0,z_1]\), its longitudinal derivative is an`
- N034 `1` — line 126: `rectangle \(S=[x_0,x_1]\times[z_0,z_1]\), its longitudinal derivative is an`
- N035 `0` — line 126: `rectangle \(S=[x_0,x_1]\times[z_0,z_1]\), its longitudinal derivative is an`
- N036 `1` — line 126: `rectangle \(S=[x_0,x_1]\times[z_0,z_1]\), its longitudinal derivative is an`
- N037 `2` — line 142: `\(\kappa=\nu\lambda^2\).  The result is an exact endpoint representation.`
- N038 `1` — line 157: `(\nu\lambda)^{-(r+1)}(\nu\lambda^2)^{-(u+1)}.`
- N039 `2` — line 157: `(\nu\lambda)^{-(r+1)}(\nu\lambda^2)^{-(u+1)}.`
- N040 `1` — line 157: `(\nu\lambda)^{-(r+1)}(\nu\lambda^2)^{-(u+1)}.`
- N041 `3` — line 159: `The resulting power of \(\lambda\) is \(n=r+2u+3\).  Both sums terminate at`
- N042 `1` — line 166: `full-multiplicity chine.  Across \(1\leq\lambda\leq20\), every discrepancy`
- N043 `5` — line 167: `falls below the scaled \(5\times10^{-10}\) threshold.`
- N044 `10` — line 167: `falls below the scaled \(5\times10^{-10}\) threshold.`
- N045 `0` — line 173: `\(z_j=0\), and submerged terms \(S\), for which \(z_j>0\).  At low Froude`
- N046 `0` — line 173: `\(z_j=0\), and submerged terms \(S\), for which \(z_j>0\).  At low Froude`
- N047 `2` — line 178: `Expanding \(|\widetilde F|^2\) before outer integration removes the`
- N048 `2` — line 181: `where \(s_{ij}=n_i+n_j-2\geq4\),`
- N049 `2` — line 192: `Expand \(|F|^2-|\widetilde F|^2\) into precisely the ordered pairs in`
- N050 `2` — line 192: `Expand \(|F|^2-|\widetilde F|^2\) into precisely the ordered pairs in`
- N051 `1` — line 194: `\(|\e^{\iu\omega\lambda}|=1\), and bound`
- N052 `2` — line 195: `\(\e^{-\nu\lambda^2(z_i+z_j)}\leq\e^{-\nu(z_i+z_j)}\) for`
- N053 `0` — line 196: `\(\lambda\geq1\).  The remaining positive integral is \(K_{s_{ij}}(0)\).`
- N054 `1935` — line 210: `\(\Ki_s\) the ``Bickley function'' and cites the 1935 paper by W.~G. Bickley`
- N055 `1` — line 216: `The substitution \(\lambda=1+t^2\) removes the square-root endpoint in`
- N056 `2` — line 216: `The substitution \(\lambda=1+t^2\) removes the square-root endpoint in`
- N057 `0` — line 219: `For \(\omega>0\), analyticity in the intervening sector permits the exact`
- N058 `4` — line 221: `\(t=\e^{\iu\pi/4}r/\sqrt\omega\), giving`
- N059 `0` — line 233: `\(0\leq\arg t\leq\pi/4\).  It is therefore analytic throughout the deformation.`
- N060 `4` — line 233: `\(0\leq\arg t\leq\pi/4\).  It is therefore analytic throughout the deformation.`
- N061 `2` — line 235: `\(|\e^{\iu\omega t^2}|=\e^{-\omega R^2\sin(2\phi)}\leq1\), while the`
- N062 `2` — line 235: `\(|\e^{\iu\omega t^2}|=\e^{-\omega R^2\sin(2\phi)}\leq1\), while the`
- N063 `2` — line 235: `\(|\e^{\iu\omega t^2}|=\e^{-\omega R^2\sin(2\phi)}\leq1\), while the`
- N064 `24` — line 240: `The implementation normally applies 24- and 48-point Gauss--Legendre rules on`
- N065 `48` — line 240: `The implementation normally applies 24- and 48-point Gauss--Legendre rules on`
- N066 `0` — line 241: `\(0\leq r\leq8\).  When \(s/|\omega|\geq2\), the continued algebraic factor is`
- N067 `48` — line 242: `narrow near the origin, so the code uses 48- and 96-point rules.  The neglected`
- N068 `96` — line 242: `narrow near the origin, so the code uses 48- and 96-point rules.  The neglected`
- N069 `1` — line 245: `using \(|1+\iu r^2/\omega|\geq1\),`
- N070 `2` — line 245: `using \(|1+\iu r^2/\omega|\geq1\),`
- N071 `2` — line 246: `\(|\sqrt{2+\iu r^2/\omega}|\geq\sqrt2\), and the standard Gaussian-tail`
- N072 `2` — line 246: `\(|\sqrt{2+\iu r^2/\omega}|\geq\sqrt2\), and the standard Gaussian-tail`
- N073 `25` — line 247: `bound.  At the minimum supported nonzero frequency \(\omega=25\), this is`
- N074 `4` — line 248: `below \(4\times10^{-30}\).  The selected coarse/fine difference estimates`
- N075 `30` — line 248: `below \(4\times10^{-30}\).  The selected coarse/fine difference estimates`
- N076 `4` — line 251: `is not certified.  A committed fixed-rule scan over \(4\leq s\leq400\) and`
- N077 `25,100,400` — line 252: `\(|\omega|\in\{25,100,400\}\), using a 192-node same-contour reference and a`
- N078 `192` — line 252: `\(|\omega|\in\{25,100,400\}\), using a 192-node same-contour reference and a`
- N079 `10` — line 253: `\(10^{-10}\) comparison floor, found a maximum relative 48-node discrepancy`
- N080 `10` — line 253: `\(10^{-10}\) comparison floor, found a maximum relative 48-node discrepancy`
- N081 `48` — line 253: `\(10^{-10}\) comparison floor, found a maximum relative 48-node discrepancy`
- N082 `1.84` — line 254: `of \(1.84\times10^{-11}\) through \(s=50\), but \(3.16\times10^{-1}\) above`
- N083 `11` — line 254: `of \(1.84\times10^{-11}\) through \(s=50\), but \(3.16\times10^{-1}\) above`
- N084 `50` — line 254: `of \(1.84\times10^{-11}\) through \(s=50\), but \(3.16\times10^{-1}\) above`
- N085 `3.16` — line 254: `of \(1.84\times10^{-11}\) through \(s=50\), but \(3.16\times10^{-1}\) above`
- N086 `1` — line 254: `of \(1.84\times10^{-11}\) through \(s=50\), but \(3.16\times10^{-1}\) above`
- N087 `24` — line 255: `that range.  The 24/48 difference covered every meaningful 48/192 discrepancy`
- N088 `48` — line 255: `that range.  The 24/48 difference covered every meaningful 48/192 discrepancy`
- N089 `48` — line 255: `that range.  The 24/48 difference covered every meaningful 48/192 discrepancy`
- N090 `192` — line 255: `that range.  The 24/48 difference covered every meaningful 48/192 discrepancy`
- N091 `2.22` — line 256: `by at least 2.22 in this scan.  This measured degree dependence motivates the`
- N092 `48` — line 257: `adaptive 48/96 rule within the validated \(s\leq98\) envelope; it is evidence,`
- N093 `96` — line 257: `adaptive 48/96 rule within the validated \(s\leq98\) envelope; it is evidence,`
- N094 `72` — line 260: `For fixed geometry, a nonzero endpoint pair costs either 72 or 144 scalar`
- N095 `144` — line 260: `For fixed geometry, a nonzero endpoint pair costs either 72 or 144 scalar`
- N096 `16` — line 263: `one longitudinal span and one vertical span.  Its 16 endpoint terms include`
- N097 `32` — line 264: `eight at the waterline.  Their full ordered-pair expansion contains 32`
- N098 `16` — line 266: `these to 16 kernel evaluations.  Each uses the 24/48 rule, for`
- N099 `24` — line 266: `these to 16 kernel evaluations.  Each uses the 24/48 rule, for`
- N100 `48` — line 266: `these to 16 kernel evaluations.  Each uses the 24/48 rule, for`
- N101 `16` — line 267: `\(16\times72=1{,}152\) scalar nodes in all; equal-position pairs use`
- N102 `1` — line 267: `\(16\times72=1{,}152\) scalar nodes in all; equal-position pairs use`
- N103 `152` — line 267: `\(16\times72=1{,}152\) scalar nodes in all; equal-position pairs use`
- N104 `0` — line 280: `its spline domain beginning exactly at \(z=0\), and has`
- N105 `0` — line 283: `\(0<|\omega_{ij}|<25\);`
- N106 `25` — line 283: `\(0<|\omega_{ij}|<25\);`
- N107 `1` — line 296: `pair \((i,j)\), let \(m_{ij}=1\) on the diagonal and 2 otherwise, and let`
- N108 `2` — line 296: `pair \((i,j)\), let \(m_{ij}=1\) on the diagonal and 2 otherwise, and let`
- N109 `2` — line 320: `\(\omega=\nu L=gL/U^2=1/\Fn^2\), which equals the implementation's`
- N110 `1` — line 320: `\(\omega=\nu L=gL/U^2=1/\Fn^2\), which equals the implementation's`
- N111 `2` — line 320: `\(\omega=\nu L=gL/U^2=1/\Fn^2\), which equals the implementation's`
- N112 `2` — line 321: `\(\nu|x_i-x_j|\).  Thus \(\omega\geq25\) corresponds to \(\Fn\leq0.2\).  The`
- N113 `1` — line 324: `\(\omega_{\rm adj}=1/(m\Fn^2)\), so the separation gate requires`
- N114 `2` — line 324: `\(\omega_{\rm adj}=1/(m\Fn^2)\), so the separation gate requires`
- N115 `0.05` — line 326: `At \(\Fn=0.05\), for example, this permits at most 16 such spans.  A retained`
- N116 `16` — line 326: `At \(\Fn=0.05\), for example, this permits at most 16 such spans.  A retained`
- N117 `360` — line 327: `multi-span endpoint map with 360 distinct nonzero-frequency pairs would use`
- N118 `360` — line 328: `\(360\times72=25{,}920\) ordinary-rule nodes, compared with 1,152 for the`
- N119 `25` — line 328: `\(360\times72=25{,}920\) ordinary-rule nodes, compared with 1,152 for the`
- N120 `920` — line 328: `\(360\times72=25{,}920\) ordinary-rule nodes, compared with 1,152 for the`
- N121 `1,152` — line 328: `\(360\times72=25{,}920\) ordinary-rule nodes, compared with 1,152 for the`
- N122 `16` — line 329: `16-pair Wigley map; this is node-count arithmetic, not a runtime benchmark.`
- N123 `10` — line 333: `default relative tolerance \(10^{-5}\), it rejects the Wigley reduction at`
- N124 `5` — line 333: `default relative tolerance \(10^{-5}\), it rejects the Wigley reduction at`
- N125 `0.08` — line 334: `\(\Fn=0.08\) and accepts it at \(\Fn=0.05\).`
- N126 `0.05` — line 334: `\(\Fn=0.08\) and accepts it at \(\Fn=0.05\).`
- N127 `4` — line 340: `classical catamaran \(4\cos^2\) interference, tolerance self-convergence, and`
- N128 `2` — line 340: `classical catamaran \(4\cos^2\) interference, tolerance self-convergence, and`
- N129 `10` — line 350: `with \(L=\SI{10}{m}\), \(B=\SI{1}{m}\), and`
- N130 `1` — line 350: `with \(L=\SI{10}{m}\), \(B=\SI{1}{m}\), and`
- N131 `0.625` — line 351: `\(T=\SI{0.625}{m}\).  The primary reference shares neither the production`
- N132 `1` — line 353: `transform on the positive real axis after \(\lambda=1+t^2\), which converts`
- N133 `2` — line 353: `transform on the positive real axis after \(\lambda=1+t^2\), which converts`
- N134 `2` — line 355: `\(2\lambda^2|F(\lambda)|^2/\sqrt{2+t^2}\).  This explicitly removes the`
- N135 `2` — line 355: `\(2\lambda^2|F(\lambda)|^2/\sqrt{2+t^2}\).  This explicitly removes the`
- N136 `2` — line 355: `\(2\lambda^2|F(\lambda)|^2/\sqrt{2+t^2}\).  This explicitly removes the`
- N137 `2` — line 355: `\(2\lambda^2|F(\lambda)|^2/\sqrt{2+t^2}\).  This explicitly removes the`
- N138 `2` — line 355: `\(2\lambda^2|F(\lambda)|^2/\sqrt{2+t^2}\).  This explicitly removes the`
- N139 `1` — line 356: `square-root singularity at \(\lambda=1\); \citet[p.~371]{Tuck1989} identifies`
- N140 `371` — line 356: `square-root singularity at \(\lambda=1\); \citet[p.~371]{Tuck1989} identifies`
- N141 `2` — line 361: `\(\nu\lambda L/2\) advances by at most \(\pi/2\).  Their endpoints are then`
- N142 `2` — line 361: `\(\nu\lambda L/2\) advances by at most \(\pi/2\).  Their endpoints are then`
- N143 `16` — line 362: `mapped to \(t\); each panel uses a locally generated 16-point Gauss--Legendre`
- N144 `4000` — line 363: `rule, the cutoff is \(\lambda_{\max}=4000\), and every nodal contribution is`
- N145 `3,182,304` — line 364: `combined by Kahan accumulation.  This requires 3,182,304 nodes at`
- N146 `0.08` — line 365: `\(\Fn=0.08\), rising to 50,916,864 at \(\Fn=0.02\).`
- N147 `50,916,864` — line 365: `\(\Fn=0.08\), rising to 50,916,864 at \(\Fn=0.02\).`
- N148 `0.02` — line 365: `\(\Fn=0.08\), rising to 50,916,864 at \(\Fn=0.02\).`
- N149 `3` — line 369: `weight is \(\sec^3\theta\,|F(\sec\theta)|^2\).  It maps the same`
- N150 `2` — line 369: `weight is \(\sec^3\theta\,|F(\sec\theta)|^2\).  It maps the same`
- N151 `8` — line 372: `map, panel-order, and cutoff convergence.  The 8/16 trend and 16/24 agreement`
- N152 `16` — line 372: `map, panel-order, and cutoff convergence.  The 8/16 trend and 16/24 agreement`
- N153 `16` — line 372: `map, panel-order, and cutoff convergence.  The 8/16 trend and 16/24 agreement`
- N154 `24` — line 372: `map, panel-order, and cutoff convergence.  The 8/16 trend and 16/24 agreement`
- N155 `16` — line 379: ```map'' compares the 16-node \(1+t^2\) and \(\sec\theta\) constructions.`
- N156 `1` — line 379: ```map'' compares the 16-node \(1+t^2\) and \(\sec\theta\) constructions.`
- N157 `2` — line 379: ```map'' compares the 16-node \(1+t^2\) and \(\sec\theta\) constructions.`
- N158 `16` — line 385: `\(\Fn\) & Nodes/map & Map 16/16 & \(t^2\) 8/16 & \(t^2\) 16/24 & \(\sec\theta\) 8/16 & \(\sec\theta\) 16/24\\`
- N159 `16` — line 385: `\(\Fn\) & Nodes/map & Map 16/16 & \(t^2\) 8/16 & \(t^2\) 16/24 & \(\sec\theta\) 8/16 & \(\sec\theta\) 16/24\\`
- N160 `2` — line 385: `\(\Fn\) & Nodes/map & Map 16/16 & \(t^2\) 8/16 & \(t^2\) 16/24 & \(\sec\theta\) 8/16 & \(\sec\theta\) 16/24\\`
- N161 `8` — line 385: `\(\Fn\) & Nodes/map & Map 16/16 & \(t^2\) 8/16 & \(t^2\) 16/24 & \(\sec\theta\) 8/16 & \(\sec\theta\) 16/24\\`
- N162 `16` — line 385: `\(\Fn\) & Nodes/map & Map 16/16 & \(t^2\) 8/16 & \(t^2\) 16/24 & \(\sec\theta\) 8/16 & \(\sec\theta\) 16/24\\`
- N163 `2` — line 385: `\(\Fn\) & Nodes/map & Map 16/16 & \(t^2\) 8/16 & \(t^2\) 16/24 & \(\sec\theta\) 8/16 & \(\sec\theta\) 16/24\\`
- N164 `16` — line 385: `\(\Fn\) & Nodes/map & Map 16/16 & \(t^2\) 8/16 & \(t^2\) 16/24 & \(\sec\theta\) 8/16 & \(\sec\theta\) 16/24\\`
- N165 `24` — line 385: `\(\Fn\) & Nodes/map & Map 16/16 & \(t^2\) 8/16 & \(t^2\) 16/24 & \(\sec\theta\) 8/16 & \(\sec\theta\) 16/24\\`
- N166 `8` — line 385: `\(\Fn\) & Nodes/map & Map 16/16 & \(t^2\) 8/16 & \(t^2\) 16/24 & \(\sec\theta\) 8/16 & \(\sec\theta\) 16/24\\`
- N167 `16` — line 385: `\(\Fn\) & Nodes/map & Map 16/16 & \(t^2\) 8/16 & \(t^2\) 16/24 & \(\sec\theta\) 8/16 & \(\sec\theta\) 16/24\\`
- N168 `16` — line 385: `\(\Fn\) & Nodes/map & Map 16/16 & \(t^2\) 8/16 & \(t^2\) 16/24 & \(\sec\theta\) 8/16 & \(\sec\theta\) 16/24\\`
- N169 `24` — line 385: `\(\Fn\) & Nodes/map & Map 16/16 & \(t^2\) 8/16 & \(t^2\) 16/24 & \(\sec\theta\) 8/16 & \(\sec\theta\) 16/24\\`
- N170 `0.08` — line 387: `0.08 & 3,182,304  & \num{1.284e-15} & \num{2.399e-10} & \num{1.426e-16} & \num{3.641e-10} & \num{0.000e0}\\`
- N171 `3,182,304` — line 387: `0.08 & 3,182,304  & \num{1.284e-15} & \num{2.399e-10} & \num{1.426e-16} & \num{3.641e-10} & \num{0.000e0}\\`
- N172 `1.284e-15` — line 387: `0.08 & 3,182,304  & \num{1.284e-15} & \num{2.399e-10} & \num{1.426e-16} & \num{3.641e-10} & \num{0.000e0}\\`
- N173 `2.399e-10` — line 387: `0.08 & 3,182,304  & \num{1.284e-15} & \num{2.399e-10} & \num{1.426e-16} & \num{3.641e-10} & \num{0.000e0}\\`
- N174 `1.426e-16` — line 387: `0.08 & 3,182,304  & \num{1.284e-15} & \num{2.399e-10} & \num{1.426e-16} & \num{3.641e-10} & \num{0.000e0}\\`
- N175 `3.641e-10` — line 387: `0.08 & 3,182,304  & \num{1.284e-15} & \num{2.399e-10} & \num{1.426e-16} & \num{3.641e-10} & \num{0.000e0}\\`
- N176 `0.000e0` — line 387: `0.08 & 3,182,304  & \num{1.284e-15} & \num{2.399e-10} & \num{1.426e-16} & \num{3.641e-10} & \num{0.000e0}\\`
- N177 `0.05` — line 388: `0.05 & 8,146,704  & \num{3.280e-16} & \num{1.849e-10} & \num{9.839e-16} & \num{1.564e-10} & \num{1.640e-16}\\`
- N178 `8,146,704` — line 388: `0.05 & 8,146,704  & \num{3.280e-16} & \num{1.849e-10} & \num{9.839e-16} & \num{1.564e-10} & \num{1.640e-16}\\`
- N179 `3.280e-16` — line 388: `0.05 & 8,146,704  & \num{3.280e-16} & \num{1.849e-10} & \num{9.839e-16} & \num{1.564e-10} & \num{1.640e-16}\\`
- N180 `1.849e-10` — line 388: `0.05 & 8,146,704  & \num{3.280e-16} & \num{1.849e-10} & \num{9.839e-16} & \num{1.564e-10} & \num{1.640e-16}\\`
- N181 `9.839e-16` — line 388: `0.05 & 8,146,704  & \num{3.280e-16} & \num{1.849e-10} & \num{9.839e-16} & \num{1.564e-10} & \num{1.640e-16}\\`
- N182 `1.564e-10` — line 388: `0.05 & 8,146,704  & \num{3.280e-16} & \num{1.849e-10} & \num{9.839e-16} & \num{1.564e-10} & \num{1.640e-16}\\`
- N183 `1.640e-16` — line 388: `0.05 & 8,146,704  & \num{3.280e-16} & \num{1.849e-10} & \num{9.839e-16} & \num{1.564e-10} & \num{1.640e-16}\\`
- N184 `0.03` — line 389: `0.03 & 22,629,712 & \num{4.664e-15} & \num{1.364e-10} & \num{2.120e-15} & \num{1.420e-10} & \num{1.484e-15}\\`
- N185 `22,629,712` — line 389: `0.03 & 22,629,712 & \num{4.664e-15} & \num{1.364e-10} & \num{2.120e-15} & \num{1.420e-10} & \num{1.484e-15}\\`
- N186 `4.664e-15` — line 389: `0.03 & 22,629,712 & \num{4.664e-15} & \num{1.364e-10} & \num{2.120e-15} & \num{1.420e-10} & \num{1.484e-15}\\`
- N187 `1.364e-10` — line 389: `0.03 & 22,629,712 & \num{4.664e-15} & \num{1.364e-10} & \num{2.120e-15} & \num{1.420e-10} & \num{1.484e-15}\\`
- N188 `2.120e-15` — line 389: `0.03 & 22,629,712 & \num{4.664e-15} & \num{1.364e-10} & \num{2.120e-15} & \num{1.420e-10} & \num{1.484e-15}\\`
- N189 `1.420e-10` — line 389: `0.03 & 22,629,712 & \num{4.664e-15} & \num{1.364e-10} & \num{2.120e-15} & \num{1.420e-10} & \num{1.484e-15}\\`
- N190 `1.484e-15` — line 389: `0.03 & 22,629,712 & \num{4.664e-15} & \num{1.364e-10} & \num{2.120e-15} & \num{1.420e-10} & \num{1.484e-15}\\`
- N191 `0.02` — line 390: `0.02 & 50,916,864 & \num{1.534e-15} & \num{8.035e-11} & \num{6.136e-15} & \num{8.284e-11} & \num{5.523e-15}\\`
- N192 `50,916,864` — line 390: `0.02 & 50,916,864 & \num{1.534e-15} & \num{8.035e-11} & \num{6.136e-15} & \num{8.284e-11} & \num{5.523e-15}\\`
- N193 `1.534e-15` — line 390: `0.02 & 50,916,864 & \num{1.534e-15} & \num{8.035e-11} & \num{6.136e-15} & \num{8.284e-11} & \num{5.523e-15}\\`
- N194 `8.035e-11` — line 390: `0.02 & 50,916,864 & \num{1.534e-15} & \num{8.035e-11} & \num{6.136e-15} & \num{8.284e-11} & \num{5.523e-15}\\`
- N195 `6.136e-15` — line 390: `0.02 & 50,916,864 & \num{1.534e-15} & \num{8.035e-11} & \num{6.136e-15} & \num{8.284e-11} & \num{5.523e-15}\\`
- N196 `8.284e-11` — line 390: `0.02 & 50,916,864 & \num{1.534e-15} & \num{8.035e-11} & \num{6.136e-15} & \num{8.284e-11} & \num{5.523e-15}\\`
- N197 `5.523e-15` — line 390: `0.02 & 50,916,864 & \num{1.534e-15} & \num{8.035e-11} & \num{6.136e-15} & \num{8.284e-11} & \num{5.523e-15}\\`
- N198 `0.02` — line 395: `At \(\Fn=0.02\), increasing \(\lambda_{\max}\) from 4000 to 8000 changes`
- N199 `4000` — line 395: `At \(\Fn=0.02\), increasing \(\lambda_{\max}\) from 4000 to 8000 changes`
- N200 `8000` — line 395: `At \(\Fn=0.02\), increasing \(\lambda_{\max}\) from 4000 to 8000 changes`
- N201 `1` — line 396: `the \(1+t^2\) and \(\sec\theta\) references by \(1.38\times10^{-15}\) and`
- N202 `2` — line 396: `the \(1+t^2\) and \(\sec\theta\) references by \(1.38\times10^{-15}\) and`
- N203 `1.38` — line 396: `the \(1+t^2\) and \(\sec\theta\) references by \(1.38\times10^{-15}\) and`
- N204 `15` — line 396: `the \(1+t^2\) and \(\sec\theta\) references by \(1.38\times10^{-15}\) and`
- N205 `1.23` — line 397: `\(1.23\times10^{-15}\), respectively.`
- N206 `15` — line 397: `\(1.23\times10^{-15}\), respectively.`
- N207 `1` — line 399: `As a published-value anchor, \citet[Table~1]{DoctorsBeck1987} report the`
- N208 `10` — line 400: `classical thin-ship value \(10^3 C_w=1.2486\) for \(B/L=0.1\),`
- N209 `3` — line 400: `classical thin-ship value \(10^3 C_w=1.2486\) for \(B/L=0.1\),`
- N210 `1.2486` — line 400: `classical thin-ship value \(10^3 C_w=1.2486\) for \(B/L=0.1\),`
- N211 `0.1` — line 400: `classical thin-ship value \(10^3 C_w=1.2486\) for \(B/L=0.1\),`
- N212 `0.0625` — line 401: `\(T/L=0.0625\), and \(\Fn=0.35\).  Our dimensional realization is the`
- N213 `0.35` — line 401: `\(T/L=0.0625\), and \(\Fn=0.35\).  Our dimensional realization is the`
- N214 `10` — line 402: `\SI{10}{m} hull above, with \(U=\Fn\sqrt{gL}\), \(\rho\) and \(g\) as stated`
- N215 `10` — line 404: `gives \(10^3 C_w=1.247922\), a relative difference of`
- N216 `3` — line 404: `gives \(10^3 C_w=1.247922\), a relative difference of`
- N217 `1.247922` — line 404: `gives \(10^3 C_w=1.247922\), a relative difference of`
- N218 `5.43` — line 405: `\(5.43\times10^{-4}\) (\(0.0543\%\)).  We classify this agreement as a known`
- N219 `4` — line 405: `\(5.43\times10^{-4}\) (\(0.0543\%\)).  We classify this agreement as a known`
- N220 `0.0543` — line 405: `\(5.43\times10^{-4}\) (\(0.0543\%\)).  We classify this agreement as a known`
- N221 `1` — line 406: `result reproduced (class 1).`
- N222 `10` — line 411: `The marcher requested \(\mathrm{rel\_tol}=10^{-8}\) with six refinements, but`
- N223 `8` — line 411: `The marcher requested \(\mathrm{rel\_tol}=10^{-8}\) with six refinements, but`
- N224 `03` — line 413: `converged.  For \(\Fn\leq0.03\), the physical omission correction is`
- N225 `6.2` — line 414: `\(6.2\times10^{-32}\) or smaller, so the reduced/reference difference measures`
- N226 `32` — line 414: `\(6.2\times10^{-32}\) or smaller, so the reduced/reference difference measures`
- N227 `0.08` — line 422: `0.08 & \num{1.945940872510e-1} & \num{1.172e-5} & \num{3.817e-5} & \num{6.918e-8} & \num{2.896e-7} & 56,032 & RC\\`
- N228 `1.945940872510e-1` — line 422: `0.08 & \num{1.945940872510e-1} & \num{1.172e-5} & \num{3.817e-5} & \num{6.918e-8} & \num{2.896e-7} & 56,032 & RC\\`
- N229 `1.172e-5` — line 422: `0.08 & \num{1.945940872510e-1} & \num{1.172e-5} & \num{3.817e-5} & \num{6.918e-8} & \num{2.896e-7} & 56,032 & RC\\`
- N230 `3.817e-5` — line 422: `0.08 & \num{1.945940872510e-1} & \num{1.172e-5} & \num{3.817e-5} & \num{6.918e-8} & \num{2.896e-7} & 56,032 & RC\\`
- N231 `6.918e-8` — line 422: `0.08 & \num{1.945940872510e-1} & \num{1.172e-5} & \num{3.817e-5} & \num{6.918e-8} & \num{2.896e-7} & 56,032 & RC\\`
- N232 `2.896e-7` — line 422: `0.08 & \num{1.945940872510e-1} & \num{1.172e-5} & \num{3.817e-5} & \num{6.918e-8} & \num{2.896e-7} & 56,032 & RC\\`
- N233 `56,032` — line 422: `0.08 & \num{1.945940872510e-1} & \num{1.172e-5} & \num{3.817e-5} & \num{6.918e-8} & \num{2.896e-7} & 56,032 & RC\\`
- N234 `0.05` — line 423: `0.05 & \num{1.057847519846e-2} & \num{1.963e-13}& \num{2.000e-8} & \num{1.493e-7} & \num{5.974e-7} & 119,936 & RC\\`
- N235 `1.057847519846e-2` — line 423: `0.05 & \num{1.057847519846e-2} & \num{1.963e-13}& \num{2.000e-8} & \num{1.493e-7} & \num{5.974e-7} & 119,936 & RC\\`
- N236 `1.963e-13` — line 423: `0.05 & \num{1.057847519846e-2} & \num{1.963e-13}& \num{2.000e-8} & \num{1.493e-7} & \num{5.974e-7} & 119,936 & RC\\`
- N237 `2.000e-8` — line 423: `0.05 & \num{1.057847519846e-2} & \num{1.963e-13}& \num{2.000e-8} & \num{1.493e-7} & \num{5.974e-7} & 119,936 & RC\\`
- N238 `1.493e-7` — line 423: `0.05 & \num{1.057847519846e-2} & \num{1.963e-13}& \num{2.000e-8} & \num{1.493e-7} & \num{5.974e-7} & 119,936 & RC\\`
- N239 `5.974e-7` — line 423: `0.05 & \num{1.057847519846e-2} & \num{1.963e-13}& \num{2.000e-8} & \num{1.493e-7} & \num{5.974e-7} & 119,936 & RC\\`
- N240 `119,936` — line 423: `0.05 & \num{1.057847519846e-2} & \num{1.963e-13}& \num{2.000e-8} & \num{1.493e-7} & \num{5.974e-7} & 119,936 & RC\\`
- N241 `0.03` — line 424: `0.03 & \num{5.113923599401e-4} & \num{1.045e-12}& \num{2.000e-8} & \num{3.322e-7} & \num{1.320e-6} & 267,856 & RC\\`
- N242 `5.113923599401e-4` — line 424: `0.03 & \num{5.113923599401e-4} & \num{1.045e-12}& \num{2.000e-8} & \num{3.322e-7} & \num{1.320e-6} & 267,856 & RC\\`
- N243 `1.045e-12` — line 424: `0.03 & \num{5.113923599401e-4} & \num{1.045e-12}& \num{2.000e-8} & \num{3.322e-7} & \num{1.320e-6} & 267,856 & RC\\`
- N244 `2.000e-8` — line 424: `0.03 & \num{5.113923599401e-4} & \num{1.045e-12}& \num{2.000e-8} & \num{3.322e-7} & \num{1.320e-6} & 267,856 & RC\\`
- N245 `3.322e-7` — line 424: `0.03 & \num{5.113923599401e-4} & \num{1.045e-12}& \num{2.000e-8} & \num{3.322e-7} & \num{1.320e-6} & 267,856 & RC\\`
- N246 `1.320e-6` — line 424: `0.03 & \num{5.113923599401e-4} & \num{1.045e-12}& \num{2.000e-8} & \num{3.322e-7} & \num{1.320e-6} & 267,856 & RC\\`
- N247 `267,856` — line 424: `0.03 & \num{5.113923599401e-4} & \num{1.045e-12}& \num{2.000e-8} & \num{3.322e-7} & \num{1.320e-6} & 267,856 & RC\\`
- N248 `0.02` — line 425: `0.02 & \num{4.417234459694e-5} & \num{1.829e-13}& \num{2.000e-8} & \num{6.334e-7} & \num{2.517e-6} & 511,504 & RC\\`
- N249 `4.417234459694e-5` — line 425: `0.02 & \num{4.417234459694e-5} & \num{1.829e-13}& \num{2.000e-8} & \num{6.334e-7} & \num{2.517e-6} & 511,504 & RC\\`
- N250 `1.829e-13` — line 425: `0.02 & \num{4.417234459694e-5} & \num{1.829e-13}& \num{2.000e-8} & \num{6.334e-7} & \num{2.517e-6} & 511,504 & RC\\`
- N251 `2.000e-8` — line 425: `0.02 & \num{4.417234459694e-5} & \num{1.829e-13}& \num{2.000e-8} & \num{6.334e-7} & \num{2.517e-6} & 511,504 & RC\\`
- N252 `6.334e-7` — line 425: `0.02 & \num{4.417234459694e-5} & \num{1.829e-13}& \num{2.000e-8} & \num{6.334e-7} & \num{2.517e-6} & 511,504 & RC\\`
- N253 `2.517e-6` — line 425: `0.02 & \num{4.417234459694e-5} & \num{1.829e-13}& \num{2.000e-8} & \num{6.334e-7} & \num{2.517e-6} & 511,504 & RC\\`
- N254 `511,504` — line 425: `0.02 & \num{4.417234459694e-5} & \num{1.829e-13}& \num{2.000e-8} & \num{6.334e-7} & \num{2.517e-6} & 511,504 & RC\\`
- N255 `0.08` — line 442: `0.08 & \num{3.815e-5} & \num{1.355e-12} & \num{9.153e-15} & \num{2.000e-8} & \num{3.817e-5}\\`
- N256 `3.815e-5` — line 442: `0.08 & \num{3.815e-5} & \num{1.355e-12} & \num{9.153e-15} & \num{2.000e-8} & \num{3.817e-5}\\`
- N257 `1.355e-12` — line 442: `0.08 & \num{3.815e-5} & \num{1.355e-12} & \num{9.153e-15} & \num{2.000e-8} & \num{3.817e-5}\\`
- N258 `9.153e-15` — line 442: `0.08 & \num{3.815e-5} & \num{1.355e-12} & \num{9.153e-15} & \num{2.000e-8} & \num{3.817e-5}\\`
- N259 `2.000e-8` — line 442: `0.08 & \num{3.815e-5} & \num{1.355e-12} & \num{9.153e-15} & \num{2.000e-8} & \num{3.817e-5}\\`
- N260 `3.817e-5` — line 442: `0.08 & \num{3.815e-5} & \num{1.355e-12} & \num{9.153e-15} & \num{2.000e-8} & \num{3.817e-5}\\`
- N261 `0.05` — line 443: `0.05 & \num{3.656e-12}& \num{7.222e-14} & \num{9.702e-15} & \num{2.000e-8} & \num{2.000e-8}\\`
- N262 `3.656e-12` — line 443: `0.05 & \num{3.656e-12}& \num{7.222e-14} & \num{9.702e-15} & \num{2.000e-8} & \num{2.000e-8}\\`
- N263 `7.222e-14` — line 443: `0.05 & \num{3.656e-12}& \num{7.222e-14} & \num{9.702e-15} & \num{2.000e-8} & \num{2.000e-8}\\`
- N264 `9.702e-15` — line 443: `0.05 & \num{3.656e-12}& \num{7.222e-14} & \num{9.702e-15} & \num{2.000e-8} & \num{2.000e-8}\\`
- N265 `2.000e-8` — line 443: `0.05 & \num{3.656e-12}& \num{7.222e-14} & \num{9.702e-15} & \num{2.000e-8} & \num{2.000e-8}\\`
- N266 `2.000e-8` — line 443: `0.05 & \num{3.656e-12}& \num{7.222e-14} & \num{9.702e-15} & \num{2.000e-8} & \num{2.000e-8}\\`
- N267 `0.03` — line 444: `0.03 & \num{6.153e-32}& \num{3.462e-12} & \num{9.292e-15} & \num{2.000e-8} & \num{2.000e-8}\\`
- N268 `6.153e-32` — line 444: `0.03 & \num{6.153e-32}& \num{3.462e-12} & \num{9.292e-15} & \num{2.000e-8} & \num{2.000e-8}\\`
- N269 `3.462e-12` — line 444: `0.03 & \num{6.153e-32}& \num{3.462e-12} & \num{9.292e-15} & \num{2.000e-8} & \num{2.000e-8}\\`
- N270 `9.292e-15` — line 444: `0.03 & \num{6.153e-32}& \num{3.462e-12} & \num{9.292e-15} & \num{2.000e-8} & \num{2.000e-8}\\`
- N271 `2.000e-8` — line 444: `0.03 & \num{6.153e-32}& \num{3.462e-12} & \num{9.292e-15} & \num{2.000e-8} & \num{2.000e-8}\\`
- N272 `2.000e-8` — line 444: `0.03 & \num{6.153e-32}& \num{3.462e-12} & \num{9.292e-15} & \num{2.000e-8} & \num{2.000e-8}\\`
- N273 `0.02` — line 445: `0.02 & \num{5.508e-70}& \num{7.620e-13} & \num{9.428e-15} & \num{2.000e-8} & \num{2.000e-8}\\`
- N274 `5.508e-70` — line 445: `0.02 & \num{5.508e-70}& \num{7.620e-13} & \num{9.428e-15} & \num{2.000e-8} & \num{2.000e-8}\\`
- N275 `7.620e-13` — line 445: `0.02 & \num{5.508e-70}& \num{7.620e-13} & \num{9.428e-15} & \num{2.000e-8} & \num{2.000e-8}\\`
- N276 `9.428e-15` — line 445: `0.02 & \num{5.508e-70}& \num{7.620e-13} & \num{9.428e-15} & \num{2.000e-8} & \num{2.000e-8}\\`
- N277 `2.000e-8` — line 445: `0.02 & \num{5.508e-70}& \num{7.620e-13} & \num{9.428e-15} & \num{2.000e-8} & \num{2.000e-8}\\`
- N278 `2.000e-8` — line 445: `0.02 & \num{5.508e-70}& \num{7.620e-13} & \num{9.428e-15} & \num{2.000e-8} & \num{2.000e-8}\\`
- N279 `0.08` — line 450: `\Cref{tab:accuracy} shows two regimes.  At \(\Fn=0.08\), submerged terms remain`
- N280 `05` — line 452: `\(\Fn\leq0.05\), reduced/reference differences range from`
- N281 `1.83` — line 453: `\(1.83\times10^{-13}\) to \(1.05\times10^{-12}\) and remain within the reduced`
- N282 `13` — line 453: `\(1.83\times10^{-13}\) to \(1.05\times10^{-12}\) and remain within the reduced`
- N283 `1.05` — line 453: `\(1.83\times10^{-13}\) to \(1.05\times10^{-12}\) and remain within the reduced`
- N284 `12` — line 453: `\(1.83\times10^{-13}\) to \(1.05\times10^{-12}\) and remain within the reduced`
- N285 `4` — line 461: `integration for \(s\in\{4,7,10,50,98,128\}\) and`
- N286 `7` — line 461: `integration for \(s\in\{4,7,10,50,98,128\}\) and`
- N287 `10` — line 461: `integration for \(s\in\{4,7,10,50,98,128\}\) and`
- N288 `50` — line 461: `integration for \(s\in\{4,7,10,50,98,128\}\) and`
- N289 `98,128` — line 461: `integration for \(s\in\{4,7,10,50,98,128\}\) and`
- N290 `25,100,400` — line 462: `\(|\omega|\in\{25,100,400\}\), with scaled discrepancy below`
- N291 `10` — line 463: `\(10^{-9}\).  The public degree envelope reaches at most \(s=98\); order 128`
- N292 `9` — line 463: `\(10^{-9}\).  The public degree envelope reaches at most \(s=98\); order 128`
- N293 `98` — line 463: `\(10^{-9}\).  The public degree envelope reaches at most \(s=98\); order 128`
- N294 `128` — line 463: `\(10^{-9}\).  The public degree envelope reaches at most \(s=98\); order 128`
- N295 `0.48` — line 468: `\begin{minipage}[t]{0.48\textwidth}`
- N296 `6` — line 472: `width=\linewidth,height=6.2cm,`
- N297 `20` — line 475: `grid=major,major grid style={black!20,densely dotted},`
- N298 `0.03` — line 477: `at={(0.03,0.03)},anchor=south west},`
- N299 `0.03` — line 477: `at={(0.03,0.03)},anchor=south west},`
- N300 `0.08` — line 481: `(0.08,1.172e-5) (0.05,1.963e-13) (0.03,1.045e-12) (0.02,1.829e-13)};`
- N301 `1.172e-5` — line 481: `(0.08,1.172e-5) (0.05,1.963e-13) (0.03,1.045e-12) (0.02,1.829e-13)};`
- N302 `0.05` — line 481: `(0.08,1.172e-5) (0.05,1.963e-13) (0.03,1.045e-12) (0.02,1.829e-13)};`
- N303 `1.963e-13` — line 481: `(0.08,1.172e-5) (0.05,1.963e-13) (0.03,1.045e-12) (0.02,1.829e-13)};`
- N304 `0.03` — line 481: `(0.08,1.172e-5) (0.05,1.963e-13) (0.03,1.045e-12) (0.02,1.829e-13)};`
- N305 `1.045e-12` — line 481: `(0.08,1.172e-5) (0.05,1.963e-13) (0.03,1.045e-12) (0.02,1.829e-13)};`
- N306 `0.02` — line 481: `(0.08,1.172e-5) (0.05,1.963e-13) (0.03,1.045e-12) (0.02,1.829e-13)};`
- N307 `1.829e-13` — line 481: `(0.08,1.172e-5) (0.05,1.963e-13) (0.03,1.045e-12) (0.02,1.829e-13)};`
- N308 `0.08` — line 484: `(0.08,6.918e-8) (0.05,1.493e-7) (0.03,3.322e-7) (0.02,6.334e-7)};`
- N309 `6.918e-8` — line 484: `(0.08,6.918e-8) (0.05,1.493e-7) (0.03,3.322e-7) (0.02,6.334e-7)};`
- N310 `0.05` — line 484: `(0.08,6.918e-8) (0.05,1.493e-7) (0.03,3.322e-7) (0.02,6.334e-7)};`
- N311 `1.493e-7` — line 484: `(0.08,6.918e-8) (0.05,1.493e-7) (0.03,3.322e-7) (0.02,6.334e-7)};`
- N312 `0.03` — line 484: `(0.08,6.918e-8) (0.05,1.493e-7) (0.03,3.322e-7) (0.02,6.334e-7)};`
- N313 `3.322e-7` — line 484: `(0.08,6.918e-8) (0.05,1.493e-7) (0.03,3.322e-7) (0.02,6.334e-7)};`
- N314 `0.02` — line 484: `(0.08,6.918e-8) (0.05,1.493e-7) (0.03,3.322e-7) (0.02,6.334e-7)};`
- N315 `6.334e-7` — line 484: `(0.08,6.918e-8) (0.05,1.493e-7) (0.03,3.322e-7) (0.02,6.334e-7)};`
- N316 `0.08` — line 487: `at (axis cs:0.08,1.172e-5) {rejected};`
- N317 `1.172e-5` — line 487: `at (axis cs:0.08,1.172e-5) {rejected};`
- N318 `0.035` — line 489: `at (axis cs:0.035,3e-5) {accepted samples};`
- N319 `3e-5` — line 489: `at (axis cs:0.035,3e-5) {accepted samples};`
- N320 `0.48` — line 493: `\begin{minipage}[t]{0.48\textwidth}`
- N321 `6` — line 497: `width=\linewidth,height=6.2cm,`
- N322 `20` — line 500: `grid=major,major grid style={black!20,densely dotted},`
- N323 `0.03` — line 502: `at={(0.03,0.97)},anchor=north west},`
- N324 `0.97` — line 502: `at={(0.03,0.97)},anchor=north west},`
- N325 `0.08` — line 506: `(0.08,1152) (0.05,1152) (0.03,1152) (0.02,1152)};`
- N326 `1152` — line 506: `(0.08,1152) (0.05,1152) (0.03,1152) (0.02,1152)};`
- N327 `0.05` — line 506: `(0.08,1152) (0.05,1152) (0.03,1152) (0.02,1152)};`
- N328 `1152` — line 506: `(0.08,1152) (0.05,1152) (0.03,1152) (0.02,1152)};`
- N329 `0.03` — line 506: `(0.08,1152) (0.05,1152) (0.03,1152) (0.02,1152)};`
- N330 `1152` — line 506: `(0.08,1152) (0.05,1152) (0.03,1152) (0.02,1152)};`
- N331 `0.02` — line 506: `(0.08,1152) (0.05,1152) (0.03,1152) (0.02,1152)};`
- N332 `1152` — line 506: `(0.08,1152) (0.05,1152) (0.03,1152) (0.02,1152)};`
- N333 `0.08` — line 509: `(0.08,56032) (0.05,119936) (0.03,267856) (0.02,511504)};`
- N334 `56032` — line 509: `(0.08,56032) (0.05,119936) (0.03,267856) (0.02,511504)};`
- N335 `0.05` — line 509: `(0.08,56032) (0.05,119936) (0.03,267856) (0.02,511504)};`
- N336 `119936` — line 509: `(0.08,56032) (0.05,119936) (0.03,267856) (0.02,511504)};`
- N337 `0.03` — line 509: `(0.08,56032) (0.05,119936) (0.03,267856) (0.02,511504)};`
- N338 `267856` — line 509: `(0.08,56032) (0.05,119936) (0.03,267856) (0.02,511504)};`
- N339 `0.02` — line 509: `(0.08,56032) (0.05,119936) (0.03,267856) (0.02,511504)};`
- N340 `511504` — line 509: `(0.08,56032) (0.05,119936) (0.03,267856) (0.02,511504)};`
- N341 `10` — line 515: `The two endpoint-aware references agree in the low \(10^{-15}\) range; the`
- N342 `15` — line 515: `The two endpoint-aware references agree in the low \(10^{-15}\) range; the`
- N343 `0.08` — line 519: `result at \(\Fn=0.08\) and accepts it at \(\Fn=0.05,0.03,0.02\).}`
- N344 `0.05` — line 519: `result at \(\Fn=0.08\) and accepts it at \(\Fn=0.05,0.03,0.02\).}`
- N345 `0.03` — line 519: `result at \(\Fn=0.08\) and accepts it at \(\Fn=0.05,0.03,0.02\).}`
- N346 `0.02` — line 519: `result at \(\Fn=0.08\) and accepts it at \(\Fn=0.05,0.03,0.02\).}`
- N347 `30` — line 525: `Timing used \texttt{cargo bench -p michell --bench wigley}, 30 samples per`
- N348 `10` — line 527: `default requested relative tolerance \(10^{-5}\).  Each case received an`
- N349 `5` — line 527: `default requested relative tolerance \(10^{-5}\).  Each case received an`
- N350 `0.1` — line 528: `unmeasured warm call.  Calls below \SI{0.1}{ms} were timed in batches of 256`
- N351 `256` — line 528: `unmeasured warm call.  Calls below \SI{0.1}{ms} were timed in batches of 256`
- N352 `256` — line 529: `and divided by 256; checksums use separate untimed calls, so timer and checksum`
- N353 `10` — line 532: `was a 10-core Apple M5 MacBook Air with 24 GB RAM, arm64 Darwin 25.5.0,`
- N354 `24` — line 532: `was a 10-core Apple M5 MacBook Air with 24 GB RAM, arm64 Darwin 25.5.0,`
- N355 `25.5` — line 532: `was a 10-core Apple M5 MacBook Air with 24 GB RAM, arm64 Darwin 25.5.0,`
- N356 `0` — line 532: `was a 10-core Apple M5 MacBook Air with 24 GB RAM, arm64 Darwin 25.5.0,`
- N357 `1.96` — line 533: `Rust/Cargo 1.96.0, and LLVM 22.1.2.  Absolute timings are host-dependent; work`
- N358 `0` — line 533: `Rust/Cargo 1.96.0, and LLVM 22.1.2.  Absolute timings are host-dependent; work`
- N359 `22.1` — line 533: `Rust/Cargo 1.96.0, and LLVM 22.1.2.  Absolute timings are host-dependent; work`
- N360 `2` — line 533: `Rust/Cargo 1.96.0, and LLVM 22.1.2.  Absolute timings are host-dependent; work`
- N361 `21` — line 544: `21 speeds, \(\Fn=0.10\ldots0.50\) & 13.059 & 12.952--13.126 & checksum \(2.553251156101\times10^5\)\\`
- N362 `0.10` — line 544: `21 speeds, \(\Fn=0.10\ldots0.50\) & 13.059 & 12.952--13.126 & checksum \(2.553251156101\times10^5\)\\`
- N363 `50` — line 544: `21 speeds, \(\Fn=0.10\ldots0.50\) & 13.059 & 12.952--13.126 & checksum \(2.553251156101\times10^5\)\\`
- N364 `13.059` — line 544: `21 speeds, \(\Fn=0.10\ldots0.50\) & 13.059 & 12.952--13.126 & checksum \(2.553251156101\times10^5\)\\`
- N365 `12.952` — line 544: `21 speeds, \(\Fn=0.10\ldots0.50\) & 13.059 & 12.952--13.126 & checksum \(2.553251156101\times10^5\)\\`
- N366 `13.126` — line 544: `21 speeds, \(\Fn=0.10\ldots0.50\) & 13.059 & 12.952--13.126 & checksum \(2.553251156101\times10^5\)\\`
- N367 `2.553251156101` — line 544: `21 speeds, \(\Fn=0.10\ldots0.50\) & 13.059 & 12.952--13.126 & checksum \(2.553251156101\times10^5\)\\`
- N368 `5` — line 544: `21 speeds, \(\Fn=0.10\ldots0.50\) & 13.059 & 12.952--13.126 & checksum \(2.553251156101\times10^5\)\\`
- N369 `0.05` — line 545: `Default API, \(\Fn=0.05\) & 0.030 & 0.030--0.030 & 1,152 nodes; batch 256\\`
- N370 `0.030` — line 545: `Default API, \(\Fn=0.05\) & 0.030 & 0.030--0.030 & 1,152 nodes; batch 256\\`
- N371 `0.030` — line 545: `Default API, \(\Fn=0.05\) & 0.030 & 0.030--0.030 & 1,152 nodes; batch 256\\`
- N372 `0.030` — line 545: `Default API, \(\Fn=0.05\) & 0.030 & 0.030--0.030 & 1,152 nodes; batch 256\\`
- N373 `1,152` — line 545: `Default API, \(\Fn=0.05\) & 0.030 & 0.030--0.030 & 1,152 nodes; batch 256\\`
- N374 `256` — line 545: `Default API, \(\Fn=0.05\) & 0.030 & 0.030--0.030 & 1,152 nodes; batch 256\\`
- N375 `0.02` — line 546: `Real-axis marcher, \(\Fn=0.02\) & 21.541 & 21.311--21.666 & 511,504 evaluations; batch 1\\`
- N376 `21.541` — line 546: `Real-axis marcher, \(\Fn=0.02\) & 21.541 & 21.311--21.666 & 511,504 evaluations; batch 1\\`
- N377 `21.311` — line 546: `Real-axis marcher, \(\Fn=0.02\) & 21.541 & 21.311--21.666 & 511,504 evaluations; batch 1\\`
- N378 `21.666` — line 546: `Real-axis marcher, \(\Fn=0.02\) & 21.541 & 21.311--21.666 & 511,504 evaluations; batch 1\\`
- N379 `511,504` — line 546: `Real-axis marcher, \(\Fn=0.02\) & 21.541 & 21.311--21.666 & 511,504 evaluations; batch 1\\`
- N380 `1` — line 546: `Real-axis marcher, \(\Fn=0.02\) & 21.541 & 21.311--21.666 & 511,504 evaluations; batch 1\\`
- N381 `0.02` — line 547: `Endpoint/NSD, \(\Fn=0.02\) & 0.030 & 0.030--0.030 & 1,152 nodes; batch 256\\`
- N382 `0.030` — line 547: `Endpoint/NSD, \(\Fn=0.02\) & 0.030 & 0.030--0.030 & 1,152 nodes; batch 256\\`
- N383 `0.030` — line 547: `Endpoint/NSD, \(\Fn=0.02\) & 0.030 & 0.030--0.030 & 1,152 nodes; batch 256\\`
- N384 `0.030` — line 547: `Endpoint/NSD, \(\Fn=0.02\) & 0.030 & 0.030--0.030 & 1,152 nodes; batch 256\\`
- N385 `1,152` — line 547: `Endpoint/NSD, \(\Fn=0.02\) & 0.030 & 0.030--0.030 & 1,152 nodes; batch 256\\`
- N386 `256` — line 547: `Endpoint/NSD, \(\Fn=0.02\) & 0.030 & 0.030--0.030 & 1,152 nodes; batch 256\\`
- N387 `0.02` — line 552: `At \(\Fn=0.02\), the frozen paired run gives a 718-fold median speedup; repeated`
- N388 `718` — line 552: `At \(\Fn=0.02\), the frozen paired run gives a 718-fold median speedup; repeated`
- N389 `608` — line 553: `protocols on this host have spanned roughly 608--865-fold, so the ratio is not`
- N390 `865` — line 553: `protocols on this host have spanned roughly 608--865-fold, so the ratio is not`
- N391 `1.83` — line 556: `\(1.83\times10^{-13}\), whereas the capped marcher's is`
- N392 `13` — line 556: `\(1.83\times10^{-13}\), whereas the capped marcher's is`
- N393 `6.33` — line 557: `\(6.33\times10^{-7}\).  Cost changes little between \(\Fn=0.05\) and 0.02, as`
- N394 `7` — line 557: `\(6.33\times10^{-7}\).  Cost changes little between \(\Fn=0.05\) and 0.02, as`
- N395 `0.05` — line 557: `\(6.33\times10^{-7}\).  Cost changes little between \(\Fn=0.05\) and 0.02, as`
- N396 `0.02` — line 557: `\(6.33\times10^{-7}\).  Cost changes little between \(\Fn=0.05\) and 0.02, as`
- N397 `21` — line 558: `\cref{eq:gaussian} predicts.  The sweep grid is the 21 speeds`
- N398 `0.10` — line 559: `\(\Fn=0.10,0.12,\ldots,0.50\); its checksum is the sum of the 21 returned`
- N399 `0.12` — line 559: `\(\Fn=0.10,0.12,\ldots,0.50\); its checksum is the sum of the 21 returned`
- N400 `0.50` — line 559: `\(\Fn=0.10,0.12,\ldots,0.50\); its checksum is the sum of the 21 returned`
- N401 `21` — line 559: `\(\Fn=0.10,0.12,\ldots,0.50\); its checksum is the sum of the 21 returned`
- N402 `1972` — line 571: `reductions from the implementation work reported here.  For the 1972 and 1988`
- N403 `1988` — line 571: `reductions from the implementation work reported here.  For the 1972 and 1988`
- N404 `0.15` — line 583: `\begin{tabular}{@{}p{0.15\textwidth}p{0.20\textwidth}p{0.20\textwidth}p{0.17\textwidth}p{0.17\textwidth}@{}}`
- N405 `0.20` — line 583: `\begin{tabular}{@{}p{0.15\textwidth}p{0.20\textwidth}p{0.20\textwidth}p{0.17\textwidth}p{0.17\textwidth}@{}}`
- N406 `0.20` — line 583: `\begin{tabular}{@{}p{0.15\textwidth}p{0.20\textwidth}p{0.20\textwidth}p{0.17\textwidth}p{0.17\textwidth}@{}}`
- N407 `0.17` — line 583: `\begin{tabular}{@{}p{0.15\textwidth}p{0.20\textwidth}p{0.20\textwidth}p{0.17\textwidth}p{0.17\textwidth}@{}}`
- N408 `0.17` — line 583: `\begin{tabular}{@{}p{0.15\textwidth}p{0.20\textwidth}p{0.20\textwidth}p{0.17\textwidth}p{0.17\textwidth}@{}}`
- N409 `39` — line 587: `Birkhoff--Kotik via \citet[eqs.~39--45]{Wehausen1973}`
- N410 `45` — line 587: `Birkhoff--Kotik via \citet[eqs.~39--45]{Wehausen1973}`
- N411 `0` — line 589: `& Geometry-independent kernel \(K(\nu u,\nu v)\); a \(Y_0\) reduction is also reported`
- N412 `3.3` — line 593: `\citet[eqs.~3.3--3.4]{Michelsen1960}`
- N413 `3.4` — line 593: `\citet[eqs.~3.3--3.4]{Michelsen1960}`
- N414 `3.21` — line 599: `\citet[eqs.~3.21--3.22, 3.38--3.40]{Michelsen1960}`
- N415 `3.22` — line 599: `\citet[eqs.~3.21--3.22, 3.38--3.40]{Michelsen1960}`
- N416 `3.38` — line 599: `\citet[eqs.~3.21--3.22, 3.38--3.40]{Michelsen1960}`
- N417 `3.40` — line 599: `\citet[eqs.~3.21--3.22, 3.38--3.40]{Michelsen1960}`
- N418 `1` — line 644: `1--16 and charges coefficient construction and endpoint cancellation;`
- N419 `16` — line 644: `1--16 and charges coefficient construction and endpoint cancellation;`
- N420 `83` — line 659: `through bow and stern waterline slopes.  \citet[pp.~83--85]{Gotman2002}`
- N421 `85` — line 659: `through bow and stern waterline slopes.  \citet[pp.~83--85]{Gotman2002}`
- N422 `6.5` — line 668: `\citet[sec.~6.5, pp.~6-10--6-12]{Lazauskas2009} gives exact`
- N423 `6` — line 668: `\citet[sec.~6.5, pp.~6-10--6-12]{Lazauskas2009} gives exact`
- N424 `10` — line 668: `\citet[sec.~6.5, pp.~6-10--6-12]{Lazauskas2009} gives exact`
- N425 `6` — line 668: `\citet[sec.~6.5, pp.~6-10--6-12]{Lazauskas2009} gives exact`
- N426 `12` — line 668: `\citet[sec.~6.5, pp.~6-10--6-12]{Lazauskas2009} gives exact`
- N427 `1` — line 701: `single hulls of degree 1--16 in each spline direction.  Its active`
- N428 `16` — line 701: `single hulls of degree 1--16 in each spline direction.  Its active`
- N429 `48` — line 705: `\item \textbf{Contour degree.}  The fixed 48-node rule was accurate through`
- N430 `400` — line 706: `\(s\simeq50\) and degraded beyond it in scans to \(s=400\); its`
- N431 `24` — line 707: `24/48 difference nevertheless remained conservative against the`
- N432 `48` — line 707: `24/48 difference nevertheless remained conservative against the`
- N433 `192` — line 708: `192-node same-contour comparator throughout that scan.  Production`
- N434 `48` — line 709: `therefore switches to 48/96 when \(s/|\omega|\geq2\), and the claimed`
- N435 `96` — line 709: `therefore switches to 48/96 when \(s/|\omega|\geq2\), and the claimed`
- N436 `98` — line 710: `envelope ends at the independently tested public maximum \(s=98\).`
- N437 `1` — line 714: `\item \textbf{Reference floor.}  The independent \(1+t^2\) and`
- N438 `2` — line 714: `\item \textbf{Reference floor.}  The independent \(1+t^2\) and`
- N439 `4000` — line 716: `\(\lambda=4000\).  At \(\Fn=0.02\), their cutoff-doubling changes are`
- N440 `0.02` — line 716: `\(\lambda=4000\).  At \(\Fn=0.02\), their cutoff-doubling changes are`
- N441 `1.38` — line 717: `\(1.38\times10^{-15}\) and \(1.23\times10^{-15}\), respectively.`
- N442 `15` — line 717: `\(1.38\times10^{-15}\) and \(1.23\times10^{-15}\), respectively.`
- N443 `1.23` — line 717: `\(1.38\times10^{-15}\) and \(1.23\times10^{-15}\), respectively.`
- N444 `15` — line 717: `\(1.38\times10^{-15}\) and \(1.23\times10^{-15}\), respectively.`
- N445 `10` — line 718: `Differences near the observed \(10^{-13}\) scale are resolution`
- N446 `13` — line 718: `Differences near the observed \(10^{-13}\) scale are resolution`
- N447 `2` — line 731: `\(\nu y\lambda\sqrt{\lambda^2-1}\) and moving saddles that`
- N448 `1` — line 731: `\(\nu y\lambda\sqrt{\lambda^2-1}\) and moving saddles that`
- N449 `0.02` — line 748: `For the standard Wigley hull at \(\Fn=0.02\), the reduced solver is more`
- N450 `700` — line 749: `accurate than the general marcher and about 700 times faster in the frozen`

### Units

- U001 `\SI{21.541}{ms}` — line 24: `\SI{21.541}{ms} to \SI{0.030}{ms} (about 700-fold), and differs from an`
- U002 `\SI{0.030}{ms}` — line 24: `\SI{21.541}{ms} to \SI{0.030}{ms} (about 700-fold), and differs from an`
- U003 `\SI{9.80665}{m.s^{-2}` — line 118: `All calculations below use \(g=\SI{9.80665}{m.s^{-2}}\), freshwater density`
- U004 `\SI{999.1}{kg.m^{-3}` — line 119: `\(\rho=\SI{999.1}{kg.m^{-3}}\), and the length Froude number`
- U005 `metres` — line 123: `speeds, and forces are in metres, metres per second, and newtons; \(\lambda\),`
- U006 `metres` — line 123: `speeds, and forces are in metres, metres per second, and newtons; \(\lambda\),`
- U007 `newtons` — line 123: `speeds, and forces are in metres, metres per second, and newtons; \(\lambda\),`
- U008 `\SI{10}{m}` — line 350: `with \(L=\SI{10}{m}\), \(B=\SI{1}{m}\), and`
- U009 `\SI{1}{m}` — line 350: `with \(L=\SI{10}{m}\), \(B=\SI{1}{m}\), and`
- U010 `\SI{0.625}{m}` — line 351: `\(T=\SI{0.625}{m}\).  The primary reference shares neither the production`
- U011 `\SI{10}{m}` — line 402: `\SI{10}{m} hull above, with \(U=\Fn\sqrt{gL}\), \(\rho\) and \(g\) as stated`
- U012 `\SI{0.1}{ms}` — line 528: `unmeasured warm call.  Calls below \SI{0.1}{ms} were timed in batches of 256`
- U013 `GB` — line 532: `was a 10-core Apple M5 MacBook Air with 24 GB RAM, arm64 Darwin 25.5.0,`
- U014 `RAM` — line 532: `was a 10-core Apple M5 MacBook Air with 24 GB RAM, arm64 Darwin 25.5.0,`
- U015 `ms` — line 542: `Case & Median (ms) & IQR (ms) & Diagnostic\\`
- U016 `ms` — line 542: `Case & Median (ms) & IQR (ms) & Diagnostic\\`

## Frozen equation manifest

### E01 `eq:michell`

Environment `align`; 311 bytes; SHA-256 `13df23d84fa2cb756cec8a82ca70d8fa7ee1fb61b5f6158669b0b1152ee163b1`.

```tex
\begin{align}
  \Rw &=\frac{4\rho g^2}{\pi U^2}
  \int_1^\infty \frac{\lambda^2}{\sqrt{\lambda^2-1}}
  \left|F(\lambda)\right|^2\,d\lambda,
  \label{eq:michell}\\
  F(\lambda) &= \iint_\Omega f_x(x,z)
  \e^{-\nu\lambda^2 z}\e^{\iu\nu\lambda x}\,dx\,dz,
  \qquad \nu=\frac{g}{U^2}.
  \label{eq:inner}
\end{align}
```

### E02 `eq:cw`

Environment `equation`; 139 bytes; SHA-256 `df87a9775412690f0654c625ac44db8408b2740103a549e93bc58a346043cace`.

```tex
\begin{equation}
 C_w=\frac{\Rw}{\tfrac12\rho U^2 S_w},\qquad
 S_w=2\iint_\Omega\sqrt{1+f_x^2+f_z^2}\,dx\,dz,
 \label{eq:cw}
\end{equation}
```

### E03 `eq:spanpoly`

Environment `equation`; 125 bytes; SHA-256 `2c75e50c177a901f96b60fa283a2f5ded4c491e11864ae8ba71f6a7c2fb4fecf`.

```tex
\begin{equation}
  f_x|_S=\sum_{a=0}^{p-1}\sum_{b=0}^{q}
  c^{S}_{ab}(x-x_0)^a(z-z_0)^b.
  \label{eq:spanpoly}
\end{equation}
```

### E04 `eq:xparts`

Environment `align`; 336 bytes; SHA-256 `84dfe41cbef9386cc80765046db7f2da604d85e8e784cb1def82f92da114fcba`.

```tex
\begin{align}
 \int_{x_a}^{x_b} p(x)\e^{\iu kx}\,dx
 &=\sum_{r=0}^{m}\frac{(-1)^r}{(\iu k)^{r+1}}
   \left[p^{(r)}(x)\e^{\iu kx}\right]_{x_a}^{x_b},
 \label{eq:xparts}\\
 \int_{z_a}^{z_b} p(z)\e^{-\kappa z}\,dz
 &=-\sum_{u=0}^{m}\frac{1}{\kappa^{u+1}}
   \left[p^{(u)}(z)\e^{-\kappa z}\right]_{z_a}^{z_b}.
 \label{eq:zparts}
\end{align}
```

### E05 `eq:endpointsum`

Environment `equation`; 164 bytes; SHA-256 `16a3d1aac21de04dfff98dfc833ee4283f05801d85e93c4fcd44023f04967456`.

```tex
\begin{equation}
 F(\lambda)=\sum_{j=1}^{M} c_j\lambda^{-n_j}
 \e^{\iu\nu\lambda x_j}\e^{-\nu\lambda^2z_j},
 \qquad n_j\geq3,
 \label{eq:endpointsum}
\end{equation}
```

### E06 `eq:waterline`

Environment `equation`; 132 bytes; SHA-256 `3538052bbee6bfe04aad761c94dd1647718a1c65324524d8b51245f0c0807171`.

```tex
\begin{equation}
 \widetilde F(\lambda)=\sum_{j\in W}c_j\lambda^{-n_j}
 \e^{\iu\nu\lambda x_j}.
 \label{eq:waterline}
\end{equation}
```

### E07 `eq:pairresistance`

Environment `equation`; 203 bytes; SHA-256 `626f5186957425ce5763c42eca0e2b2c5b33dc962d8414c310bd1c02d82f66b9`.

```tex
\begin{equation}
 \widetilde R_{\mathrm w}=C
 \sum_{i,j\in W}\Re\!\left\{
 c_i\overline{c_j}K_{s_{ij}}(\omega_{ij})\right\},
 \quad
 C=\frac{4\rho g^2}{\pi U^2},
 \label{eq:pairresistance}
\end{equation}
```

### E08 `eq:kernel`

Environment `equation`; 149 bytes; SHA-256 `b3ff75e31edbc0f7bff254fd71b0785d32b15c209696cb1414c05188009061f1`.

```tex
\begin{equation}
 K_s(\omega)=\int_1^\infty
 \frac{\e^{\iu\omega\lambda}}{
 \lambda^s\sqrt{\lambda^2-1}}\,d\lambda.
 \label{eq:kernel}
\end{equation}
```

### E09 `eq:omitbound`

Environment `equation`; 174 bytes; SHA-256 `1543d07538d7d3002f5c2ccdf6d8112cdfce7d94ad5a3b8c9c58cb8f37370532`.

```tex
\begin{equation}
 \left|\Rw-\widetilde R_{\mathrm w}\right|
 \leq C\sum_{(i,j)\in\mathcal P_S}
 |c_i c_j|\e^{-\nu(z_i+z_j)}K_{s_{ij}}(0).
 \label{eq:omitbound}
\end{equation}
```

### E10 `eq:kzero`

Environment `equation`; 137 bytes; SHA-256 `f459864e3b623a7f333babffd56a9a93f280db3bec1fa5f111bf120a5711c1c8`.

```tex
\begin{equation}
 K_s(0)=\int_0^\infty\cosh^{-s}t\,dt
 = \frac{\sqrt\pi\,\Gamma(s/2)}{2\Gamma((s+1)/2)},
 \label{eq:kzero}
\end{equation}
```

### E11 `eq:bickley`

Environment `equation`; 137 bytes; SHA-256 `7de646e06e97e3bc88317a4268e09369e70284517c416cf2d64296fd2c2798c9`.

```tex
\begin{equation}
 K_s(\omega)=\int_0^\infty
 \e^{\iu\omega\cosh t}\cosh^{-s}t\,dt
 =\Ki_s(-\iu\omega),
 \label{eq:bickley}
\end{equation}
```

### E12 `eq:tform`

Environment `equation`; 148 bytes; SHA-256 `1aa049d5f87e5bd5d6e59bda7294bbf56392e8dd5381606334a0bf5d05cdc40d`.

```tex
\begin{equation}
 K_s(\omega)=2\e^{\iu\omega}\int_0^\infty
 \frac{\e^{\iu\omega t^2}}{
 (1+t^2)^s\sqrt{2+t^2}}\,dt.
 \label{eq:tform}
\end{equation}
```

### E13 `eq:gaussian`

Environment `equation`; 193 bytes; SHA-256 `f9b374b0680b02c71c6b176ee22fd6030552c45f469c0dbc9d34e0b36715c575`.

```tex
\begin{equation}
 K_s(\omega)=\frac{2\e^{\iu(\omega+\pi/4)}}{\sqrt\omega}
 \int_0^\infty \frac{\e^{-r^2}\,dr}{
 (1+\iu r^2/\omega)^s\sqrt{2+\iu r^2/\omega}}.
 \label{eq:gaussian}
\end{equation}
```

### E14 `eq:gausstail`

Environment `equation`; 127 bytes; SHA-256 `56e140dee2e00b2a937f0093b7ff462e149faa7f0e92f63b4683a769941084f6`.

```tex
\begin{equation}
 |K_s-K_s^{(8)}|\leq
 \frac{\e^{-64}}{8\sqrt{2\omega}},
 \qquad \omega>0,
 \label{eq:gausstail}
\end{equation}
```

### E15 `eq:coeffround`

Environment `equation`; 132 bytes; SHA-256 `ea2af592dc39e1e79206cd238a7aeae9d78e8327ca484552bf1832440b6fbe71`.

```tex
\begin{equation}
 \delta_i=\gamma_{2a_i+2}\sum_k|q_{ik}|,\qquad
 \gamma_n=\frac{n\,u}{1-n\,u},
 \label{eq:coeffround}
\end{equation}
```

### E16 `eq:bo`

Environment `align`; 327 bytes; SHA-256 `79fd8fbe097b6dba92898dda665bafc0094db63a3605ecf365eae592300b9e10`.

```tex
\begin{align}
 B_o &= C\!\sum_{(i,j)\in\mathcal P_S}
 |c_i c_j|\e^{-\nu(z_i+z_j)}K_{s_{ij}}(0),
 \label{eq:bo}\\
 B_q &= C\!\sum_{i\leq j}m_{ij}|c_i c_j|
 |K^{f}_{ij}-K^{c}_{ij}|,
 \label{eq:bq}\\
 B_r &= C\!\sum_{i,j}
 (\delta_i|c_j|+|c_i|\delta_j+\delta_i\delta_j)
 \e^{-\nu(z_i+z_j)}K_{s_{ij}}(0).
 \label{eq:br}
\end{align}
```

### E17 `eq:floor`

Environment `equation`; 99 bytes; SHA-256 `f530441d1aeeba6cc57b840089d9f74998fd496acd7a156fcfd2d7a9444ca229`.

```tex
\begin{equation}
 \epsilon_c=\max(2\times10^{-8},\,4r_{\rm span}),
 \label{eq:floor}
\end{equation}
```

### E18 `eq:dispatcherror`

Environment `equation`; 124 bytes; SHA-256 `4774dc686fdb71b17cca1240daeee56566f5e80f965c8a9e7e847c703976bfb4`.

```tex
\begin{equation}
 \widehat\epsilon=\epsilon_c+
 \frac{B}{|\widetilde R_{\rm w}|-B}.
 \label{eq:dispatcherror}
\end{equation}
```

### E19 `eq:spangate`

Environment `equation`; 77 bytes; SHA-256 `8384ee52c6a44fbe52e5dcf4930f08f0eb1a2b02241ad47eed049933db514986`.

```tex
\begin{equation}
 m\leq\frac{1}{25\Fn^2}.
 \label{eq:spangate}
\end{equation}
```

### E20 `eq:wigley`

Environment `equation`; 174 bytes; SHA-256 `1c1c1356176c45f69e0d7d93e15f7734080e94b8cea9ae314dfc2e713e0289fc`.

```tex
\begin{equation}
 f(x,z)=\frac{B}{2}\left(1-\frac{4x^2}{L^2}\right)
 \left(1-\frac{z^2}{T^2}\right),
 \quad |x|\leq L/2,\quad 0\leq z\leq T,
 \label{eq:wigley}
\end{equation}
```

## Citation-support manifest

- C01 `Michell1898`, `Tuck1989` — Michell's 1898 thin-ship theory expresses wave resistance as a one-dimensional integral of the squared Fourier--Laplace transform of the longitudinal hull slope \citep{Michell1898,Tuck1989}.
- C02 `Lazauskas2009`, `DambrinePierreRousseaux2016` — Designers therefore use it for preliminary design, multihull studies, and mathematical shape optimization \citep{Lazauskas2009,DambrinePierreRousseaux2016}.
- C03 `Michelsen1960` — \citet{Michelsen1960}, building on Birkhoff and Kotik's transformation, splits the calculation into a hull function and a speed-dependent Michell function, derives a convergent special-function series for polynomial hull functions, and proposes tabulation for systematic design.
- C04 `Michelsen1963` — Michelsen next treated polynomial centerline singularity distributions \citep{Michelsen1963} and high- and low-speed asymptotic approximations \citep{Michelsen1966}.
- C05 `Michelsen1966` — Michelsen next treated polynomial centerline singularity distributions \citep{Michelsen1963} and high- and low-speed asymptotic approximations \citep{Michelsen1966}.
- C06 `Michelsen1972` — \citet{Michelsen1972} later uses Gegenbauer source distributions and orthogonality to reduce the integral to a finite double sum.
- C07 `SendagortaGrases1988` — The verified abstract of \citet{SendagortaGrases1988} likewise describes rapidly convergent series, products of shape functions, and shape-independent velocity functions for tabulation and computer-aided design.
- C08 `Tuck1989`, `Lazauskas2009` — Tuck used Filon-like longitudinal integration, while Lazauskas used exact piecewise-quadratic formulas and a fixed angular rule; Lazauskas reports that very low Froude number remains exceptional \citep{Tuck1989,Lazauskas2009}.
- C09 `Wehausen1973` — \citet{Wehausen1973} surveys classical thin-ship theory and its limits.
- C10 `KellerAhluwalia1976` — Keller and Ahluwalia show that bow and stern waterline data control the small-Froude wave field and resistance \citep{KellerAhluwalia1976}.
- C11 `Gotman2002` — Gotman gives finite endpoint-derivative sums and products that separate bow, stern, and bow--stern interaction terms \citep[pp.~83--85]{Gotman2002}.
- C12 `HuybrechsVandewalle2006` — Huybrechs and Vandewalle develop numerical steepest descent \citep{HuybrechsVandewalle2006}, and Motygin applies it to the Kelvin Green-function integral in linear ship-wave theory \citep{Motygin2017}.
- C13 `Motygin2017` — Huybrechs and Vandewalle develop numerical steepest descent \citep{HuybrechsVandewalle2006}, and Motygin applies it to the Kelvin Green-function integral in linear ship-wave theory \citep{Motygin2017}.
- C14 `Michelsen1972` — The accessible publisher and index records establish the statements above for \citet{Michelsen1972} and \citet{SendagortaGrases1988}; their full texts remain due-diligence items.
- C15 `SendagortaGrases1988` — The accessible publisher and index records establish the statements above for \citet{Michelsen1972} and \citet{SendagortaGrases1988}; their full texts remain due-diligence items.
- C16 `BickleyNayler1935`, `DLMF1043` — With \(\lambda=\cosh t\), \cref{eq:kernel} becomes \begin{equation} K_s(\omega)=\int_0^\infty \e^{\iu\omega\cosh t}\cosh^{-s}t\,dt =\Ki_s(-\iu\omega), \label{eq:bickley} \end{equation} where the final equality denotes analytic continuation of the Bickley function \citep{BickleyNayler1935,DLMF1043}.\footnote{DLMF calls \(\Ki_s\) the ``Bickley function'' and cites the 1935 paper by W.~G.
- C17 `RuffaToni2026` — The compound spelling ``Bickley--Naylor'' occurs in later literature, including \citet{RuffaToni2026}; it is not DLMF terminology.} This identity supplies recurrences and other possible implementations.
- C18 `GibbsEtAl2020` — The same numerical-asymptotic principle yields frequency-independent work in high-frequency scattering \citep{GibbsEtAl2020}; here, the quadratic phase gives the contour in closed form.
- C19 `Tuck1989` — This explicitly removes the square-root singularity at \(\lambda=1\); \citet[p.~371]{Tuck1989} identifies that singularity and recommends its removal before quadrature.
- C20 `DoctorsBeck1987` — As a published-value anchor, \citet[Table~1]{DoctorsBeck1987} report the classical thin-ship value \(10^3 C_w=1.2486\) for \(B/L=0.1\), \(T/L=0.0625\), and \(\Fn=0.35\).
- C21 `Wehausen1973` — \midrule Birkhoff--Kotik via \citet[eqs.~39--45]{Wehausen1973} & Hull autocorrelation \(M(u,v)\), with equivalent \((P,Q)\) transforms and separable elementary ships & Geometry-independent kernel \(K(\nu u,\nu v)\); a \(Y_0\) reduction is also reported & Absolute convergence justifies the change of integration order & Establishes separation of hull data from a reusable kernel
- C22 `Michelsen1960` — \addlinespace \citet[eqs.~3.3--3.4]{Michelsen1960} & Polynomial hull function \(H(\xi,\zeta)\) after the Birkhoff--Kotik transformation & Michell function \(C(s,t)\) contains speed dependence & Convergence conditions accompany the transformation & Direct historical antecedent of shape/kernel separation
- C23 `Michelsen1960` — \addlinespace \citet[eqs.~3.21--3.22, 3.38--3.40]{Michelsen1960} & Monomial expansion of \(H\); each coefficient contributes through a closed series expression & Confluent-hypergeometric, Bessel, and Struve functions; speed-and-degree terms can be tabulated & Appendix proves series convergence; computer tabulation is proposed & Antecedent of analytic polynomial reduction for design use
- C24 `Michelsen1963` — \addlinespace \citet{Michelsen1963}, program record & Polynomial centerline singularity distributions & Direct evaluation of their wave-resistance contribution & The verified seminar program establishes topic and venue; detailed controls await source access & Connects the dissertation's polynomial program to centerplane distributions
- C25 `Michelsen1966` — \addlinespace \citet{Michelsen1966}, bibliographic record & High- and low-speed asymptotic approximations to Michell's integral & Speed-limit formulas rather than a general finite-speed kernel & The verified journal record establishes scope; detailed remainder control awaits source access & Extends the lineage explicitly into both speed limits
- C26 `Michelsen1972` — \addlinespace \citet{Michelsen1972}, verified record & Gegenbauer expansion of the centerplane Havelock-source distribution & Orthogonality reduces resistance to a finite double sum depending on source coefficients, \(\Fn\), and length--draft ratio & The accessible abstract states the reduction and identifies its governing inputs & Published JSR continuation of the basis-reduction program
- C27 `SendagortaGrases1988` — \addlinespace \citet{SendagortaGrases1988}, verified record & Products of integral shape functions & Rapidly convergent series with linear, shape-independent velocity functions suitable for tabulation & The abstract reports that few terms give suitable accuracy; detailed controls await source access & Antecedent of separated, table-driven computer-aided design
- C28 `KellerAhluwalia1976` — \citet{KellerAhluwalia1976} express leading small-Froude resistance and waves through bow and stern waterline slopes.
- C29 `Gotman2002` — \citet[pp.~83--85]{Gotman2002} derives finite endpoint-derivative sums, forms their products in resistance, and separates bow, stern, and bow--stern contributions; this is closer to the present pair structure than a generic integration-by-parts citation.
- C30 `Lazauskas2009` — In particular, \citet[sec.~6.5, pp.~6-10--6-12]{Lazauskas2009} gives exact piecewise-quadratic hull integrals, compares equally spaced Simpson and cosine-spaced trapezoidal angular rules, and reports that very low Froude number remains exceptional.
- C31 `RuffaToni2026` — Ruffa and Toni give finite Bessel--Struve module representations for integer-order Bickley functions \citep{RuffaToni2026}; their stability on the imaginary axis remains to be tested as a possible kernel backend.
- C32 `HuybrechsVandewalle2006` — Contour deformation is well established for oscillatory quadrature \citep{HuybrechsVandewalle2006}.
- C33 `Motygin2017` — \citet{Motygin2017} applies steepest descent and Clenshaw--Curtis quadrature in ship-wave theory, but to the oscillatory part of the Kelvin Green function rather than Michell's resistance integral.
- C34 `GibbsHewettHuybrechs2024` — Automated methods now handle general phases and coalescing saddles \citep{GibbsHewettHuybrechs2024}.
- C35 `TrinhChapman2015` — For free-surface flow past singular obstacles, low-Froude waves can be exponentially small, beyond all algebraic orders, and controlled by Stokes phenomena; that different problem is not resolved by better quadrature \citep{TrinhChapman2015}.

## Claim-strength map

Each line paraphrases one frozen claim or qualification before editing.
The final audit maps each identifier to its revised section and confirms unchanged strength.

- K01 [Abstract] Low-Froude Michell integration is increasingly oscillatory and has a difficult algebraic tail. — After: pending register pass.
- K02 [Abstract] Michelsen developed analytic polynomial and orthogonal-basis reductions in the cited 1960 and 1972 works. — After: pending register pass.
- K03 [Abstract] The present implementation extends that lineage to tensor-product B-spline hulls within the stated degree range. — After: pending register pass.
- K04 [Abstract] Endpoint integration, bounded submerged-term omission, and contour-evaluated Bickley kernels produce speed-independent reduced work. — After: pending register pass.
- K05 [Abstract] The combined estimate controls selection; direct real-axis integration is used when the requested tolerance is not met. — After: pending register pass.
- K06 [Abstract] The stated Wigley work, timing, agreement, geometry scope, and non-certified contour qualification all hold together. — After: pending register pass.
- K07 [Introduction] Michell theory is inexpensive, preserves bow-stern interference, and remains useful in the cited design applications. — After: pending register pass.
- K08 [Introduction] At low speed, computational speed is useful only with a trustworthy error estimate. — After: pending register pass.
- K09 [Introduction] Michelsen's dissertation separated hull and speed functions and proposed convergent, tabulated polynomial reductions. — After: pending register pass.
- K10 [Introduction] Michelsen's later records cover polynomial centerline distributions, speed-limit asymptotics, and a finite Gegenbauer double sum. — After: pending register pass.
- K11 [Introduction] The verified Sendagorta-Grases abstract establishes separated, rapidly convergent shape and velocity functions for design. — After: pending register pass.
- K12 [Introduction] The present contribution is implementation engineering with measured error accounting, not method priority. — After: pending register pass.
- K13 [Introduction] Falling Froude number drives increasing real-axis oscillation and makes adaptive truncation expensive and delicate. — After: pending register pass.
- K14 [Introduction] Tuck and Lazauskas used the stated inner and angular treatments, while Lazauskas identified very low Froude number as exceptional. — After: pending register pass.
- K15 [Introduction] The cited pre-fix quiet-window estimate under-covered measured error by the stated factor. — After: pending register pass.
- K16 [Introduction] The hardened direct calculation adds phase-based windows, refinement, and tail estimates, but all four reported comparators hit the refinement cap. — After: pending register pass.
- K17 [Introduction] Wehausen, Keller-Ahluwalia, Gotman, Huybrechs-Vandewalle, and Motygin supply the stated theoretical ingredients. — After: pending register pass.
- K18 [Introduction] The four stated implementation deltas are exact decomposition, omission accounting, contour kernels, and tolerance-based selection. — After: pending register pass.
- K19 [Introduction] The inaccessible Michelsen and Sendagorta-Grases full texts remain due-diligence items; claims rely only on verified records. — After: pending register pass.
- K20 [Michell resistance] The coordinates, Michell normalization, physical constants, Froude convention, resistance coefficient, wetted area, and units are as defined. — After: pending register pass.
- K21 [Michell resistance] Each nonzero B-spline knot rectangle has an exact polynomial longitudinal derivative. — After: pending register pass.
- K22 [Michell resistance] The method obtains coefficients from exact derivatives, drops zero-length repeated-knot intervals, and represents a full-multiplicity chine exactly. — After: pending register pass.
- K23 [Endpoint reduction] Repeated integration by parts terminates for finite polynomial degree and yields the stated exact endpoint representation. — After: pending register pass.
- K24 [Endpoint reduction] Endpoint coefficients are finite derivative combinations and equal endpoint-power triples can be combined exactly. — After: pending register pass.
- K25 [Endpoint reduction] Independent moment tests meet the stated scaled discrepancy over the stated lambda interval for Wigley and chine geometries. — After: pending register pass.
- K26 [Error bound] Submerged endpoint waves are exponentially suppressed at low Froude number, permitting the stated waterline reduction. — After: pending register pass.
- K27 [Error bound] Pair expansion produces the stated reusable kernel representation. — After: pending register pass.
- K28 [Error bound] The submerged-pair proposition is an absolute resistance bound proved by ordered-pair expansion and the triangle inequality. — After: pending register pass.
- K29 [Error bound] The zero-frequency kernel has the stated analytic form and stable even/odd recurrence. — After: pending register pass.
- K30 [Contour evaluation] The pair kernel is an analytically continued Bickley function with the historical naming qualification in the footnote. — After: pending register pass.
- K31 [Contour evaluation] The endpoint substitution and exact contour rotation replace oscillation with Gaussian decay for nonzero frequency. — After: pending register pass.
- K32 [Contour evaluation] The deformation has no intervening poles or branch points, and its closing arc vanishes for the stated kernel orders. — After: pending register pass.
- K33 [Contour evaluation] The implementation uses the stated ordinary and stiff Gauss-Legendre rules and has the stated analytic Gaussian-tail bound. — After: pending register pass.
- K34 [Contour evaluation] The coarse/fine difference is empirical; only the omission and contour-tail bounds have the stated rigorous status. — After: pending register pass.
- K35 [Contour evaluation] The fixed-rule scan has the stated range, errors, coverage factor, and degree-dependent interpretation. — After: pending register pass.
- K36 [Contour evaluation] Nonzero-pair work is frequency independent, with the stated Wigley endpoint, pair, evaluation, and node counts. — After: pending register pass.
- K37 [Method selection] The reduced calculation is considered only for the stated physical, spline, endpoint-spacing, and error conditions. — After: pending register pass.
- K38 [Method selection] The coefficient, omission, contour, rounding, construction-floor, denominator, and refusal calculations are exactly those stated. — After: pending register pass.
- K39 [Method selection] The reverse-triangle denominator is necessary; the total estimate combines analytic bounds with empirical components and is not an interval certificate. — After: pending register pass.
- K40 [Method selection] The frequency condition implies the stated Froude and equal-span restrictions. — After: pending register pass.
- K41 [Method selection] The multi-span node arithmetic, whole-hull refusal, and default Wigley decisions hold only as qualified in the text. — After: pending register pass.
- K42 [Method selection] The implementation and validation tests have the stated language, dependency, and property-test coverage. — After: pending register pass.
- K43 [Numerical results] The primary Wigley reference uses the stated dimensions, endpoint regularization, phase panels, order, cutoff, node counts, and Kahan accumulation. — After: pending register pass.
- K44 [Numerical results] The independently constructed secant reference has the stated weight and separate construction. — After: pending register pass.
- K45 [Numerical results] The map, order, and cutoff studies support only the reported digits and do not constitute interval bounds. — After: pending register pass.
- K46 [Numerical results] The Doctors-Beck comparison uses the stated nondimensionalization and reproduces the published value within the stated difference. — After: pending register pass.
- K47 [Numerical results] Every direct-integration row reached RefinementCap, and the lowest-Froude reduced differences measure reference quadrature rather than omitted physics. — After: pending register pass.
- K48 [Numerical results] The estimate decomposition uses the guarded denominator and sums the stated components. — After: pending register pass.
- K49 [Numerical results] The reduced method is rejected at the highest tabulated Froude number and agrees within its estimates below it; the capped direct error grows as speed falls. — After: pending register pass.
- K50 [Numerical results] The endpoint-aware references, rather than the capped direct calculation, are the comparison standard. — After: pending register pass.
- K51 [Numerical results] Kernel tests cover the stated orders and frequencies; the claimed degree range ends at the stated maximum. — After: pending register pass.
- K52 [Numerical results] The timing protocol, batching, alternation, hardware, toolchain, deterministic work counts, and host dependence are exactly qualified as stated. — After: pending register pass.
- K53 [Numerical results] The measured speed ratio is not portable and compares equal requested tolerance at unequal observed accuracy. — After: pending register pass.
- K54 [Numerical results] The speed grid and literal checksums are defined as untimed sums, and the selected path follows the acceptance condition. — After: pending register pass.
- K55 [Earlier work] The lineage table distinguishes inspected equations from verified records and leaves inaccessible controls unresolved. — After: pending register pass.
- K56 [Earlier work] The table attributes the stated basis, kernel, convergence, and design-use facts to each historical source. — After: pending register pass.
- K57 [Earlier work] The five implementation differences remain engineering deltas without a priority or head-to-head performance claim. — After: pending register pass.
- K58 [Earlier work] Keller-Ahluwalia and Gotman establish the stated endpoint physics and finite derivative-product structure. — After: pending register pass.
- K59 [Earlier work] Tuck and Lazauskas establish the stated piecewise-polynomial practices, but the accessible thesis does not establish every Michlet internal. — After: pending register pass.
- K60 [Earlier work] Ruffa-Toni offer a possible Bickley backend whose imaginary-axis stability remains untested. — After: pending register pass.
- K61 [Earlier work] The cited contour literature covers numerical steepest descent, Kelvin-wave integration, and more general phases; multihull phases remain future work here. — After: pending register pass.
- K62 [Discussion] The evidence is confined to upright symmetric monohulls in linear deep-water thin-ship theory. — After: pending register pass.
- K63 [Discussion] The geometry limitation retains the exact degree, spacing, and excluded-configuration scope. — After: pending register pass.
- K64 [Discussion] The contour limitation retains the measured fixed-rule behavior, adaptive rule, and maximum-order envelope. — After: pending register pass.
- K65 [Discussion] The analytic bounds and empirical or non-interval-certified error components remain distinguished. — After: pending register pass.
- K66 [Discussion] The independent references retain their summation, cutoff, convergence, and resolution-floor qualifications. — After: pending register pass.
- K67 [Discussion] Michell theory retains its linear, inviscid, slender, deep-water, fixed-attitude limitations and does not resolve exponentially small nonlinear wave phenomena. — After: pending register pass.
- K68 [Discussion] Multihull and shallow-endpoint phases require different contours; the listed certification, recurrence, finite-depth, and differentiation extensions remain future work. — After: pending register pass.
- K69 [Conclusions] Endpoint reduction, pair kernels, contour rotation, and the submerged-pair bound give geometry-controlled accepted work. — After: pending register pass.
- K70 [Conclusions] The stated Wigley advantage is confined to the supported geometry and does not justify claims beyond the present error control. — After: pending register pass.
- K71 [Reproducibility] The repository contents, placeholder DOI process, immutable-tag identification, and no-move rule remain unchanged. — After: pending register pass.
- K72 [Acknowledgments] The funding, AI assistance, author review, independent checks, and author responsibility disclosure remain unchanged. — After: pending register pass.

## Build and visual audit

The final audit records the Tectonic diagnostics, PDF page count, and rendered-page inspection.

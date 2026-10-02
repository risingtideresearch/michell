#!/usr/bin/env node
// Golden outputs from propopt's web cores, for this crate's port to match.
//
//     node crates/propeller/tests/golden/make.js ../propopt
//
// Runs propcore.js's sweep and motorcore.js's motorScatter on a few
// configurations and writes each answer, as JSON, beside this script. The
// motor database is the one vendored in crates/propeller/data, so the port
// and the original read the same records.

const fs = require('fs');
const path = require('path');
const propopt = process.argv[2] || path.join(__dirname, '..', '..', '..', '..', '..', 'propopt');
const web = path.join(propopt, 'web');
global.MOTOR_DB = require(path.join(__dirname, '..', '..', 'data', 'motors.js')).MOTOR_DB;
eval(fs.readFileSync(path.join(web, 'propcore.js'), 'utf8'));
eval(fs.readFileSync(path.join(web, 'motorcore.js'), 'utf8'));

const KN = 0.514444, IN = 0.0254;
const CASES = {
  // The web page's defaults.
  web_default: {
    opts: { V_s: 8 * KN, T: 1000, D_min: 0.04, D_max: 16 * IN, Z: [2, 3, 4, 5, 6, 7],
            shafts: 1, depth: 0.3, kellerK: 0.2, wake: 0, thrustDeduction: 0,
            cavitation: true, reCorrect: true, strictEAR: false },
    powerCap: 1.6,
  },
  // The defaults with the second operating point on.
  top_speed: {
    opts: { V_s: 8 * KN, T: 1000, D_min: 0.04, D_max: 16 * IN, Z: [2, 3, 4, 5, 6, 7],
            shafts: 1, depth: 0.3, kellerK: 0.2, wake: 0, thrustDeduction: 0,
            cavitation: true, reCorrect: true, strictEAR: false,
            V_top: 12 * KN, T_top: 2250 },
    powerCap: 1.6,
  },
  // A small twin: two shafts, some wake, a lighter cavitation margin.
  small_twin: {
    opts: { V_s: 4.0, T: 230, D_min: 0.04, D_max: 0.30, Z: [3, 4], shafts: 2,
            depth: 0.25, kellerK: 0.1, wake: 0.1, thrustDeduction: 0.05,
            cavitation: true, reCorrect: true, strictEAR: true },
    powerCap: 1.6,
  },
};

// The curve's points, less the circular perZ links: each point keeps its
// blade counts' answers as plain summaries.
function point(p) {
  if (!p.ok) return { rpm: p.rpm, ok: false };
  const out = {};
  for (const [k, v] of Object.entries(p)) {
    if (k === 'perZ') {
      out.perZ = v.map((z) => ({ Z: z.Z, D: z.D, PD: z.PD, EAR: z.EAR, eta0: z.eta0,
                                 P_shaft: z.P_shaft }));
    } else out[k] = v;
  }
  return out;
}

function dot(d) {
  return { motor: d.motor.id, rpm: d.rpm, P_shaft: d.P_shaft, P_elec: d.P_elec, eta: d.eta,
           etaDrive: d.etaDrive, ratio: d.ratio, motorRpm: d.motorRpm, motorQ: d.motorQ,
           V_bus: d.V_bus, I_dc: d.I_dc, I_arms: d.I_arms, load: d.load,
           extrapolated: d.extrapolated, top: d.top };
}

for (const [name, c] of Object.entries(CASES)) {
  const cfg = makeConfig(c.opts);
  const res = sweep(cfg, { powerCap: c.powerCap });
  const mo = { rating: 'cont', gearEta: 0.97, // motorcore's GEAR_ETA_DEFAULT
               allWindings: false,
               controllerEta: MOTOR_DB.controller_eta || 1, rankBy: 'power',
               checkTop: !!cfg.top };
  const sc = motorScatter(res.points, mo);
  const all = motorScatter(res.points, Object.assign({}, mo, { allWindings: true }));
  const out = {
    opts: c.opts, powerCap: c.powerCap, config: cfg,
    sweep: {
      feasible: res.feasible, rpmLo: res.rpmLo, rpmHi: res.rpmHi,
      feasLo: res.feasLo, feasHi: res.feasHi,
      best: res.best ? point(res.best) : null,
      bestUnconstrained: res.bestUnconstrained ? point(res.bestUnconstrained) : null,
      topOkCount: res.topOkCount, topFeasible: res.topFeasible,
      points: res.points.map(point),
    },
    motors: { dots: sc.dots.map(dot), unreachable: sc.unreachable.map((m) => m.id),
              all: all.dots.map(dot) },
  };
  const file = path.join(__dirname, name + '.json');
  fs.writeFileSync(file, JSON.stringify(out, null, 1) + '\n');
  console.log(`${name}: best ${res.best ? res.best.P_shaft.toFixed(1) + ' W at ' + res.best.rpm.toFixed(0) + ' rpm, Z ' + res.best.Z : 'none'}; ${sc.dots.length} motor families, ${sc.unreachable.length} unreachable`);
}

# Vendored motor database

`motors.js` is copied unchanged from
[propopt](https://github.com/risingtideresearch/propopt) `web/motors.js`, at
commit `636720c` (2026-09-19). It's generated there by
`tools/make_motor_db.py`, from the manufacturers' datasheets under
`data/sources/` plus hand specs in `tools/motor_specs.py` and the map
digitisations in `data/labels/`. Don't edit it here. To update it, regenerate
it in propopt and copy it over, then rerun this crate's golden tests.

The file is a JS literal, `const MOTOR_DB = {…};`, whose object is plain JSON.
`propeller::motor::Database` reads that object straight out of it.

It holds 159 motor records (EMRAX, Oswald, Zapi, Golden Motor, Flipsky, Beyond,
YASA, Dana, Maytech) plus two digitised efficiency maps. Each record's `tier`
says how far its efficiency is to be trusted:

| tier | efficiency from |
|---|---|
| `A` | the manufacturer's published map, digitised |
| `A-` | the frame's published map, corrected to this winding's copper loss |
| `A*` | the published map read by hand |
| `B` | a loss model fitted to the datasheet (and any published points) |
| `C` | none published |

# Running boatmath-web on our own server

boatmath-web keeps its hulls, cases, studies and results in one
directory (a SQLite database and a folder of gzipped blobs) and works through
its queue on the machine it runs on. It has no login of its own: it listens
on loopback only, and Tailscale decides who can reach it.

```
tailnet ──https──▶ tailscale serve ──http──▶ 127.0.0.1:8080 boatmath-web ──▶ /var/lib/boatmath
```

Behind `tailscale serve` each request carries the tailnet user
(`Tailscale-User-Login`, `Tailscale-User-Name`), and uploads and studies are
labelled with it. boatmath-web believes those headers only when it is bound
to loopback, where nothing but the local proxy can reach it.

## Build

On the server, from a git checkout (the build stamps results with the last
commit to touch the solver, so build from a clean tree):

```sh
git clone https://github.com/risingtideresearch/michell.git && cd michell
cargo build --release -p boatmath-web
sudo install -m 755 target/release/boatmath-web /usr/local/bin/boatmath-web
```

The solver uses every core; a machine with many is worth it. Nothing else is
needed at run time (SQLite is built in). The pages load three.js from
jsDelivr, so browsers need to reach it.

## Install the service

```sh
sudo useradd --system --home /var/lib/boatmath --shell /usr/sbin/nologin boatmath
sudo install -m 644 deploy/boatmath-web.service /etc/systemd/system/
sudo systemctl daemon-reload
sudo systemctl enable --now boatmath-web
journalctl -u boatmath-web -f        # "serving http://127.0.0.1:8080/ (identity from tailscale serve's headers)"
```

## Put it on the tailnet

```sh
sudo tailscale serve --bg --https=443 http://127.0.0.1:8080
tailscale serve status              # https://<machine>.<tailnet>.ts.net/
```

Who may reach it is the tailnet's policy: grant the team access to this
machine's port 443 (and nothing else it serves) in the ACLs. Use
`tailscale serve`, not `funnel`: funnel would put it on the internet.

A device tagged rather than owned by a user has no user login; its requests
are labelled with whatever name the page sends.

## Upgrade

```sh
git pull && cargo build --release -p boatmath-web
sudo install -m 755 target/release/boatmath-web /usr/local/bin/boatmath-web
sudo systemctl restart boatmath-web
```

The queue survives a restart: a study cut off part-way goes back on the queue
and starts again. If the solver's code changed, the results already saved are
marked stale (the queue page offers to re-run them); they are kept until
they are.

## Back up

Everything is in `/var/lib/boatmath`. The blobs never change once written, so
copy the database consistently and the blobs as files:

```sh
sqlite3 /var/lib/boatmath/boatmath.db ".backup '/backups/boatmath.db'"
rsync -a /var/lib/boatmath/blobs/ /backups/blobs/
```

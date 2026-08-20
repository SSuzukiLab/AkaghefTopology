# Topology.jl static site

Dependency-free Cloudflare Pages projection of the Julia-generated SVG artifacts.
The Julia package and its artifacts remain canonical; `dist/` is disposable output.

## Build

```sh
python3 site/build.py
```

The builder copies the three explicitly selected SVGs in `site/build.py` from
`julia/artifacts/julia_plots/` into `site/dist/` and creates a small gallery.
Add artifacts deliberately rather than publishing the full workbench by default.

## Preview

```sh
python3 -m http.server 8765 --directory site/dist
```

## Deploy

```sh
npx wrangler pages deploy site/dist --project-name akaghef-topology
```

No Node package, JavaScript framework, database, or Julia server is required by the
published site.

# Evo Genesis and Sticks

A small static site hosting two browser games:

- EvoGenesis (simulation sandbox)
- Stick Gladiators (physics arena)

This repo contains a hub at `evo-genesis-and-sticks/index.html` that links to each game. For now the game pages redirect to the existing single-file HTMLs in the project root (so you can keep developing them in place). Later, you can copy the full HTML content into the `games/*/index.html` files if you prefer a clean structure.

## Local development
Just open the HTML files in a browser. No build step is required.

- EvoGenesis: `EvoGenesis.html`
- Stick Gladiators: `Stick Gladiators.html`

Hub page: `evo-genesis-and-sticks/index.html`

## GitHub Pages
A workflow `.github/workflows/pages.yml` deploys the `evo-genesis-and-sticks` folder as the site. After it runs on `main`, set the Pages source to "GitHub Actions" in repository settings if not already.

Site entry: `https://<your-username>.github.io/<repo>/`

Replace the placeholder social image at `evo-genesis-and-sticks/assets/og-image.png`.

## License
MIT © 2025 Santiago Garcia

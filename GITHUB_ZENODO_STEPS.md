# How to publish this code on GitHub and get a Zenodo DOI

Claude prepared this repository but cannot log in to your GitHub/Zenodo accounts,
so the final publishing steps are done by you. It takes ~10 minutes.

## 0. First, edit the placeholders
Fill in author name(s) in: `LICENSE`, `CITATION.cff`, `.zenodo.json`.

## 1. Put the code on GitHub
Option A — website (no command line):
1. Sign in at https://github.com (create a free account if needed).
2. Click **New repository** → name e.g. `icca-spp1-tam-mycaf-niche` → **Public** →
   do NOT add a README (this repo already has one) → **Create repository**.
3. On the new repo page click **uploading an existing file** → drag in ALL files/folders
   from this repository → **Commit changes**.

Option B — git command line (if installed):
```bash
cd repo_iCCA_SPP1_niche
git init && git add . && git commit -m "Initial release: iCCA SPP1+ TAM–myCAF niche analysis code"
git branch -M main
git remote add origin https://github.com/<your-username>/<repo-name>.git
git push -u origin main
```

## 2. Mint a DOI with Zenodo (GitHub integration)
1. Go to https://zenodo.org and **Log in with GitHub** (authorize when asked — this is you granting access, not Claude).
2. Top-right menu → **GitHub**. You will see your repositories.
3. Flip the toggle **ON** next to your new repository. (Zenodo now watches it for releases.)
4. Back on GitHub → your repo → **Releases** → **Create a new release** →
   tag `v1.0.0`, title `v1.0.0`, → **Publish release**.
5. Zenodo automatically archives that release and issues a **DOI** within a minute.
   Find it at https://zenodo.org → **Upload** (or the DOI badge on your repo).

## 3. Use the DOI
- Add to the manuscript **Code Availability**: “Analysis code is available at
  https://github.com/<user>/<repo> and archived at Zenodo (DOI: 10.5281/zenodo.XXXXXXX).”
- Zenodo gives a **concept DOI** (always latest) and a **version DOI** (this release) —
  cite the concept DOI in the paper.

## For peer review (keep anonymous)
- Use an anonymized view: https://anonymous.4open.science (upload the repo), or
- Share the repo privately and note “code available to reviewers on request”, or
- Zenodo lets you keep a record **restricted** and share a secret reviewer link.

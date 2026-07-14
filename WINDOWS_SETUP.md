# Finalizing traceAsthma on Windows

This is the concrete, ordered checklist to close the remaining gaps from
`vignette("clinical-deployment-notes")` and get the package onto GitHub +
r-universe. Everything here runs on your Windows desktop, in RStudio.

## 1. One-time setup

**Install R and RStudio** if you don't already have them:
- R: https://cran.r-project.org/bin/windows/base/
- RStudio: https://posit.co/download/rstudio-desktop/

**Install Rtools** (required to compile several dependencies from source on
Windows — glmnet, xgboost, arrow, and others have compiled components):
- https://cran.r-project.org/bin/windows/Rtools/
- Pick the Rtools version matching your R version (e.g. R 4.4.x → Rtools44)
- After installing, restart RStudio

**Install Git for Windows** if you don't have it:
- https://git-scm.com/download/win
- During setup, the defaults are fine; just make sure "Git from the command
  line and also from 3rd-party software" is selected

## 2. Unpack the package

Unzip `traceAsthma_source.zip` somewhere on your machine, e.g.
`C:\Users\<you>\Documents\traceAsthma\`. Open RStudio, then:

```r
setwd("C:/Users/<you>/Documents/traceAsthma")
```

(Use forward slashes even on Windows — R accepts them.)

## 3. Run the finalization script

```r
source("finalize_package.R")
```

This installs every dependency, regenerates documentation, runs the real
test suite, builds the real vignettes, and runs a full `R CMD check --as-cran`.
Expect 15–30 minutes the first time. Read the console output at the end —
it will tell you plainly whether the check is clean.

**If something fails**, the error message will name the specific function
and file. Common first-run issues on Windows:
- A Bioconductor package fails to install → re-run just that line,
  e.g. `BiocManager::install("minfi", update = FALSE, ask = FALSE)`,
  and read the actual error (sometimes it's a transient CRAN/Bioconductor
  mirror timeout — just retry)
- `pkgbuild::has_rtools()` returns `FALSE` → Rtools isn't on PATH; restart
  RStudio after installing it, or run `Sys.which("make")` to check

## 4. Push to GitHub

Create the empty repository first at https://github.com/new (name:
`traceAsthma`, owner: your account `Eskezeia`, leave it empty — no README,
no license, no .gitignore, since your local folder already has all of those).

Then, in RStudio's **Terminal** tab (not the R console):

```bash
cd "C:/Users/<you>/Documents/traceAsthma"
git init
git add .
git commit -m "Initial commit: traceAsthma v0.3.0"
git branch -M main
git remote add origin https://github.com/Eskezeia/traceAsthma.git
git push -u origin main
```

You'll be prompted to authenticate — the easiest path is a GitHub Personal
Access Token used as your password (GitHub → Settings → Developer settings
→ Personal access tokens → generate one with `repo` scope), since GitHub
retired plain password authentication for git operations.

## 5. Confirm CI runs

Once pushed, go to `https://github.com/Eskezeia/traceAsthma/actions` — the
`R-CMD-check.yaml` workflow (already included in `.github/workflows/`)
should start automatically and run `R CMD check --as-cran` on Windows,
macOS, and two Linux R versions. This is the actual multi-platform
verification that wasn't possible in the sandboxed environment this
package was authored in — green checkmarks here are the real confirmation.

## 6. Register on r-universe (makes it `install.packages()`-able)

Once GitHub Actions is green:
1. Go to https://github.com/r-universe-dev and follow the "get started"
   instructions, or more directly: create a repo named exactly
   `Eskezeia.r-universe.dev` containing a `packages.json` listing
   `traceAsthma` (r-universe's docs walk through this — it's a few lines
   of JSON, not a review process).
2. Once live, anyone can run:
   ```r
   install.packages("traceAsthma", repos = "https://eskezeia.r-universe.dev")
   ```

## 7. Only after all of the above: consider CRAN

If steps 1–6 are clean and you want CRAN specifically:
```r
devtools::release()
```
walks through CRAN's own pre-submission checklist interactively (it will
ask about `cran-comments.md`, re-run checks, and open the actual submission
form). Reserve a name check first: search
https://cran.r-project.org/web/packages/available_packages_by_name.html
for `traceAsthma` to confirm it's still unclaimed before you invest time here.

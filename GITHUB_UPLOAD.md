# Upload this project to GitHub

1. Create a new empty repository at https://github.com/new. Choose a repository name such as `ecoli-wgs-virulence-pipeline`. Do not initialize it with a README because this folder already contains one.
2. From this folder, run:

```bash
git init
git branch -M main
git add README.md GITHUB_UPLOAD.md REFERENCE_NOTES.md LICENSE-REVIEW.txt .gitattributes .gitignore environment.yml run_ecoli_pipeline.sh build_presence_absence.py make_stx_matrix.py markers.fna docs tests assembly

git status --short
git diff --cached --stat
```

3. Inspect the staged files. Ensure that no customer genomes, private data, credentials, or unauthorized material is staged.
4. Commit and connect the remote (replace `YOUR_USERNAME` with your GitHub username):

```bash
git commit -m "Initial release of E. coli WGS pipeline"
git remote add origin https://github.com/YOUR_USERNAME/ecoli-wgs-virulence-pipeline.git
git push -u origin main
```

If `origin` already exists, inspect it with `git remote -v` instead of adding it again.

**Before making the repository public:** review `REFERENCE_NOTES.md` and `LICENSE-REVIEW.txt`. The marker sequence provenance and redistribution rights still need to be confirmed. Do not claim the workflow has been end-to-end validated until you have run it with known test assemblies and checked the expected results.

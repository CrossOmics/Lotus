# Deploy Documentation to GitHub Pages

This guide explains how to deploy the Lotus documentation to GitHub Pages.

## Automatic Deployment (Recommended)

The repository includes a GitHub Actions workflow that automatically builds and deploys the documentation when you push changes to the `main` branch.

### Setup Steps

1. **Enable GitHub Pages in your repository:**
   - Go to your GitHub repository
   - Click **Settings** → **Pages**
   - Under **Source**, select **GitHub Actions**
   - Save the settings

2. **Push your code:**
   ```bash
   git add .
   git commit -m "Add documentation"
   git push origin main
   ```

3. **Wait for deployment:**
   - Go to **Actions** tab in your GitHub repository
   - Wait for the "Build and Deploy Documentation" workflow to complete
   - Your documentation will be available at:
     `https://<your-username>.github.io/Lotus/`

## Manual Deployment

If you prefer to deploy manually:

1. **Build the documentation:**
   ```bash
   cd docs
   make html
   ```

2. **Create a `gh-pages` branch:**
   ```bash
   git checkout -b gh-pages
   git rm -rf .
   cp -r docs/_build/html/* .
   git add .
   git commit -m "Deploy documentation"
   git push origin gh-pages
   ```

3. **Enable GitHub Pages:**
   - Go to **Settings** → **Pages**
   - Select **gh-pages** branch as source
   - Save

## Alternative: Read the Docs

You can also use Read the Docs for hosting:

1. Sign up at https://readthedocs.org/
2. Import your GitHub repository
3. Configure build settings:
   - Python version: 3.10
   - Requirements file: `requirements.txt`
   - Documentation type: Sphinx
   - Configuration file: `docs/conf.py`
4. Build and publish

## Documentation URL

After deployment, your documentation will be available at:
- GitHub Pages: `https://<your-username>.github.io/Lotus/`
- Read the Docs: `https://lotus.readthedocs.io/`


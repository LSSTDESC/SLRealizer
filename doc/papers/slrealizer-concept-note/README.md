
### _SLRealizer Concept DESC Note_
# LSST Catalog-level Realization of Gravitationally-lensed Quasars

*Jenny Kim, Mike Baumer, Phil Marshall, Steve Kahn, Rahul Biswas*

The scale of the LSST dataset will be such that, when considering the
problem of finding lensed quasars, we should anticipate extracting as
much information out of the the catalogs as possible before turning to
the pixel-level data. In this work we explore the use of simple, low
multiplicity Gaussian mixture models for realizing gravitational lens
systems in LSST catalog space, to enable both large-scale sample
simulation and direct model inference.

## Editing this Note

Fork and/or clone the project repo, and then
edit the primary file, `main.tex`, and the figures and tables in their respective sub-folders.

## Building this Note
```
make main
```
You can also `make clean` to remove the intermediate products.

<!--
## Automatic PDF Sharing

If this project is in a public GitHub repo, you can use the `.travis.yml` file in this folder to cause [travis-ci](http://travis-ci.org) to compile your paper into a PDF in the base repo at GitHub every time you push a commit to the master branch. The paper should appear as:

**https://github.com/DarkEnergyScienceCollaboration/SLRealizer/tree/pdf/desc-0000-slrealizer-concept.pdf**

To enable this service, you need to follow these steps:

1. Turn on travis continuous integration, by [toggling your repo on your travis profile](https://travis-ci.org/profile). If you don't see your repo listed, you may not have permission to do this: in this case, [contact an admin via the issues](https://github.com/DarkEnergyScienceCollaboration/SLRealizer/issues/new?body=@DarkEnergyScienceCollaboration/admin).
2. Get a [GitHub "personal access token"](https://github.com/settings/tokens). Choose the "repo" option.
3. Set the `GITHUB_API_KEY` environment variable with the value of this token at your repo's [travis settings page](https://travis-ci.org/DarkEnergyScienceCollaboration/SLRealizer/settings).
4. Copy the `.travis.yml` file in this folder to the top level of your repo (or merge its contents with your existing `.travis.yml` file).
Edit the final `git push` command with your GitHub username.  
Commit and push to trigger your travis build, but note that the PDF will only be deployed if the master branch is updated.

-->

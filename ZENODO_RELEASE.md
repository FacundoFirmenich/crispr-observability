# Archival citation and GitHub to Zenodo release handoff

The [Zenodo v0.8.1 record](https://doi.org/10.5281/zenodo.24121127) contains the versioned source ZIP. The native reproduction release uses a separate release identifier and links to this archived source version for attribution.

Version 0.7.1 is a packaging and citation update to the scientifically checked
0.7.0 release. It adds the 235 historical evidence files to the Git tree itself,
so a repository snapshot includes them rather than relying on separate release
attachments. No scientific count, sequence or algorithm is changed. The versioned
COL17A1 companion retains its original 0.7.0 input manifest and execution receipts.

The 0.7.1 release is published for the author's GitHub-to-Zenodo workflow.
Activate or synchronize FacundoFirmenich/crispr-observability in the Zenodo
GitHub panel and select the current release. No manual ZIP upload is requested.
The integration handles new release events once enabled; if the pre-existing
release is not offered for archiving, a subsequent release event is required.

Manuscript, formal supplement and their editable sources remain excluded until
the author confirms bioRxiv acceptance and authorizes adding them. Cite an archive DOI only for the version deposited under that identifier.

Official instructions:
https://help.zenodo.org/docs/github/enable-repository/
https://help.zenodo.org/docs/github/archive-software/github-upload/

The repository has component-specific licensing. Original code, third-party
interfaces and original-study data retain their respective explicit terms;
the release does not replace them with a blanket manuscript or data license.

One 42 MB source reference is stored in verified parts to satisfy the GitHub API transport limit. All original bytes are retained. Run `python extensions/historical_scientific_evidence/restore_large_files.py --check` to verify them without writing, or omit `--check` to reconstruct the original file before a historical raw-source replay.

# Initial decisions and unresolved items

- Supported expression workflow: strand-matched interval GeneMaskScorer mean;
  historical variant-based pathway retained for an ablation, not exposed by launcher.
- Human/mouse, sequence length and mask width are configuration fields.
- Client pinned to 0.7.0 as in upstream; model_version remains independently explicit.
- Preserve representative and all-replicate mappings; never treat them as identical.
- Both raw and qnorm truth analyses belong in the sensitivity study, with preprocessing verified first.
- Initial paper scope: expression, MPRA, ATAC. eQTL optional.
- No authors, affiliations, quantitative scores, dataset sizes, or performance ranking inferred.
- Need exact historical configs/results to substantiate strand/version/split effects.
- Need desired GitHub owner/repository and visibility before creating a hosted remote.

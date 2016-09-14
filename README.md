# bioentities

Resources for grounding biological entities from text and describing their
hierarchical relationships.

## Gene prefixes

The file ```gene_prefixes.csv``` enumerates prefixes and suffixes frequently
appended to named entities. Some of these represent subtleties of experimental
context (for example, that a protein of interest was tagged with a fluorescent
protein in an experiment) that can safely be ignored when determining the logic
of a sentence. However, others carry essential meaning: for example, a sentence
describing the effect of 'AKT shRNA' on a downstream target has the opposite
meaning of a sentence involving 'AKT', because 'AKT shRNA' represents
inhibition of AKT by genetic silencing.

The patterns included in this file were found by manually reviewing 70,000
named entities extracted by the REACH parser from a corpus of roughly 32,000
papers focused on growth factor signaling.

**Important note: the prefixes/suffixes may be applied additively, for example
```Myr-Flag-Akt1```, indicating myristoylated, FLAG-tagged AKT1; or
```GFP-KRAS-G12V```, indicating GFP-tagged KRAS with a G12V mutation.**

The file contains three columns:

1. A case-sensitive pattern, e.g., ```mEGFP-{Gene name}```, where ```{Gene name}``` represents a protein/gene name.
2. A category, described below.
3. Notes: spelling out acronyms, etc.

The category of the prefix/suffix determines whether it can be stripped off
with minimal effect on the meaning, or whether it carries meaning that needs to
be incorporated by a parser. The categories are as follows:

* ```experimental context```. Protein tags, gene delivery techniques, etc. **Can
  generally be ignored.**

* ```species```. Prefixes denoting human, mouse, primate, or mammalian versions
  of a gene. **In most use cases can be ignored.**

* ```generic descriptor```. Additional words extracted by the entity recognizer
  that might designate that an entity is a "protein", a "protease",
  "transcription factor", etc. **In most use cases can be ignored.**

* ```mrna grounding```. In most cases, entities can be grounded to proteins; in
  the case of ```{Gene name} mRNA```, the entity **must be explicitly grounded
  as an mRNA.**

* ```protein state```. Designate activation state, post-translational
  modification, cellular localization, etc. **Must be captured by the
  parser.**

* ```inhibition```. Designate protein forms or interventions that represent an
  inhibition of the protein, that is, a loss-of-function experiment.  Have the
  effect of switching the polarity of the extracted mechanism. For example, the
  sentence "DUSP6 silencing leads to MAPK1 phosphorylation" indicates that DUSP6
  **inhibits** MAPK1 phosphorylation. **Must be captured by the parser.**


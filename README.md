# Bioentities

*Bioentities* is a collection of resources for grounding biological entities
from text and describing their hierarchical relationships. Resources were
developed by manual curation for use by natural language processing and
biological modeling teams in the [DARPA Big
Mechanism](http://www.darpa.mil/program/big-mechanism) and [Communicating with
Computers](http://www.darpa.mil/program/communicating-with-computers) programs.
The repository contains the following files:

* ```relations.csv```. Defines membership of specific genes/proteins in
  families and protein complexes. For example, ```PIK3CA isa PIK3C```, where
PIK3C represents the class of catalytic subunits of PI3K; and ```PIK3C partof
PI3K```, where PI3K represents a named complex consisting of a catalytic and
regulatory subunit.

* ```entities.csv```. A registry of the families and complexes defined in the
  Bioentities namespace.

* ```grounding_map.csv```. Explicit mapping of text strings to identifiers in
  biological databases.

* ```gene_prefixes.csv```. Patterns of prefixes and suffixes on named entities.

* ```check_references.py```. A script to check the integrity and completeness
  of the cross-references among the various files.

## Entities and Relations

*Bioentities* contains resources for defining the relationships between
genes/proteins and their membership in families and named complexes. Entities
defined within the Bioentities namespace are listed in the ```entities.csv```
file. Cross-referencing the entries among the various files maintains
consistency and prevents errors.

Relationships are defined in ```relations.csv``` as a triples using two
relationships:

* ```isa```, denoting membership in a family;

* ```partof```, denoting membership in a protein complex.

These two relationships can be combined to capture complex hierarchical
relationships, including sub-families (families within families) and complexes
consisting of families of related subunits (e.g., PI3K, NF-kB).

The ```relations.csv``` file consists of five columns: (1) the namespace for
the subject (e.g., ```HGNC``` for gene names, ```UP``` for Uniprot, or ```BE```
for the Bioentities namespace), (2) the identifier for the subject, (3) the
relationship (```isa``` or ```partof```), (4) the namespace for the object, and
(5) the identifier for the object.

## Grounding Map

Using mechanisms extracted from text mining to explain biological datasets
requires that the entities in text are correctly grounded to the canonical
names and IDs of genes, proteins, and chemicals. The problem is that simple
lookups based on string matching often fail, particularly for protein families
and named complexes, which appear frequently in text but lack corresponding
entries in databases.

The grounding map addresses this by providing explicit grounding for frequently
encountered entities in the biological literature. The text strings were drawn
from a corpus of roughly 32,000 papers focused on growth factor signaling in
cancer.

Entities are grounded to the following databases:

* Genes/proteins: [Uniprot](http://www.uniprot.org)

* Chemicals: [PubChem](https://pubchem.ncbi.nlm.nih.gov/),
  [CHEBI](https://www.ebi.ac.uk/chebi/), and [HMDB](http://www.hmdb.ca/) (for
  metabolites)

* Biological processes: [GO](http://geneontology.org/) and
  [MeSH](http://www.ncbi.nlm.nih.gov/mesh)

* Protein families and named complexes: grounded to entities defined within
  the Bioentities repository in the ```entities.csv``` and ```relations.csv```
  files, and to identifiers in [PFAM](http://pfam.xfam.org/)
  and [Interpro](https://www.ebi.ac.uk/interpro/) when possible.

**Note: Some text strings in the map have no grounding.** This was originally
used to identify entities that represent parsing errors and that should not be
included in downstream output. For example, "MAP" a degenerate extraction that
could signify many entities, including MAP kinase, MAP kinase inhibitor, MAP
kinase kinase, etc. However, these empty entries could be used differently
depending on the downstream application.

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

## Contributing

Contributions are welcome! If making additions or revisions to the CSV files
take care to handle quotations and newlines correctly. This allows diffs to be
handled correctly so changes can be reviewed. Please submit updates via pull
requests on Github.

The CSV files in the Bioentities repo are set up to be edited natively using
Microsoft Excel. The CSV files in the repo have Windows line terminators
('\r\n'), and are not ragged (i.e., missing entries in a row are padded out
with empty strings to reach the full width of the longest row).

To preserve correct newlines, take the following steps:

1. If saving from Excel (Windows or Mac OS X), save to the "Windows Comma
   Separated (.csv)" format.

2. If reading (or writing) the files using a Python script, use the following
   set of csv format parameters::

    csvreader = csv.reader(f, delimiter=',', quotechar='"',
                           quoting=csv.QUOTE_MINIMAL, lineterminator='\r\n')

3. If editing the files on Linux, post-process files using ```unix2dos``` or a
   similar program.
